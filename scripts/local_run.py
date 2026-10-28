#!/usr/bin/env python3
"""Execute an existing Python pipeline and retain its local validation record.

No environment creation, package installation, download or old-run deletion.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import platform
import re
import subprocess
import sys
import time
import zipfile

def utc_now():
    return datetime.now(timezone.utc).isoformat()

def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()

def write_json(path, payload):
    Path(path).write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')

def contained(path, parent):
    try:
        Path(path).resolve().relative_to(Path(parent).resolve())
        return True
    except ValueError:
        return False

def verify_package(package_root):
    package_root = Path(package_root)
    manifest_path = package_root / 'package_manifest.json'
    manifest = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    entries = manifest.get('files')
    if not isinstance(entries, list) or not entries:
        raise ValueError('The package manifest has no file entries.')
    seen = set()
    total_bytes = 0
    for entry in entries:
        relative = entry.get('relative_path')
        if not isinstance(relative, str) or not relative or '\\' in relative:
            raise ValueError('Manifest paths must be nonempty relative POSIX paths.')
        parsed = PurePosixPath(relative)
        if parsed.is_absolute() or '..' in parsed.parts or ':' in relative or parsed.as_posix() != relative:
            raise ValueError('Unsafe manifest path: ' + relative)
        if relative.casefold() in seen:
            raise ValueError('Duplicate or case-colliding manifest path: ' + relative)
        seen.add(relative.casefold())
        target = package_root.joinpath(*parsed.parts)
        if not contained(target, package_root) or not target.is_file() or target.is_symlink():
            raise ValueError('Missing or unsafe package file: ' + relative)
        expected = entry.get('sha256', '')
        if not isinstance(expected, str) or not re.fullmatch(r'[0-9a-fA-F]{64}', expected):
            raise ValueError('Invalid SHA256 for ' + relative)
        size = target.stat().st_size
        if size != entry.get('bytes'):
            raise ValueError('Size mismatch for ' + relative)
        if sha256(target) != expected.lower():
            raise ValueError('SHA256 mismatch for ' + relative)
        total_bytes += size
    required = {'scripts/run_pipeline.py', 'scripts/local_run.py', 'run.ps1'}
    if not required.issubset(seen):
        raise ValueError('The package manifest omits a required entry point.')
    return {'status': 'PASS', 'manifest_sha256': sha256(manifest_path),
            'verified_files': len(entries), 'verified_file_bytes': total_bytes}

EXCLUDED_DIRECTORY_NAMES = {'input', 'inputs', 'data', 'cache', 'tmp', 'temp', '.cache', 'workspace', 'work', 'working', 'reference', '__pycache__'}

def make_result_bundle(run_root, destination, record):
    """Archive generated reports and diagnostics, excluding source data/caches."""
    run_root, destination = Path(run_root), Path(destination)
    if destination.exists():
        raise FileExistsError('Refusing to replace an existing result archive: ' + str(destination))
    destination.parent.mkdir(parents=True, exist_ok=True)
    selected, excluded = [], []
    for path in sorted(run_root.rglob('*')):
        if not path.is_file() or path.is_symlink():
            continue
        relative = path.relative_to(run_root)
        reason = None
        if any(part.casefold() in EXCLUDED_DIRECTORY_NAMES for part in relative.parts):
            reason = 'input_or_reproducible_cache_directory'
        elif path.suffix.lower() in {'.npz', '.npy', '.pyc', '.zip'}:
            reason = 'array_cache_or_nested_archive'
        elif path.name.startswith('GSE') and path.suffix.lower() == '.gz':
            reason = 'original_expression_matrix'
        elif path.name.endswith('.csv.gz') and path.stat().st_size > 5 * 1024 * 1024:
            reason = 'large_generated_table_preserved_in_run_directory'
        if reason:
            excluded.append({'path': relative.as_posix(), 'bytes': path.stat().st_size, 'reason': reason})
        else:
            selected.append((path, {'path': relative.as_posix(), 'bytes': path.stat().st_size, 'sha256': sha256(path)}))
    inventory = {'created_utc': utc_now(), 'run_status': record.get('status'),
                 'included_files': [item for _, item in selected], 'excluded_files': excluded,
                 'original_input_data_included': False, 'generated_arrays_included': False,
                 'generated_compressed_table_limit_bytes': 5 * 1024 * 1024}
    partial = destination.with_suffix(destination.suffix + '.partial')
    if partial.exists():
        raise FileExistsError('Refusing to replace an existing partial archive: ' + str(partial))
    with zipfile.ZipFile(partial, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for path, item in selected:
            archive.write(path, item['path'])
        archive.writestr('bundle_inventory.json', json.dumps(inventory, indent=2, ensure_ascii=False) + '\n')
    with zipfile.ZipFile(partial) as archive:
        if archive.testzip() is not None:
            raise ValueError('Result archive CRC verification failed.')
        for _, item in selected:
            if hashlib.sha256(archive.read(item['path'])).hexdigest() != item['sha256']:
                raise ValueError('Result archive content verification failed: ' + item['path'])
    partial.replace(destination)
    return {'path': str(destination), 'bytes': destination.stat().st_size,
            'sha256': sha256(destination), 'included_files': len(selected) + 1,
            'excluded_files': len(excluded), 'CRC_and_content_hashes_verified': True}

def run_local(package_root, project_root, run_root, result_directory):
    package_root = Path(package_root).resolve()
    project_root = Path(project_root).resolve()
    run_root = Path(run_root).resolve()
    result_directory = Path(result_directory).resolve()
    if not contained(run_root, project_root) or run_root == project_root:
        raise ValueError('The run directory must be a child of the project directory.')
    run_root.mkdir(parents=True, exist_ok=True)
    log_dir = run_root / 'logs'
    log_dir.mkdir(exist_ok=True)
    output = run_root / 'results'
    status_path = run_root / 'local_execution.json'
    if status_path.exists() or output.exists():
        raise FileExistsError('This run already has an execution record or results; use a fresh run directory.')
    start = time.perf_counter()
    record = {'status': 'FAIL', 'started_utc': utc_now(), 'python': platform.python_version(),
              'python_executable': sys.executable, 'platform': platform.platform(),
              'native_windows_execution': platform.system() == 'Windows',
              'powershell_launcher_confirmed': platform.system() == 'Windows' and os.environ.get('STRESS_PIPELINE_LAUNCHER') == 'powershell',
              'project_root': str(project_root), 'package_root': str(package_root),
              'run_root': str(run_root), 'environment_creation': False,
              'package_installation': False, 'network_access_requested': False,
              'pipeline_exit_code': None, 'pipeline_validation_status': None, 'error': None}
    destination = result_directory / ('RESULT_STRESS_COMBINATION_VALIDATION_' + run_root.name + '.zip')
    exit_code = 1
    try:
        verified = verify_package(package_root)
        record['package_verification'] = verified
        write_json(log_dir / 'package_verification.json', verified)
        print('PACKAGE_SHA256_VERIFIED_FILES=' + str(verified['verified_files']), flush=True)
        environment = os.environ.copy()
        for key, folder in [('TEMP', 'tmp'), ('TMP', 'tmp'), ('MPLCONFIGDIR', 'cache/matplotlib'), ('XDG_CACHE_HOME', 'cache')]:
            directory = run_root / folder
            directory.mkdir(parents=True, exist_ok=True)
            environment[key] = str(directory)
        for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
            environment[key] = '1'
        environment['PYTHONDONTWRITEBYTECODE'] = '1'
        environment['PYTHONUNBUFFERED'] = '1'
        environment['PYTHONIOENCODING'] = 'utf-8'
        command = [sys.executable, '-B', '-u', str(package_root / 'scripts/run_pipeline.py'), '--output', str(output)]
        record['pipeline_command'] = command
        record['cache_and_temporary_paths'] = {key: environment[key] for key in ('TEMP', 'TMP', 'MPLCONFIGDIR', 'XDG_CACHE_HOME')}
        write_json(log_dir / 'local_invocation.json', record)
        print('START_FULL_PIPELINE', flush=True)
        with (log_dir / 'pipeline.log').open('w', encoding='utf-8', newline='\n') as log:
            with subprocess.Popen(command, cwd=package_root, env=environment,
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  text=True, encoding='utf-8', errors='replace', bufsize=1) as process:
                assert process.stdout is not None
                for line in process.stdout:
                    log.write(line)
                    log.flush()
                    print(line.rstrip('\n'), flush=True)
                record['pipeline_exit_code'] = process.wait()
        if record['pipeline_exit_code'] != 0:
            raise RuntimeError('Pipeline exited with code ' + str(record['pipeline_exit_code']))
        validation_path = output / 'validation.json'
        if not validation_path.is_file():
            raise RuntimeError('Pipeline did not produce results/validation.json.')
        validation = json.loads(validation_path.read_text(encoding='utf-8-sig'))
        record['pipeline_validation_status'] = validation.get('status')
        record['pipeline_validation_sha256'] = sha256(validation_path)
        if validation.get('status') != 'PASS':
            raise RuntimeError('Pipeline validation did not report PASS.')
        record['status'] = 'PASS'
        exit_code = 0
    except Exception as error:
        record['error'] = type(error).__name__ + ': ' + str(error)
        print('FAILED=' + record['error'], flush=True)
    finally:
        record['completed_utc'] = utc_now()
        record['runtime_seconds'] = time.perf_counter() - start
        write_json(status_path, record)
    try:
        receipt = make_result_bundle(run_root, destination, record)
        write_json(run_root / 'bundle_receipt.json', receipt)
        print('UPLOAD_THIS_ZIP=' + str(destination), flush=True)
        print('RESULT_ZIP_SHA256=' + receipt['sha256'], flush=True)
    except Exception as error:
        exit_code = 1
        record['bundle_error'] = type(error).__name__ + ': ' + str(error)
        record['status'] = 'FAIL'
        write_json(status_path, record)
        print('BUNDLE_FAILED=' + record['bundle_error'], flush=True)
    print('VALIDATION_STATUS=' + record['status'], flush=True)
    print('NATIVE_WINDOWS_EXECUTION=' + str(record['native_windows_execution']).upper(), flush=True)
    print('EXIT_CODE=' + str(exit_code), flush=True)
    return exit_code

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', required=True, type=Path)
    parser.add_argument('--run-root', required=True, type=Path)
    parser.add_argument('--result-directory', required=True, type=Path)
    args = parser.parse_args()
    try:
        return run_local(Path(__file__).resolve().parents[1], args.project_root,
                         args.run_root, args.result_directory)
    except Exception as error:
        print('LAUNCH_FAILED=' + type(error).__name__ + ': ' + str(error), flush=True)
        return 1

if __name__ == '__main__':
    raise SystemExit(main())
