"""Local launcher tests; PowerShell checks are explicitly static only."""
import contextlib
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('local_run', ROOT / 'scripts/local_run.py')
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)

class LocalRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='stress_local_runner_')
        self.base = Path(self.temporary.name)
        self.package = self.base / 'package with spaces'
        self.package.mkdir()
        (self.package / 'scripts').mkdir()
        (self.package / 'data').mkdir()
        shutil.copy2(ROOT / 'scripts/local_run.py', self.package / 'scripts/local_run.py')
        shutil.copy2(ROOT / 'run.ps1', self.package / 'run.ps1')
        (self.package / 'data/pinned.txt').write_text('original', encoding='utf-8')
        self.project = self.base / 'project with spaces'
        self.run = self.project / 'runs' / 'sample run'
        self.downloads = self.base / 'Downloads'
        self.make_pipeline('PASS', 0)

    def tearDown(self):
        self.temporary.cleanup()

    def make_pipeline(self, status='PASS', exit_code=0, write_validation=True):
        script = (
            'import argparse,json,os\nfrom pathlib import Path\n'
            'p=argparse.ArgumentParser();p.add_argument("--output",required=True);a=p.parse_args()\n'
            'out=Path(a.output);out.mkdir(parents=True)\n'
            '(out/"tables").mkdir();(out/"tables"/"values.csv").write_text("x,y\\n1,2\\n")\n'
            '(out/"analysis"/"working").mkdir(parents=True)\n'
            '(out/"analysis"/"working"/"raw.txt").write_text("do not archive")\n'
            '(out/"binary.npz").write_bytes(b"array")\n'
            'print("mock pipeline executed",flush=True)\n')
        if write_validation:
            script += '(out/"validation.json").write_text(json.dumps({"status":' + repr(status) + ',"temp":os.environ.get("TEMP"),"mpl":os.environ.get("MPLCONFIGDIR")}))\n'
        script += 'raise SystemExit(' + str(exit_code) + ')\n'
        (self.package / 'scripts/run_pipeline.py').write_text(script, encoding='utf-8')
        self.refresh_manifest()

    def refresh_manifest(self):
        items = []
        for path in sorted(self.package.rglob('*')):
            if path.is_file() and path.name != 'package_manifest.json':
                items.append({'relative_path': path.relative_to(self.package).as_posix(), 'sha256': RUNNER.sha256(path), 'bytes': path.stat().st_size})
        RUNNER.write_json(self.package / 'package_manifest.json', {'files': items})

    def execute(self):
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            result = RUNNER.run_local(self.package, self.project, self.run, self.downloads)
        self.output_text = stream.getvalue()
        return result

    def bundle(self):
        return self.downloads / ('RESULT_STRESS_COMBINATION_VALIDATION_' + self.run.name + '.zip')

    def test_success_with_space_paths_and_cache_routing(self):
        prior_temp = os.environ.get('TEMP')
        self.assertEqual(self.execute(), 0)
        record = json.loads((self.run / 'local_execution.json').read_text())
        self.assertEqual(record['status'], 'PASS')
        self.assertEqual(record['pipeline_exit_code'], 0)
        self.assertEqual(record['pipeline_validation_status'], 'PASS')
        self.assertEqual(os.environ.get('TEMP'), prior_temp)
        validation = json.loads((self.run / 'results/validation.json').read_text())
        self.assertTrue(RUNNER.contained(validation['temp'], self.run))
        self.assertTrue(RUNNER.contained(validation['mpl'], self.run))
        self.assertEqual((self.package / 'data/pinned.txt').read_text(), 'original')
        with zipfile.ZipFile(self.bundle()) as archive:
            names = archive.namelist()
            self.assertIn('results/tables/values.csv', names)
            self.assertIn('local_execution.json', names)
            self.assertNotIn('results/analysis/working/raw.txt', names)
            self.assertNotIn('results/binary.npz', names)
            self.assertIsNone(archive.testzip())
            inventory = json.loads(archive.read('bundle_inventory.json'))
            for item in inventory['included_files']:
                self.assertEqual(hashlib.sha256(archive.read(item['path'])).hexdigest(), item['sha256'])

    def test_nonzero_pipeline_exit_preserves_failure_bundle(self):
        self.make_pipeline('PASS', 7)
        self.assertEqual(self.execute(), 1)
        record = json.loads((self.run / 'local_execution.json').read_text())
        self.assertEqual(record['pipeline_exit_code'], 7)
        self.assertEqual(record['status'], 'FAIL')
        self.assertTrue(self.bundle().is_file())

    def test_exit_zero_but_failed_validation_is_failure(self):
        self.make_pipeline('FAIL', 0)
        self.assertEqual(self.execute(), 1)
        record = json.loads((self.run / 'local_execution.json').read_text())
        self.assertEqual(record['pipeline_validation_status'], 'FAIL')
        self.assertTrue(self.bundle().is_file())

    def test_missing_validation_is_failure(self):
        self.make_pipeline(write_validation=False)
        self.assertEqual(self.execute(), 1)
        self.assertIn('did not produce', self.output_text)
        self.assertTrue(self.bundle().is_file())

    def test_corrupt_input_stops_before_pipeline(self):
        (self.package / 'data/pinned.txt').write_text('tampered', encoding='utf-8')
        self.assertEqual(self.execute(), 1)
        self.assertFalse((self.run / 'results').exists())
        self.assertIn('SHA256 mismatch', self.output_text)
        self.assertTrue(self.bundle().is_file())

    def test_unsafe_manifest_path_rejected(self):
        RUNNER.write_json(self.package / 'package_manifest.json', {'files': [{'relative_path': '../outside', 'bytes': 0, 'sha256': '0' * 64}]})
        with self.assertRaisesRegex(ValueError, 'Unsafe manifest path'):
            RUNNER.verify_package(self.package)

    def test_large_generated_table_is_preserved_but_not_bundled(self):
        self.run.mkdir(parents=True)
        large = self.run / 'large_genes.csv.gz'
        large.write_bytes(b'x' * (5 * 1024 * 1024 + 1))
        small = self.run / 'summary.csv'
        small.write_text('x\n1\n')
        receipt = RUNNER.make_result_bundle(self.run, self.bundle(), {'status': 'PASS'})
        self.assertEqual(receipt['excluded_files'], 1)
        self.assertTrue(large.is_file())
        with zipfile.ZipFile(self.bundle()) as archive:
            self.assertIn('summary.csv', archive.namelist())
            self.assertNotIn('large_genes.csv.gz', archive.namelist())

    def test_case_colliding_manifest_paths_rejected(self):
        path = self.package / 'package_manifest.json'
        manifest = json.loads(path.read_text())
        duplicate = dict(manifest['files'][0])
        duplicate['relative_path'] = duplicate['relative_path'].upper()
        manifest['files'].append(duplicate)
        RUNNER.write_json(path, manifest)
        with self.assertRaisesRegex(ValueError, 'case-colliding'):
            RUNNER.verify_package(self.package)

    def test_existing_run_is_never_overwritten(self):
        self.run.mkdir(parents=True)
        (self.run / 'results').mkdir()
        sentinel = self.run / 'results/keep.txt'
        sentinel.write_text('keep')
        with self.assertRaises(FileExistsError):
            self.execute()
        self.assertEqual(sentinel.read_text(), 'keep')

    def test_existing_result_archive_is_never_overwritten(self):
        self.downloads.mkdir()
        self.bundle().write_bytes(b'keep previous archive')
        self.assertEqual(self.execute(), 1)
        self.assertEqual(self.bundle().read_bytes(), b'keep previous archive')
        record = json.loads((self.run / 'local_execution.json').read_text())
        self.assertEqual(record['status'], 'FAIL')
        self.assertIn('Refusing to replace', record['bundle_error'])

    def test_run_path_cannot_escape_project(self):
        with self.assertRaisesRegex(ValueError, 'child of the project'):
            RUNNER.run_local(self.package, self.project, self.base / 'outside', self.downloads)

    def test_powershell_contract_static_only(self):
        text = (ROOT / 'run.ps1').read_text(encoding='utf-8')
        bootstrap = (ROOT / 'scripts/bootstrap_windows.ps1').read_text(encoding='utf-8')
        for marker in ('STRESS_PIPELINE_LAUNCHER', 'savedEnvironment', 'PYTHONDONTWRITEBYTECODE', 'MPLCONFIGDIR', 'ConvertTo-Json', 'Compress-Archive'):
            self.assertIn(marker, text)
        for unsupported in ('??', '&&', '||', 'ForEach-Object -Parallel'):
            self.assertNotIn(unsupported, text)
        for destructive in ('Remove-Item', 'pip install', 'conda create', '-m venv', 'Invoke-WebRequest', 'Invoke-RestMethod'):
            self.assertNotIn(destructive, text)
        self.assertIn("'anaconda3\\python.exe'", text)
        self.assertNotIn('19292', text)
        self.assertIn('STRESS_COMBINATION_LOCAL_v1*.zip', bootstrap)
        self.assertIn('stress-combination-prediction\\run.ps1', bootstrap)

if __name__ == '__main__':
    unittest.main()
