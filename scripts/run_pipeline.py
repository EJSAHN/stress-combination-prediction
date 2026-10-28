#!/usr/bin/env python3
"""Run the offline release with an existing 64-bit Python environment.

This entry point installs nothing. It verifies the release before and after
executing analysis, figures, and workbook stages, and preserves their logs.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib
import json
import os
from pathlib import Path, PurePosixPath, PureWindowsPath
import platform
import re
import shutil
import struct
import subprocess
import sys
import time
import traceback
import uuid

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = ("numpy", "pandas", "matplotlib", "openpyxl", "PIL")
REQUIRED_ENTRY_POINTS = {
    "scripts/run_pipeline.py", "scripts/run_analysis.py",
    "scripts/make_figures.py", "scripts/make_workbook.py",
}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".writing")
    temporary.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def verify_package(package_root):
    """Fail closed on unsafe paths, missing entries, or changed file content."""
    root = Path(package_root).resolve()
    manifest_path = root / "package_manifest.json"
    if not manifest_path.is_file() or manifest_path.is_symlink():
        raise ValueError("A regular package_manifest.json is required; no stage was authorized.")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    entries = manifest.get("files") if isinstance(manifest, dict) else None
    if not isinstance(entries, list) or not entries:
        raise ValueError("Package manifest must contain a nonempty files list.")
    seen, total_bytes = set(), 0
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("Each package manifest entry must be an object.")
        relative = entry.get("relative_path")
        if not isinstance(relative, str) or not relative or "\\" in relative or "\x00" in relative:
            raise ValueError("Manifest paths must be nonempty relative POSIX paths.")
        posix, windows = PurePosixPath(relative), PureWindowsPath(relative)
        if (posix.is_absolute() or windows.drive or windows.root or ":" in relative
                or any(part in {"", ".", ".."} for part in relative.split("/"))):
            raise ValueError("Unsafe manifest path: " + relative)
        normalized = posix.as_posix()
        if normalized.casefold() in seen:
            raise ValueError("Duplicate or case-colliding manifest path: " + relative)
        seen.add(normalized.casefold())
        target = root.joinpath(*posix.parts)
        try:
            target.resolve().relative_to(root)
        except ValueError as exc:
            raise ValueError("Manifest path escapes the package: " + relative) from exc
        # Reject directory symlinks as well as file symlinks for a portable archive.
        current = root
        for part in posix.parts:
            current = current / part
            if current.is_symlink():
                raise ValueError("Symlink is not a portable package input: " + relative)
        if not target.is_file():
            raise ValueError("Missing package file: " + relative)
        expected_hash = entry.get("sha256")
        expected_bytes = entry.get("bytes")
        if not isinstance(expected_hash, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", expected_hash):
            raise ValueError("Invalid SHA256 for " + relative)
        if type(expected_bytes) is not int or expected_bytes < 0:
            raise ValueError("Invalid byte count for " + relative)
        size = target.stat().st_size
        if size != expected_bytes:
            raise ValueError("Package byte count mismatch: " + relative)
        if sha256(target) != expected_hash.lower():
            raise ValueError("Package SHA256 mismatch: " + relative)
        total_bytes += size
    if not REQUIRED_ENTRY_POINTS.issubset(seen):
        missing = sorted(REQUIRED_ENTRY_POINTS - seen)
        raise ValueError("Manifest omits required stage entry points: " + ", ".join(missing))
    return {
        "status": "PASS", "manifest_sha256": sha256(manifest_path),
        "verified_files": len(entries), "verified_file_bytes": total_bytes,
    }


def prepare_environment(output):
    environment = os.environ.copy()
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        environment[key] = "1"
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    environment["PYTHONUTF8"] = "1"
    environment["PYTHONUNBUFFERED"] = "1"
    temporary = output / "temp"
    cache = output / "cache" / "matplotlib"
    temporary.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    for key in ("TMP", "TEMP", "TMPDIR"):
        environment[key] = str(temporary)
    environment["MPLCONFIGDIR"] = str(cache)
    # Imports in this process receive the same cache and threading settings.
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "PYTHONDONTWRITEBYTECODE", "PYTHONUTF8",
                "PYTHONUNBUFFERED", "TMP", "TEMP", "TMPDIR", "MPLCONFIGDIR"):
        os.environ[key] = environment[key]
    sys.dont_write_bytecode = True
    return environment


def dependency_versions():
    if struct.calcsize("P") * 8 != 64:
        raise RuntimeError("A pre-existing 64-bit Python interpreter is required.")
    versions = {}
    for name in DEPENDENCIES:
        try:
            module = importlib.import_module(name)
        except Exception as exc:
            raise RuntimeError(
                "Existing Python cannot import " + name + "; no installation was attempted."
            ) from exc
        versions[name] = str(getattr(module, "__version__", "unknown"))
    return versions


def run_stage(name, arguments, report_path, root, output, environment, records):
    command = [sys.executable, "-u", *map(str, arguments)]
    stage = {
        "name": name, "status": "FAIL", "command": command,
        "started_utc": utc_now(), "return_code": None,
        "log": str((output / "logs" / (name + ".log")).relative_to(output)),
        "validation_report": str(report_path.relative_to(output)),
    }
    records.append(stage)
    started = time.perf_counter()
    print("START_STAGE=" + name, flush=True)
    try:
        with (output / stage["log"]).open("w", encoding="utf-8", newline="\n") as log:
            log.write("COMMAND_JSON=" + json.dumps(command) + "\n")
            log.flush()
            process = subprocess.Popen(
                command, cwd=str(root), env=environment, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
                bufsize=1,
            )
            try:
                assert process.stdout is not None
                for line in process.stdout:
                    log.write(line)
                    log.flush()
                    print(line, end="", flush=True)
                stage["return_code"] = process.wait()
            except BaseException:
                if process.poll() is None:
                    process.terminate()
                    try:
                        process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
                stage["return_code"] = process.returncode
                raise
            finally:
                if process.stdout is not None:
                    process.stdout.close()
        if stage["return_code"] != 0:
            raise RuntimeError(name + " failed with exit code " + str(stage["return_code"]))
        if not report_path.is_file():
            raise RuntimeError(name + " did not produce its required validation report.")
        report = json.loads(report_path.read_text(encoding="utf-8-sig"))
        if not isinstance(report, dict) or report.get("status") != "PASS":
            raise RuntimeError(name + " validation report is not PASS.")
        stage["validation_report_sha256"] = sha256(report_path)
        stage["validation_status"] = report["status"]
        stage["status"] = "PASS"
        print("PASS_STAGE=" + name, flush=True)
    except BaseException as exc:
        stage["error"] = str(exc)
        raise
    finally:
        stage["runtime_seconds"] = round(time.perf_counter() - started, 6)
        stage["finished_utc"] = utc_now()


def failure_output(requested):
    """Preserve an occupied output; put rejection diagnostics beside it."""
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    path = requested.with_name(requested.name + "_rejected_" + stamp + "_" + uuid.uuid4().hex[:8])
    return path


def run_pipeline(package_root, requested_output):
    root = Path(package_root).resolve()
    requested = Path(requested_output).expanduser().resolve()
    output = requested
    started = time.perf_counter()
    report = {
        "status": "FAIL", "started_utc": utc_now(), "package_root": str(root),
        "requested_output": str(requested), "output": str(output),
        "python_executable": sys.executable, "python_version": platform.python_version(),
        "platform": platform.platform(), "interpreter_bits": struct.calcsize("P") * 8,
        "native_windows_launcher_executed": (
            os.name == "nt" and os.environ.get("STRESS_PIPELINE_LAUNCHER") == "powershell"
        ),
        "environment_creation": False, "package_installation": False,
        "network_access_requested": False,
        "source_inputs": "Bundled author-processed expression matrices; no FASTQ processing.",
        "stages": [], "errors": [],
    }
    error_details = []
    package_before = None
    try:
        if requested.exists() and (not requested.is_dir() or any(requested.iterdir())):
            output = failure_output(requested)
            report["output"] = str(output)
            report["requested_output_preserved"] = True
            output.mkdir(parents=True, exist_ok=False)
            raise FileExistsError("Refusing nonempty or non-directory --output: " + str(requested))
        output.mkdir(parents=True, exist_ok=True)
        (output / "logs").mkdir(exist_ok=True)
        package_before = verify_package(root)
        report["package_before"] = package_before
        write_json(output / "logs" / "package_verification_before.json", package_before)
        shutil.copyfile(root / "package_manifest.json", output / "logs" / "package_manifest_before.json")
        environment = prepare_environment(output)
        report["dependency_versions"] = dependency_versions()
        for folder in ("analysis", "figures", "tables"):
            (output / folder).mkdir(exist_ok=True)
        stages = [
            ("analysis", ["scripts/run_analysis.py", "--output", output / "analysis"],
             output / "analysis" / "validation" / "analysis_validation.json"),
            ("figures", ["scripts/make_figures.py", "--results", output / "analysis", "--output", output / "figures"],
             output / "figures" / "figure_validation.json"),
            ("workbook", ["scripts/make_workbook.py", "--results", output / "analysis", "--output", output / "tables" / "Results.xlsx", "--figure-sources", output / "figures" / "source_data"],
             output / "tables" / "workbook_validation.json"),
        ]
        for name, arguments, validation_path in stages:
            run_stage(name, arguments, validation_path, root, output, environment, report["stages"])
    except BaseException as exc:
        report["errors"].append(type(exc).__name__ + ": " + str(exc))
        error_details.append(traceback.format_exc())
    finally:
        if package_before is not None:
            try:
                after = verify_package(root)
                report["package_after"] = after
                write_json(output / "logs" / "package_verification_after.json", after)
                shutil.copyfile(root / "package_manifest.json", output / "logs" / "package_manifest_after.json")
                if after["manifest_sha256"] != package_before["manifest_sha256"]:
                    raise RuntimeError("The package manifest changed during execution.")
                report["package_unchanged"] = True
            except BaseException as exc:
                report["package_unchanged"] = False
                report["package_after"] = {
                    "status": "FAIL", "error": type(exc).__name__ + ": " + str(exc),
                }
                try:
                    if (root / "package_manifest.json").is_file():
                        report["package_after"]["manifest_sha256"] = sha256(root / "package_manifest.json")
                        shutil.copyfile(root / "package_manifest.json", output / "logs" / "package_manifest_after.json")
                    write_json(output / "logs" / "package_verification_after.json", report["package_after"])
                except OSError:
                    pass
                report["errors"].append("Post-run verification: " + type(exc).__name__ + ": " + str(exc))
                error_details.append(traceback.format_exc())
        if not report["errors"] and len(report["stages"]) == 3 and all(
                stage["status"] == "PASS" for stage in report["stages"]):
            report["status"] = "PASS"
        report["finished_utc"] = utc_now()
        report["runtime_seconds"] = round(time.perf_counter() - started, 6)
        try:
            (output / "logs").mkdir(parents=True, exist_ok=True)
            (output / "logs" / "error.log").write_text(
                "\n".join(error_details) if error_details else "No errors.\n", encoding="utf-8",
            )
            write_json(output / "validation.json", report)
        except OSError as exc:
            # A filesystem that refuses writes cannot retain a local report.
            report["status"] = "FAIL"
            report["errors"].append("Unable to write diagnostics: " + str(exc))
        print(json.dumps({
            "status": report["status"], "validation_json": str(output / "validation.json"),
            "runtime_seconds": report["runtime_seconds"], "errors": report["errors"],
        }, ensure_ascii=True), flush=True)
    return 0 if report["status"] == "PASS" else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True,
                        help="Fresh output directory. Existing nonempty output is preserved and refused.")
    args = parser.parse_args()
    return run_pipeline(PACKAGE_ROOT, args.output)


if __name__ == "__main__":
    raise SystemExit(main())
