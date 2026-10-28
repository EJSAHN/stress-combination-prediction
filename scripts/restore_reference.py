#!/usr/bin/env python3
"""Restore byte-identical reference assets from bundled verified parts.

Uses only the standard library. No network, installation or replacement of
existing mismatched files is permitted.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import tempfile

BLOCK_BYTES = 1024 * 1024


def safe_path(root: Path, value: str) -> Path:
    if not isinstance(value, str) or not value or "\\" in value or ":" in value:
        raise ValueError("Expected a canonical repository-relative path")
    relative = PurePosixPath(value)
    if not relative.parts or relative.is_absolute() or ".." in relative.parts or relative.as_posix() != value:
        raise ValueError("Unsafe asset path: " + value)
    path = root.joinpath(*relative.parts)
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError("Asset path escapes the repository: " + value)
    current = path
    while current != root:
        if current.is_symlink():
            raise ValueError("Symbolic links are not accepted for asset paths: " + value)
        current = current.parent
    return path


def file_record(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(BLOCK_BYTES), b""):
            size += len(block)
            digest.update(block)
    return size, digest.hexdigest()


def validate_record(record: dict) -> None:
    if not isinstance(record, dict):
        raise ValueError("Each asset and part must be an object")
    if type(record.get("bytes")) is not int or record["bytes"] < 0:
        raise ValueError("Invalid recorded byte count")
    if not isinstance(record.get("sha256"), str) or not re.fullmatch(r"[0-9a-f]{64}", record["sha256"]):
        raise ValueError("Invalid recorded SHA256")


def verify_file(path: Path, record: dict) -> None:
    if not path.is_file():
        raise FileNotFoundError("Required file is missing: " + str(path))
    if file_record(path) != (record["bytes"], record["sha256"]):
        raise ValueError("Byte count or SHA256 mismatch: " + str(path))


def load_assets(root: Path) -> list[dict]:
    manifest = json.loads(safe_path(root, "config/chunked_assets.json").read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1 or not isinstance(manifest.get("assets"), list) or not manifest["assets"]:
        raise ValueError("Unsupported or empty chunk manifest")
    seen: set[str] = set()
    assets = manifest["assets"]
    for asset in assets:
        validate_record(asset)
        safe_path(root, asset.get("path"))
        if asset["path"] in seen:
            raise ValueError("Duplicate asset path")
        seen.add(asset["path"])
        if not isinstance(asset.get("parts"), list) or not asset["parts"]:
            raise ValueError("Asset has no parts")
        for part in asset["parts"]:
            validate_record(part)
            safe_path(root, part.get("path"))
            if part["path"] in seen:
                raise ValueError("Duplicate or overlapping asset/part path")
            seen.add(part["path"])
        if sum(p["bytes"] for p in asset["parts"]) != asset["bytes"]:
            raise ValueError("Part byte counts do not sum to the asset size")
    return assets


def restore_asset(root: Path, asset: dict, check_only: bool = False) -> dict:
    target = safe_path(root, asset["path"])
    for part in asset["parts"]:
        verify_file(safe_path(root, part["path"]), part)
    if target.exists():
        try:
            verify_file(target, asset)
        except (ValueError, FileNotFoundError) as error:
            raise FileExistsError("Refusing to overwrite an existing mismatched asset: " + str(target)) from error
        return {"path": asset["path"], "status": "ALREADY_VERIFIED", "bytes": asset["bytes"], "sha256": asset["sha256"]}
    if check_only:
        raise FileNotFoundError("Asset is not restored. Run without --check-only: " + asset["path"])
    target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=target.name + ".", suffix=".partial", dir=target.parent)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, "wb") as output:
            combined = hashlib.sha256()
            written = 0
            for part in asset["parts"]:
                part_hash = hashlib.sha256()
                part_size = 0
                with safe_path(root, part["path"]).open("rb") as source:
                    for block in iter(lambda: source.read(BLOCK_BYTES), b""):
                        output.write(block)
                        combined.update(block)
                        part_hash.update(block)
                        written += len(block)
                        part_size += len(block)
                if (part_size, part_hash.hexdigest()) != (part["bytes"], part["sha256"]):
                    raise ValueError("Part changed during assembly: " + part["path"])
            output.flush()
            os.fsync(output.fileno())
        if (written, combined.hexdigest()) != (asset["bytes"], asset["sha256"]):
            raise ValueError("Assembled asset byte count or SHA256 mismatch")
        # The handle is closed before publication. A same-directory hard link
        # publishes atomically without overwriting a concurrently created file.
        # Windows NTFS supports this operation.
        try:
            os.link(temporary, target)
        except FileExistsError:
            verify_file(target, asset)
        verify_file(target, asset)
        return {"path": asset["path"], "status": "RESTORED", "bytes": written, "sha256": combined.hexdigest()}
    finally:
        temporary.unlink(missing_ok=True)


def restore_assets(root: Path, check_only: bool = False) -> list[dict]:
    root = Path(root).resolve()
    return [restore_asset(root, asset, check_only=check_only) for asset in load_assets(root)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1], help="Repository root")
    parser.add_argument("--check-only", action="store_true", help="Verify parts and restored assets without writing")
    args = parser.parse_args()
    records = restore_assets(args.root, check_only=args.check_only)
    print(json.dumps({"status": "PASS", "network_access": False, "installation": False, "assets": records}, indent=2))


if __name__ == "__main__":
    main()
