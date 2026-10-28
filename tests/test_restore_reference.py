"""Restoration checks: corruption, safe paths, preservation and real bytes."""
from pathlib import Path
import hashlib
import importlib.util
import json
import shutil
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("restore_reference", ROOT / "scripts/restore_reference.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class RestoreReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="reference_restoration_")
        self.root = Path(self.temp.name)
        self.content = b"original scientific reference bytes\x00\xff" * 13
        self.target = self.root / "provenance/frozen/reference.npz"
        parts = []
        for number, data in enumerate((self.content[:100], self.content[100:])):
            relative = "provenance/chunks/part-%02d.bin" % number
            path = self.root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
            parts.append({"path": relative, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
        self.asset = {"path": self.target.relative_to(self.root).as_posix(), "bytes": len(self.content),
                      "sha256": hashlib.sha256(self.content).hexdigest(), "parts": parts}
        self.write_manifest()

    def tearDown(self):
        self.temp.cleanup()

    def write_manifest(self):
        p = self.root / "config/chunked_assets.json"
        p.parent.mkdir(exist_ok=True)
        p.write_text(json.dumps({"schema_version": 1, "assets": [self.asset]}))

    def test_exact_restore_idempotence_and_check_only(self):
        self.assertEqual(MODULE.restore_assets(self.root)[0]["status"], "RESTORED")
        self.assertEqual(self.target.read_bytes(), self.content)
        stamp = self.target.stat().st_mtime_ns
        self.assertEqual(MODULE.restore_assets(self.root)[0]["status"], "ALREADY_VERIFIED")
        MODULE.restore_assets(self.root, check_only=True)
        self.assertEqual(self.target.stat().st_mtime_ns, stamp)
        self.assertEqual(list(self.target.parent.glob("*.partial")), [])

    def test_corrupted_part_rejected_without_output(self):
        (self.root / self.asset["parts"][0]["path"]).write_bytes(b"broken")
        with self.assertRaises(ValueError): MODULE.restore_assets(self.root)
        self.assertFalse(self.target.exists())

    def test_missing_part_rejected_without_output(self):
        (self.root / self.asset["parts"][1]["path"]).unlink()
        with self.assertRaises(FileNotFoundError): MODULE.restore_assets(self.root)
        self.assertFalse(self.target.exists())

    def test_existing_mismatched_file_is_preserved(self):
        self.target.parent.mkdir(parents=True)
        self.target.write_bytes(b"keep this file")
        with self.assertRaises(FileExistsError): MODULE.restore_assets(self.root)
        self.assertEqual(self.target.read_bytes(), b"keep this file")

    def test_bad_whole_hash_cleans_partial(self):
        self.asset["sha256"] = "0" * 64
        self.write_manifest()
        with self.assertRaises(ValueError): MODULE.restore_assets(self.root)
        self.assertFalse(self.target.exists())
        self.assertEqual(list(self.target.parent.glob("*.partial")), [])

    def test_check_only_does_not_create_missing_asset(self):
        with self.assertRaises(FileNotFoundError): MODULE.restore_assets(self.root, check_only=True)
        self.assertFalse(self.target.exists())

    def test_traversal_is_rejected(self):
        self.asset["path"] = "../outside.npz"
        self.write_manifest()
        with self.assertRaises(ValueError): MODULE.restore_assets(self.root)

    def test_publish_failure_cleans_partial(self):
        with mock.patch.object(MODULE.os, "link", side_effect=OSError("simulated publication failure")):
            with self.assertRaises(OSError): MODULE.restore_assets(self.root)
        self.assertFalse(self.target.exists())
        self.assertEqual(list(self.target.parent.glob("*.partial")), [])

    def test_production_parts_restore_exact_reference(self):
        config = json.loads((ROOT / "config/chunked_assets.json").read_text())
        with tempfile.TemporaryDirectory(prefix="production_reference_restore_") as temp:
            root = Path(temp)
            (root / "config").mkdir()
            (root / "config/chunked_assets.json").write_text(json.dumps(config))
            for asset in config["assets"]:
                for part in asset["parts"]:
                    dest = root / part["path"]
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(ROOT / part["path"], dest)
            results = MODULE.restore_assets(root)
            for asset, result in zip(config["assets"], results):
                self.assertEqual(result["status"], "RESTORED")
                self.assertEqual(MODULE.file_record(root / asset["path"]), (asset["bytes"], asset["sha256"]))
                if (ROOT / asset["path"]).exists():
                    self.assertEqual((root / asset["path"]).read_bytes(), (ROOT / asset["path"]).read_bytes())


if __name__ == "__main__":
    unittest.main()
