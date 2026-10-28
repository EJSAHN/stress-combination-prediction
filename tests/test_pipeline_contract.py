"""Offline orchestration contracts using tiny synthetic stages, not plant data."""
from __future__ import annotations
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

RELEASE = Path(__file__).resolve().parents[1]
ORCHESTRATOR = RELEASE / 'scripts' / 'run_pipeline.py'
spec = importlib.util.spec_from_file_location('release_orchestrator', ORCHESTRATOR)
pipeline = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pipeline)


class PipelineContractTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='stress_pipeline_contract_')
        self.base = Path(self.temporary.name)
        self.package = self.base / 'package'
        (self.package / 'scripts').mkdir(parents=True)
        shutil.copyfile(ORCHESTRATOR, self.package / 'scripts' / 'run_pipeline.py')
        self.output = self.base / 'output'
        self.stage_status = {}
        self.mutate_package = False

    def tearDown(self):
        self.temporary.cleanup()

    def build(self):
        for stage, filename in (('analysis', 'run_analysis.py'),
                                ('figures', 'make_figures.py'),
                                ('workbook', 'make_workbook.py')):
            code = """import argparse,json,os
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--results');p.add_argument('--figure-sources');a=p.parse_args()
out=Path(a.output)
stage=STAGE_LITERAL
if stage=='workbook':
 out.parent.mkdir(parents=True,exist_ok=True);out.write_bytes(b'synthetic workbook');out=out.parent
else:out.mkdir(parents=True,exist_ok=True)
report=out/'validation'/'analysis_validation.json' if stage=='analysis' else out/('figure_validation.json' if stage=='figures' else 'workbook_validation.json')
report.parent.mkdir(parents=True,exist_ok=True)
report.write_text(json.dumps({'status':STATUS_LITERAL}),encoding='utf-8')
(out/('executed_'+stage+'.json')).write_text(json.dumps({'stage':stage,'threads':os.environ.get('OPENBLAS_NUM_THREADS'),'temp':os.environ.get('TEMP')}),encoding='utf-8')
print('synthetic '+stage+' stage executed',flush=True)
""".replace('STAGE_LITERAL', repr(stage)).replace('STATUS_LITERAL', repr(self.stage_status.get(stage, 'PASS')))
            if stage == 'analysis' and self.mutate_package:
                code += "(Path(__file__).resolve().parents[1]/'payload.txt').write_text('changed input',encoding='utf-8')\n"
            (self.package / 'scripts' / filename).write_text(code, encoding='utf-8')
        (self.package / 'payload.txt').write_text('original input', encoding='utf-8')
        files = []
        for file in sorted(self.package.rglob('*')):
            if file.is_file() and file.name != 'package_manifest.json':
                files.append({'relative_path': file.relative_to(self.package).as_posix(),
                              'sha256': hashlib.sha256(file.read_bytes()).hexdigest(), 'bytes': file.stat().st_size})
        (self.package / 'package_manifest.json').write_text(json.dumps({'files': files}), encoding='utf-8')

    def execute(self):
        environment = os.environ.copy()
        environment['STRESS_PIPELINE_LAUNCHER'] = 'powershell'
        environment['PYTHONDONTWRITEBYTECODE'] = '1'
        return subprocess.run([sys.executable, str(self.package / 'scripts' / 'run_pipeline.py'),
                               '--output', str(self.output)], capture_output=True, text=True,
                              env=environment, timeout=90)

    def test_full_synthetic_pipeline_preserves_provenance(self):
        self.build()
        process = self.execute()
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
        report = json.loads((self.output / 'validation.json').read_text())
        self.assertEqual(report['status'], 'PASS')
        self.assertEqual(len(report['stages']), 3)
        self.assertTrue(report['package_unchanged'])
        self.assertEqual(report['package_before']['manifest_sha256'], report['package_after']['manifest_sha256'])
        self.assertEqual(report['native_windows_launcher_executed'], os.name == 'nt')
        for package in pipeline.DEPENDENCIES:
            self.assertIn(package, report['dependency_versions'])
        marker = json.loads((self.output / 'analysis' / 'executed_analysis.json').read_text())
        self.assertEqual(marker['threads'], '1')
        self.assertEqual(Path(marker['temp']), self.output / 'temp')

    def test_missing_manifest_fails_before_any_stage(self):
        self.build()
        (self.package / 'package_manifest.json').unlink()
        process = self.execute()
        self.assertEqual(process.returncode, 1)
        report = json.loads((self.output / 'validation.json').read_text())
        self.assertEqual(report['status'], 'FAIL')
        self.assertEqual(report['stages'], [])
        self.assertTrue((self.output / 'logs' / 'error.log').is_file())

    def test_manifest_paths_and_input_tampering_are_rejected(self):
        self.build()
        manifest_path = self.package / 'package_manifest.json'
        original = json.loads(manifest_path.read_text())
        for path in ('../outside.txt', '/absolute', 'C:/absolute', 'folder\\file', './payload.txt'):
            altered = json.loads(json.dumps(original))
            altered['files'][0]['relative_path'] = path
            manifest_path.write_text(json.dumps(altered))
            with self.subTest(path=path), self.assertRaises(ValueError):
                pipeline.verify_package(self.package)
        manifest_path.write_text(json.dumps(original))
        (self.package / 'payload.txt').write_text('replaced input')
        with self.assertRaisesRegex(ValueError, 'mismatch'):
            pipeline.verify_package(self.package)

    def test_nonempty_output_is_preserved(self):
        self.build()
        self.output.mkdir()
        sentinel = self.output / 'existing_result.txt'
        sentinel.write_text('keep this result')
        process = self.execute()
        self.assertEqual(process.returncode, 1)
        self.assertEqual(list(self.output.iterdir()), [sentinel])
        self.assertEqual(sentinel.read_text(), 'keep this result')
        rejected = list(self.base.glob('output_rejected_*'))
        self.assertEqual(len(rejected), 1)
        self.assertEqual(json.loads((rejected[0] / 'validation.json').read_text())['status'], 'FAIL')

    def test_nonpass_stage_report_stops_later_stages(self):
        self.stage_status['analysis'] = 'FAIL'
        self.build()
        process = self.execute()
        self.assertEqual(process.returncode, 1)
        report = json.loads((self.output / 'validation.json').read_text())
        self.assertEqual([stage['name'] for stage in report['stages']], ['analysis'])
        self.assertTrue(report['package_unchanged'])
        self.assertFalse((self.output / 'figures' / 'executed_figures.json').exists())

    def test_postrun_package_change_prevents_pass(self):
        self.mutate_package = True
        self.build()
        process = self.execute()
        self.assertEqual(process.returncode, 1)
        report = json.loads((self.output / 'validation.json').read_text())
        self.assertEqual(len(report['stages']), 3)
        self.assertEqual(report['status'], 'FAIL')
        self.assertFalse(report['package_unchanged'])
        self.assertEqual(report['package_after']['status'], 'FAIL')


if __name__ == '__main__':
    unittest.main()
