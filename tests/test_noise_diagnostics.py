"""Meaningful boundary and fail-closed tests for diagnostics-only correction."""
import csv
import importlib.util
import math
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('noise_diagnostics', ROOT / 'src/stress_prediction/noise_diagnostics.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class NoiseDiagnosticTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.folder = Path(self.temporary.name)
        mean = np.array([0., 1e-16, 1., 2., 3.])
        variance = np.array([0., 1e-32, 1., 1., 1.])
        self.sources = {(g, c): {'mean': mean.copy(), 'variance': variance.copy(),
                        'features': 5, 'maximum_identity_error': 1e-16,
                        'arithmetic_roundoff_bound': 1e-14, 'source_minimum': 0.,
                        'source_negative_fraction': 0., 'source_mean_noise_trace': float(variance.mean()/3)}
                        for g in MODULE.GENOTYPES for c in MODULE.CONDITIONS[:4]}
        latent, probability = [], []
        for g in MODULE.GENOTYPES:
            for null in ('single_mean_span', 'joint_rank3'):
                for c in MODULE.CONDITIONS:
                    identity = null == 'single_mean_span' and c in MODULE.CONDITIONS[:4]
                    boundary = identity and c != 'CT'
                    raw = mean.copy()
                    if boundary:
                        raw[0] = -1e-16
                    latent.append({'genotype': g, 'latent_null': null, 'condition': c,
                        'mean_displacement_fit_mse': '1e-34' if boundary else ('0.0' if identity else '0.1'),
                        'sample_mean_noise_trace_per_gene': repr(float(variance.mean()/3)),
                        'latent_log_mean_min': '-1e-16' if boundary else ('0.0' if identity else '-0.3'),
                        'latent_log_mean_negative_fraction': '0.2' if boundary or not identity else '0.0',
                        'latent_single_rank': '3', 'unchanged_note': 'Keep, unchanged'})
                    for scale in (.5, 1., 2.):
                        probability.append({'genotype': g, 'latent_null': null, 'condition': c,
                            'noise_sd_scale': repr(scale),
                            'expected_negative_condition_mean_log_fraction': repr(MODULE._negative_probability(raw, scale*np.sqrt(variance/3))) if identity else '0.27',
                            'method': 'Unchanged method'})
        self.write(MODULE.FILES[0], latent)
        self.write(MODULE.FILES[1], probability)
        (self.folder/'simulation_draws.csv.gz').write_bytes(b'unchanged simulation sentinel')

    def tearDown(self):
        self.temporary.cleanup()

    def write(self, name, rows):
        with (self.folder/name).open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)

    def run_correction(self):
        with patch.object(MODULE, '_source_identity', return_value=(self.sources, [])):
            return MODULE.canonicalize_noise_diagnostics(ROOT, self.folder)

    def test_exact_identity_correction_preserves_other_rows_and_simulation(self):
        before = {name: MODULE._read_csv(self.folder/name)[1] for name in MODULE.FILES}
        report = self.run_correction()
        self.assertEqual(report['status'], 'PASS')
        self.assertFalse(report['simulation_draws_or_predictions_modified'])
        for name in MODULE.FILES:
            after = MODULE._read_csv(self.folder/name)[1]
            self.assertEqual(len(after), len(before[name]))
            for old, new in zip(before[name], after):
                allowed = old['latent_null'] == 'single_mean_span' and old['condition'] in MODULE.CONDITIONS[:4]
                if not allowed:
                    self.assertEqual(old, new)
                if allowed and name == MODULE.FILES[0]:
                    self.assertEqual(float(new['latent_log_mean_negative_fraction']), 0.)
                    self.assertEqual(float(new['latent_log_mean_min']), 0.)
                self.assertEqual(old.keys(), new.keys())
            self.assertEqual(MODULE._read_csv(self.folder/'raw_diagnostics'/name)[1], before[name])
        self.assertEqual((self.folder/'simulation_draws.csv.gz').read_bytes(), b'unchanged simulation sentinel')

    def test_tiny_positive_sd_probability_is_not_forced_to_zero(self):
        mean = np.array([0., 1e-16])
        sd = np.array([0., 1e-16])
        correct = .25 * math.erfc(1 / math.sqrt(2))
        self.assertAlmostEqual(MODULE._negative_probability(mean, sd), correct)
        self.assertGreater(correct, 0)
        unstable = MODULE._negative_probability(np.array([-1e-16, -1e-16]), sd)
        self.assertGreater(unstable, .5)

    def test_substantive_projection_error_fails_without_changing_files(self):
        fields, rows = MODULE._read_csv(self.folder/MODULE.FILES[0])
        rows[1]['mean_displacement_fit_mse'] = '0.001'
        self.write(MODULE.FILES[0], rows)
        hashes = [MODULE._sha(self.folder/n) for n in MODULE.FILES]
        with self.assertRaisesRegex(ValueError, 'MSE'):
            self.run_correction()
        self.assertEqual(hashes, [MODULE._sha(self.folder/n) for n in MODULE.FILES])
        self.assertFalse((self.folder/'raw_diagnostics').exists())

    def test_excessive_negative_count_is_rejected(self):
        _, rows = MODULE._read_csv(self.folder/MODULE.FILES[0])
        rows[1]['latent_log_mean_negative_fraction'] = '0.8'
        self.write(MODULE.FILES[0], rows)
        with self.assertRaisesRegex(ValueError, 'negative count'):
            self.run_correction()

    def test_unjustified_gaussian_probability_is_rejected(self):
        _, rows = MODULE._read_csv(self.folder/MODULE.FILES[1])
        rows[3]['expected_negative_condition_mean_log_fraction'] = '0.99'
        self.write(MODULE.FILES[1], rows)
        with self.assertRaisesRegex(ValueError, 'negative probability'):
            self.run_correction()
        self.assertFalse((self.folder/'raw_diagnostics').exists())

    def test_idempotence_and_postcorrection_tampering(self):
        first = self.run_correction()
        self.assertEqual(first, self.run_correction())
        with (self.folder/MODULE.FILES[0]).open('a') as f:
            f.write('tampered\n')
        with self.assertRaisesRegex(ValueError, 'changed after correction'):
            self.run_correction()

    def test_immutable_provenance_directory_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'immutable provenance'):
            MODULE.canonicalize_noise_diagnostics(ROOT, ROOT/'provenance/reference_tables/noise')

    def test_summary_hashes_are_updated_and_original_is_preserved(self):
        summary = {'code_sha256': 'unchanged historical code hash',
                   'outputs': {name: MODULE._sha(self.folder/name) for name in MODULE.FILES}}
        path = self.folder/'RUN_SUMMARY.json'
        path.write_text(json.dumps(summary))
        original = path.read_bytes()
        report = self.run_correction()
        current = json.loads(path.read_text())
        self.assertEqual(current['code_sha256'], summary['code_sha256'])
        self.assertEqual((self.folder/'raw_diagnostics/RUN_SUMMARY.json').read_bytes(), original)
        for name in MODULE.FILES:
            self.assertEqual(current['outputs'][name], MODULE._sha(self.folder/name))
        self.assertEqual(current['diagnostic_postprocessing']['amendment_report'], 'PORTABILITY_DIAGNOSTIC_AMENDMENT.json')
        self.assertEqual(report['run_summary']['corrected_sha256'], MODULE._sha(path))
        self.assertEqual(report, self.run_correction())

    def test_mismatched_original_summary_is_rejected_before_correction(self):
        path = self.folder/'RUN_SUMMARY.json'
        path.write_text(json.dumps({'outputs': {name: 'wrong' for name in MODULE.FILES}}))
        with self.assertRaisesRegex(ValueError, 'run-summary diagnostic hash mismatch'):
            self.run_correction()
        self.assertFalse((self.folder/'raw_diagnostics').exists())


if __name__ == '__main__':
    unittest.main()
