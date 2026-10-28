#!/usr/bin/env python3
"""Independent primary prediction algebra from deposited processed matrices.

Uses small Gram-matrix projection rather than the production SVD basis.
No production analysis modules, fitted arrays or old outcome tables are used
to calculate predictions. This is an outcome-known audit, not a new test.
"""
from __future__ import annotations
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from audit_release_tables import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args(); repo = args.repository.resolve(); out = args.output.resolve()
    historical = repo/'provenance/historical'
    sources = json.loads((repo/'config/data_sources.json').read_text())['sources']
    source_checks = []
    for source in sources:
        if source['accession'] == 'GSE147962':
            continue
        path = repo/source['relative_path']
        assert digest(path) == source['sha256']
        source_checks.append({'file': source['relative_path'], 'sha256': source['sha256']})
    designs = [
        ('external_prediction', historical/'stress_combinations_code6/fallback_metadata/sample_inventory.csv',
         'genotype_metadata', 'library_name', 'gene', ['S', 'EL', 'HS']),
        ('citrus_prediction', historical/'stress_combinations_code8/input/GSE203331_design_verified.csv',
         'genotype_code', 'expression_column', 'Feature_GID', ['WS', 'HL', 'HS'])]
    checks = []; design_checks = []
    for stage, design_path, genotype_col, sample_col, id_col, factors in designs:
        design = pd.read_csv(design_path)
        design_checks.append({'file': str(design_path.relative_to(repo)), 'sha256': digest(design_path)})
        summary = pd.read_csv(out/'tables'/stage/'model_summary.csv')
        for genotype, sample in design.groupby(genotype_col, sort=False):
            if stage == 'external_prediction':
                suffix = 'Col0' if genotype == 'WT' else 'RAP2.3'
                path = historical/f'stress_combinations_code6/input/GSE330879_{suffix}_fpkm.txt.gz'
            else:
                path = repo/next(s['relative_path'] for s in sources if s['accession'] == 'GSE203331')
            columns = sample[sample_col].tolist()
            raw = pd.read_csv(path, sep='\t', usecols=[id_col, *columns], dtype={id_col: str})
            raw = raw.loc[~raw[id_col].duplicated(keep=False)].reset_index(drop=True)
            assert raw[id_col].notna().all() and raw[id_col].is_unique
            numerical = raw[columns].to_numpy(float).T
            assert np.isfinite(numerical).all() and np.all(numerical >= 0)
            bits = sample[factors].to_numpy(int); orders = bits.sum(axis=1)
            assert sorted(orders.tolist()) == [0]*3+[1]*9+[2]*9+[3]*3
            target_ids = np.flatnonzero(orders == 3)
            filtering_ids = np.flatnonzero(orders <= 1)
            assert not set(target_ids) & set(filtering_ids)
            for threshold in (.5, 1., 2.):
                keep = numerical[filtering_ids].mean(axis=0) >= threshold
                values = np.log2(1. + numerical[:, keep])
                c = values[orders == 0].mean(axis=0)
                singles = np.array([values[(orders == 1) & (bits[:, j] == 1)].mean(axis=0)
                                    for j in range(3)])
                pair_bits = np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1]])
                pairs = np.array([values[np.all(bits == b, axis=1)].mean(axis=0) for b in pair_bits])
                e = singles - c
                centroid = pairs.mean(axis=0)
                # The small 3x3 Gram inverse yields the same fitted subspace
                # while avoiding import/reuse of the production SVD routine.
                gram = e @ e.T
                assert np.linalg.matrix_rank(gram) == 3
                within = ((centroid-c) @ e.T) @ np.linalg.pinv(gram) @ e
                predictions = {'control': c, 'full_pair_centroid': centroid,
                               'single_additive': singles.sum(axis=0)-2*c,
                               'second_order_factorial': pairs.sum(axis=0)-singles.sum(axis=0)+c,
                               ('within_only' if stage == 'external_prediction' else 'within_V3'): c+within}
                if stage == 'citrus_prediction':
                    x = np.concatenate([e.T * b for b in pair_bits], axis=0)
                    response = (pairs-c).reshape(-1)
                    beta = np.linalg.lstsq(x, response, rcond=None)[0]
                    additive = x.sum(axis=1)
                    alpha = (additive @ response)/(additive @ additive)
                    predictions['single_mean'] = singles.mean(axis=0)
                    predictions['pair_trained_scalar'] = c + alpha*e.sum(axis=0)
                    predictions['pair_trained_factor_OLS'] = c + beta@e
                target = values[target_ids]
                for model, prediction in predictions.items():
                    observed = summary[(summary.genotype == genotype) &
                                       (summary.threshold_fpkm == threshold) &
                                       (summary.model == model)]
                    assert len(observed) == 1
                    observed = observed.iloc[0]
                    assert int(observed.genes) == int(keep.sum())
                    expected = float(observed['mse' if stage == 'external_prediction' else 'mean_library_mse'])
                    calculated = float(np.mean((target-prediction)**2))
                    assert np.isclose(calculated, expected, atol=1e-10, rtol=1e-9), (
                        genotype, threshold, model, calculated, expected)
                    checks.append({'stage': stage, 'genotype': genotype, 'threshold_fpkm': threshold,
                                   'model': model, 'features': int(keep.sum()), 'target_pools': len(target),
                                   'independently_calculated_mse': calculated, 'released_mse': expected,
                                   'absolute_difference': abs(calculated-expected)})
    assert len(checks) == 78
    report = {'status': 'PASS_INDEPENDENT_PRIMARY_PREDICTION_ALGEBRA',
              'completed_utc': datetime.now(timezone.utc).isoformat(),
              'production_analysis_modules_or_fitted_arrays_used': False,
              'scope': 'Outcome-known independent recomputation from three deposited FPKM files, four genotype applications and three fixed thresholds. Independent primary control, centroid, additive, second-order and mean-single-span predictions, plus citrus single-mean and pair-trained scalar/factor fits. Does not cover raw-read preprocessing or establish new blind biological validation.',
              'mse_comparisons': len(checks),
              'maximum_absolute_difference': max(v['absolute_difference'] for v in checks),
              'filter_uses_control_and_singles_only': True,
              'target_numerics_used_for_prediction_or_filtering': False,
              'absolute_tolerance': 1e-10, 'relative_tolerance': 1e-9,
              'audit_code_sha256': digest(__file__), 'sources': source_checks,
              'designs': design_checks, 'comparisons': checks}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('sources', 'designs', 'comparisons')}, indent=2))


if __name__ == '__main__':
    main()
