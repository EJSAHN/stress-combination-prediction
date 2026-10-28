#!/usr/bin/env python3
"""Independently recompute every numerical figure trace from released tables.

This release audit imports neither the analysis engine nor the plotting code.
Optional original records provide an additional untouched-source comparison.
It verifies numerical provenance and exports, not new biological replication.
"""
from __future__ import annotations
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from audit_release_tables import compare_csv, digest


ATOL = 1e-10
RTOL = 1e-9


def close(actual, expected, label):
    assert np.isclose(float(actual), float(expected), atol=ATOL, rtol=RTOL), (
        label, actual, expected)


def original_for(rel, mappings, archive_root, repo):
    key = 'tables/' + rel
    if key in mappings:
        item = mappings[key]
        expected = (archive_root / item['original_relative_path'] if archive_root
                    else repo / item['reference'])
        assert digest(expected) == item['sha256'], ('reference hash', rel)
        return expected
    if rel == 'noise/simulation_draws.csv.gz' and archive_root:
        return archive_root / 'stress_combinations_code7/results/simulation/simulation_draws.csv.gz'
    assert rel == 'noise/simulation_draws.csv.gz', ('unmapped figure source', rel)
    return None


def audit_noise(draws, summary):
    checks = []
    for _, row in summary.iterrows():
        selected = draws[(draws.genotype == row.genotype) &
                         (draws.latent_null == row.latent_null) &
                         (draws['lambda'] == row['lambda']) &
                         (draws.noise_sd_scale == row.noise_sd_scale)]
        values = selected[row.score].to_numpy(float)
        assert len(values) == int(row.draws) == 1000
        assert selected['draw'].nunique() == 1000
        close(row['mean'], values.mean(), 'noise mean')
        close(row['sd'], values.std(ddof=1), 'noise SD')
        for field, q in [('q025', .025), ('q50', .5), ('q95', .95),
                         ('q975', .975), ('q99', .99)]:
            close(row[field], np.quantile(values, q), 'noise ' + field)
        if pd.notna(row.n_at_least_observed):
            count = int(np.count_nonzero(values >= row.observed_statistic))
            close(row.n_at_least_observed, count, 'conditional exceedance count')
            close(row.conditional_exceedance_fraction, count / len(values),
                  'conditional exceedance fraction')
        checks.append({'genotype': row.genotype, 'latent_null': row.latent_null,
                       'diagonal_fraction': float(row['lambda']),
                       'noise_sd_scale': float(row.noise_sd_scale), 'score': row.score})
    assert len(checks) == 180
    return checks


def trace_value(record, table):
    ids = [int(x) for x in str(record.source_rows).split(';') if x != '']
    fields = record.source_fields.split(';')
    operation = record.operation
    if operation.startswith('histogram count;'):
        return float(len(ids)), ids
    assert ids and all(0 <= i < len(table) for i in ids)
    selected = table.iloc[ids]
    if operation == 'identity':
        field = fields[1] if fields[0] == 'requested_fraction' else fields[0]
        value = float(selected.iloc[0][field])
        if fields[0] == 'requested_fraction':
            close(record.x, selected.iloc[0][fields[0]], 'requested fraction x')
        if 'q025' in fields and 'q975' in fields:
            close(record.lower, selected.iloc[0].q025, 'lower partition/simulation range')
            close(record.upper, selected.iloc[0].q975, 'upper partition/simulation range')
        if 'n_features' in fields:
            close(record.n, selected.iloc[0].n_features, 'group feature count')
    elif operation in ('numerator_row / control_row',
                        'library_model_MSE / same_library_control_MSE',
                        'log10(numerator_row / control_row)'):
        assert len(ids) == 2 and len(fields) == 1
        value = float(selected.iloc[0][fields[0]]) / float(selected.iloc[1][fields[0]])
        for key in ('scenario', 'genotype', 'sample', 'library'):
            if key in selected:
                assert selected.iloc[0][key] == selected.iloc[1][key], ('ratio pairing', key)
        if operation.startswith('log10'):
            value = float(np.log10(value))
        if 'denominator' in record and record.denominator != '':
            close(record.denominator, selected.iloc[1][fields[0]], 'denominator')
    elif operation == 'within_span_MSE - centroid_MSE':
        assert len(ids) == 2 and len(fields) == 1
        value = float(selected.iloc[0][fields[0]]) - float(selected.iloc[1][fields[0]])
    elif operation == 'model_MSE / original_centroid_MSE':
        assert len(ids) == 1 and len(fields) == 2
        value = float(selected.iloc[0][fields[0]]) / float(selected.iloc[0][fields[1]])
    elif operation == 'net_mass / itself':
        net = float(selected.iloc[0].net_mass)
        assert net != 0
        value = net / net
    elif operation.replace(' ', '') == '(V11_MSE-centroid_MSE)/(V3_MSE-centroid_MSE)':
        assert len(ids) == 3 and len(fields) == 1
        a, b, c = (float(x) for x in selected[fields[0]])
        value = (a-b)/(c-b)
    else:
        raise AssertionError(('unsupported numerical operation', operation))
    return value, ids


def audit_histograms(traces, tables):
    hist = traces[traces.operation.str.startswith('histogram count;')]
    checks = []
    for (figure, panel, series), bins in hist.groupby(['figure', 'panel', 'series']):
        table = tables[bins.iloc[0].source_file]
        field = bins.iloc[0].source_fields
        bins = bins.sort_values('x', key=lambda x: pd.to_numeric(x))
        all_ids = [int(i) for v in bins.source_rows for i in str(v).split(';') if i != '']
        assert len(all_ids) == len(set(all_ids)) == 1000
        data = table.iloc[all_ids]
        for key in ['genotype', 'latent_null', 'lambda', 'noise_sd_scale']:
            assert data[key].nunique() == 1, ('mixed histogram scenario', key)
        assert data.iloc[0].genotype == series
        assert data.iloc[0].latent_null == 'joint_rank3'
        assert float(data.iloc[0]['lambda']) == 0.
        assert float(data.iloc[0].noise_sd_scale) == 1.
        edges = np.array([float(v) for v in bins.lower] + [float(bins.iloc[-1].upper)])
        assert len(edges) == 31 and np.all(np.diff(edges) > 0)
        for index, (_, b) in enumerate(bins.iterrows()):
            close(b.upper, edges[index+1], 'contiguous histogram upper edge')
            close(b.x, (edges[index]+edges[index+1])/2, 'histogram bin centre')
            got = set(int(i) for i in str(b.source_rows).split(';') if i != '')
            # Text round trips of the extrema can move them by < 1e-15.
            lo = edges[index] - (1e-12 if index == 0 else 0.)
            hi = edges[index+1] + (1e-12 if index == len(bins)-1 else 0.)
            condition = (data[field] >= lo) & (
                data[field] <= hi if index == len(bins)-1 else data[field] < hi)
            expected = set(data.index[condition])
            assert got == expected, ('histogram membership', figure, panel, index)
            close(b.value, len(expected), 'histogram count')
        checks.append({'figure': figure, 'panel': panel, 'genotype': series,
                       'bins': len(bins), 'unique_draws': len(all_ids)})
    assert len(checks) == 2
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path, required=True, help='Analysis output containing tables/')
    parser.add_argument('--figures', type=Path, required=True, help='Final figure output directory')
    parser.add_argument('--archive-root', type=Path)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    repo = args.repository.resolve(); output = args.output.resolve(); figures = args.figures.resolve()
    mapping_path = repo/'config/table_reference_map.json'
    mappings = {r['output']: r for r in json.loads(mapping_path.read_text())['tables']}
    manifest_path = figures/'figure_validation.json'
    manifest = json.loads(manifest_path.read_text())
    assert manifest['figure_count'] == 8 and manifest['main_figures'] == 4
    assert manifest['plot_code_sha256'] == digest(repo/'scripts/make_figures.py'), 'plot code changed since rendering'
    tables = {}; input_checks = []; figure_checks = []; all_traces = []
    for rel, info in manifest['input_files'].items():
        path = output/'tables'/rel
        assert digest(path) == info['sha256'], ('figure input hash', rel)
        table = pd.read_csv(path)
        assert len(table) == info['rows']
        tables[rel] = table
        ref = original_for(rel, mappings, args.archive_root, repo)
        check = {'source_file': rel, 'output_sha256': digest(path), 'rows': len(table)}
        if ref is not None:
            check.update(compare_csv(path, ref))
            check['original_sha256'] = digest(ref)
        input_checks.append(check)
    for item in manifest['figures']:
        source_path = figures/item['source_csv']
        assert digest(source_path) == item['source_sha256']
        traces = pd.read_csv(source_path, dtype=str, keep_default_na=False)
        assert len(traces) == item['source_rows'] and len(traces) > 0
        assert traces.figure.eq(item['figure']).all()
        for index, record in traces.iterrows():
            assert record.source_file in tables
            expected, ids = trace_value(record, tables[record.source_file])
            close(record.value, expected, (item['figure'], index, 'declared value'))
            assert record.value_axis in ('x', 'y', 'color', 'annotation')
            if record.value_axis in ('x', 'y'):
                close(record[record.value_axis], expected,
                      (item['figure'], index, 'plotted coordinate'))
        all_traces.append(traces)
        for kind in ('pdf', 'png'):
            path = figures/item[kind]
            assert digest(path) == item[kind+'_sha256']
        png = Image.open(figures/item['png'])
        assert list(png.size) == item['png_pixels']
        dpi = png.info.get('dpi', (0, 0))
        assert min(dpi) >= 299.9
        pdf = (figures/item['pdf']).read_bytes()
        assert pdf.startswith(b'%PDF-') and b'%%EOF' in pdf[-1024:]
        assert b'/Subtype /Image' not in pdf
        figure_checks.append({'figure': item['figure'], 'trace_rows': len(traces),
                              'source_sha256': digest(source_path),
                              'pdf_sha256': digest(figures/item['pdf']),
                              'png_sha256': digest(figures/item['png']),
                              'png_dpi': list(dpi), 'raster_image_objects_in_pdf': False})
    combined = pd.concat(all_traces, ignore_index=True).fillna('')
    noise_checks = audit_noise(tables['noise/simulation_draws.csv.gz'], tables['noise/simulation_summary.csv'])
    histogram_checks = audit_histograms(combined, tables)
    report = {
        'status': 'PASS_INDEPENDENT_FIGURE_SOURCE_AUDIT',
        'completed_utc': datetime.now(timezone.utc).isoformat(),
        'analysis_or_plot_code_imported': False,
        'untouched_original_archive_checked': args.archive_root is not None,
        'scope': 'Every declared numerical figure trace recomputed from its source cells; source tables compared to pinned original records; histogram membership and all conditional-noise summary moments/quantiles independently recomputed. This is a numerical provenance/export audit, not new biological validation or native Windows execution.',
        'absolute_tolerance': ATOL, 'relative_tolerance': RTOL,
        'figure_count': len(figure_checks), 'trace_rows': len(combined),
        'input_table_count': len(input_checks), 'conditional_noise_summary_rows': len(noise_checks),
        'conditional_noise_draw_rows': len(tables['noise/simulation_draws.csv.gz']),
        'figure_manifest_sha256': digest(manifest_path),
        'figure_code_sha256': digest(repo/'scripts/make_figures.py'),
        'audit_code_sha256': digest(__file__),
        'table_audit_helper_sha256': digest(repo/'scripts/audit_release_tables.py'),
        'mapping_sha256': digest(mapping_path),
        'inputs': input_checks, 'figures': figure_checks, 'histograms': histogram_checks}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('inputs', 'figures', 'histograms')}, indent=2))


if __name__ == '__main__':
    main()
