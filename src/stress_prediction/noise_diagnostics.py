"""Correct only the exact-identity rows of two working-scale diagnostics.

The audited simulation engine and its predictions/draws are unchanged. The
three defining single-response rows project to themselves mathematically; their
absolute latent means therefore equal the observed nonnegative means. Passing
those means through SVD again makes strict sign/CDF diagnostics nonportable at
zero. This adapter retains raw CSVs and uses the exact identity for diagnostics.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

GENOTYPES = ('WT', 'RAP2.3')
CONDITIONS = ('CT', 'S', 'EL', 'HS', 'S+EL', 'S+HS', 'EL+HS', 'S+EL+HS')
FILES = ('latent_mean_diagnostics.csv', 'negative_working_scale_probabilities.csv')
IDENTITY_ERROR_LIMIT = 5e-10
VERSION = 'exact-single-span-diagnostics-v1'


def _sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def _check(condition, message):
    if not condition:
        raise ValueError(message)


def _negative_probability(mean, sd):
    """Marginal Gaussian P(X<0), including a deterministic SD=0 atom."""
    mean, sd = np.asarray(mean, float), np.asarray(sd, float)
    _check(mean.shape == sd.shape and np.isfinite(mean).all()
           and np.isfinite(sd).all() and (sd >= 0).all(), 'Invalid Gaussian inputs')
    out = np.zeros_like(mean)
    positive = sd > 0
    ratio = mean[positive] / (np.sqrt(2) * sd[positive])
    out[positive] = np.fromiter((0.5 * math.erfc(float(z)) for z in ratio),
                               dtype=float, count=int(positive.sum()))
    out[~positive] = (mean[~positive] < 0).astype(float)
    return float(out.mean())


def _source_identity(repo):
    base = repo / 'provenance/historical/stress_combinations_code6'
    pins = {x['relative_path']: x for x in json.loads(
        (repo / 'config/analysis_assets.json').read_text(encoding='utf-8'))['files']}

    def pinned(path):
        relative = path.relative_to(repo).as_posix()
        _check(relative in pins, 'Missing source pin: ' + relative)
        item = pins[relative]
        _check(path.is_file() and path.stat().st_size == item['bytes']
               and _sha(path) == item['sha256'], 'Source pin failed: ' + relative)
        return item['sha256']

    manifest_path = base / 'frozen_verified/FREEZE_MANIFEST.json'
    pinned(manifest_path)
    fm = json.loads(manifest_path.read_text(encoding='utf-8'))
    design_path = base / 'fallback_metadata/sample_inventory.csv'
    _check(pinned(design_path) == fm['design_sha256'], 'Design hash mismatch')
    design = pd.read_csv(design_path)
    records, provenance = {}, []
    for genotype in GENOTYPES:
        scenario = next(x for x in fm['scenarios']
                        if x['genotype'] == genotype and x['threshold_fpkm'] == 1)
        asset = base / 'frozen_verified' / scenario['assets']
        _check(pinned(asset) == scenario['assets_sha256'], 'Fixed panel hash mismatch')
        with np.load(asset, allow_pickle=False) as archive:
            genes = archive['genes'].astype(str)
        source = base / 'input' / Path(fm['source_hashes'][genotype]['path']).name
        _check(pinned(source) == fm['source_hashes'][genotype]['sha256'], 'Expression hash mismatch')
        dd = design[design.genotype_metadata == genotype]
        names = dd.library_name.tolist()
        _check(len(names) == 24 and len(set(names)) == 24, 'Expected 24 distinct RNA pools')
        frame = pd.read_csv(source, sep='\t', usecols=['gene', *names],
                            dtype={'gene': str}, keep_default_na=False)
        duplicated = frame.gene.duplicated(keep=False)
        _check(frame.loc[duplicated, 'gene'].value_counts().to_dict()
               == {'AT1G79790': 2, 'AT2G04725': 2}, 'Unexpected duplicate gene IDs')
        frame = frame.loc[~duplicated].set_index('gene')
        _check(frame.index.is_unique and len(set(genes)) == len(genes)
               and set(genes).issubset(frame.index), 'Invalid fixed feature panel')
        frame = frame.loc[genes]
        arrays = []
        for condition in CONDITIONS:
            columns = dd.loc[dd.condition == condition, 'library_name'].tolist()
            _check(len(columns) == 3, 'Exactly three pools per condition are required')
            raw = frame[columns].to_numpy(dtype=float).T
            _check(np.isfinite(raw).all() and (raw >= 0).all(), 'Invalid source FPKM')
            arrays.append(np.log2(1 + raw))
        x = np.stack(arrays)
        means = x.mean(axis=1)
        residuals = x - means[:, None, :]
        variances = residuals.var(axis=1, ddof=1)
        effects = means - means[0]
        _, singular, q = np.linalg.svd(effects[1:4], full_matrices=False)
        rank_limit = max(effects[1:4].shape) * np.finfo(float).eps * singular[0]
        _check(int(np.sum(singular > rank_limit)) == 3, 'Single-response identity needs rank three')
        _check(np.allclose(q @ q.T, np.eye(3), atol=IDENTITY_ERROR_LIMIT, rtol=0),
               'Single-span basis is not orthonormal')
        reconstructed = means[0] + (effects @ q.T) @ q
        for index, condition in enumerate(CONDITIONS[:4]):
            mean = means[index]
            error = float(np.max(np.abs(reconstructed[index] - mean)))
            _check(error <= IDENTITY_ERROR_LIMIT, 'Source projection identity failed')
            _check((mean >= 0).all(), 'Observed means must be nonnegative')
            roundoff = float(64 * np.finfo(float).eps *
                             (1 + np.max(np.abs(means[0])) + np.max(np.abs(mean))
                              + np.max(np.abs(effects[index]))))
            records[(genotype, condition)] = {
                'mean': mean, 'variance': variances[index], 'features': len(genes),
                'maximum_identity_error': error, 'arithmetic_roundoff_bound': roundoff,
                'source_minimum': float(mean.min()),
                'source_negative_fraction': float(np.mean(mean < 0)),
                'source_mean_noise_trace': float(variances[index].mean() / 3),
            }
        provenance.append({'genotype': genotype, 'source': source.relative_to(repo).as_posix(),
                           'source_sha256': _sha(source), 'fixed_panel_sha256': _sha(asset),
                           'gene_count': len(genes), 'single_span_singular_values': singular.tolist(),
                           'maximum_identity_error': max(records[(genotype, c)]['maximum_identity_error']
                                                         for c in CONDITIONS[:4])})
    return records, provenance


def _legacy_identity_bound(row, source):
    """Bound each old reconstruction error by sqrt(G*MSE), plus roundoff."""
    mse = float(row['mean_displacement_fit_mse'])
    _check(math.isfinite(mse) and 0 <= mse <= IDENTITY_ERROR_LIMIT ** 2,
           'Legacy identity-row MSE is not a roundoff-scale quantity')
    bound = math.sqrt(source['features'] * mse) + source['arithmetic_roundoff_bound']
    _check(bound <= IDENTITY_ERROR_LIMIT, 'Legacy per-feature identity bound is too large')
    raw_min = float(row['latent_log_mean_min'])
    fraction = float(row['latent_log_mean_negative_fraction'])
    _check(math.isfinite(raw_min) and math.isfinite(fraction) and 0 <= fraction <= 1,
           'Invalid legacy sign diagnostics')
    _check(abs(raw_min - source['source_minimum']) <= bound,
           'Legacy minimum is inconsistent with the exact source identity')
    count = fraction * source['features']
    _check(abs(count - round(count)) < 1e-6, 'Legacy negative fraction is not a feature count')
    _check(round(count) <= int(np.count_nonzero(source['mean'] <= bound)),
           'Legacy negative count exceeds all features that could cross zero')
    noise = float(row['sample_mean_noise_trace_per_gene'])
    _check(np.isclose(noise, source['source_mean_noise_trace'], atol=5e-10, rtol=5e-10),
           'Legacy source noise trace differs')
    _check(int(row['latent_single_rank']) == 3, 'Legacy single span does not have rank three')
    return bound


def _read_csv(path):
    with path.open('r', encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream)
        return list(reader.fieldnames), list(reader)


def _write_csv(path, fields, rows):
    temporary = path.with_name(path.name + '.portability-writing')
    _check(not temporary.exists(), 'Refusing to overwrite an incomplete diagnostic write')
    with temporary.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def canonicalize_noise_diagnostics(repo, noise_dir):
    """Amend two diagnostic CSVs, keeping raw backups and a cell-level record.

    This never reads or writes simulation draws, prediction arrays or scores.
    All other diagnostic rows and field strings are retained exactly.
    """
    repo, noise_dir = Path(repo).resolve(), Path(noise_dir).resolve()
    _check(not noise_dir.is_relative_to(repo / 'provenance'),
           'Refusing to modify immutable provenance or reference directories')
    report_path = noise_dir / 'PORTABILITY_DIAGNOSTIC_AMENDMENT.json'
    backup = noise_dir / 'raw_diagnostics'
    if report_path.exists():
        old = json.loads(report_path.read_text(encoding='utf-8'))
        _check(old.get('definition') == VERSION and old.get('status') == 'PASS',
               'Unexpected existing diagnostic amendment')
        for item in old['files']:
            _check(_sha(noise_dir / item['filename']) == item['corrected_sha256']
                   and _sha(backup / item['filename']) == item['raw_sha256'],
                   'Amended or raw diagnostics changed after correction')
        if old.get('run_summary') is not None:
            item = old['run_summary']
            _check(_sha(noise_dir / 'RUN_SUMMARY.json') == item['corrected_sha256']
                   and _sha(backup / 'RUN_SUMMARY.json') == item['raw_sha256'],
                   'Amended or raw run summary changed after correction')
        return old
    _check(not backup.exists(), 'Existing raw backup without an amendment; preserve it and use a fresh run')
    sources, provenance = _source_identity(repo)
    headers, rows = {}, {}
    for name in FILES:
        headers[name], rows[name] = _read_csv(noise_dir / name)
    summary_path = noise_dir / 'RUN_SUMMARY.json'
    summary = None
    if summary_path.exists():
        summary = json.loads(summary_path.read_text(encoding='utf-8'))
        _check(isinstance(summary.get('outputs'), dict), 'Missing noise output hash map')
        for name in FILES:
            _check(summary['outputs'].get(name) == _sha(noise_dir / name),
                   'Original run-summary diagnostic hash mismatch: ' + name)
    latent, probabilities = rows[FILES[0]], rows[FILES[1]]
    expected_latent = {(g, n, c) for g in GENOTYPES for n in ('single_mean_span', 'joint_rank3') for c in CONDITIONS}
    actual_latent = [(r['genotype'], r['latent_null'], r['condition']) for r in latent]
    _check(len(actual_latent) == len(expected_latent) and set(actual_latent) == expected_latent,
           'Unexpected latent diagnostic rows')
    expected_probability = {(g, n, c, s) for g in GENOTYPES for n in ('single_mean_span', 'joint_rank3')
                            for c in CONDITIONS for s in (0.5, 1.0, 2.0)}
    actual_probability = [(r['genotype'], r['latent_null'], r['condition'], float(r['noise_sd_scale'])) for r in probabilities]
    _check(len(actual_probability) == len(expected_probability)
           and set(actual_probability) == expected_probability, 'Unexpected probability diagnostic rows')
    changes, checks, bounds = [], [], {}

    def replace_cell(name, index, row, field, value):
        before, after = row[field], repr(float(value))
        if before != after:
            changes.append({'file': name, 'csv_row_number': index + 2,
                            'genotype': row['genotype'], 'latent_null': row['latent_null'],
                            'condition': row['condition'], 'noise_sd_scale': row.get('noise_sd_scale'),
                            'field': field, 'raw_text': before, 'corrected_text': after})
            row[field] = after

    for index, row in enumerate(latent):
        key = (row['genotype'], row['condition'])
        if row['latent_null'] != 'single_mean_span' or key not in sources:
            continue
        source = sources[key]
        bound = _legacy_identity_bound(row, source)
        bounds[key] = bound
        checks.append({'genotype': key[0], 'condition': key[1], 'features': source['features'],
                       'maximum_source_projection_error': source['maximum_identity_error'],
                       'legacy_identity_error_bound': bound,
                       'canonical_minimum': source['source_minimum'],
                       'canonical_negative_fraction': source['source_negative_fraction']})
        replace_cell(FILES[0], index, row, 'latent_log_mean_min', source['source_minimum'])
        replace_cell(FILES[0], index, row, 'latent_log_mean_negative_fraction', source['source_negative_fraction'])
    for index, row in enumerate(probabilities):
        key = (row['genotype'], row['condition'])
        if row['latent_null'] != 'single_mean_span' or key not in sources:
            continue
        source = sources[key]
        scale = float(row['noise_sd_scale'])
        sd = scale * np.sqrt(source['variance'] / 3)
        canonical = _negative_probability(source['mean'], sd)
        low = _negative_probability(source['mean'] + bounds[key], sd)
        high = _negative_probability(source['mean'] - bounds[key], sd)
        raw = float(row['expected_negative_condition_mean_log_fraction'])
        _check(math.isfinite(raw) and low - 5e-10 <= raw <= high + 5e-10,
               'Legacy negative probability is outside its projection-error bound')
        replace_cell(FILES[1], index, row, 'expected_negative_condition_mean_log_fraction', canonical)
    # Validate everything before modifying either result table.
    backup.mkdir()
    file_records = []
    for name in FILES:
        path = noise_dir / name
        raw_hash = _sha(path)
        shutil.copy2(path, backup / name)
        _write_csv(path, headers[name], rows[name])
        file_records.append({'filename': name, 'raw_sha256': raw_hash,
                             'corrected_sha256': _sha(path)})
    summary_record = None
    if summary is not None:
        raw_summary_hash = _sha(summary_path)
        shutil.copy2(summary_path, backup / summary_path.name)
        for item in file_records:
            summary['outputs'][item['filename']] = item['corrected_sha256']
        summary['diagnostic_postprocessing'] = {
            'definition': VERSION, 'amendment_report': report_path.name,
            'raw_summary_sha256': raw_summary_hash, 'adapter_sha256': _sha(Path(__file__)),
            'changed_outputs': list(FILES), 'simulation_draws_or_predictions_modified': False,
            'historical_engine_code_sha256_preserved': summary.get('code_sha256'),
        }
        temporary_summary = summary_path.with_suffix('.json.portability-writing')
        _check(not temporary_summary.exists(), 'An incomplete summary write already exists')
        temporary_summary.write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')
        temporary_summary.replace(summary_path)
        summary_record = {'filename': summary_path.name, 'raw_sha256': raw_summary_hash,
                          'corrected_sha256': _sha(summary_path)}
    report = {'status': 'PASS', 'definition': VERSION,
              'scope': 'single_mean_span CT/S/EL/HS diagnostic rows only',
              'identity_error_limit': IDENTITY_ERROR_LIMIT,
              'explanation': 'A full-rank single span reproduces its defining rows exactly; use observed absolute means for strict sign and Gaussian CDF diagnostics.',
              'simulation_draws_or_predictions_modified': False,
              'historical_engine_or_reference_modified': False,
              'source_provenance': provenance, 'identity_checks': checks,
              'files': file_records, 'run_summary': summary_record, 'changed_cells': changes,
              'changed_cell_count': len(changes), 'adapter_sha256': _sha(Path(__file__))}
    temporary = report_path.with_suffix('.json.writing')
    temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(report_path)
    return report
