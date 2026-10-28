#!/usr/bin/env python3
"""Frozen pair-to-triple prediction in GSE330879 author-deposited FPKM.

Freeze parses 21 non-target columns per genotype; evaluate parses the triple
columns only after all six genotype/threshold scenarios have been frozen.
No downloads, installations, or inferred matching of replicate labels.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import tempfile
import numpy as np
import pandas as pd

GENOTYPES = ('WT', 'RAP2.3')
SINGLES = ('S', 'EL', 'HS')
PAIRS = ('EL+HS', 'S+EL', 'S+HS')
TARGET = 'S+EL+HS'
THRESHOLDS = (1.0, 0.5, 2.0)
MODELS = ('control', 'full_pair_centroid', 'within_only', 'outside_only',
          'single_additive', 'second_order_factorial')
WEIGHTS = np.array([[1/3, 1/3, 1/3, 0], [0, 0, 0, 1.]])

def now():
    return datetime.now(timezone.utc).isoformat()

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    return h.hexdigest()

def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')

def key_for(genotype, threshold):
    return genotype.replace('.', '_') + '_fpkm' + str(threshold).replace('.', 'p')

def checked_protocol(path):
    p = json.loads(path.read_text())
    assert p['source'] == 'GSE330879'
    assert p['primary_genotype'] == 'WT' and p['secondary_genotype'] == 'RAP2.3'
    assert set(p['single_conditions']) == set(SINGLES)
    assert set(p['donor_conditions']) == set(PAIRS)
    assert p['target_condition'] == TARGET
    assert set(p['thresholds']) == set(THRESHOLDS) and p['primary_threshold'] == 1
    assert p['models'] == list(MODELS)
    assert p['primary_support_criterion']['numerical_zero_tolerance'] == 1e-12
    return p

def checked_design(path):
    d = pd.read_csv(path, dtype={'library_name': str, 'condition': str})
    assert len(d) == 48 and d.library_name.is_unique
    assert set(d.genotype_metadata) == set(GENOTYPES)
    assert set(d.condition) == {'CT', *SINGLES, *PAIRS, TARGET}
    assert d.groupby(['genotype_metadata', 'condition']).size().eq(3).all()
    return d

def header_mapping(path, group_design):
    """Read exactly one header line, never a numerical row."""
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8-sig') as h:
        header = h.readline().rstrip('\r\n').split('\t')
    assert header[0] == 'gene', header[0]
    assert len(header) == 25 and len(set(header)) == 25
    assert set(header[1:]) == set(group_design.library_name)
    return header[0], {c: c for c in group_design.library_name}

def identifier_inventory(path, gene_column):
    identifiers = pd.read_csv(path, sep='\t', usecols=[gene_column],
                              dtype={gene_column: str}, keep_default_na=False)[gene_column]
    assert identifiers.notna().all() and identifiers.ne('').all()
    duplicated = identifiers.duplicated(keep=False)
    digest = hashlib.sha256(('\n'.join(identifiers.tolist())+'\n').encode()).hexdigest()
    return {'policy': 'exclude_all_rows_with_duplicated_gene_identifier',
            'source_rows': len(identifiers), 'retained_unique_identifier_rows': int((~duplicated).sum()),
            'excluded_duplicate_rows': int(duplicated.sum()),
            'source_identifier_sha256': digest,
            'excluded_identifier_counts': identifiers[duplicated].value_counts().sort_index().to_dict()}

def validate_inventory(inventory, expected):
    assert inventory['source_rows'] == expected['source_gene_rows_per_matrix']
    assert inventory['excluded_duplicate_rows'] == expected['excluded_rows_per_matrix']
    assert inventory['excluded_identifier_counts'] == {k: expected['rows_per_duplicated_id'] for k in expected['duplicate_gene_ids']}

def checked_amendment(path, protocol_path):
    amendment = json.loads(path.read_text())
    assert amendment['original_protocol_sha256'] == sha256(protocol_path)
    assert amendment['amendment_id'] == 'CODE6_INPUT_AMENDMENT_01'
    return amendment

def read_allowed(path, gene_column, columns):
    """ID-only exclusion precedes numeric QC; only permitted columns are parsed."""
    t = pd.read_csv(path, sep='\t', usecols=[gene_column, *columns],
                    dtype={gene_column: str}, keep_default_na=False)
    assert t[gene_column].notna().all() and t[gene_column].ne('').all()
    t = t.loc[~t[gene_column].duplicated(keep=False)]
    assert t[gene_column].is_unique
    values = t[columns].to_numpy(dtype=float).T
    assert np.isfinite(values).all() and (values >= 0).all()
    return t[gene_column].to_numpy(str), values

def span(single_means, control):
    e = single_means - control
    _, s, vh = np.linalg.svd(e, full_matrices=False)
    tol = max(e.shape) * np.finfo(float).eps * s[0] if s[0] else 0
    r = int(np.sum(s > tol))
    q = vh[:r]
    assert r <= 3 and np.allclose(q @ q.T, np.eye(r), atol=1e-10)
    return q, s, r, tol

def project(x, q):
    return (x @ q.T) @ q

def make_predictions(groups, q):
    c = groups['CT'].mean(axis=0)
    singles = np.stack([groups[k].mean(axis=0) for k in SINGLES])
    pairs = np.stack([groups[k].mean(axis=0) for k in PAIRS])
    donor = pairs.mean(axis=0)
    within = project(donor-c, q)
    outside = donor-c-within
    p = np.stack([c, donor, c+within, c+outside, singles.sum(axis=0)-2*c,
                  pairs.sum(axis=0)-singles.sum(axis=0)+c])
    assert np.allclose(p[1], p[2]+p[3]-p[0], atol=2e-10)
    e2 = p[5]-c
    assert np.allclose(e2-project(e2, q), 3*outside, atol=2e-10)
    return p

def freeze_one(path, design, genotype, threshold, expected_identifier_inventory=None):
    gene_col, mapping = header_mapping(path, design)
    permitted = design[design.condition != TARGET]
    columns = [mapping[x] for x in permitted.library_name]
    targets = [mapping[x] for x in design.loc[design.condition == TARGET, 'library_name']]
    assert set(columns).isdisjoint(targets)
    identifier_qc = identifier_inventory(path, gene_col)
    if expected_identifier_inventory is not None:
        validate_inventory(identifier_qc, expected_identifier_inventory)
    genes, raw = read_allowed(path, gene_col, columns)
    cond = permitted.condition.to_numpy(str)
    lower = np.isin(cond, ['CT', *SINGLES])
    singles_mask = np.isin(cond, SINGLES)
    assert lower.sum() == 12 and singles_mask.sum() == 9 and len(columns) == 21
    keep = raw[lower].mean(axis=0) >= threshold
    dk = raw[singles_mask].mean(axis=0) >= threshold
    assert keep.sum() >= 10 and dk.sum() >= 10
    v = np.log2(1+raw[:, keep])
    groups = {c: v[cond == c] for c in ('CT', *SINGLES, *PAIRS)}
    sm = np.stack([groups[c].mean(axis=0) for c in SINGLES])
    q, s, r, tol = span(sm, groups['CT'].mean(axis=0))
    arrays = {'genes': genes[keep], 'predictions': make_predictions(groups, q),
              'basis': q, 'singular_values': s, 'model_names': np.array(MODELS),
              'target_columns': np.array(targets), 'diagnostic_genes': genes[dk]}
    dv = np.log2(1+raw[:, dk])
    dg = {c: dv[cond == c] for c in ('CT', *SINGLES, *PAIRS)}
    arrays['diagnostic_controls'] = dg['CT']
    arrays['diagnostic_pair_groups'] = np.stack([dg[c] for c in PAIRS])
    dsm = np.stack([dg[c].mean(axis=0) for c in SINGLES])
    cc = [mapping[x] for x in design.loc[design.condition == 'CT', 'library_name']]
    folds = []
    for k in range(3):
        qk, sk, rk, tk = span(dsm, dg['CT'][k])
        arrays[f'diagnostic_basis_{k}'] = qk
        arrays[f'diagnostic_singular_{k}'] = sk
        folds.append({'fold': k+1, 'basis_control': cc[k],
                      'reference_controls': [cc[j] for j in range(3) if j != k],
                      'rank': rk, 'rank_tolerance': tk})
    meta = {'genotype': genotype, 'threshold_fpkm': threshold,
            'primary_threshold': threshold == 1, 'primary_genotype': genotype == 'WT',
            'gene_column': gene_col, 'header_mapping': mapping,
            'allowed_numeric_columns': columns, 'excluded_target_columns': targets,
            'source_gene_rows': identifier_qc['source_rows'], 'identifier_qc': identifier_qc,
            'eligible_unique_identifier_rows': len(genes), 'prediction_gene_count': int(keep.sum()),
            'diagnostic_gene_count': int(dk.sum()), 'span_rank': r, 'rank_tolerance': tol,
            'singular_values': s.tolist(), 'diagnostic_folds': folds}
    return arrays, meta

def save_verified_npz(path, arrays):
    """Atomically persist and round-trip every frozen array before hashing."""
    temporary = path.with_suffix(path.suffix + '.partial')
    with temporary.open('wb') as handle:
        np.savez_compressed(handle, **arrays)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    with np.load(path, allow_pickle=False) as saved:
        assert set(saved.files) == set(arrays)
        for name, expected in arrays.items():
            assert np.array_equal(saved[name], expected), name
    return sha256(path)

def fresh(path):
    assert not path.exists() or not any(path.iterdir()), 'Use a fresh output directory'
    path.mkdir(parents=True, exist_ok=True)

def freeze(args):
    fresh(args.output)
    checked_protocol(args.protocol)
    amendment = checked_amendment(args.amendment, args.protocol)
    design = checked_design(args.design)
    sources = {'WT': args.wt, 'RAP2.3': args.mutant}
    m = {'status': 'FREEZE_IN_PROGRESS', 'started_utc': now(),
         'source': 'GSE330879 author-deposited FPKM',
         'source_hashes': {g: {'path': str(p.resolve()), 'sha256': sha256(p)} for g,p in sources.items()},
         'protocol_sha256': sha256(args.protocol), 'amendment_sha256': sha256(args.amendment),
         'design_sha256': sha256(args.design),
         'freeze_code_sha256': sha256(__file__), 'target_numerical_values_parsed': False,
         'sample_pairing_inferred': False, 'scenarios': []}
    for g in GENOTYPES:
        for t in THRESHOLDS:
            arrays, item = freeze_one(sources[g], design[design.genotype_metadata == g], g, t, amendment['observed_identifier_metadata'])
            key = key_for(g, t)
            name = key+'.npz'
            asset_digest = save_verified_npz(args.output/name, arrays)
            item.update({'key': key, 'assets': name, 'assets_sha256': asset_digest,
                         'atomic_write_round_trip_verified': True})
            m['scenarios'].append(item)
            print(f'FROZEN {key}: {item["prediction_gene_count"]} genes, rank {item["span_rank"]}', flush=True)
    # Recheck all files only after every scenario has been written.
    for item in m['scenarios']:
        path = args.output/item['assets']
        with np.load(path, allow_pickle=False) as saved:
            for name in saved.files:
                saved[name]  # Read every compressed member and verify its CRC.
        assert sha256(path) == item['assets_sha256']
    m.update({'all_asset_integrity_rechecked_at_completion': True,
              'status': 'ALL_SIX_SCENARIOS_FROZEN_BEFORE_TARGET_PARSE', 'completed_utc': now(),
              'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__})
    write_json(args.output/'FREEZE_MANIFEST.json', m)
    print('FREEZE_COMPLETE; no triple-stress numerical values parsed', flush=True)

def corrected_geometry(groups, controls):
    """Symmetrized distinct-reference Gram; sample-covariance noise corrected.

    groups = three pair conditions and triple, each with three RNA pools.
    controls = two different CT pools, neither used to fit the projector.
    The donor mean uses weights 1/3; the triple mean uses weight 1.
    """
    n = groups.shape[-1]
    means = WEIGHTS @ groups.mean(axis=1)
    sampling = groups.var(axis=1, ddof=1).mean(axis=1)/groups.shape[1]
    correction = (WEIGHTS*sampling) @ WEIGHTS.T
    a, b = means-controls[0], means-controls[1]
    gram = (a @ b.T + b @ a.T)/(2*n)-correction
    shared = means-controls.mean(axis=0)
    shared_gram = shared @ shared.T/n-correction
    inflation = np.mean((controls[0]-controls[1])**2)/4
    assert np.allclose(shared_gram-gram, inflation, atol=2e-10)
    distance = gram[0,0]+gram[1,1]-2*gram[0,1]
    direct = np.mean((means[0]-means[1])**2)-np.sum((WEIGHTS[0]-WEIGHTS[1])**2*sampling)
    assert np.isclose(distance, direct, atol=2e-10)
    return {'donor_squared_effect': float(gram[0,0]),
            'target_squared_effect': float(gram[1,1]),
            'donor_target_inner_product': float(gram[0,1]),
            'population_proximity_margin': float(2*gram[0,1]-gram[0,0]),
            'corrected_squared_distance': float(distance),
            'shared_reference_inner_product': float(shared_gram[0,1]),
            'shared_reference_inflation': float(inflation),
            'donor_mean_sampling_variance': float(correction[0,0]),
            'target_mean_sampling_variance': float(correction[1,1])}, gram

def align_target(genes, values, wanted):
    index = pd.Index(genes).get_indexer(wanted)
    assert (index >= 0).all(), 'Frozen feature missing from target matrix'
    return np.log2(1+values[:, index])

def evaluate(args):
    fresh(args.output)
    protocol = checked_protocol(args.protocol)
    amendment = checked_amendment(args.amendment, args.protocol)
    fm = args.frozen/'FREEZE_MANIFEST.json'
    manifest = json.loads(fm.read_text())
    assert manifest['status'] == 'ALL_SIX_SCENARIOS_FROZEN_BEFORE_TARGET_PARSE'
    assert not manifest['target_numerical_values_parsed']
    assert manifest['protocol_sha256'] == sha256(args.protocol)
    assert manifest['amendment_sha256'] == sha256(args.amendment)
    assert manifest['design_sha256'] == sha256(args.design)
    assert manifest['freeze_code_sha256'] == sha256(__file__), 'Code changed after freeze'
    sources = {'WT': args.wt, 'RAP2.3': args.mutant}
    for g,p in sources.items():
        assert manifest['source_hashes'][g]['sha256'] == sha256(p)
        inventory = identifier_inventory(p, 'gene')
        validate_inventory(inventory, amendment['observed_identifier_metadata'])
        expected_qc = next(x['identifier_qc'] for x in manifest['scenarios'] if x['genotype'] == g)
        assert inventory == expected_qc
    for item in manifest['scenarios']:
        assert item['assets_sha256'] == sha256(args.frozen/item['assets'])
    write_json(args.output/'TARGET_OPENING_RECORD.json', {
        'opening_started_utc': now(), 'freeze_manifest_sha256': sha256(fm),
        'protocol_sha256': sha256(args.protocol), 'amendment_sha256': sha256(args.amendment),
        'evaluation_code_sha256': sha256(__file__),
        'source_hashes': manifest['source_hashes'],
        'statement': 'All predictions, panels, models and reference-split bases are frozen; target-only columns are parsed next.'})
    targets = {}
    for g in GENOTYPES:
        item = next(s for s in manifest['scenarios'] if s['genotype'] == g)
        targets[g] = read_allowed(sources[g], item['gene_column'], item['excluded_target_columns'])
    errors, components, geometry, bases = [], [], [], []
    for item in manifest['scenarios']:
        g, t = item['genotype'], item['threshold_fpkm']
        info = {'scenario': item['key'], 'genotype': g, 'threshold_fpkm': t,
                'primary_threshold': t == 1, 'primary_genotype': g == 'WT'}
        with np.load(args.frozen/item['assets'], allow_pickle=False) as a:
            assert a['model_names'].tolist() == list(MODELS)
            genes, raw = targets[g]
            y = align_target(genes, raw, a['genes'])
            pred, q = a['predictions'], a['basis']
            c = pred[0]
            donor = pred[1]-c
            di = project(donor, q)
            do = donor-di
            for j, target_col in enumerate(a['target_columns']):
                meta = {**info, 'target_library': str(target_col), 'genes': len(a['genes']), 'span_rank': len(q)}
                target = y[j]-c
                ti = project(target, q)
                to = target-ti
                mse0 = np.mean(target**2)
                gain = mse0-np.mean((y[j]-pred[1])**2)
                gi = 2*np.mean(ti*di)-np.mean(di**2)
                go = 2*np.mean(to*do)-np.mean(do**2)
                assert np.isclose(gain, gi+go, atol=2e-10)
                for mi, model in enumerate(MODELS):
                    e = pred[mi]-y[j]
                    mse, mae = np.mean(e**2), np.mean(np.abs(e))
                    errors.append({**meta, 'model': model, 'mse': float(mse), 'mae': float(mae),
                                   'mse_gain_over_control': float(mse0-mse),
                                   'mse_gain_percent_over_control': float((mse0-mse)/mse0*100) if mse0 else 0.0})
                for name, gain_c, d, target_c in [('raw', gain, donor, target),
                    ('within_single_span', gi, di, ti), ('outside_single_span', go, do, to)]:
                    components.append({**meta, 'component': name, 'mse_gain_over_control': float(gain_c),
                        'donor_squared_effect': float(np.mean(d**2)), 'target_squared_effect': float(np.mean(target_c**2)),
                        'donor_target_inner_product': float(np.mean(d*target_c)), 'control_total_mse': float(mse0)})
            ydiag = align_target(genes, raw, a['diagnostic_genes'])
            groups = np.concatenate([a['diagnostic_pair_groups'], ydiag[None,:,:]], axis=0)
            controls = a['diagnostic_controls']
            for k in range(3):
                qk = a[f'diagnostic_basis_{k}']
                refs = controls[[j for j in range(3) if j != k]]
                pg, pr = project(groups, qk), project(refs, qk)
                f = item['diagnostic_folds'][k]
                meta = {**info, 'fold': k+1, 'genes': len(a['diagnostic_genes']), 'span_rank': len(qk),
                        'basis_control': f['basis_control'], 'reference_control_a': f['reference_controls'][0],
                        'reference_control_b': f['reference_controls'][1]}
                bases.append(meta)
                grams = {}
                for name, obs, ref in [('raw', groups, refs), ('within_single_span', pg, pr),
                                       ('outside_single_span', groups-pg, refs-pr)]:
                    row, gram = corrected_geometry(obs, ref)
                    geometry.append({**meta, 'component': name, **row})
                    grams[name] = gram
                assert np.allclose(grams['raw'], grams['within_single_span']+grams['outside_single_span'], atol=2e-10)
    ef, cf, gf = map(pd.DataFrame, (errors, components, geometry))
    ef.to_csv(args.output/'library_errors.csv', index=False)
    cf.to_csv(args.output/'library_gain_components.csv', index=False)
    gf.to_csv(args.output/'disjoint_control_geometry.csv', index=False)
    pd.DataFrame(bases).to_csv(args.output/'diagnostic_basis_inventory.csv', index=False)
    keys = ['scenario', 'genotype', 'threshold_fpkm', 'primary_threshold', 'primary_genotype', 'genes', 'span_rank']
    ms = ef.groupby([*keys, 'model'], sort=False)[['mse','mae','mse_gain_over_control','mse_gain_percent_over_control']].mean().reset_index()
    # Summary percentages use the ratio of mean MSEs, not the mean of library ratios.
    controls_by_scenario = ms[ms.model == 'control'].set_index('scenario')['mse']
    ms['mse_gain_percent_over_control'] = [100*gain/controls_by_scenario.loc[scenario] if controls_by_scenario.loc[scenario] else 0.0 for scenario,gain in zip(ms.scenario,ms.mse_gain_over_control)]
    ms.to_csv(args.output/'model_summary.csv', index=False)
    cs = cf.groupby([*keys, 'component'], sort=False)[['mse_gain_over_control','donor_squared_effect','target_squared_effect','donor_target_inner_product','control_total_mse']].agg(['mean','min','max'])
    cs.columns = ['_'.join(c) for c in cs.columns]
    cs.reset_index().to_csv(args.output/'component_summary.csv', index=False)
    outcomes = []
    tol = protocol['primary_support_criterion']['numerical_zero_tolerance']
    for g in GENOTYPES:
        for t in THRESHOLDS:
            key = key_for(g,t)
            c, gg = cf[cf.scenario == key], gf[gf.scenario == key]
            full_gain = float(c.loc[c.component == 'raw', 'mse_gain_over_control'].mean())
            outside = c.loc[c.component == 'outside_single_span', 'mse_gain_over_control']
            margin = gg.loc[gg.component == 'outside_single_span', 'population_proximity_margin']
            criteria = {'full_centroid_mean_gain_positive': full_gain > tol,
                        'outside_mean_gain_positive': float(outside.mean()) > tol,
                        'outside_gain_positive_in_all_three_target_libraries': bool((outside > tol).all()),
                        'outside_margin_positive_in_all_three_disjoint_control_folds': bool((margin > tol).all())}
            outcomes.append({'genotype': g, 'threshold_fpkm': t, 'primary_endpoint': g == 'WT' and t == 1,
                'classification': 'DESCRIPTIVE_SUPPORT_FOR_NARROWED_THREE_FACTOR_HYPOTHESIS' if all(criteria.values()) else 'PREDECLARED_SUPPORT_CRITERION_NOT_MET',
                'criteria': criteria, 'full_centroid_mean_mse_gain': full_gain,
                'outside_mean_mse_gain': float(outside.mean()), 'outside_library_mse_gains': outside.tolist(),
                'outside_disjoint_control_margins': margin.tolist()})
    summary = {'status': 'EXTERNAL_THREE_FACTOR_EVALUATION_COMPLETE', 'completed_utc': now(),
        'experiment': 'GSE330879', 'primary_genotype': 'WT', 'secondary_genotype': 'RAP2.3',
        'primary_filter_fpkm': 1, 'outcomes': outcomes, 'formal_significance_test': False,
        'folds_are_independent_replicates': False, 'freeze_manifest_sha256': sha256(fm),
        'protocol_sha256': sha256(args.protocol), 'amendment_sha256': sha256(args.amendment),
        'identifier_inventories': {g: next(x['identifier_qc'] for x in manifest['scenarios'] if x['genotype'] == g) for g in GENOTYPES},
        'code_sha256': sha256(__file__),
        'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__,
        'scope': 'Separate experiment, same species/research group; generalized pair-to-triple test, not exact replication of the discovery four-factor selected-donor finding.',
        'limits': 'Three RNA pools per condition; source-wide preprocessing; noisy fitted single-stress span; outside-span is not pure higher-order interaction; no causal mechanism or crop generality inferred.',
        'output_hashes': {p.name: sha256(p) for p in args.output.iterdir() if p.is_file()}}
    write_json(args.output/'RUN_SUMMARY.json', summary)
    print(json.dumps(summary, indent=2))

def selftest():
    """Target exclusion, known projection recovery, and Gram identities."""
    rng = np.random.default_rng(34071)
    with tempfile.TemporaryDirectory() as td:
        directory = Path(td)
        conditions = ['CT','EL','HS','EL+HS','S','S+EL','S+HS',TARGET]
        records = []
        table = {'gene': [f'g{i}' for i in range(70)]}
        for rep in range(1,4):
            for j,condition in enumerate(conditions,1):
                library = f'C{rep}{j}'
                records.append({'library_name': library, 'condition': condition})
                table[library] = rng.lognormal(1,0.8,70)
        design = pd.DataFrame(records)
        clean_path = directory/'clean.tsv'
        pd.DataFrame(table).to_csv(clean_path, sep='\t', index=False)
        clean, _ = freeze_one(clean_path, design, 'WT', 1)
        for poison in ('UNREADABLE_TARGET', np.nan, 1e100):
            bad = pd.DataFrame(table)
            for column in design.loc[design.condition == TARGET, 'library_name']:
                bad[column] = poison
            path = directory/'poison.tsv'
            bad.to_csv(path, sep='\t', index=False)
            poisoned, _ = freeze_one(path, design, 'WT', 1)
            assert clean.keys() == poisoned.keys()
            for k in clean:
                assert np.array_equal(clean[k],poisoned[k]), k
        duplicate = pd.DataFrame(table)
        duplicate.loc[0, 'gene'] = 'g1'
        duplicate_path = directory/'duplicate.tsv'
        duplicate.to_csv(duplicate_path, sep='\t', index=False)
        ids, _ = read_allowed(duplicate_path, 'gene', design.loc[design.condition != TARGET, 'library_name'].tolist())
        iq = identifier_inventory(duplicate_path, 'gene')
        assert len(ids) == 68 and 'g1' not in ids and iq['excluded_duplicate_rows'] == 2
        observations, refs = rng.normal(size=(4,3,70)), rng.normal(size=(2,70))
        _, gram = corrected_geometry(observations,refs)
        q,_,_,_ = span(rng.normal(size=(3,70)),np.zeros(70))
        _, inside = corrected_geometry(project(observations,q),project(refs,q))
        _, outside = corrected_geometry(observations-project(observations,q),refs-project(refs,q))
        assert np.allclose(gram,inside+outside,atol=1e-12)
        groups = {'CT': np.zeros((3,70))}
        for i,c in enumerate(SINGLES):
            v = np.zeros((3,70)); v[:,i] = 1
            groups[c] = v
        for c in PAIRS:
            v = np.zeros((3,70)); v[:,:4] = 1
            groups[c] = v
        q,_,_,_ = span(np.stack([groups[c].mean(axis=0) for c in SINGLES]), np.zeros(70))
        p = make_predictions(groups,q)
        assert np.max(np.abs(p[3,:3])) < 1e-12 and p[3,3] == 1
        assert np.allclose(p[1],p[2]+p[3]-p[0])
    print('SELFTEST_PASS: target text/NaN/extreme poisoning leaves all frozen arrays unchanged; all duplicate-ID occurrences excluded; projection and covariance identities hold')

def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command',required=True)
    sub.add_parser('selftest')
    for command in ('freeze','evaluate'):
        parser = sub.add_parser(command)
        for arg in ('wt','mutant','design','protocol','amendment','output'):
            parser.add_argument('--'+arg,type=Path,required=True)
        if command == 'evaluate':
            parser.add_argument('--frozen',type=Path,required=True)
    args = p.parse_args()
    if args.command == 'selftest': selftest()
    elif args.command == 'freeze': freeze(args)
    else: evaluate(args)

if __name__ == '__main__':
    main()
