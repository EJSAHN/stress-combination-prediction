#!/usr/bin/env python3
"""Post-result sensitivity to expanded observed single-stress expression spans.

No target-driven rank selection, gene bootstrap, P-values or true-rank claim.
CODE6 targets were already opened. This does not constitute new validation.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
import pandas as pd

SINGLES = ('S', 'EL', 'HS')
PAIRS = ('EL+HS', 'S+EL', 'S+HS')
TARGET = 'S+EL+HS'
CONDITIONS = ('CT', *SINGLES, *PAIRS, TARGET)
WEIGHTS = np.array([[1/3, 1/3, 1/3, 0], [0, 0, 0, 1.]])
TOL = 1e-12

def now():
    return datetime.now(timezone.utc).isoformat()

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()

def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')

def fresh(path):
    if path.exists() and any(path.iterdir()):
        raise ValueError('Use a fresh output directory; existing results are preserved.')
    path.mkdir(parents=True, exist_ok=True)

def basis(rows, maximum_rank):
    _, s, vt = np.linalg.svd(rows, full_matrices=False)
    tolerance = max(rows.shape) * np.finfo(float).eps * s[0] if s[0] else 0.0
    rank = int(np.sum(s > tolerance))
    assert rank <= maximum_rank, (rank, maximum_rank)
    q = vt[:rank]
    assert np.allclose(q @ q.T, np.eye(rank), atol=1e-10)
    return q, s, float(tolerance)

def project(x, q):
    return (x @ q.T) @ q

def nested_check(small, large):
    residual = small - project(small, large)
    norm = float(np.linalg.norm(residual, ord='fro'))
    assert norm < 1e-8, norm
    cosines = np.linalg.svd(small @ large.T, compute_uv=False)
    angles = np.arccos(np.clip(cosines, -1.0, 1.0))
    return {'nested_projection_residual_frobenius': norm,
            'minimum_principal_angle_cosine': float(cosines.min()) if len(cosines) else 1.0,
            'maximum_principal_angle_radians': float(angles.max()) if len(angles) else 0.0}

def span_bundle(groups, reference=None, include_control_residuals=True):
    c = groups['CT'].mean(axis=0) if reference is None else reference
    means = np.stack([groups[k].mean(axis=0) for k in SINGLES])
    individual = np.concatenate([groups[k] for k in SINGLES], axis=0)
    q3, s3, t3 = basis(means-c, 3)
    q9, s9, t9 = basis(individual-c, 9)
    bases = {'V3': q3, 'V9': q9}
    singular = {'V3': s3, 'V9': s9}
    tolerances = {'V3': t3, 'V9': t9}
    checks = {'V3_in_V9': nested_check(q3, q9)}
    if include_control_residuals:
        assert reference is None, 'Expanded controls are forbidden for disjoint-reference folds.'
        q11, s11, t11 = basis(np.concatenate([individual-c, groups['CT']-c]), 11)
        bases['V11'], singular['V11'], tolerances['V11'] = q11, s11, t11
        checks['V9_in_V11'] = nested_check(q9, q11)
    return bases, singular, tolerances, checks

def component_vectors(vector, bases):
    a, b, c = [project(vector, bases[k]) for k in ('V3', 'V9', 'V11')]
    parts = {'within_V3': a, 'V9_minus_V3': b-a,
             'V11_minus_V9': c-b, 'outside_V11': vector-c}
    assert np.allclose(sum(parts.values()), vector, atol=2e-10)
    keys = list(parts)
    for i, k in enumerate(keys):
        for l in keys[i+1:]:
            assert abs(float(np.mean(parts[k]*parts[l]))) < 2e-10
    return parts

def corrected_geometry(groups, references):
    """The two reference libraries are distinct and excluded from basis fitting."""
    n_genes = groups.shape[-1]
    state_means = WEIGHTS @ groups.mean(axis=1)
    variance_of_means = groups.var(axis=1, ddof=1).mean(axis=1)/groups.shape[1]
    correction = (WEIGHTS*variance_of_means) @ WEIGHTS.T
    a = state_means-references[0]
    b = state_means-references[1]
    gram = (a @ b.T + b @ a.T)/(2*n_genes)-correction
    distance = gram[0,0]+gram[1,1]-2*gram[0,1]
    direct = np.mean((state_means[0]-state_means[1])**2)-np.sum((WEIGHTS[0]-WEIGHTS[1])**2*variance_of_means)
    assert np.isclose(distance, direct, atol=2e-10)
    return {'donor_squared_effect': float(gram[0,0]),
            'target_squared_effect': float(gram[1,1]),
            'donor_target_inner_product': float(gram[0,1]),
            'population_proximity_margin': float(2*gram[0,1]-gram[0,0]),
            'corrected_squared_distance': float(distance),
            'donor_mean_sampling_variance': float(correction[0,0]),
            'target_mean_sampling_variance': float(correction[1,1])}, gram

def load_source(path, design, expected_ids):
    columns = design.library_name.tolist()
    frame = pd.read_csv(path, sep='\t', dtype={'gene': str}, keep_default_na=False)
    assert frame.columns[0] == 'gene' and set(frame.columns[1:]) == set(columns)
    assert len(frame.columns) == 25
    duplicated = frame.gene.duplicated(keep=False)
    counts = frame.loc[duplicated, 'gene'].value_counts().sort_index().to_dict()
    assert counts == expected_ids['excluded_identifier_counts']
    assert int(duplicated.sum()) == expected_ids['excluded_duplicate_rows']
    assert len(frame) == expected_ids['source_rows']
    idhash = hashlib.sha256(('\n'.join(frame.gene.tolist())+'\n').encode()).hexdigest()
    assert idhash == expected_ids['source_identifier_sha256']
    frame = frame.loc[~duplicated]
    assert frame.gene.ne('').all() and frame.gene.is_unique
    raw = frame[columns].to_numpy(float).T
    assert np.isfinite(raw).all() and (raw >= 0).all()
    return frame.gene.to_numpy(str), raw, design.condition.to_numpy(str)

def panel_groups(genes, raw, conditions, wanted):
    idx = pd.Index(genes).get_indexer(wanted)
    assert (idx >= 0).all()
    y = np.log2(1+raw[:, idx])
    groups = {k: y[conditions == k] for k in CONDITIONS}
    assert all(v.shape[0] == 3 for v in groups.values())
    return groups

def summarize(frame, by, values, destination):
    frame.groupby(by, sort=False)[values].mean().reset_index().to_csv(destination, index=False)

def run(args):
    fresh(args.output)
    protocol = json.loads(args.protocol.read_text())
    assert protocol['source'] == 'GSE330879'
    assert protocol['primary_genotype'] == 'WT' and protocol['secondary_genotype'] == 'RAP2.3'
    assert protocol['primary_threshold_fpkm'] == 1 and protocol['numerical_zero_tolerance'] == TOL
    root = args.code6_root
    manifest_path = root/'frozen_verified/FREEZE_MANIFEST.json'
    manifest = json.loads(manifest_path.read_text())
    previous = json.loads((root/'results/RUN_SUMMARY.json').read_text())
    design_path = root/'fallback_metadata/sample_inventory.csv'
    design = pd.read_csv(design_path)
    assert len(design) == 48 and design.library_name.is_unique
    assert sha256(design_path) == manifest['design_sha256']
    assert sha256(root/'CODE6_PROTOCOL.json') == manifest['protocol_sha256']
    assert sha256(root/'INPUT_AMENDMENT_01.json') == manifest['amendment_sha256']
    assert sha256(root/'code/external_three_factor.py') == manifest['freeze_code_sha256']
    assert sha256(manifest_path) == previous['freeze_manifest_sha256'] == protocol['code6_freeze_manifest_sha256']
    assert manifest['protocol_sha256'] == protocol['code6_protocol_sha256']
    assert manifest['amendment_sha256'] == protocol['code6_amendment_sha256']
    for name in ('library_errors.csv', 'disjoint_control_geometry.csv'):
        assert sha256(root/'results'/name) == previous['output_hashes'][name]
    previous_errors = pd.read_csv(root/'results/library_errors.csv')
    previous_geometry = pd.read_csv(root/'results/disjoint_control_geometry.csv')
    sources = {}
    inputs = {'code6_freeze_manifest': sha256(manifest_path),
              'code6_run_summary': sha256(root/'results/RUN_SUMMARY.json'),
              'code6_code': manifest['freeze_code_sha256'], 'code7_protocol': sha256(args.protocol)}
    for genotype, item in manifest['source_hashes'].items():
        path = root/'input'/Path(item['path']).name
        assert sha256(path) == item['sha256']
        inputs[str(path)] = item['sha256']
        idqc = previous['identifier_inventories'][genotype]
        sources[genotype] = load_source(path, design[design.genotype_metadata == genotype], idqc)
    errors, components, geometry, ranks, inventory, crosschecks = [], [], [], [], [], []
    for item in manifest['scenarios']:
        genotype = item['genotype']
        threshold = item['threshold_fpkm']
        asset = root/'frozen_verified'/item['assets']
        assert sha256(asset) == item['assets_sha256']
        inputs[str(asset)] = item['assets_sha256']
        genes, raw, conditions = sources[genotype]
        common = {'scenario': item['key'], 'genotype': genotype,
                  'threshold_fpkm': threshold, 'primary_endpoint': genotype == 'WT' and threshold == 1.0}
        with np.load(asset, allow_pickle=False) as fixed:
            main_genes = fixed['genes']
            dgenes = fixed['diagnostic_genes']
            expected = genes[raw[np.isin(conditions, ['CT', *SINGLES])].mean(axis=0) >= threshold]
            expected_d = genes[raw[np.isin(conditions, SINGLES)].mean(axis=0) >= threshold]
            assert np.array_equal(main_genes, expected) and np.array_equal(dgenes, expected_d)
            groups = panel_groups(genes, raw, conditions, main_genes)
            bases, singular, tolerance, nesting = span_bundle(groups)
            c = groups['CT'].mean(axis=0)
            donor = np.stack([groups[k].mean(axis=0) for k in PAIRS]).mean(axis=0)
            d = donor-c
            dparts = component_vectors(d, bases)
            predictions = {'control': c, 'full_pair_centroid': donor}
            for name,q in bases.items():
                predictions['within_'+name] = c+project(d,q)
            assert np.allclose(predictions['control'], fixed['predictions'][0], atol=2e-10)
            assert np.allclose(predictions['full_pair_centroid'], fixed['predictions'][1], atol=2e-10)
            assert np.allclose(predictions['within_V3'], fixed['predictions'][2], atol=2e-10)
            assert np.allclose(project(fixed['basis'], bases['V3']), fixed['basis'], atol=2e-10)
            info = {**common, 'genes':len(main_genes), 'rank_V3':len(bases['V3']), 'rank_V9':len(bases['V9']), 'rank_V11':len(bases['V11'])}
            for name,q in bases.items():
                ranks.append({**common, 'panel':'main', 'span':name, 'rank':len(q), 'genes':len(main_genes), 'svd_tolerance':tolerance[name], 'smallest_retained_singular_value':float(singular[name][len(q)-1]) if len(q) else 0.0})
            inventory.append({**info, 'nesting':nesting, 'basis_libraries':'CT and nine single-stress libraries only; no pair or triple libraries'})
            for j, target in enumerate(groups[TARGET]):
                col = str(fixed['target_columns'][j])
                m = {**info, 'target_library':col}
                y = target-c
                yparts = component_vectors(y,bases)
                mse0 = float(np.mean(y*y))
                fullgain = float(np.mean(2*d*y-d*d))
                gain_sum = 0.0
                outside_gain = float(np.mean(2*(d-project(d,bases['V3']))*(y-project(y,bases['V3']))-(d-project(d,bases['V3']))**2))
                for name, dp in dparts.items():
                    yp = yparts[name]
                    gain = float(np.mean(2*dp*yp-dp*dp))
                    gain_sum += gain
                    components.append({**m, 'component':name, 'mse_gain_over_control':gain,
                        'donor_squared_effect':float(np.mean(dp*dp)), 'target_squared_effect':float(np.mean(yp*yp)),
                        'donor_target_inner_product':float(np.mean(dp*yp))})
                assert abs(gain_sum-fullgain) < 2e-10
                outside_parts = sum(float(np.mean(2*dparts[k]*yparts[k]-dparts[k]**2)) for k in ('V9_minus_V3','V11_minus_V9','outside_V11'))
                assert abs(outside_gain-outside_parts) < 2e-10
                crosschecks.append({**m, 'full_gain_identity_error':gain_sum-fullgain, 'outside_gain_identity_error':outside_parts-outside_gain})
                for name,pred in predictions.items():
                    mse = float(np.mean((target-pred)**2))
                    mae = float(np.mean(np.abs(target-pred)))
                    errors.append({**m, 'model':name, 'mse':mse, 'mae':mae, 'mse_gain_over_control':mse0-mse})
                    oldname = {'control':'control','full_pair_centroid':'full_pair_centroid','within_V3':'within_only'}.get(name)
                    if oldname:
                        old = previous_errors[(previous_errors.scenario == item['key']) & (previous_errors.target_library == col) & (previous_errors.model == oldname)].iloc[0]
                        assert abs(old.mse-mse) < 2e-10 and abs(old.mae-mae) < 2e-10
            dg = panel_groups(genes,raw,conditions,dgenes)
            pair_target = np.stack([dg[k] for k in (*PAIRS,TARGET)])
            controls = dg['CT']
            for k in range(3):
                db, ds, dt, checks = span_bundle(dg, reference=controls[k], include_control_residuals=False)
                refs = controls[[i for i in range(3) if i != k]]
                for name,q in db.items():
                    ranks.append({**common, 'panel':'disjoint', 'fold':k+1, 'span':name, 'rank':len(q), 'genes':len(dgenes), 'svd_tolerance':dt[name], 'smallest_retained_singular_value':float(ds[name][len(q)-1]) if len(q) else 0.0})
                    inside, refinside = project(pair_target,q),project(refs,q)
                    grams = {}
                    for component, data, ref in [('raw',pair_target,refs),('within',inside,refinside),('outside',pair_target-inside,refs-refinside)]:
                        row,gram = corrected_geometry(data,ref)
                        grams[component] = gram
                        geometry.append({**common, 'genes':len(dgenes), 'fold':k+1, 'span':name, 'rank':len(q), 'component':component,
                            'basis_control':item['diagnostic_folds'][k]['basis_control'],
                            'reference_control_a':item['diagnostic_folds'][k]['reference_controls'][0],
                            'reference_control_b':item['diagnostic_folds'][k]['reference_controls'][1], **row})
                        if name == 'V3':
                            oldcomponent = {'raw':'raw','within':'within_single_span','outside':'outside_single_span'}[component]
                            old = previous_geometry[(previous_geometry.scenario == item['key']) & (previous_geometry.fold == k+1) & (previous_geometry.component == oldcomponent)].iloc[0]
                            for value in row:
                                assert abs(row[value]-old[value]) < 2e-10, value
                    assert np.allclose(grams['raw'],grams['within']+grams['outside'],atol=2e-10)
        print(f'EMPIRICAL {item["key"]}: fixed panels and nested-span identities verified',flush=True)
    ef,cf,gf,rf = map(pd.DataFrame,(errors,components,geometry,ranks))
    ef.to_csv(args.output/'library_errors.csv',index=False)
    cf.to_csv(args.output/'library_gain_components.csv',index=False)
    gf.to_csv(args.output/'disjoint_control_geometry.csv',index=False)
    rf.to_csv(args.output/'span_inventory.csv',index=False)
    pd.DataFrame(crosschecks).to_csv(args.output/'algebra_checks.csv',index=False)
    by=['scenario','genotype','threshold_fpkm','primary_endpoint','genes','rank_V3','rank_V9','rank_V11']
    summarize(ef,[*by,'model'],['mse','mae','mse_gain_over_control'],args.output/'model_summary.csv')
    summarize(cf,[*by,'component'],['mse_gain_over_control','donor_squared_effect','target_squared_effect','donor_target_inner_product'],args.output/'component_summary.csv')
    outcomes=[]
    for item in manifest['scenarios']:
        key=item['key']; comp=cf[cf.scenario==key]; geom=gf[gf.scenario==key]
        outside=comp.loc[comp.component=='outside_V11','mse_gain_over_control']
        margin=geom.loc[(geom.span=='V9')&(geom.component=='outside'),'population_proximity_margin']
        original=comp.loc[comp.component!='within_V3','mse_gain_over_control'].sum()/3
        outcomes.append({'scenario':key,'genotype':item['genotype'],'threshold_fpkm':item['threshold_fpkm'],
            'primary_endpoint':item['genotype']=='WT' and item['threshold_fpkm']==1,
            'original_outside_V3_mean_gain':float(original),'outside_V11_mean_gain':float(outside.mean()),
            'outside_V11_library_gains':outside.tolist(),'outside_V9_disjoint_control_margins':margin.tolist(),
            'outside_V11_positive_in_all_libraries':bool((outside>TOL).all()),
            'outside_V9_disjoint_margin_positive_in_all_folds':bool((margin>TOL).all()),
            'empirical_status':'SURVIVES_OBSERVED_VARIATION_EXPANSION' if (outside>TOL).all() and (margin>TOL).all() else 'SENSITIVE_TO_OBSERVED_VARIATION_EXPANSION'})
    write_json(args.output/'BASIS_NESTING_RECORD.json',inventory)
    summary={'status':'POST_RESULT_EMPIRICAL_SPAN_DIAGNOSTIC_COMPLETE','completed_utc':now(),
        'targets_previously_opened':True,'target_used_to_choose_basis':False,'formal_significance_test':False,
        'observed_span_not_true_biological_subspace':True,'outcomes':outcomes,'input_hashes':inputs,
        'code_sha256':sha256(__file__),'runtime':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__},
        'output_hashes':{p.name:sha256(p) for p in args.output.iterdir() if p.is_file()},
        'interpretation':'Observed single-library/control variability enlarges a descriptive nuisance span; it neither identifies pure noise nor proves true biological rank. Disjoint diagnostics use V3/V9 only and different fixed gene panels.'}
    write_json(args.output/'RUN_SUMMARY.json',summary)
    print(json.dumps(summary,indent=2))

def selftest():
    rng=np.random.default_rng(472901)
    genes=47
    groups={k:rng.normal(size=(3,genes)) for k in CONDITIONS}
    bases,_,_,checks=span_bundle(groups)
    assert [len(bases[k]) for k in ('V3','V9','V11')]==[3,9,11]
    d=rng.normal(size=genes);y=rng.normal(size=genes)
    dp,yp=component_vectors(d,bases),component_vectors(y,bases)
    gains={k:np.mean(2*dp[k]*yp[k]-dp[k]**2) for k in dp}
    assert np.allclose(sum(gains.values()),np.mean(2*d*y-d*d),atol=1e-11)
    assert np.allclose(sum(v for k,v in gains.items() if k!='within_V3'),np.mean(2*(d-project(d,bases['V3']))*(y-project(y,bases['V3']))-(d-project(d,bases['V3']))**2),atol=1e-11)
    altered={k:v.copy() for k,v in groups.items()}
    for k in (*PAIRS,TARGET):altered[k]=rng.normal(100,30,size=(3,genes))
    changed,_,_,_=span_bundle(altered)
    for k in bases:assert np.array_equal(bases[k],changed[k])
    for k in range(3):
        db,_,_,_=span_bundle(groups,reference=groups['CT'][k],include_control_residuals=False)
        assert set(db)=={'V3','V9'}
        for q in db.values():
            obs=np.stack([groups[x] for x in (*PAIRS,TARGET)])
            refs=groups['CT'][[i for i in range(3) if i!=k]]
            _,raw=corrected_geometry(obs,refs)
            _,inside=corrected_geometry(project(obs,q),project(refs,q))
            _,outside=corrected_geometry(obs-project(obs,q),refs-project(refs,q))
            assert np.allclose(raw,inside+outside,atol=1e-11)
    # Degenerate constant data must yield rank zero, not an artificial rank.
    zero={k:np.ones((3,genes)) for k in CONDITIONS}
    z,_,_,_=span_bundle(zero)
    assert all(len(q)==0 for q in z.values())
    print('SELFTEST_PASS: nested ranks, gain partition, donor/target exclusion, disjoint additivity and rank-zero edge case')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    sub.add_parser('selftest')
    runparser=sub.add_parser('run')
    for name in ('code6-root','protocol','output'):runparser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    selftest() if args.command=='selftest' else run(args)

if __name__=='__main__':main()
