#!/usr/bin/env python3
"""Exploratory target-excluded donor selection; known data, not a new holdout."""
import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
from pathlib import Path
import platform
import numpy as np
import pandas as pd

FACTORS = ['salt', 'paraquat', 'high_light', 'heat']
MODELS = ['control_mean', 'single_effect_sum', 'pairwise', 'all_donors_mean',
          'all_donors_half_shrink', 'common_low1', 'common_low2', 'raw_low1',
          'raw_low2', 'residual_low1', 'residual_low2', 'common_low2_residual_only']
RULES = ['common_low1', 'common_low2', 'raw_low1', 'raw_low2', 'residual_low1', 'residual_low2']
RANKS = [2, 1]

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()

def require(condition, message):
    if not condition:
        raise ValueError(message)

def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')

def source_check(args):
    protocol = json.loads(args.protocol.read_text())
    require(protocol['models'] == MODELS and protocol['factors'] == FACTORS, 'Changed model/factor order')
    for path in [args.data, args.design]:
        require(sha256(path) == protocol['input_hashes'][path.name], 'Source hash mismatch: ' + path.name)
    design = pd.read_csv(args.design)
    require(len(design) == 48 and design.matrix_column.is_unique, 'Invalid sample design')
    require(design.groupby('condition').size().eq(3).all(), 'Expected three RNA pools per cell')
    require(design[FACTORS].sum(axis=1).eq(design.order).all(), 'Inconsistent factor order')
    require(set(design.loc[design.order >= 3, 'condition']) == set(protocol['targets']), 'Changed targets')
    return protocol, design

def immediate_donors(design, target):
    rows = design.loc[design.condition.eq(target)]
    require(len(rows) == 3, 'Invalid target')
    bits = rows[FACTORS].iloc[0].to_numpy(int)
    cells = design.drop_duplicates('condition')
    candidates = cells.loc[(cells.order == bits.sum()-1) &
                           (cells[FACTORS].to_numpy(int) <= bits).all(axis=1), 'condition'].tolist()
    require(len(candidates) == bits.sum() and target not in candidates, 'Invalid immediate subsets')
    return sorted(candidates), bits

def load_task(data, design, target):
    """Parse only current-task allowed columns; target values are never loaded."""
    donors, bits = immediate_donors(design, target)
    used = design.loc[design.order.le(2) | design.condition.isin(donors)].copy()
    columns = used.matrix_column.tolist()
    target_columns = design.loc[design.condition.eq(target), 'matrix_column'].tolist()
    require(not set(columns).intersection(target_columns), 'Target leaked into parser')
    require(len(columns) == (45 if bits.sum() == 4 else 33), 'Unexpected information budget')
    table = pd.read_csv(data, sep='\t', usecols=['genes'] + columns)
    require(set(table.columns) == set(['genes'] + columns), 'Parser returned forbidden columns')
    require(table.genes.is_unique, 'Duplicate source genes')
    return table, used, donors, bits, target_columns

def fit_lower_order(table, used):
    counts = table[used.matrix_column].to_numpy(float).T
    require(np.isfinite(counts).all() and (counts >= 0).all(), 'Invalid normalized counts')
    lower = used.order.le(2).to_numpy()
    require(lower.sum() == 33, 'Invalid lower-order training size')
    keep = counts[lower].mean(axis=0) >= 10
    genes = table.loc[keep, 'genes'].to_numpy(str)
    expression = np.log2(1 + counts[:, keep])
    factor_array = used[FACTORS].to_numpy(int)
    keys = sorted(set(map(tuple, factor_array[lower])))
    means = {key: expression[(factor_array == key).all(axis=1)].mean(axis=0) for key in keys}
    control = means[(0, 0, 0, 0)]
    _, singular, basis = np.linalg.svd(np.stack([means[key]-control for key in keys]), full_matrices=False)
    for row in basis:
        if row[np.argmax(np.abs(row))] < 0:
            row *= -1
    require(len(genes) == 21495, 'Unexpected training filter size')
    return genes, expression, means, control, basis[:2], singular

def factorial_predictions(means, bits):
    """Exact degree<=2 extrapolation, including intercept and all pair contrasts."""
    control = means[(0, 0, 0, 0)]
    present = np.flatnonzero(bits)
    singles, pairs = [], []
    for i in present:
        key = np.zeros(4, dtype=int); key[i] = 1
        singles.append(means[tuple(key)])
    for i, j in itertools.combinations(present, 2):
        key = np.zeros(4, dtype=int); key[[i, j]] = 1
        pairs.append(means[tuple(key)])
    m = len(present)
    return (sum(singles)-(m-1)*control,
            sum(pairs)-(m-2)*sum(singles)+((m-1)*(m-2)/2)*control)

def donor_scores(donor_means, control, basis):
    names = sorted(donor_means)
    effects = np.stack([donor_means[n]-control for n in names])
    common = (effects @ basis.T) @ basis
    residual = effects-common
    scores = {n: {'raw': float(np.mean(e**2)), 'common': float(np.mean(c**2)),
                  'residual': float(np.mean(r**2))}
              for n, e, c, r in zip(names, effects, common, residual)}
    require(all(np.isclose(x['raw'], x['common']+x['residual'], atol=1e-10) for x in scores.values()),
            'Norm decomposition failed')
    return scores

def select_donors(scores, rule):
    component, number = rule.split('_low')
    return sorted(scores, key=lambda name: (scores[name][component], name))[:int(number)]

def make_predictions(donor_means, control, basis, single, pairwise):
    scores = donor_scores(donor_means, control, basis)
    all_mean = np.mean(list(donor_means.values()), axis=0)
    predictions = {'control_mean': control, 'single_effect_sum': single, 'pairwise': pairwise,
                   'all_donors_mean': all_mean, 'all_donors_half_shrink': control+0.5*(all_mean-control)}
    selected = {}
    for rule in RULES:
        selected[rule] = select_donors(scores, rule)
        predictions[rule] = np.mean([donor_means[n] for n in selected[rule]], axis=0)
    selected['common_low2_residual_only'] = selected['common_low2']
    effect = predictions['common_low2']-control
    predictions['common_low2_residual_only'] = control+effect-(effect @ basis.T) @ basis
    return np.stack([predictions[name] for name in MODELS]), scores, selected

def synthetic_test():
    rng = np.random.default_rng(51028)
    beta = rng.normal(size=(11, 9))
    def design(bits):
        return np.array([1, *bits, *[bits[i]*bits[j] for i,j in itertools.combinations(range(4),2)]])
    cells = list(itertools.product((0,1), repeat=4))
    means = {b: design(b) @ beta for b in cells if sum(b)<=2}
    error = 0.0
    for b in cells:
        if sum(b)>=3:
            _, pred = factorial_predictions(means, np.array(b))
            error = max(error, float(np.max(np.abs(pred-design(b) @ beta))))
    require(error < 1e-10, 'Synthetic degree-two extrapolation failed')
    tie_scores = {'z': {'raw': 1.0}, 'a': {'raw': 1.0}, 'm': {'raw': 2.0}}
    require(select_donors(tie_scores, 'raw_low2') == ['a','z'], 'Tie break failed')
    return {'pairwise_max_error': error, 'lexicographic_tie_break': 'PASS'}

def freeze(args):
    protocol, design = source_check(args)
    require(not args.output.exists() or not any(args.output.iterdir()), 'Use a fresh freeze directory')
    args.output.mkdir(parents=True, exist_ok=True)
    scores_rows, selected_rows, stability_rows, access, target_predictions = [], [], [], [], []
    shared_genes = shared_basis = shared_control = None
    for target in protocol['targets']:
        table, used, donors, bits, target_columns = load_task(args.data, design, target)
        genes, expression, means, control, basis, singular = fit_lower_order(table, used)
        if shared_genes is None:
            shared_genes, shared_basis, shared_control = genes, basis, control
        else:
            require(np.array_equal(genes,shared_genes) and np.allclose(basis,shared_basis) and
                    np.array_equal(control,shared_control), 'Lower-order assets vary across tasks')
        groups = {n: expression[used.condition.eq(n).to_numpy()] for n in donors}
        donor_means = {n: x.mean(axis=0) for n,x in groups.items()}
        single, pairwise = factorial_predictions(means,bits)
        rank_predictions = []
        for rank in RANKS:
            q = basis[:rank]
            predictions,scores,selected = make_predictions(donor_means,control,q,single,pairwise)
            rank_predictions.append(predictions)
            for n in donors:
                scores_rows.append({'target':target,'rank':rank,'donor':n,**scores[n]})
            for model in MODELS:
                chosen = selected.get(model, donors if model.startswith('all_donors') else [])
                selected_rows.append({'target':target,'rank':rank,'model':model,
                    'selected_donors_ranked':';'.join(chosen),'selected_donors_membership':';'.join(sorted(chosen))})
            for omitted_donor in donors:
                pool_names = used.loc[used.condition.eq(omitted_donor),'matrix_column'].tolist()
                for omitted_index,pool in enumerate(pool_names):
                    alternate = dict(donor_means)
                    alternate[omitted_donor] = np.delete(groups[omitted_donor],omitted_index,axis=0).mean(axis=0)
                    alternate_scores = donor_scores(alternate,control,q)
                    for rule in RULES:
                        alt_selected = select_donors(alternate_scores,rule)
                        stability_rows.append({'target':target,'rank':rank,'rule':rule,
                            'omitted_donor':omitted_donor,'omitted_library':pool,
                            'full_selected':';'.join(sorted(selected[rule])),
                            'omission_selected':';'.join(sorted(alt_selected)),
                            'same_selection_as_full':sorted(selected[rule])==sorted(alt_selected)})
        target_predictions.append(np.stack(rank_predictions))
        access.append({'target':target,'target_factor_bits':bits.tolist(),'candidate_donors':donors,
            'allowed_numerical_columns':used.matrix_column.tolist(),'target_columns_not_parsed':target_columns,
            'other_high_order_columns_not_parsed':design.loc[(design.order>=3)&~design.condition.isin(donors)&
                                                          ~design.condition.eq(target),'matrix_column'].tolist()})
    numeric = np.stack(target_predictions)
    require(numeric.shape==(5,2,12,21495) and np.isfinite(numeric).all(),'Invalid prediction arrays')
    np.savez_compressed(args.output/'predictions.npz',predictions=numeric,target_names=np.array(protocol['targets']),
        ranks=np.array(RANKS),model_names=np.array(MODELS),genes=shared_genes,basis_rank2=shared_basis,control_mean=shared_control)
    pd.DataFrame(scores_rows).to_csv(args.output/'donor_scores.csv',index=False)
    pd.DataFrame(selected_rows).to_csv(args.output/'selected_donors.csv',index=False)
    stability = pd.DataFrame(stability_rows)
    stability.to_csv(args.output/'donor_omission_stability.csv',index=False)
    stability.groupby(['target','rank','rule'],sort=False).same_selection_as_full.agg(
        comparisons='count',unchanged='sum',unchanged_fraction='mean').reset_index().to_csv(
        args.output/'donor_omission_summary.csv',index=False)
    write_json(args.output/'ALLOWED_INPUTS.json',access)
    manifest = {'status':'PREDICTIONS_SAVED_BEFORE_SEPARATE_EVALUATION_KNOWN_DATA',
        'created_utc':datetime.now(timezone.utc).isoformat(),'protocol_sha256':sha256(args.protocol),
        'code_sha256':sha256(__file__),'input_hashes':protocol['input_hashes'],
        'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,
        'target_count':5,'ranks':RANKS,'models':MODELS,'genes':len(shared_genes),
        'donor_omission_rows':len(stability),'synthetic':synthetic_test(),
        'output_hashes':{p.name:sha256(p) for p in args.output.iterdir() if p.is_file()}}
    write_json(args.output/'FREEZE_MANIFEST.json',manifest)
    print(json.dumps({'status':manifest['status'],'targets':5,'models':12,'genes':len(shared_genes),'stability_rows':len(stability)}))

def evaluate(args):
    protocol,design = source_check(args)
    manifest = json.loads((args.frozen/'FREEZE_MANIFEST.json').read_text())
    require(manifest['protocol_sha256']==sha256(args.protocol),'Protocol changed after generation')
    require(manifest['code_sha256']==sha256(__file__),'Code changed after generation')
    for name,expected in manifest['output_hashes'].items():
        require(sha256(args.frozen/name)==expected,'Frozen file changed: '+name)
    require(not args.output.exists() or not any(args.output.iterdir()),'Use a fresh evaluation directory')
    args.output.mkdir(parents=True,exist_ok=True)
    saved = np.load(args.frozen/'predictions.npz',allow_pickle=False)
    require(saved['target_names'].tolist()==protocol['targets'] and saved['model_names'].tolist()==MODELS,'Unexpected frozen order')
    genes = saved['genes'].astype(str)
    columns = design.loc[design.condition.isin(protocol['targets']),'matrix_column'].tolist()
    table = pd.read_csv(args.data,sep='\t',usecols=['genes']+columns).set_index('genes')
    require(set(genes).issubset(table.index),'Missing target genes')
    table = table.loc[genes]
    metrics = []
    for ti,target in enumerate(protocol['targets']):
        pool_names = design.loc[design.condition.eq(target),'matrix_column'].tolist()
        observed = np.log2(1+table[pool_names].to_numpy(float).T)
        require(np.isfinite(observed).all(),'Invalid target expression')
        for ri,rank in enumerate(RANKS):
            q = saved['basis_rank2'][:rank]
            predictions = saved['predictions'][ti,ri]
            for pool,actual in zip(pool_names,observed):
                errors = predictions-actual
                common = (errors @ q.T) @ q
                residual = errors-common
                mse = np.mean(errors**2,axis=1)
                common_mse,residual_mse = np.mean(common**2,axis=1),np.mean(residual**2,axis=1)
                require(np.allclose(mse,common_mse+residual_mse,atol=1e-10),'Error decomposition failed')
                for mi,model in enumerate(MODELS):
                    metrics.append({'target':target,'rank':rank,'model':model,'library':pool,
                        'MSE':float(mse[mi]),'MAE':float(np.mean(np.abs(errors[mi]))),
                        'common_MSE':float(common_mse[mi]),'residual_MSE':float(residual_mse[mi])})
    frame = pd.DataFrame(metrics)
    frame.to_csv(args.output/'library_errors.csv',index=False)
    summary = frame.groupby(['target','rank','model'],sort=False)[['MSE','MAE','common_MSE','residual_MSE']].mean().reset_index()
    summary.to_csv(args.output/'target_model_errors.csv',index=False)
    paired = []
    for (target,rank,pool),subset in frame.groupby(['target','rank','library'],sort=False):
        lookup = subset.set_index('model')
        for comparator in ['raw_low2','all_donors_mean','control_mean','all_donors_half_shrink']:
            paired.append({'target':target,'rank':rank,'library':pool,'proposed':'common_low2','comparator':comparator,
                'MSE_difference_proposed_minus_comparator':float(lookup.loc['common_low2','MSE']-lookup.loc[comparator,'MSE']),
                'MAE_difference_proposed_minus_comparator':float(lookup.loc['common_low2','MAE']-lookup.loc[comparator,'MAE'])})
    pd.DataFrame(paired).to_csv(args.output/'paired_primary_comparisons.csv',index=False)
    run = {'status':'EXPLORATORY_TARGET_EXCLUDED_DONOR_EVALUATION_COMPLETE',
        'created_utc':datetime.now(timezone.utc).isoformat(),'library_rows':len(frame),'target_model_rows':len(summary),
        'target_count':5,'protocol_sha256':sha256(args.protocol),'code_sha256':sha256(__file__),
        'freeze_manifest_sha256':sha256(args.frozen/'FREEZE_MANIFEST.json'),'input_hashes':protocol['input_hashes'],
        'output_hashes':{p.name:sha256(p) for p in args.output.iterdir() if p.is_file()},
        'limits':'All values known previously. Code-level target exclusion, not a new holdout. No tuning or pooled inference; information and donor-basis relation differ across task order.'}
    write_json(args.output/'RUN_SUMMARY.json',run)
    print(json.dumps(run,indent=2))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['freeze','evaluate','self-test'])
    for name in ['data','design','protocol','output','frozen']:
        parser.add_argument('--'+name,type=Path)
    args = parser.parse_args()
    if args.command=='self-test':
        print(json.dumps(synthetic_test()))
    else:
        require(all(getattr(args,n) is not None for n in ['data','design','protocol','output']),'Missing input arguments')
        if args.command=='freeze':
            freeze(args)
        else:
            require(args.frozen is not None,'Missing --frozen directory')
            evaluate(args)

if __name__=='__main__':
    main()
