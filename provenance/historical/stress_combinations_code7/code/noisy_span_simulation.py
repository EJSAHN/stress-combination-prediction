#!/usr/bin/env python3
"""Conditional, post-result sensitivity of a fitted stress-expression span.

This is not a calibrated test of a composite biological null. Gene panels are
fixed from CODE6. No synthetic expression matrices are saved or interpreted as
raw counts. Common Gaussian condition-mean errors are used across SD scales
and the two declared latent-mean constructions.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import time
import numpy as np
import pandas as pd

GENOTYPES = ('WT', 'RAP2.3')
CONDITIONS = ('CT', 'S', 'EL', 'HS', 'S+EL', 'S+HS', 'EL+HS', 'S+EL+HS')
NULLS = ('single_mean_span', 'joint_rank3')
LAMBDAS = (0.0, 0.5, 1.0)
SD_SCALES = (0.5, 1.0, 2.0)
N_DRAWS = 1000
BASE_SEED = 20261028


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for b in iter(lambda: stream.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def full_span(rows, reference_features=None):
    _, singular, q = np.linalg.svd(rows, full_matrices=False)
    tolerance = max(rows.shape[0], rows.shape[1] if reference_features is None else reference_features)*np.finfo(float).eps*singular[0] if singular[0] else 0.0
    rank = int(np.sum(singular > tolerance))
    return q[:rank], singular, rank, tolerance


def outside_gain(donor, target, q, gene_count):
    """Batch-safe outside gain, with row-orthonormal q and original G scale."""
    dr = donor-(donor@q.T)@q
    yr = target-(target@q.T)@q
    return (2*np.sum(dr*yr, axis=-1)-np.sum(dr*dr, axis=-1))/gene_count


def refitted_gain(means, gene_count):
    """Compute full-rank three-single span gains using small Gram matrices.

For poorly conditioned Gram matrices, use direct SVD on the actual rows,
which avoids squaring the condition number. No gene-selected truncation.
"""
    if means.ndim == 2:
        means = means[None]
    e = means[:,1:4]-means[:,0,None,:]
    d = means[:,4:7].mean(axis=1)-means[:,0]
    t = means[:,7]-means[:,0]
    gram = np.einsum('big,bjg->bij', e, e)/gene_count
    gd = np.einsum('big,bg->bi',e,d)/gene_count
    gt = np.einsum('big,bg->bi',e,t)/gene_count
    eigen = np.linalg.eigvalsh(gram)
    regular = eigen[:,0] > eigen[:,-1]*1e-10
    output = np.empty(len(means))
    ranks = np.empty(len(means),dtype=int)
    if regular.any():
        coeff = np.linalg.solve(gram[regular],gd[regular,:,None])[...,0]
        within_gain = 2*np.sum(coeff*gt[regular],axis=1)-np.sum(coeff*gd[regular],axis=1)
        raw_gain = (2*np.sum(d[regular]*t[regular],axis=1)-np.sum(d[regular]**2,axis=1))/gene_count
        output[regular] = raw_gain-within_gain
        ranks[regular] = 3
    for k in np.flatnonzero(~regular):
        q,_,r,_=full_span(e[k],gene_count)
        output[k]=outside_gain(d[k],t[k],q,gene_count)
        ranks[k]=r
    return output,ranks,int((~regular).sum())


def latent_means(observed_means, null_name):
    effects=observed_means-observed_means[0]
    if null_name=='single_mean_span':
        q,s,r,tol=full_span(effects[1:4])
        assert r==3
    else:
        _,s,vh=np.linalg.svd(effects[1:],full_matrices=False)
        q=vh[:3]
        tol=max(effects[1:].shape)*np.finfo(float).eps*s[0]
        assert (s>tol).sum()>=3
    latent=(effects@q.T)@q
    qsingle,ss,r,tol=full_span(latent[1:4])
    assert r==3, 'Joint null does not have three spanning latent single effects'
    assert np.max(np.abs(latent-(latent@qsingle.T)@qsingle))<1e-9
    assert np.max(np.abs(qsingle-(qsingle@q.T)@q))<1e-9
    z,ranks,_=refitted_gain(latent,len(observed_means[0]))
    oracle=outside_gain(latent[4:7].mean(axis=0),latent[7],q,len(observed_means[0]))
    assert abs(float(z[0]))<1e-9 and abs(float(oracle))<1e-9
    return latent,q,{'latent_single_rank':r,'latent_single_singular_values':ss.tolist(),
                    'latent_single_condition_number':float(ss[0]/ss[-1]),
                    'latent_basis_singular_values':s.tolist(),
                    'zero_noise_refitted_gain':float(z[0]),'zero_noise_fixed_gain':float(oracle)}


def covariance_trace_outside(residuals, variances, q, lam, gene_count):
    total=variances.sum(axis=1)
    empirical_inside=np.sum((residuals@q.T)**2,axis=(1,2))/2
    diagonal_inside=variances@np.sum(q*q,axis=0)
    return (total-((1-lam)*empirical_inside+lam*diagonal_inside))/gene_count


def analytic_oracles(residuals, variances, q, lam, sd_scale, gene_count):
    v=covariance_trace_outside(residuals,variances,q,lam,gene_count)/3
    donor=v[4:7].sum()/9
    return {'fixed_basis_estimated_control':float(sd_scale**2*(v[0]-donor)),
            'fixed_basis_true_control':float(-sd_scale**2*donor),
            'shared_control_increment':float(sd_scale**2*v[0])}


def read_group(code6, genotype, frozen_manifest, design):
    item=next(s for s in frozen_manifest['scenarios'] if s['genotype']==genotype and s['threshold_fpkm']==1)
    source=code6/'input'/Path(frozen_manifest['source_hashes'][genotype]['path']).name
    assert sha256(source)==frozen_manifest['source_hashes'][genotype]['sha256']
    asset=code6/'frozen_verified'/item['assets']
    assert sha256(asset)==item['assets_sha256']
    with np.load(asset,allow_pickle=False) as a:
        wanted=a['genes'].copy()
        frozen_q=a['basis'].copy()
    dd=design[design.genotype_metadata==genotype]
    names=dd.library_name.tolist()
    frame=pd.read_csv(source,sep='\t',usecols=['gene',*names],dtype={'gene':str},keep_default_na=False)
    dup=frame.gene.duplicated(keep=False)
    assert frame.loc[dup,'gene'].value_counts().to_dict()=={'AT1G79790':2,'AT2G04725':2}
    frame=frame.loc[~dup].set_index('gene')
    assert frame.index.is_unique and set(wanted).issubset(frame.index)
    frame=frame.loc[wanted]
    raw=frame[names].to_numpy(dtype=float)
    assert np.isfinite(raw).all() and (raw>=0).all()
    x=np.stack([np.log2(1+frame[dd.loc[dd.condition==c,'library_name'].tolist()].to_numpy(dtype=float).T)
                for c in CONDITIONS])
    assert x.shape==(8,3,len(wanted))
    means=x.mean(axis=1)
    residuals=x-means[:,None,:]
    variances=residuals.var(axis=1,ddof=1)
    observed=float(outside_gain(means[4:7].mean(axis=0)-means[0],means[7]-means[0],frozen_q,len(wanted)))
    recomputed,ranks,_=refitted_gain(means,len(wanted))
    assert ranks[0]==3 and np.isclose(recomputed[0],observed,atol=1e-10)
    return x,means,residuals,variances,observed,{'source':str(source),'source_sha256':sha256(source),
        'frozen_asset_sha256':sha256(asset),'gene_count':len(wanted),
        'gene_panel_sha256':hashlib.sha256(('\n'.join(wanted)+'\n').encode()).hexdigest(),
        'observed_mean_outside_gain':observed}


def conditional_summary(values, observed):
    q=np.quantile(values,[.025,.5,.95,.975,.99])
    count=int(np.count_nonzero(values>=observed))
    fraction=count/len(values)
    return {'draws':len(values),'mean':float(values.mean()),'sd':float(values.std(ddof=1)),
            'q025':float(q[0]),'q50':float(q[1]),'q95':float(q[2]),'q975':float(q[3]),'q99':float(q[4]),
            'observed_statistic':float(observed),'n_at_least_observed':count,
            'conditional_exceedance_fraction':float(fraction),
            'monte_carlo_standard_error':float(np.sqrt(fraction*(1-fraction)/len(values))),
            'calibrated_biological_p_value':False}


def run(args):
    started=time.monotonic()
    assert args.protocol.is_file()
    protocol=json.loads(args.protocol.read_text())
    if args.output.exists() and any(args.output.iterdir()):
        raise RuntimeError('Use a fresh simulation output directory')
    args.output.mkdir(parents=True,exist_ok=True)
    code6=args.code6_root
    fm_path=code6/'frozen_verified/FREEZE_MANIFEST.json'
    fm=json.loads(fm_path.read_text())
    summary6=json.loads((code6/'results/RUN_SUMMARY.json').read_text())
    design=pd.read_csv(code6/'fallback_metadata/sample_inventory.csv')
    configuration={'genotypes':list(GENOTYPES),'nulls':list(NULLS),'covariance_diagonal_fractions':list(LAMBDAS),
                   'noise_sd_scales':list(SD_SCALES),'draws_per_cell':N_DRAWS,'base_seed':BASE_SEED,
                   'batch_size':args.batch_size,'fixed_filter_fpkm':1.0}
    config=protocol['simulation']
    assert config['conditions']==list(CONDITIONS)
    assert config['latent_nulls']==list(NULLS)
    assert config['covariance_mixture_lambda']==list(LAMBDAS)
    assert config['noise_sd_scales']==list(SD_SCALES)
    assert config['draws_per_scenario']==N_DRAWS and config['base_seed']==BASE_SEED
    execution={'status':'CONDITIONAL_SIMULATION_STARTED','started_utc':now(),'protocol_sha256':sha256(args.protocol),
               'code_sha256':sha256(__file__),'code6_freeze_manifest_sha256':sha256(fm_path),
               'code6_summary_sha256':sha256(code6/'results/RUN_SUMMARY.json'),'configuration':configuration,
               'stage':'Post-result conditional sensitivity; all CODE6 target outcomes were already seen',
               'covariance_model':'(1-lambda)*sample covariance + lambda*diagonal sample variance; n=3 per cell',
               'fixed_panel':True,'formal_null_test':False,'synthetic_count_data':False,'source_groups':{},
               'rng_draw_order':'For each genotype/lambda seed, batches in draw order: standard_normal(batch,8,3), then for lambda>0 standard_normal(batch,8,G). Same errors reused across both latent nulls and all SD scales.',
               'seed_rule':'20261028 + 10000*genotype_index + 100*lambda_index; genotype and lambda order are in configuration',
               'compact_lambda0':'Exact observed-library row-space coordinates for empirical-only covariance, preserving dot products and original G normalization'}
    write_json(args.output/'EXECUTION_STARTED.json',execution)
    all_draws=[]; all_summaries=[]; latent_records=[]; noise_records=[]; negative_records=[]; fallback_count=0
    for gi,g in enumerate(GENOTYPES):
        x,means,residuals,variances,observed,meta=read_group(code6,g,fm,design)
        expected=next(o['outside_mean_mse_gain'] for o in summary6['outcomes'] if o['genotype']==g and o['threshold_fpkm']==1)
        assert np.isclose(observed,expected,atol=1e-10)
        n_gene=x.shape[-1]
        execution['source_groups'][g]=meta
        latent={}
        for null_name in NULLS:
            mu,q,detail=latent_means(means,null_name)
            latent[null_name]=(mu,q)
            original_effects=means-means[0]
            for si,c in enumerate(CONDITIONS):
                mse=float(np.mean((original_effects[si]-mu[si])**2))
                mean_noise=float(variances[si].mean()/3)
                latent_records.append({'genotype':g,'latent_null':null_name,'condition':c,
                    'mean_displacement_fit_mse':mse,'sample_mean_noise_trace_per_gene':mean_noise,
                    'fit_mse_over_mean_noise_trace':mse/mean_noise if mean_noise>0 else None,
                    'latent_log_mean_min':float((means[0]+mu[si]).min()),
                    'latent_log_mean_negative_fraction':float(np.mean(means[0]+mu[si]<0)),**detail})
            for scale in SD_SCALES:
                absolute_mu=means[0]+mu
                sd=scale*np.sqrt(variances/3)
                prob=np.zeros_like(absolute_mu)
                valid=sd>0
                ratio=absolute_mu[valid]/(np.sqrt(2)*sd[valid])
                prob[valid]=np.fromiter((0.5*math.erfc(float(z)) for z in ratio),dtype=float,count=len(ratio))
                prob[~valid]=(absolute_mu[~valid]<0).astype(float)
                for si,c in enumerate(CONDITIONS):
                    negative_records.append({'genotype':g,'latent_null':null_name,'noise_sd_scale':scale,'condition':c,
                        'expected_negative_condition_mean_log_fraction':float(prob[si].mean()),
                        'method':'Exact marginal Gaussian probability; same for all covariance mixtures because diagonal variances are preserved'})
            for lam in LAMBDAS:
                for scale in SD_SCALES:
                    analytic=analytic_oracles(residuals,variances,q,lam,scale,n_gene)
                    noise_records.append({'genotype':g,'latent_null':null_name,'lambda':lam,'noise_sd_scale':scale,**analytic})
        # Empirical-covariance draws lie exactly in the span of all observed rows.
        q_all,_,_,_=full_span(x.reshape(24,n_gene))
        residuals_small=residuals@q_all.T
        latent_small={name:(mu@q_all.T,q@q_all.T) for name,(mu,q) in latent.items()}
        assert max(float(np.max(np.abs(q@q.T-np.eye(3)))) for _,q in latent_small.values())<1e-9
        for li,lam in enumerate(LAMBDAS):
            seed=BASE_SEED+10000*gi+100*li
            rng=np.random.default_rng(seed)
            accum={(name,scale):[] for name in NULLS for scale in SD_SCALES}
            for first in range(0,N_DRAWS,args.batch_size):
                count=min(args.batch_size,N_DRAWS-first)
                z=rng.standard_normal((count,8,3))
                if lam==0:
                    eps=np.einsum('bsi,sig->bsg',z,residuals_small)/np.sqrt(6)
                    model=latent_small
                else:
                    eps=np.einsum('bsi,sig->bsg',z,residuals)*np.sqrt((1-lam)/6)
                    eps+=rng.standard_normal((count,8,n_gene))*np.sqrt(lam*variances/3)[None,:,:]
                    model=latent
                for null_name,(mu,q) in model.items():
                    for scale in SD_SCALES:
                        errors=scale*eps
                        simulated=mu[None,:,:]+errors
                        fitted,ranks,fallback=refitted_gain(simulated,n_gene)
                        fallback_count+=fallback
                        # All latent true effects lie in q; use noise terms directly for oracles.
                        ed=errors[:,4:7].mean(axis=1)
                        et=errors[:,7]
                        estimated=outside_gain(ed-errors[:,0],et-errors[:,0],q,n_gene)
                        truth=outside_gain(ed,et,q,n_gene)
                        block=np.column_stack([fitted,estimated,truth,fitted-estimated,estimated-truth,ranks])
                        assert np.isfinite(block).all()
                        accum[(null_name,scale)].append(block)
            for (null_name,scale),blocks in accum.items():
                values=np.concatenate(blocks)
                for j,row in enumerate(values):
                    all_draws.append({'genotype':g,'latent_null':null_name,'lambda':lam,'noise_sd_scale':scale,
                        'draw':j+1,'seed':seed,'refitted_basis_estimated_control':float(row[0]),
                        'fixed_basis_estimated_control':float(row[1]),'fixed_basis_true_control':float(row[2]),
                        'basis_refit_increment':float(row[3]),'shared_control_increment':float(row[4]),'refitted_rank':int(row[5])})
                for j,score in enumerate(['refitted_basis_estimated_control','fixed_basis_estimated_control','fixed_basis_true_control',
                                          'basis_refit_increment','shared_control_increment']):
                    row={'genotype':g,'primary_genotype':g=='WT','latent_null':null_name,'lambda':lam,
                         'noise_sd_scale':scale,'score':score,'seed':seed,
                         **conditional_summary(values[:,j],observed)}
                    # Only score 1 is directly commensurable with the observed CODE6 statistic.
                    row['observed_comparison_same_statistic']=(j==0)
                    if j!=0:
                        for k in ['observed_statistic','n_at_least_observed','conditional_exceedance_fraction','monte_carlo_standard_error']:
                            row[k]=None
                    all_summaries.append(row)
                print(f'{g} {null_name} lambda={lam:g} SD={scale:g}: simulated {N_DRAWS} means',flush=True)
    draws=pd.DataFrame(all_draws); summaries=pd.DataFrame(all_summaries)
    draws.to_csv(args.output/'simulation_draws.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    summaries.to_csv(args.output/'simulation_summary.csv',index=False)
    pd.DataFrame(latent_records).to_csv(args.output/'latent_mean_diagnostics.csv',index=False)
    analytic=pd.DataFrame(noise_records)
    analytic.to_csv(args.output/'analytic_oracle_expectations.csv',index=False)
    pd.DataFrame(negative_records).to_csv(args.output/'negative_working_scale_probabilities.csv',index=False)
    oracle_checks=[]
    for _,row in analytic.iterrows():
        selected=draws[(draws.genotype==row.genotype)&(draws.latent_null==row.latent_null)&
                       (draws['lambda']==row['lambda'])&(draws.noise_sd_scale==row.noise_sd_scale)]
        for score in ['fixed_basis_estimated_control','fixed_basis_true_control','shared_control_increment']:
            scores=selected[score]
            se=float(scores.std(ddof=1)/np.sqrt(len(scores)))
            delta=float(scores.mean()-row[score])
            oracle_checks.append({'genotype':row.genotype,'latent_null':row.latent_null,'lambda':float(row['lambda']),
                'noise_sd_scale':float(row.noise_sd_scale),'score':score,'analytic_mean':float(row[score]),
                'simulated_mean':float(scores.mean()),'monte_carlo_se_of_mean':se,
                'difference_in_mc_se':delta/se if se>0 else None})
    pd.DataFrame(oracle_checks).to_csv(args.output/'oracle_expectation_checks.csv',index=False)
    nominal=summaries[(summaries.genotype=='WT')&(summaries.noise_sd_scale==1)&
                      (summaries.score=='refitted_basis_estimated_control')]
    execution.update({'status':'CONDITIONAL_SIMULATION_COMPLETE','completed_utc':now(),
        'runtime_seconds':time.monotonic()-started,'simulation_cells':36,'draw_rows':len(draws),
        'summary_rows':len(summaries),'direct_svd_fallback_draw_evaluations':fallback_count,
        'wt_observed_exceeds_all_six_nominal_conditional_q95':bool((nominal.observed_statistic>nominal.q95).all()),
        'interpretation':'Conditional T-only sensitivity. Does not calibrate disjoint-control margins or establish a biological P-value; covariance and latent orientation remain uncertain.',
        'working_scale':'Gaussian log-expression errors are not truncated and may produce negative working-scale values; these are not synthetic counts.',
        'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,
        'outputs':{p.name:sha256(p) for p in args.output.iterdir() if p.is_file()}})
    write_json(args.output/'RUN_SUMMARY.json',execution)
    print(json.dumps({k:execution[k] for k in ['status','runtime_seconds','simulation_cells','draw_rows','wt_observed_exceeds_all_six_nominal_conditional_q95']},indent=2))


def selftest():
    rng=np.random.default_rng(9241)
    g=43
    q,_,_,_=full_span(rng.normal(size=(3,g)))
    coefficients=rng.normal(size=(8,3));coefficients[0]=0
    latent=coefficients@q
    _,_,r,_=full_span(latent[1:4]);assert r==3
    values,ranks,_=refitted_gain(latent,g)
    assert abs(values[0])<1e-12
    means=rng.normal(size=(13,8,g))
    fast,ranks,_=refitted_gain(means,g)
    for i,x in enumerate(means):
        qq,_,rr,_=full_span(x[1:4]-x[0]);assert rr==ranks[i]
        direct=outside_gain(x[4:7].mean(axis=0)-x[0],x[7]-x[0],qq,g)
        assert np.isclose(fast[i],direct,atol=1e-12)
        fullpred=x[4:7].mean(axis=0)
        targetlibs=x[7]+rng.normal(size=(3,g))
        targetlibs-=targetlibs.mean(axis=0)-x[7]
        dr=fullpred-x[0];dr-=dr@qq.T@qq
        perlib=np.mean((targetlibs-x[0])**2,axis=1)-np.mean((targetlibs-(x[0]+dr))**2,axis=1)
        assert np.isclose(perlib.mean(),fast[i],atol=1e-12)
    # Analytic oracle means versus independent Gaussian covariance construction.
    residuals=rng.normal(size=(8,3,g));residuals-=residuals.mean(axis=1,keepdims=True)
    variances=residuals.var(axis=1,ddof=1)
    lam=.5; b=30000
    z=rng.normal(size=(b,8,3))
    eps=np.einsum('bsi,sig->bsg',z,residuals)*np.sqrt((1-lam)/6)
    eps+=rng.normal(size=(b,8,g))*np.sqrt(lam*variances/3)
    ed=eps[:,4:7].mean(axis=1);et=eps[:,7]
    oracle1=outside_gain(ed-eps[:,0],et-eps[:,0],q,g)
    oracle2=outside_gain(ed,et,q,g)
    analytic=analytic_oracles(residuals,variances,q,lam,1,g)
    for score,a in [('fixed_basis_estimated_control',oracle1),('fixed_basis_true_control',oracle2)]:
        assert abs(a.mean()-analytic[score])<5*a.std(ddof=1)/np.sqrt(b)
    # Exact compact-coordinate preservation for empirical covariance.
    observed=rng.normal(size=(24,g))
    qc,_,_,_=full_span(observed)
    noise=rng.normal(size=(8,24))@observed
    qq,_,_,_=full_span(observed[:3])
    a=outside_gain(noise[4:7].mean(axis=0)-noise[0],noise[7]-noise[0],qq,g)
    small=noise@qc.T
    b=outside_gain(small[4:7].mean(axis=0)-small[0],small[7]-small[0],qq@qc.T,g)
    assert np.isclose(a,b,atol=1e-10)
    print('SELFTEST_PASS: direct-SVD/Gram identity, mean-gain identity, zero-noise null, analytic oracle expectations, compact-coordinate identity')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    sub.add_parser('selftest')
    p=sub.add_parser('run')
    p.add_argument('--code6-root',type=Path,required=True)
    p.add_argument('--protocol',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--batch-size',type=int,default=16)
    args=parser.parse_args()
    if args.command=='selftest':selftest()
    else:
        assert 1<=args.batch_size<=64
        run(args)

if __name__=='__main__': main()
