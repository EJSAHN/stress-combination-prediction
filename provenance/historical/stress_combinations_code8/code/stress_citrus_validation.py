#!/usr/bin/env python3
"""Frozen citrus pair-to-triple prediction. No installation or network use.
Freeze/poison commands never convert numerical target values. Evaluate verifies
frozen assets and a real-data poison record before target parsing.
"""
from __future__ import annotations
import argparse, gzip, hashlib, json, os, platform, tempfile, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
GENOTYPES=('Car','Cleo')
SINGLES=('WS','HL','HS')
PAIRS=('WS+HL','WS+HS','HL+HS')
TARGET='WS+HL+HS'
CONDITIONS=('CT',*SINGLES,*PAIRS,TARGET)
THRESHOLDS=(1.0,0.5,2.0)
MODELS=('control','full_pair_centroid','within_V3','outside_V3','single_additive',
        'second_order_factorial','single_mean','pair_trained_scalar','pair_trained_factor_OLS','within_V9','within_V11')
COMPARATORS=('control','within_V3','single_additive','second_order_factorial','single_mean','pair_trained_scalar','pair_trained_factor_OLS')
COMPONENTS=('within_V3','added_V9','added_CT','outside_V11')
EPSILON=1e-12

def now():return datetime.now(timezone.utc).isoformat()
def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()
def write_json(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def fresh(path):
    if path.exists() and any(path.iterdir()):raise ValueError('Use a fresh directory: '+str(path))
    path.mkdir(parents=True,exist_ok=True)
def key_for(g,t):return g+'_fpkm'+str(t).replace('.','p')

def checked_protocol(path,amendment):
    p=json.loads(path.read_text());a=json.loads(amendment.read_text())
    assert p['protocol_id']=='CODE8_CITRUS_FIXED_PAIR_TO_TRIPLE_V1' and p['source']=='GSE203331'
    assert p['genotypes']==list(GENOTYPES) and p['singles']==list(SINGLES) and p['pairs']==list(PAIRS)
    assert p['target']==TARGET and p['thresholds']==list(THRESHOLDS)
    assert p['models']==list(MODELS) and p['ordinary_primary_comparators']==list(COMPARATORS)
    assert p['primary_rule']['epsilon']==EPSILON
    assert a['amendment_id']=='CODE8_INPUT_AMENDMENT01' and a['original_protocol_sha256']==sha256(path)
    assert a['expression_values_parsed'] is False and a['replacements']['gene_column']=='Feature_GID'
    return p,a

def checked_design(path):
    d=pd.read_csv(path,dtype=str,keep_default_na=False)
    required=['expression_column','genotype_code','condition','replicate_label']
    assert all(k in d for k in required) and len(d)==48 and d.expression_column.is_unique
    assert set(d.genotype_code)==set(GENOTYPES) and set(d.condition)==set(CONDITIONS)
    assert d.expression_column.str.endswith('_FPKM').all()
    assert d.groupby(['genotype_code','condition']).size().eq(3).all()
    for _,r in d.groupby(['genotype_code','condition']):assert set(r.replicate_label)=={'1','2','3'}
    assert not (d[required]=='').any().any()
    return d.sort_values(['genotype_code','condition','replicate_label'],kind='stable').reset_index(drop=True)

def header_info(path,design,gene_column):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt',encoding='utf-8-sig') as f:header=f.readline().rstrip('\r\n').split('\t')
    assert len(header)==len(set(header)) and gene_column in header
    assert gene_column not in set(design.expression_column) and set(design.expression_column)<=set(header)
    assert set(c for c in header if c.endswith('_FPKM'))==set(design.expression_column)
    return {'gene_column':gene_column,'header':header,
            'ignored_columns':[h for h in header if h not in set(design.expression_column) and h!=gene_column]}

def identifier_inventory(path,gene_column):
    ids=pd.read_csv(path,sep='\t',usecols=[gene_column],dtype=str,keep_default_na=False)[gene_column]
    assert ids.ne('').all() and ids.notna().all(),'Empty gene identifier'
    dup=ids.duplicated(keep=False)
    return {'source_rows':len(ids),'eligible_unique_rows':int((~dup).sum()),'excluded_duplicate_rows':int(dup.sum()),
            'duplicate_identifier_counts':ids[dup].value_counts().sort_index().to_dict(),
            'identifier_order_sha256':hashlib.sha256(('\n'.join(ids)+'\n').encode()).hexdigest()}

def read_allowed(path,gene_column,columns):
    """Only explicitly selected expression fields are numerically converted."""
    df=pd.read_csv(path,sep='\t',usecols=[gene_column,*columns],dtype={gene_column:str},keep_default_na=False)
    assert df[gene_column].ne('').all() and df[gene_column].notna().all()
    df=df.loc[~df[gene_column].duplicated(keep=False)]
    values=df[columns].to_numpy(dtype=np.float64).T
    assert np.isfinite(values).all(),'Nonfinite permitted value'
    assert (values>=0).all(),'Negative permitted value'
    return df[gene_column].to_numpy(str),values

def span(rows):
    _,s,vh=np.linalg.svd(rows,full_matrices=False)
    tol=max(rows.shape)*np.finfo(np.float64).eps*s[0] if s[0] else 0.0
    rank=int(np.sum(s>tol));q=vh[:rank]
    assert np.allclose(q@q.T,np.eye(rank),atol=2e-10)
    return q,s,rank,float(tol)
def project(x,q):return (x@q.T)@q

def make_assets(groups):
    c=groups['CT'].mean(axis=0)
    singles=np.stack([groups[k].mean(axis=0) for k in SINGLES])-c
    pairs=np.stack([groups[k].mean(axis=0) for k in PAIRS])-c
    donor=pairs.mean(axis=0)
    rows9=np.concatenate([groups[k]-c for k in SINGLES])
    rows11=np.concatenate([rows9,groups['CT']-c])
    qs={};spectra={};ranks=[]
    for name,rows in [('V3',singles),('V9',rows9),('V11',rows11)]:
        q,s,rank,tol=span(rows);qs[name]=q;spectra[name]=s
        ranks.append({'span':name,'rank':rank,'tolerance':tol,'singular_values':s.tolist()})
    assert len(qs['V3'])<=3 and len(qs['V9'])<=9 and len(qs['V11'])<=11
    assert np.allclose(qs['V3'],project(qs['V3'],qs['V9']),atol=2e-9)
    assert np.allclose(qs['V9'],project(qs['V9'],qs['V11']),atol=2e-9)
    inside=project(donor,qs['V3'])
    X=np.stack([np.stack([singles[i] if s in pair.split('+') else np.zeros_like(c)
                         for i,s in enumerate(SINGLES)],axis=1) for pair in PAIRS]).reshape(-1,3)
    z=pairs.reshape(-1);additive_pairs=X.sum(axis=1)
    numerator=float(additive_pairs@z);denominator=float(additive_pairs@additive_pairs)
    alpha=numerator/denominator if denominator else 0.0
    beta,_,rank,sv=np.linalg.lstsq(X,z,rcond=None)
    pred=np.stack([c,c+donor,c+inside,c+donor-inside,c+singles.sum(axis=0),
                   c+pairs.sum(axis=0)-singles.sum(axis=0),c+singles.mean(axis=0),
                   c+alpha*singles.sum(axis=0),c+beta@singles,
                   c+project(donor,qs['V9']),c+project(donor,qs['V11'])])
    assert np.allclose(pred[1],pred[2]+pred[3]-c,atol=2e-10)
    second=pred[5]-c
    assert np.allclose(second-project(second,qs['V3']),3*(donor-inside),atol=2e-9)
    scalar_mse=float(np.mean((alpha*additive_pairs-z)**2));ols_mse=float(np.mean((X@beta-z)**2))
    assert ols_mse<=scalar_mse+1e-10
    fitting={'alpha':alpha,'alpha_numerator':numerator,'alpha_denominator':denominator,
             'alpha_zero_denominator':denominator==0,'beta':beta.tolist(),'beta_factors':list(SINGLES),
             'ols_design_rank':int(rank),'ols_design_singular_values':sv.tolist(),
             'ols_condition_number':float(sv[0]/sv[-1]) if sv[-1] else None,
             'pair_training_gene_contexts':int(len(z)),
             'identity_pair_training_mse':float(np.mean((additive_pairs-z)**2)),
             'scalar_pair_training_mse':scalar_mse,'factor_OLS_pair_training_mse':ols_mse,
             'coefficient_constraints':'unconstrained, no intercept, no target tuning'}
    arrays={'predictions':pred,'model_names':np.array(MODELS),'alpha':np.array(alpha),'beta':beta,
            'lower_group_names':np.array(('CT',*SINGLES,*PAIRS)),
            'lower_log_groups':np.stack([groups[k] for k in ('CT',*SINGLES,*PAIRS)])}
    for name,q in qs.items():arrays['basis_'+name]=q;arrays['singular_'+name]=spectra[name]
    return arrays,ranks,fitting

def freeze_one(source,design,genotype,threshold,gene_column):
    group=design[design.genotype_code==genotype];permitted=group[group.condition!=TARGET]
    columns=permitted.expression_column.tolist();targets=group.loc[group.condition==TARGET,'expression_column'].tolist()
    assert len(columns)==21 and len(targets)==3 and set(columns).isdisjoint(targets)
    genes,raw=read_allowed(source,gene_column,columns);cond=permitted.condition.to_numpy(str)
    lower=np.isin(cond,('CT',*SINGLES));assert lower.sum()==12
    keep=raw[lower].mean(axis=0)>=threshold;assert keep.sum()>=12
    log=np.log2(1+raw[:,keep]);groups={k:log[cond==k] for k in ('CT',*SINGLES,*PAIRS)}
    assert all(v.shape[0]==3 for v in groups.values())
    arrays,ranks,fitting=make_assets(groups)
    arrays.update({'genes':genes[keep],'target_columns':np.array(targets),'allowed_columns':np.array(columns)})
    item={'key':key_for(genotype,threshold),'genotype':genotype,'threshold_fpkm':threshold,
          'primary_threshold':threshold==1,'gene_column':gene_column,'genes':int(keep.sum()),
          'eligible_unique_rows':len(genes),'allowed_numeric_columns':columns,'excluded_target_columns':targets,
          'spans':ranks,'pair_training':fitting}
    return arrays,item

def save_verified_npz(path,arrays):
    tmp=path.with_suffix(path.suffix+'.partial')
    with tmp.open('wb') as f:np.savez_compressed(f,**arrays);f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)
    with np.load(path,allow_pickle=False) as a:
        assert set(a.files)==set(arrays)
        for k in a.files:assert np.array_equal(a[k],arrays[k]),k
    return sha256(path)

def freeze(args):
    start=time.perf_counter();fresh(args.output);checked_protocol(args.protocol,args.amendment)
    design=checked_design(args.design);hi=header_info(args.source,design,args.gene_column)
    inventory=identifier_inventory(args.source,hi['gene_column'])
    m={'status':'FREEZE_IN_PROGRESS','started_utc':now(),'source_sha256':sha256(args.source),
       'source_path':str(args.source.resolve()),'protocol_sha256':sha256(args.protocol),
       'amendment_sha256':sha256(args.amendment),'design_sha256':sha256(args.design),
       'code_sha256':sha256(__file__),'header_info':hi,'identifier_inventory':inventory,
       'target_numerical_values_parsed':False,'scenarios':[]}
    for g in GENOTYPES:
        for t in THRESHOLDS:
            arrays,item=freeze_one(args.source,design,g,t,hi['gene_column']);name=item['key']+'.npz'
            item.update({'assets':name,'assets_sha256':save_verified_npz(args.output/name,arrays),'atomic_roundtrip_verified':True})
            m['scenarios'].append(item)
            print(f'FROZEN {item["key"]}: {item["genes"]} genes; ranks '+','.join(str(v['rank']) for v in item['spans']),flush=True)
    for item in m['scenarios']:
        p=args.output/item['assets'];assert sha256(p)==item['assets_sha256']
        with np.load(p,allow_pickle=False) as a:
            for k in a.files:a[k]
    m.update({'status':'ALL_SIX_SCENARIOS_FROZEN_BEFORE_TARGET_PARSE','completed_utc':now(),
              'runtime_seconds':time.perf_counter()-start,'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__})
    write_json(args.output/'FREEZE_MANIFEST.json',m)
    coefficients=[];spans=[]
    for item in m['scenarios']:
        base={k:item[k] for k in ('key','genotype','threshold_fpkm','genes')};fit=item['pair_training']
        coefficients.append({**base,**{k:v for k,v in fit.items() if not isinstance(v,list)},
                             **{'beta_'+n:fit['beta'][i] for i,n in enumerate(SINGLES)}})
        for s in item['spans']:spans.append({**base,**{k:v for k,v in s.items() if k!='singular_values'},'singular_values_json':json.dumps(s['singular_values'])})
    pd.DataFrame(coefficients).to_csv(args.output/'training_coefficients.csv',index=False)
    pd.DataFrame(spans).to_csv(args.output/'span_inventory.csv',index=False)
    print('FREEZE_COMPLETE; numerical target columns not parsed',flush=True)

def verify_frozen(args):
    checked_protocol(args.protocol,args.amendment);design=checked_design(args.design)
    m=json.loads((args.frozen/'FREEZE_MANIFEST.json').read_text())
    assert m['status']=='ALL_SIX_SCENARIOS_FROZEN_BEFORE_TARGET_PARSE' and m['target_numerical_values_parsed'] is False
    for key,path in [('source',args.source),('protocol',args.protocol),('amendment',args.amendment),('design',args.design),('code',Path(__file__))]:
        assert m[key+'_sha256']==sha256(path),key+' changed after freeze'
    assert m['identifier_inventory']==identifier_inventory(args.source,m['header_info']['gene_column'])
    assert len(m['scenarios'])==6 and {s['key'] for s in m['scenarios']}=={key_for(g,t) for g in GENOTYPES for t in THRESHOLDS}
    for item in m['scenarios']:
        p=args.frozen/item['assets'];assert sha256(p)==item['assets_sha256']
        with np.load(p,allow_pickle=False) as a:
            for k in a.files:a[k]
    return design,m

def poison_check(args):
    """Overwrite target string fields without converting/inspecting their values."""
    fresh(args.output);design,m=verify_frozen(args)
    targets=set(design.loc[design.condition==TARGET,'expression_column']);opener=gzip.open if str(args.source).endswith('.gz') else open
    with tempfile.TemporaryDirectory(prefix='poison_',dir=args.output) as tmp:
        poisoned=Path(tmp)/'poisoned.tsv.gz'
        with opener(args.source,'rt',encoding='utf-8-sig') as f,gzip.open(poisoned,'wt') as out:
            line=f.readline();header=line.rstrip('\r\n').split('\t');indices=[i for i,k in enumerate(header) if k in targets]
            assert len(indices)==6;out.write(line);nrows=0
            for line in f:
                fields=line.rstrip('\r\n').split('\t');assert len(fields)==len(header)
                for i in indices:fields[i]='NaN'
                out.write('\t'.join(fields)+'\n');nrows+=1
        checks=[]
        for item in m['scenarios']:
            arrays,_=freeze_one(poisoned,design,item['genotype'],item['threshold_fpkm'],item['gene_column'])
            with np.load(args.frozen/item['assets'],allow_pickle=False) as a:
                assert set(a.files)==set(arrays)
                for k in a.files:assert np.array_equal(a[k],arrays[k]),(item['key'],k)
            checks.append({'scenario':item['key'],'all_frozen_arrays_bitwise_identical':True,'arrays_checked':len(arrays)})
    record={'status':'PASS_REAL_TARGET_POISON_ALL_SIX_SCENARIOS','completed_utc':now(),
            'source_sha256':m['source_sha256'],'protocol_sha256':m['protocol_sha256'],'amendment_sha256':m['amendment_sha256'],
            'code_sha256':sha256(__file__),'freeze_manifest_sha256':sha256(args.frozen/'FREEZE_MANIFEST.json'),
            'target_columns_replaced':sorted(targets),'replacement':'NaN','original_target_values_numerically_parsed':False,
            'source_gene_rows':nrows,'checks':checks}
    write_json(args.output/'POISON_CHECK.json',record);print(json.dumps(record,indent=2))

def align_target(ids,raw,wanted):
    idx=pd.Index(ids).get_indexer(wanted);assert (idx>=0).all(),'Frozen gene missing'
    return np.log2(1+raw[:,idx])
def direction_metrics(d,t):
    dn=float(np.linalg.norm(d));tn=float(np.linalg.norm(t))
    return (float(d@t/(dn*tn)) if dn and tn else None,dn/tn if tn else None)

def evaluate(args):
    start=time.perf_counter();fresh(args.output);design,m=verify_frozen(args)
    poison=json.loads(args.poison_report.read_text())
    assert poison['status']=='PASS_REAL_TARGET_POISON_ALL_SIX_SCENARIOS'
    assert poison['freeze_manifest_sha256']==sha256(args.frozen/'FREEZE_MANIFEST.json')
    assert poison['source_sha256']==sha256(args.source) and poison['code_sha256']==sha256(__file__)
    assert poison['amendment_sha256']==sha256(args.amendment)
    write_json(args.output/'TARGET_OPENING_RECORD.json',{
        'opening_started_utc':now(),'freeze_manifest_sha256':sha256(args.frozen/'FREEZE_MANIFEST.json'),
        'source_sha256':m['source_sha256'],'protocol_sha256':m['protocol_sha256'],'amendment_sha256':m['amendment_sha256'],
        'code_sha256':sha256(__file__),'poison_record_sha256':sha256(args.poison_report),
        'statement':'All six model sets and panels frozen and poisoned successfully; numerical targets parsed next.'})
    targets={}
    for g in GENOTYPES:
        item=next(s for s in m['scenarios'] if s['genotype']==g)
        targets[g]=read_allowed(args.source,item['gene_column'],item['excluded_target_columns'])
    errors=[];summaries=[];comparisons=[];components=[];oracles=[];checks=[]
    for item in m['scenarios']:
        info={'scenario':item['key'],'genotype':item['genotype'],'threshold_fpkm':item['threshold_fpkm'],
              'primary_threshold':item['primary_threshold'],'genes':item['genes']}
        with np.load(args.frozen/item['assets'],allow_pickle=False) as a:
            assert a['model_names'].tolist()==list(MODELS)
            y=align_target(*targets[item['genotype']],a['genes']);ym=y.mean(axis=0)
            pred=a['predictions'];c=pred[0];d=pred[1]-c
            q3,q9,q11=[a['basis_'+k] for k in ('V3','V9','V11')]
            d3,d9,d11=[project(d,q) for q in (q3,q9,q11)]
            pieces=[d3,d9-d3,d11-d9,d-d11];assert np.allclose(sum(pieces),d,atol=2e-10)
            gram=np.stack(pieces)@np.stack(pieces).T/len(c)
            offdiag=float(np.max(np.abs(gram-np.diag(np.diag(gram)))));assert offdiag<2e-9
            mse=np.mean((pred[:,None,:]-y[None,:,:])**2,axis=2)
            mae=np.mean(np.abs(pred[:,None,:]-y[None,:,:]),axis=2)
            noise=float(np.mean((y-ym)**2))
            oracle=c+project(ym-c,q3);oracle_mse=np.mean((y-oracle)**2,axis=1)
            oracle_mean=float(oracle_mse.mean());oracle_ym=float(np.mean((ym-oracle)**2))
            assert np.isclose(oracle_mean,oracle_ym+noise,atol=2e-10)
            for k,model in enumerate(MODELS):
                mean_mse=float(mse[k].mean());cosine,amplitude=direction_metrics(pred[k]-c,ym-c)
                ym_mse=float(np.mean((pred[k]-ym)**2));assert np.isclose(mean_mse,ym_mse+noise,atol=2e-10)
                summaries.append({**info,'model':model,'mean_library_mse':mean_mse,'mean_library_mae':float(mae[k].mean()),
                    'mse_to_target_mean':ym_mse,'target_library_variation_mse':noise,
                    'mean_mse_gain_over_control':float(mse[0].mean()-mean_mse),
                    'percent_mse_reduction_vs_control':100*(float(mse[0].mean())-mean_mse)/float(mse[0].mean()) if mse[0].mean() else None,
                    'direction_cosine_to_target_mean':cosine,'amplitude_ratio_to_target_mean':amplitude})
                for j,lib in enumerate(a['target_columns']):
                    errors.append({**info,'target_library':str(lib),'model':model,'mse':float(mse[k,j]),
                                   'mae':float(mae[k,j]),'mse_gain_over_control':float(mse[0,j]-mse[k,j])})
                if model!='full_pair_centroid':
                    diff=mse[k]-mse[1]
                    comparisons.append({**info,'comparator':model,'ordinary_primary_comparator':model in COMPARATORS,
                        'mean_comparator_mse_minus_centroid_mse':float(diff.mean()),
                        'centroid_better_mean':bool(diff.mean()>EPSILON),
                        'centroid_better_all_three_libraries':bool((diff>EPSILON).all()),
                        'minimum_library_advantage':float(diff.min()),
                        **{'library_'+str(j+1)+'_advantage':float(diff[j]) for j in range(3)}})
            for j,lib in enumerate(a['target_columns']):
                oracles.append({**info,'target_library':str(lib),'model':'TARGET_FITTED_ORACLE_NOT_PREDICTION',
                    'oracle_library_mse':float(oracle_mse[j]),'oracle_mean_library_mse':oracle_mean,'oracle_mse_to_target_mean':oracle_ym,
                    'centroid_library_mse':float(mse[1,j]),'oracle_minus_centroid_library_mse':float(oracle_mse[j]-mse[1,j]),
                    'oracle_minus_centroid_mean_mse':oracle_mean-float(mse[1].mean()),
                    'guaranteed_lower_bound_scope':'mean-library MSE in observed affine V3; not per-library optimum'})
            indices=[MODELS.index(k) for k in ('control','within_V3','single_additive','single_mean','pair_trained_scalar','pair_trained_factor_OLS')]
            bound=float(np.min(mse[indices].mean(axis=1)-oracle_mean));assert bound>=-2e-10
            max_error=0.0
            for j,lib in enumerate(a['target_columns']):
                target=y[j]-c;gain=float(mse[0,j]-mse[1,j]);parts=[]
                for name,piece in zip(COMPONENTS,pieces):
                    inner=float(np.mean(target*piece));power=float(np.mean(piece**2));value=2*inner-power;parts.append(value)
                    components.append({**info,'target_library':str(lib),'component':name,'mse_gain_contribution':value,
                        'donor_squared_effect':power,'donor_target_inner_product':inner,'total_centroid_gain_vs_control':gain})
                max_error=max(max_error,abs(sum(parts)-gain));assert abs(sum(parts)-gain)<2e-9
                for k,offset in [('within_V3',1),('within_V9',2),('within_V11',3)]:
                    assert np.isclose(sum(parts[offset:]),mse[MODELS.index(k),j]-mse[1,j],atol=2e-9)
            checks.append({**info,'max_component_offdiagonal_inner_product':offdiag,'max_gain_partition_error':max_error,
                'minimum_in_V3_prediction_mse_minus_oracle_mse':bound,'all_algebra_checks_pass':True})
    frames={'library_errors':pd.DataFrame(errors),'model_summary':pd.DataFrame(summaries),
            'comparator_summary':pd.DataFrame(comparisons),'library_gain_components':pd.DataFrame(components),
            'target_fitted_oracle_diagnostic':pd.DataFrame(oracles),'algebra_checks':pd.DataFrame(checks)}
    ckeys=['scenario','genotype','threshold_fpkm','primary_threshold','genes','component']
    frames['component_summary']=frames['library_gain_components'].groupby(ckeys,sort=False)[
        ['mse_gain_contribution','donor_squared_effect','donor_target_inner_product','total_centroid_gain_vs_control']].mean().reset_index()
    for name,df in frames.items():df.to_csv(args.output/(name+'.csv'),index=False)
    outcomes=[]
    for item in m['scenarios']:
        cf=frames['comparator_summary'];sel=cf[(cf.scenario==item['key'])&cf.ordinary_primary_comparator]
        assert len(sel)==len(COMPARATORS)
        ct=sel[sel.comparator=='control'].iloc[0];pm=bool(sel.centroid_better_mean.all());pc=bool(ct.centroid_better_all_three_libraries)
        outcomes.append({'scenario':item['key'],'genotype':item['genotype'],'threshold_fpkm':item['threshold_fpkm'],
            'primary_threshold':item['primary_threshold'],'genes':item['genes'],
            'centroid_beats_all_seven_comparators_on_mean_mse':pm,'centroid_beats_control_in_all_three_target_libraries':pc,
            'predeclared_per_genotype_pass':pm and pc,'failed_mean_comparators':sel.loc[~sel.centroid_better_mean,'comparator'].tolist()})
    count=sum(o['predeclared_per_genotype_pass'] for o in outcomes if o['primary_threshold'])
    decision={2:'DESCRIPTIVE_TWO_GENOTYPE_CROP_TRANSPORT_SUCCESS',1:'HETEROGENEOUS_ONE_GENOTYPE_SUPPORT',0:'PREDECLARED_CROP_TRANSPORT_CRITERION_NOT_MET'}[count]
    summary={'status':'CODE8_EVALUATION_COMPLETE','decision':decision,'completed_utc':now(),'source':'GSE203331',
        'primary_threshold_fpkm':1.0,'outcomes':outcomes,'runtime_seconds':time.perf_counter()-start,
        'source_sha256':m['source_sha256'],'protocol_sha256':m['protocol_sha256'],'amendment_sha256':m['amendment_sha256'],
        'freeze_manifest_sha256':sha256(args.frozen/'FREEZE_MANIFEST.json'),'code_sha256':sha256(__file__),
        'poison_report_sha256':sha256(args.poison_report),'formal_significance_test':False,
        'genotypes_are_independent_experiments':False,'target_fitted_oracle_is_prediction':False,
        'CODE7_noisy_span_limitation_resolved':False,'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,
        'output_hashes':{p.name:sha256(p) for p in sorted(args.output.iterdir()) if p.is_file()}}
    write_json(args.output/'RUN_SUMMARY.json',summary);print(json.dumps(summary,indent=2))

def selftest():
    rng=np.random.default_rng(283901)
    with tempfile.TemporaryDirectory() as td:
        root=Path(td);records=[];table={'Feature_GID':[f'g{i}' for i in range(80)]}
        for g in GENOTYPES:
            for cond in CONDITIONS:
                for rep in range(1,4):
                    col=g+'_'+cond+'_'+str(rep)+'_FPKM'
                    records.append({'expression_column':col,'genotype_code':g,'condition':cond,'replicate_label':str(rep)})
                    table[col]=rng.lognormal(1.2,0.5,80)
        design=pd.DataFrame(records);clean=root/'clean.tsv';pd.DataFrame(table).to_csv(clean,sep='\t',index=False)
        for g in GENOTYPES:
            original,_=freeze_one(clean,design,g,1.0,'Feature_GID')
            for poison in ('UNREADABLE_TARGET',np.nan,1e100):
                bad=pd.DataFrame(table)
                for col in design.loc[design.condition==TARGET,'expression_column']:bad[col]=poison
                path=root/'poison.tsv';bad.to_csv(path,sep='\t',index=False)
                poisoned,_=freeze_one(path,design,g,1.0,'Feature_GID')
                assert original.keys()==poisoned.keys()
                for k in original:assert np.array_equal(original[k],poisoned[k]),k
        duplicate=pd.DataFrame(table);duplicate.loc[0,'Feature_GID']='g1';path=root/'duplicate.tsv'
        duplicate.to_csv(path,sep='\t',index=False)
        ids,_=read_allowed(path,'Feature_GID',[records[0]['expression_column']]);assert len(ids)==78 and 'g1' not in ids
        assert identifier_inventory(path,'Feature_GID')['excluded_duplicate_rows']==2
    size=25;singles=np.zeros((3,size));singles[:,:3]=np.eye(3);beta=np.array([.2,1.4,-.5])
    groups={'CT':np.zeros((3,size))}
    for i,s in enumerate(SINGLES):groups[s]=np.tile(singles[i],(3,1))
    for pair in PAIRS:groups[pair]=np.tile(sum(beta[i]*singles[i] for i,s in enumerate(SINGLES) if s in pair.split('+')),(3,1))
    a,_,fit=make_assets(groups)
    assert np.allclose(a['beta'],beta,atol=1e-12) and np.allclose(a['predictions'][8],beta@singles,atol=1e-12)
    assert fit['ols_design_rank']==3 and fit['factor_OLS_pair_training_mse']<1e-24
    zero={k:np.zeros((3,size)) for k in ('CT',*SINGLES,*PAIRS)};zero['WS+HL'][:,3]=1
    a,_,fit=make_assets(zero)
    assert fit['alpha_zero_denominator'] and np.array_equal(a['beta'],np.zeros(3)) and a['basis_V3'].shape==(0,size)
    groups={k:rng.normal(size=(3,70)) for k in ('CT',*SINGLES,*PAIRS)};a,_,_=make_assets(groups)
    c=a['predictions'][0];d=a['predictions'][1]-c;q3,q9,q11=[a['basis_'+k] for k in ('V3','V9','V11')]
    d3,d9,d11=[project(d,q) for q in (q3,q9,q11)];pieces=[d3,d9-d3,d11-d9,d-d11];y=rng.normal(size=(3,70))
    gain=np.mean((y-c)**2,axis=1)-np.mean((y-c-d)**2,axis=1)
    comp=np.array([2*np.mean((y-c)*v,axis=1)-np.mean(v*v) for v in pieces]);assert np.allclose(gain,comp.sum(axis=0),atol=1e-12)
    oracle=c+project(y.mean(axis=0)-c,q3)
    for k in (0,2,4,6,7,8):assert np.mean((a['predictions'][k]-y)**2)>=np.mean((oracle-y)**2)-1e-12
    print('SELFTEST_PASS: target text/NaN/extreme invariance; duplicate exclusion; factor OLS recovery; rank-zero baselines; nested gain partition; affine mean-MSE oracle bound')

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True);sub.add_parser('selftest')
    for command in ('freeze','poison-check','evaluate'):
        q=sub.add_parser(command)
        for arg in ('source','design','protocol','amendment','output'):q.add_argument('--'+arg,type=Path,required=True)
        q.add_argument('--gene-column',default='Feature_GID')
        if command!='freeze':q.add_argument('--frozen',type=Path,required=True)
        if command=='evaluate':q.add_argument('--poison-report',type=Path,required=True)
    args=p.parse_args()
    if args.command=='selftest':selftest()
    elif args.command=='freeze':freeze(args)
    elif args.command=='poison-check':poison_check(args)
    else:evaluate(args)
if __name__=='__main__':main()
