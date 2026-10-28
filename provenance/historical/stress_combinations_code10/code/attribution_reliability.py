#!/usr/bin/env python3
"""Post-result repeatability and disjoint-pool refit diagnostics.
Original predictions stay immutable. Genes are not biological replicates.
"""
from __future__ import annotations
import argparse,hashlib,json,math,platform,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
GENOTYPES=('Car','Cleo')
SINGLES=('WS','HL','HS')
PAIRS=('WS+HL','WS+HS','HL+HS')
TARGET='WS+HL+HS'
CONDITIONS=('CT',*SINGLES,*PAIRS,TARGET)
COMPARATORS=('second_order_factorial','pair_trained_factor_OLS','within_V3')
MODEL_NAMES=('full_pair_centroid',*COMPARATORS)
PAIRINGS=(('Car1','Car2'),('Cleo1','Cleo2'),('Car1','Cleo2'),('Car2','Cleo1'),('Car1','Cleo1'),('Car2','Cleo2'))
METRICS=('Pearson','Spearman','top10pct_Jaccard')
def now():return datetime.now(timezone.utc).isoformat()
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def write_json(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def fresh(path):
 if path.exists() and any(path.iterdir()):raise ValueError('Use a fresh output directory: '+str(path))
 path.mkdir(parents=True,exist_ok=True)
def centered(x):return x-x.mean(axis=-1,keepdims=True)
def ordinary_correlation(a,b):
 ac=centered(np.asarray(a));bc=centered(np.asarray(b));va=float(np.mean(ac*ac));vb=float(np.mean(bc*bc))
 if va==0 or vb==0:return None,'zero_centered_variance_both' if va==vb==0 else ('zero_centered_variance_left' if va==0 else 'zero_centered_variance_right')
 return float(np.mean(ac*bc)/np.sqrt(va*vb)),''
def pair_association(a,b):
 pearson,pr=ordinary_correlation(a,b)
 ar=pd.Series(a).rank(method='average').to_numpy();br=pd.Series(b).rank(method='average').to_numpy()
 spearman,sr=ordinary_correlation(ar,br)
 return {'Pearson':pearson,'Pearson_unavailable_reason':pr,'Spearman':spearman,'Spearman_unavailable_reason':sr}
def gain_vectors(b,d,y):
 delta=(y-b)**2-(y-d)**2
 assert np.allclose(delta,(b-d)*(b+d-2*y),atol=2e-10,rtol=2e-10)
 mean=delta.mean(axis=0)
 assert np.allclose(mean,(b-d)*(b+d-2*y.mean(axis=0)),atol=2e-10,rtol=2e-10)
 assert np.allclose(delta-mean,-2*(b-d)*(y-y.mean(axis=0)),atol=2e-10,rtol=2e-10)
 return delta

def moment_reliability(delta):
 n,g=delta.shape;assert n>1 and g>0
 z=centered(delta);m=z.mean(axis=0);residual=z-m
 observed=float(np.mean(m*m));noise=float(np.sum(residual*residual)/(n*(n-1)*g));stable=observed-noise
 gram=z@z.T/g;offdiag=float((gram.sum()-np.trace(gram))/(n*(n-1)))
 assert np.isclose(stable,offdiag,atol=2e-10,rtol=2e-10)
 assert noise>=0 and stable<=observed+2e-10
 reliability=stable/observed if observed>0 else None
 return {'n_target_pools':n,'n_common_features':g,'observed_variance':observed,
  'pool_variance_of_mean':noise,'stable_variance':stable,'offdiagonal_stable_variance':offdiag,
  'identity_absolute_error':abs(stable-offdiag),'conditional_reliability':reliability,
  'reliability_unavailable_reason':'' if observed>0 else 'zero_observed_centered_variance',
  'stable_variance_nonpositive':stable<=0,'reliability_below_zero':reliability is not None and reliability<0,
  'reliability_above_one':reliability is not None and reliability>1},gram,m

def adjusted_correlation(ma,mb,va,vb):
 cov=float(np.mean(ma*mb))
 if va<=0 or vb<=0:
  reason='nonpositive_stable_variance_both' if va<=0 and vb<=0 else ('nonpositive_stable_variance_Car' if va<=0 else 'nonpositive_stable_variance_Cleo')
  return cov,None,reason,False
 result=cov/np.sqrt(va*vb)
 return cov,float(result),'',bool(abs(result)>1)

def fit_lower(groups):
 """The fitting interface excludes target observations."""
 assert set(groups)==set(('CT',*SINGLES,*PAIRS))
 c=groups['CT'].mean(axis=0)
 singles=np.stack([groups[k].mean(axis=0) for k in SINGLES])-c
 pairs=np.stack([groups[k].mean(axis=0) for k in PAIRS])-c;donor=pairs.mean(axis=0)
 _,sv,vh=np.linalg.svd(singles,full_matrices=False)
 tol=max(singles.shape)*np.finfo(float).eps*sv[0] if sv[0] else 0.
 rank=int((sv>tol).sum());q=vh[:rank];inside=(donor@q.T)@q
 X=np.stack([np.stack([singles[i] if s in pair.split('+') else np.zeros_like(c)
   for i,s in enumerate(SINGLES)],axis=1) for pair in PAIRS]).reshape(-1,3)
 beta,_,olsrank,olssv=np.linalg.lstsq(X,pairs.reshape(-1),rcond=None)
 predictions=np.stack([c+donor,c+pairs.sum(axis=0)-singles.sum(axis=0),c+beta@singles,c+inside])
 info={'single_span_rank':rank,'single_span_tolerance':float(tol),'single_span_rank_deficient':rank<3,
  'single_singular_1':float(sv[0]),'single_singular_2':float(sv[1]),'single_singular_3':float(sv[2]),
  'OLS_rank':int(olsrank),'OLS_rank_deficient':int(olsrank)<3,
  'OLS_condition_number':float(olssv[0]/olssv[-1]) if olssv[-1]>0 else None,
  'beta_WS':float(beta[0]),'beta_HL':float(beta[1]),'beta_HS':float(beta[2]),
  'fit_failed':False,'fit_failure_reason':''}
 return predictions,info

def top_indices(values,genes,k):
 assert 0<k<=len(genes) and len(values)==len(genes)
 return set(np.lexsort((genes,-values))[:k].tolist())
def side_indices(choice,side):
 choice=np.asarray(choice,int);assert choice.shape==(8,) and np.isin(choice,[0,1,2]).all()
 return [[int(i)] if side==1 else [j for j in range(3) if j!=i] for i in choice]

def check_inputs(root,protocol_path):
 protocol=json.loads(protocol_path.read_text())
 assert protocol['protocol_id']=='CODE10_CITRUS_ATTRIBUTION_REPEATABILITY_V1'
 assert protocol['comparators']==list(COMPARATORS) and protocol['condition_order']==list(CONDITIONS)
 assert protocol['part_B']['draws']==256 and protocol['part_B']['top_k']==1547
 for record in protocol['parent_pins'].values():assert sha(root.parent/record['relative_path'])==record['sha256'],record['relative_path']
 for section,path_key,hash_key in [('common_features','path','sha256'),('part_B','partition_path','partition_sha256')]:
  assert sha(root/protocol[section][path_key])==protocol[section][hash_key]
 common=np.array((root/protocol['common_features']['path']).read_text().splitlines())
 assert len(common)==protocol['common_features']['n']==15466 and common.tolist()==sorted(set(common))
 design=pd.read_csv(root.parent/protocol['parent_pins']['design']['relative_path'],dtype=str,keep_default_na=False)
 design=design.sort_values(['genotype_code','condition','replicate_label'],kind='stable')
 partition=pd.read_csv(root/protocol['part_B']['partition_path'],dtype=str,keep_default_na=False)
 assert len(partition)==2*256*8
 choices={g:{} for g in GENOTYPES}
 for g in GENOTYPES:
  rng=np.random.Generator(np.random.PCG64(protocol['part_B']['seeds'][g]));codes=rng.choice(3**8,size=256,replace=False)
  assert len(set(codes.tolist()))==256
  for draw,code in enumerate(codes,start=1):
   rows=partition[(partition.genotype==g)&(partition.draw_id==str(draw))]
   assert len(rows)==8 and set(rows.condition)==set(CONDITIONS)
   choice=[]
   for ci,cond in enumerate(CONDITIONS):
    row=rows[rows.condition==cond].iloc[0];expected=int((code//(3**ci))%3)
    assert int(row.partition_code)==int(code) and int(row.side1_index)==expected
    cols=design[(design.genotype_code==g)&(design.condition==cond)].expression_column.tolist();assert len(cols)==3
    assert row.side1_column==cols[expected]
    complement=[i for i in range(3) if i!=expected]
    assert row.side2_indices==';'.join(map(str,complement)) and row.side2_columns==';'.join(cols[i] for i in complement)
    choice.append(expected)
   choices[g][draw]=np.array(choice,dtype=int)
 return protocol,common,design,partition,choices

def load_data(root,protocol,common,design):
 path=root.parent/protocol['parent_pins']['source']['relative_path']
 columns=design.expression_column.tolist()
 df=pd.read_csv(path,sep='\t',usecols=['Feature_GID',*columns],dtype={'Feature_GID':str},keep_default_na=False)
 assert df.Feature_GID.is_unique and df.Feature_GID.ne('').all()
 raw=df[columns].to_numpy(float).T
 assert np.isfinite(raw).all() and (raw>=0).all()
 log=np.log2(1+raw);bycol={col:log[i] for i,col in enumerate(columns)}
 data={};frozen={};panel_ids={};common_indices={};checks=[]
 for g in GENOTYPES:
  with np.load(root.parent/protocol['parent_pins']['frozen_'+g+'_fpkm1p0']['relative_path'],allow_pickle=False) as a:
   genes=a['genes'];assert len(genes)==protocol['full_fit_panels'][g]
   indices=pd.Index(df.Feature_GID).get_indexer(genes);assert (indices>=0).all()
   data[g]=np.stack([np.stack([bycol[col][indices] for col in design[(design.genotype_code==g)&(design.condition==cond)].expression_column]) for cond in CONDITIONS])
   frozen[g]={name:a['predictions'][a['model_names'].tolist().index(name)].copy() for name in MODEL_NAMES}
   panel_ids[g]=genes.copy();common_indices[g]=pd.Index(genes).get_indexer(common);assert (common_indices[g]>=0).all()
   pred,fit=fit_lower({cond:data[g][i] for i,cond in enumerate(CONDITIONS[:-1])})
   for mi,name in enumerate(MODEL_NAMES):
    discrepancy=float(np.max(np.abs(pred[mi]-frozen[g][name])));assert discrepancy<2e-10
    checks.append({'genotype':g,'model':name,'n_full_panel_features':len(genes),
     'max_absolute_prediction_difference':discrepancy,'same_as_frozen_CODE8':True,
     'single_span_rank':fit['single_span_rank'],'OLS_rank':fit['OLS_rank']})
 assert sorted(set(panel_ids['Car'])&set(panel_ids['Cleo']))==common.tolist()
 return data,frozen,panel_ids,common_indices,checks

def part_A(data,frozen,common_indices,common,prior):
 moments=[];poolpairs=[];grams=[];cross=[];means={};rawmeans={};stats={}
 for g in GENOTYPES:
  ix=common_indices[g];y=data[g][-1][:,ix];d=frozen[g]['full_pair_centroid'][ix]
  for comparator in COMPARATORS:
   b=frozen[g][comparator][ix];delta=gain_vectors(b,d,y)
   stat,gram,m=moment_reliability(delta);means[(g,comparator)]=m;rawmeans[(g,comparator)]=delta.mean(axis=0);stats[(g,comparator)]=stat
   moments.append({'genotype':g,'comparator':comparator,**stat})
   for j in range(3):
    for k in range(3):grams.append({'genotype':g,'comparator':comparator,'pool1':j+1,'pool2':k+1,'centered_inner_product':float(gram[j,k])})
    for k in range(j+1,3):poolpairs.append({'genotype':g,'comparator':comparator,'pool1':j+1,'pool2':k+1,**pair_association(delta[j],delta[k])})
 for comparator in COMPARATORS:
  sa=stats[('Car',comparator)];sb=stats[('Cleo',comparator)]
  ma=means[('Car',comparator)];mb=means[('Cleo',comparator)]
  cov,adj,reason,flag=adjusted_correlation(ma,mb,sa['stable_variance'],sb['stable_variance'])
  raw=pair_association(rawmeans[('Car',comparator)],rawmeans[('Cleo',comparator)]);old=prior[prior.comparator==comparator].iloc[0]
  assert np.isclose(raw['Pearson'],old.pearson_gain,atol=2e-11)
  assert np.isclose(raw['Spearman'],old.spearman_gain,atol=2e-11)
  cross.append({'comparator':comparator,'n_common_features':len(common),'cross_genotype_covariance':cov,
   'raw_Pearson':raw['Pearson'],'raw_Spearman':raw['Spearman'],'adjusted_Pearson':adj,
   'raw_Pearson_unavailable_reason':raw['Pearson_unavailable_reason'],'raw_Spearman_unavailable_reason':raw['Spearman_unavailable_reason'],
   'adjusted_Pearson_unavailable_reason':reason,'adjusted_Pearson_outside_unit_interval':flag,
   'Car_observed_variance':sa['observed_variance'],'Cleo_observed_variance':sb['observed_variance'],
   'Car_stable_variance':sa['stable_variance'],'Cleo_stable_variance':sb['stable_variance'],
   'CODE9_Pearson_reproduced':True,'CODE9_Spearman_reproduced':True})
 return {'part_A_moments':pd.DataFrame(moments),'part_A_pool_pairs':pd.DataFrame(poolpairs),
         'part_A_gram':pd.DataFrame(grams),'part_A_cross_genotype':pd.DataFrame(cross)}

def summarize(frame,groupcols,metriccols):
 rows=[]
 for keys,group in frame.groupby(groupcols,sort=False,dropna=False):
  if not isinstance(keys,tuple):keys=(keys,)
  base=dict(zip(groupcols,keys))
  for metric in metriccols:
   values=pd.to_numeric(group[metric],errors='coerce').to_numpy(float);finite=values[np.isfinite(values)]
   if len(finite):
    q=np.quantile(finite,[0,.025,.5,.975,1],method='linear');stats=dict(zip(('min','q025','median','q975','max'),map(float,q)))
   else:stats={k:None for k in ('min','q025','median','q975','max')}
   rows.append({**base,'metric':metric,'n_total':len(values),'n_finite':len(finite),'n_unavailable':len(values)-len(finite),
    **stats,'summary_type':'overlapping_partition_sensitivity_NOT_CONFIDENCE_INTERVAL'})
 return pd.DataFrame(rows)

def part_B(data,panel_ids,common_indices,common,choices,protocol):
 pairrows=[];contrastrows=[];fitrows=[];scorerows=[]
 k=protocol['part_B']['top_k']
 for draw in range(1,257):
  prepared={}
  for genotype in GENOTYPES:
   for side in (1,2):
    label=genotype+str(side);indices=side_indices(choices[genotype][draw],side)
    groups={cond:data[genotype][ci,indices[ci],:] for ci,cond in enumerate(CONDITIONS)}
    low={cond:groups[cond] for cond in CONDITIONS[:-1]};ypools=groups[TARGET]
    base={'draw_id':draw,'genotype':genotype,'side':side,'side_label':label,
          'n_target_pools':len(ypools),'n_fit_features':len(panel_ids[genotype]),'n_common_features':len(common)}
    try:
     predictions,fit=fit_lower(low);fitrows.append({**base,**fit})
    except np.linalg.LinAlgError as exc:
     fitrows.append({**base,'fit_failed':True,'fit_failure_reason':str(exc)})
     for comparator in COMPARATORS:
      prepared[(label,comparator)]=None
      scorerows.append({**base,'comparator':comparator,'fit_failed':True,'fit_failure_reason':str(exc)})
     continue
    ix=common_indices[genotype];target=ypools[:,ix];d=predictions[0,ix]
    for mi,comparator in enumerate(COMPARATORS,start=1):
     b=predictions[mi,ix];delta=gain_vectors(b,d,target);gain=delta.mean(axis=0)
     bl=(target-b)**2;dl=(target-d)**2
     avg_gain=float(gain.mean());assert np.isclose(avg_gain,float(bl.mean()-dl.mean()),atol=2e-10)
     full_b=predictions[mi];full_d=predictions[0]
     scorerows.append({**base,'comparator':comparator,'fit_failed':False,'fit_failure_reason':'',
      'mean_gain_common':avg_gain,'baseline_mse_common':float(bl.mean()),'centroid_mse_common':float(dl.mean()),
      'baseline_mse_target_mean_common':float(np.mean((b-target.mean(axis=0))**2)),
      'centroid_mse_target_mean_common':float(np.mean((d-target.mean(axis=0))**2)),
      'target_pool_variation_mse_common':float(np.mean((target-target.mean(axis=0))**2)),
      'baseline_mse_full_panel':float(np.mean((ypools-full_b)**2)),
      'centroid_mse_full_panel':float(np.mean((ypools-full_d)**2)),
      'mean_gain_full_panel':float(np.mean((ypools-full_b)**2-(ypools-full_d)**2)),
      'fraction_positive_gain_common':float(np.mean(gain>0))})
     ranks=pd.Series(gain).rank(method='average').to_numpy()
     prepared[(label,comparator)]={'gain':gain,'rank':ranks,'top':top_indices(gain,common,k)}
  for comparator in COMPARATORS:
   bypair={}
   for pi,(left,right) in enumerate(PAIRINGS):
    pairing=left+'_'+right;la=prepared[(left,comparator)];ra=prepared[(right,comparator)]
    if la is None or ra is None:
     metric={'Pearson':None,'Spearman':None,'top10pct_Jaccard':None,
      'Pearson_unavailable_reason':'fit_failure','Spearman_unavailable_reason':'fit_failure','Jaccard_unavailable_reason':'fit_failure',
      'top_intersection':None,'top_union':None}
    else:
     p,pr=ordinary_correlation(la['gain'],ra['gain']);s,sr=ordinary_correlation(la['rank'],ra['rank'])
     intersection=len(la['top']&ra['top']);union=len(la['top']|ra['top'])
     assert len(la['top'])==len(ra['top'])==k and union>0
     metric={'Pearson':p,'Spearman':s,'top10pct_Jaccard':intersection/union,
      'Pearson_unavailable_reason':pr,'Spearman_unavailable_reason':sr,'Jaccard_unavailable_reason':'',
      'top_intersection':intersection,'top_union':union}
    row={'draw_id':draw,'comparator':comparator,'pairing':pairing,
     'pairing_role':'primary_size_matched' if pi<4 else 'secondary_same_size',
     'left_side':left,'right_side':right,'top_k_each':k,**metric}
    pairrows.append(row);bypair[pairing]=row
   for metric in METRICS:
    needed=[bypair[x][metric] for x in ('Car1_Car2','Cleo1_Cleo2','Car1_Cleo2','Car2_Cleo1')]
    if any(v is None for v in needed):within=cross=difference=None;reason='one_or_more_required_pair_metrics_unavailable'
    else:within=float(np.mean(needed[:2]));cross=float(np.mean(needed[2:]));difference=within-cross;reason=''
    contrastrows.append({'draw_id':draw,'comparator':comparator,'metric':metric,
     'within_mean':within,'cross_mean':cross,'difference':difference,'unavailable_reason':reason})
  if draw%32==0:print('PARTITION_DRAW_COMPLETED='+str(draw),flush=True)
 pairs=pd.DataFrame(pairrows);contrasts=pd.DataFrame(contrastrows);fits=pd.DataFrame(fitrows);scores=pd.DataFrame(scorerows)
 pairsummary=summarize(pairs,['comparator','pairing','pairing_role'],list(METRICS))
 # metric is the association metric, component identifies within/cross/difference.
 cs=summarize(contrasts.rename(columns={'metric':'association_metric'}),['comparator','association_metric'],['within_mean','cross_mean','difference'])
 cs=cs.rename(columns={'metric':'component','association_metric':'metric'})
 scoresummary=summarize(scores,['genotype','side','comparator'],['mean_gain_common','baseline_mse_common','centroid_mse_common','baseline_mse_target_mean_common','centroid_mse_target_mean_common','target_pool_variation_mse_common','mean_gain_full_panel'])
 fitsummary=summarize(fits,['genotype','side'],['single_span_rank','OLS_rank','OLS_condition_number','beta_WS','beta_HL','beta_HS'])
 return {'part_B_pair_metrics':pairs,'part_B_matched_contrasts':contrasts,'part_B_fit_details':fits,
  'part_B_side_scores':scores,'part_B_pair_summary':pairsummary,'part_B_contrast_summary':cs,
  'part_B_side_score_summary':scoresummary,'part_B_fit_summary':fitsummary}

def run(args):
 start=time.perf_counter();root=args.protocol.resolve().parent;fresh(args.output)
 protocol,common,design,partitions,choices=check_inputs(root,args.protocol)
 write_json(args.output/'RUN_STARTED.json',{'started_utc':now(),'code_sha256':sha(__file__),
  'protocol_sha256':sha(args.protocol),'partition_sha256':protocol['part_B']['partition_sha256'],
  'common_features_sha256':protocol['common_features']['sha256'],'parent_pins':protocol['parent_pins'],
  'scope':'Post-result diagnostic; CODE8/CODE9 outcomes already known; no new validation.'})
 data,frozen,panels,common_indices,fullchecks=load_data(root,protocol,common,design)
 prior=pd.read_csv(root.parent/protocol['parent_pins']['code9_correlations']['relative_path'])
 frames=part_A(data,frozen,common_indices,common,prior)
 frames['full_pool_reconstruction_checks']=pd.DataFrame(fullchecks)
 for name,frame in frames.items():frame.to_csv(args.output/(name+'.csv'),index=False)
 print('PART_A_COMPLETED',flush=True)
 frames.update(part_B(data,panels,common_indices,common,choices,protocol))
 for name,frame in frames.items():frame.to_csv(args.output/(name+'.csv'),index=False)
 # Recheck immutable sources, the fixed partition table, and common universe.
 check_inputs(root,args.protocol)
 summary={'status':'CODE10_POSTRESULT_RELIABILITY_COMPLETE','completed_utc':now(),
  'runtime_seconds':time.perf_counter()-start,'code_sha256':sha(__file__),'protocol_sha256':sha(args.protocol),
  'partition_sha256':protocol['part_B']['partition_sha256'],'common_features_sha256':protocol['common_features']['sha256'],
  'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,
  'parent_pins_unchanged':True,'new_validation':False,'confidence_intervals_computed':False,
  'n_draws':256,'n_side_fits':len(frames['part_B_fit_details']),
  'n_failed_side_fits':int(frames['part_B_fit_details'].fit_failed.sum()),
  'n_rank_deficient_single_spans':int(frames['part_B_fit_details'].single_span_rank_deficient.sum()),
  'n_rank_deficient_OLS_fits':int(frames['part_B_fit_details'].OLS_rank_deficient.sum()),
  'rows':{name:len(frame) for name,frame in frames.items()},
  'output_hashes':{p.name:sha(p) for p in sorted(args.output.iterdir()) if p.is_file()}}
 write_json(args.output/'RUN_SUMMARY.json',summary);print(json.dumps(summary,indent=2))

def selftest():
 rng=np.random.default_rng(20261030);b=rng.normal(size=80);d=rng.normal(size=80);y=rng.normal(size=(3,80))
 delta=gain_vectors(b,d,y);mean_gain=(b-y.mean(axis=0))**2-(d-y.mean(axis=0))**2
 assert np.allclose(delta.mean(axis=0),mean_gain,atol=2e-14)
 pooled_mse=float(np.mean((y-d)**2));mean_mse=float(np.mean((y.mean(axis=0)-d)**2));variation=float(np.mean((y-y.mean(axis=0))**2))
 assert variation>0 and np.isclose(pooled_mse,mean_mse+variation)
 info,gram,m=moment_reliability(delta);assert np.isclose(info['stable_variance'],sum(gram[i,j] for i in range(3) for j in range(3) if i!=j)/6)
 z=rng.normal(size=80);negative=np.stack([z,-z,np.zeros_like(z)]);info,_,_=moment_reliability(negative)
 assert info['stable_variance']<0 and info['conditional_reliability'] is None
 zz,_,_=moment_reliability(np.zeros((3,80)));assert zz['stable_variance_nonpositive'] and zz['conditional_reliability'] is None
 _,adj,reason,flag=adjusted_correlation(centered(z),centered(z),-1,1);assert adj is None and reason and not flag
 _,adj,reason,flag=adjusted_correlation(centered(z),centered(z),.01,.01);assert adj>1 and flag and not reason
 assert ordinary_correlation(np.ones(10),np.arange(10))[0] is None
 assert top_indices(np.ones(3),np.array(['z','a','m']),2)=={1,2}
 assert top_indices(np.array([-3.,-2.,-1.]),np.array(['a','b','c']),2)=={1,2}
 codes=np.random.Generator(np.random.PCG64(100)).choice(3**8,256,replace=False);assert len(set(codes))==256
 for code in codes:
  choice=np.array([(int(code)//3**i)%3 for i in range(8)])
  one=side_indices(choice,1);two=side_indices(choice,2)
  assert all(len(a)==1 and len(b)==2 and set(a).isdisjoint(b) and sorted(a+b)==[0,1,2] for a,b in zip(one,two))
 v=rng.normal(size=30);sm=np.stack([v,2*v,np.zeros_like(v)]);beta=np.array([.6,1.2,0.]);groups={'CT':np.zeros((2,30))}
 for i,name in enumerate(SINGLES):groups[name]=np.tile(sm[i],(2,1))
 for pair in PAIRS:groups[pair]=np.tile(sum(beta[i]*sm[i] for i,name in enumerate(SINGLES) if name in pair.split('+')),(2,1))
 pred,fit=fit_lower(groups);assert fit['single_span_rank']==1 and fit['OLS_rank']==2
 assert np.allclose([fit['beta_WS'],fit['beta_HL'],fit['beta_HS']],beta,atol=2e-12)
 assert np.allclose(pred[0],pred[3],atol=2e-12)
 zero={k:np.zeros((1,30)) for k in ('CT',*SINGLES,*PAIRS)};zero['WS+HL'][:]=v
 pred,fit=fit_lower(zero);assert fit['single_span_rank']==0 and fit['OLS_rank']==0 and np.allclose(pred[3],0)
 fake=pd.DataFrame({'group':['g']*4,'x':[0.,1.,2.,np.nan]});s=summarize(fake,['group'],['x']).iloc[0]
 assert s.n_finite==3 and s.n_unavailable==1 and s['median']==1 and np.isclose(s.q025,.05)
 print('SELFTEST_PASS: gain/target-variance identities; centered offdiagonal moments; negative/zero stable variance; unclipped adjusted ratio; exact top-k ties; unique complete partitions; disjoint sides; rank-deficient span/OLS; quantile/unavailable handling')

def main():
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True);sub.add_parser('selftest')
 q=sub.add_parser('run');q.add_argument('--protocol',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
 args=p.parse_args()
 if args.command=='selftest':selftest()
 else:run(args)
if __name__=='__main__':main()
