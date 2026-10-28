#!/usr/bin/env python3
"""Post-result gene attribution of frozen citrus predictions.

No model refitting or significance tests. All CODE8 outcomes were already known.
The score-removal and target-fitted direction diagnostics are explanatory only.
"""
from __future__ import annotations
import argparse,gzip,hashlib,json,math,platform,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd

PRIMARY_COMPARATORS=('second_order_factorial','pair_trained_factor_OLS','within_V3')
COMPARATORS=(*PRIMARY_COMPARATORS,'within_V11')
ORACLE_MODELS=('full_pair_centroid','second_order_factorial','pair_trained_factor_OLS')
EPS=1e-12
FRACTIONS=(.01,.05,.10)


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
def ratio(a,b):return float(a/b) if abs(b)>EPS else None

def gain_arrays(b,d,y):
 baseline_loss=(y-b)**2;centroid_loss=(y-d)**2
 delta=baseline_loss-centroid_loss;mean=delta.mean(axis=0)
 identity=(b-d)*(b+d-2*y.mean(axis=0))
 assert np.allclose(mean,identity,rtol=2e-10,atol=2e-11)
 return delta,baseline_loss,centroid_loss

def sign_summary(delta):
 g=delta.mean(axis=0);n=len(g)
 quantiles=(0.,.01,.05,.25,.5,.75,.95,.99,1.)
 qnames=('q00','q01','q05','q25','q50','q75','q95','q99','q100')
 if n==0:return {'n_features':0,'mean_gain':None,'median_gain':None,**{k+'_gain':None for k in qnames},'n_mean_positive':0,'n_mean_negative':0,'n_mean_zero_tolerance':0,'n_positive_all3':0,'fraction_mean_positive':None,'fraction_positive_all3':None,'positive_mass':0.,'negative_mass':0.,'negative_signed_mass':0.,'net_mass':0.}
 q=np.quantile(g,quantiles);pos=g>EPS;neg=g < -EPS
 return {'n_features':n,'mean_gain':float(g.mean()),'median_gain':float(np.median(g)),
  **{k+'_gain':float(v) for k,v in zip(qnames,q)},
  'n_mean_positive':int(pos.sum()),'n_mean_negative':int(neg.sum()),'n_mean_zero_tolerance':int((~pos&~neg).sum()),
  'n_positive_all3':int((delta>EPS).all(axis=0).sum()),'fraction_mean_positive':float(pos.mean()),
  'fraction_positive_all3':float((delta>EPS).all(axis=0).mean()),
  'positive_mass':float(g[g>0].sum()),'negative_mass':float(-g[g<0].sum()),'negative_signed_mass':float(g[g<0].sum()),'net_mass':float(g.sum())}

def removal_analysis(genes,delta,baseline_loss,centroid_loss):
 g=delta.mean(axis=0);n=len(g);order=np.lexsort((genes,-g));ranked=g[order]
 total=float(g.sum());positive=float(g[g>0].sum());cumulative=np.cumsum(ranked)
 curves=pd.DataFrame({'removed_rank':np.arange(1,n+1),'Feature_GID':genes[order],
  'removed_gene_mean_gain':ranked,'cumulative_removed_net_gain':cumulative,
  'remaining_net_gain':total-cumulative})
 rows=[]
 for fraction in FRACTIONS:
  k=math.ceil(fraction*n);keep=order[k:];removed=order[:k]
  if not len(keep):raise ValueError('Removal would leave no features')
  perlib=delta[:,keep].mean(axis=1)
  rows.append({'requested_fraction':fraction,'removed_count':k,'removed_fraction':k/n,
   'remaining_count':len(keep),'removed_net_mass':float(g[removed].sum()),
   'retained_net_mass':float(g[keep].sum()),'retained_mean_gain':float(g[keep].mean()),
   'retained_over_original_net_mass':ratio(g[keep].sum(),total),
   'positive_mass_share_removed':ratio(np.maximum(g[removed],0).sum(),positive),
   'retained_baseline_mse':float(baseline_loss[:,keep].mean()),'retained_centroid_mse':float(centroid_loss[:,keep].mean()),
   'retained_advantage_all3libraries':bool((perlib>EPS).all()),
   **{f'library_{j+1}_retained_gain':float(perlib[j]) for j in range(3)}})
 p=ranked[ranked>EPS]
 if total<=EPS:critical=0
 else:
  candidates=np.flatnonzero(total-np.cumsum(p)<=EPS)
  if len(candidates):critical=int(candidates[0]+1)
  elif total-float(p.sum())<=EPS:critical=len(p)
  else:critical=None
 critical_info={'minimum_top_gain_features_to_erase_net_advantage':critical,'minimum_removal_fraction':critical/n if critical is not None else None,
  'critical_prefix_status':'NO_POSITIVE_INITIAL_NET_GAIN' if total<=EPS else ('FOUND' if critical is not None else 'UNAVAILABLE_AFTER_EPSILON_POSITIVE_PREFIX'),
  'remaining_after_all_epsilon_positive':float(total-p.sum()),
  'n_positive_mass_features':len(p),'original_net_gain':total,'original_positive_mass':positive,
  'critical_is_outcome_dependent_score_diagnostic':True}
 return rows,curves,critical_info

def direction_oracle(pred,c,y):
 u=pred-c;t=y.mean(axis=0)-c;den=float(u@u);alpha=float(u@t)/den if den else 0.
 fitted=c+alpha*u
 original=(y-pred)**2;oracle=(y-fitted)**2
 drop=float(original.mean()-oracle.mean());identity=float((alpha-1)**2*np.mean(u*u))
 assert np.isclose(drop,identity,atol=2e-10,rtol=1e-10)
 assert drop>=-2e-10
 return {'alpha_target_fitted':alpha,'zero_direction_norm':den==0,'original_mean_mse':float(original.mean()),
  'target_fitted_mean_mse':float(oracle.mean()),'rescaling_recoverable_mse':drop,
  'rescaling_identity_mse':identity,'direction_norm_squared_per_feature':float(np.mean(u*u)),
  'oracle_label':'TARGET_FITTED_GLOBAL_RESCALE_DIAGNOSTIC_NOT_PREDICTION',
  'rescaled_mse_to_target_mean':float(np.mean((fitted-y.mean(axis=0))**2)),
  'direction_cosine_to_target_mean':float(u@t/np.sqrt(den*float(t@t))) if den and float(t@t) else None,
  'direction_cosine_unavailable_zero_norm':not bool(den and float(t@t)),
  **{f'library_{j+1}_oracle_mse':float(oracle[j].mean()) for j in range(3)}}

def correlation(a,b,rank=False):
 if rank:a=pd.Series(a).rank(method='average').to_numpy();b=pd.Series(b).rank(method='average').to_numpy()
 if len(a)<2 or np.std(a)==0 or np.std(b)==0:return None
 return float(np.corrcoef(a,b)[0,1])

def read_target(source,columns):
 t=pd.read_csv(source,sep='\t',usecols=['Feature_GID',*columns],dtype={'Feature_GID':str},keep_default_na=False)
 assert t.Feature_GID.is_unique and t.Feature_GID.ne('').all()
 values=t[columns].to_numpy(float).T
 assert np.isfinite(values).all() and (values>=0).all()
 return t.Feature_GID.to_numpy(str),np.log2(1+values)

def check_parent(root,protocol):
 m=json.loads((root/'frozen/FREEZE_MANIFEST.json').read_text());s=json.loads((root/'results/RUN_SUMMARY.json').read_text())
 paths={'protocol':root/'CODE8_PROTOCOL.json','amendment':root/'INPUT_AMENDMENT01.json',
  'freeze_manifest':root/'frozen/FREEZE_MANIFEST.json','result_summary':root/'results/RUN_SUMMARY.json',
  'source':root/'input/GSE203331_Expression_Profile.GCF_000317415.1_Csi_valencia_1.0.gene.txt.gz'}
 for key,record in protocol['parent_pins'].items():
  path=root.parent/record['relative_path'];assert sha(path)==record['sha256'],key
 assert s['decision']=='DESCRIPTIVE_TWO_GENOTYPE_CROP_TRANSPORT_SUCCESS'
 assert sha(root/'code/stress_citrus_validation.py')==m['code_sha256']
 for item in m['scenarios']:assert sha(root/'frozen'/item['assets'])==item['assets_sha256']
 for filename,digest in s['output_hashes'].items():assert sha(root/'results'/filename)==digest
 return m,paths['source']

def bool_column(series):
 vals=series.astype(str).str.lower();assert vals.isin(['true','false','1','0']).all()
 return vals.isin(['true','1']).to_numpy()

def read_annotation(directory,annex_path):
 annex=json.loads(annex_path.read_text())
 for filename in ('feature_annotation.csv.gz','memberships.csv.gz','categories.csv'):
  assert sha(directory/filename)==annex['output_hashes'][filename],filename
 a=pd.read_csv(directory/'feature_annotation.csv.gz',dtype=str,keep_default_na=False)
 assert a.Feature_GID.is_unique
 for col in ('EntrezID','Gene_Symbol','Product','ambiguous_entrez','GOcovered','BPcovered'):assert col in a,col
 for col in ('ambiguous_entrez','GOcovered','BPcovered'):a[col]=bool_column(a[col])
 a['valid_entrez']=a.EntrezID.str.fullmatch(r'[0-9]+')
 a=a.set_index('Feature_GID',drop=False)
 memberships=pd.read_csv(directory/'memberships.csv.gz',dtype=str,keep_default_na=False)
 categories=pd.read_csv(directory/'categories.csv',dtype=str,keep_default_na=False)
 assert categories.category_id.is_unique and len(categories)==6
 assert not memberships.duplicated(['Feature_GID','category_id']).any()
 assert set(memberships.Feature_GID)<=set(a.index)
 assert set(memberships.category_id)<=set(categories.category_id)
 sets={str(row.category_id):set(memberships.loc[memberships.category_id==row.category_id,'Feature_GID']) for row in categories.itertuples()}
 names=dict(zip(categories.category_id,categories.category_name))
 a['categories']=memberships.groupby('Feature_GID',sort=True).category_id.apply(lambda v:';'.join(sorted(v))).reindex(a.index).fillna('')
 return a,sets,names,annex

def group_masks(genes,annotation,sets,names):
 a=annotation.loc[genes];eligible=(~a.ambiguous_entrez&a.valid_entrez).to_numpy()
 covered=a.GOcovered.to_numpy();bp=a.BPcovered.to_numpy()
 anygroup=np.zeros(len(genes),dtype=bool);masks=[]
 for category,ids in sets.items():
  mask=np.array([g in ids for g in genes])&eligible;anygroup|=mask
  assert not np.any(mask&~bp)
  masks.append((category,names[category],'functional_category_overlaps_allowed',mask))
 masks += [('BP_annotated_other','BP annotated outside the six categories','annotation_partition',eligible&bp&~anygroup),
           ('GO_without_BP','GO annotated without biological-process annotation','annotation_partition',eligible&covered&~bp),
           ('GO_unannotated','No admissible GO annotation','annotation_partition',eligible&~covered),
           ('BP_annotated_background','All eligible BP-annotated features','background',eligible&bp),
           ('eligible_background','All valid nonambiguous eligible features','background',eligible)]
 return masks,eligible,covered,bp

def functional_rows(info,genes,delta,annotation,sets,names):
 masks,eligible,covered,bp=group_masks(genes,annotation,sets,names)
 overall=sign_summary(delta);bg=sign_summary(delta[:,eligible&bp]);rows=[]
 for category,name,role,mask in masks:
  stats=sign_summary(delta[:,mask]);n=int(mask.sum())
  rows.append({**info,'category_id':category,'category_name':name,'group_role':role,
   'small_group_n_below20':n<20,**stats,
   'net_gain_share_of_original_scenario':ratio(stats['net_mass'],overall['net_mass']),
   'positive_mass_share_of_original_scenario':ratio(stats['positive_mass'],overall['positive_mass']),
   'fraction_positive_minus_BP_background':None if n==0 or bg['n_features']==0 else stats['fraction_mean_positive']-bg['fraction_mean_positive'],
   'n_original_scenario_features':len(genes),'n_functionally_ineligible_excluded':int((~eligible).sum()),
   'n_BP_eligible_background':bg['n_features']})
 return rows

def exemplar_rows(info,genes,delta,baseline_loss,centroid_loss,annotation):
 mean=delta.mean(axis=0);orderpos=np.lexsort((genes,-mean));orderneg=np.lexsort((genes,mean))
 selected=[('largest_positive_gain',i) for i in orderpos[mean[orderpos]>EPS][:20]]+[
  ('largest_negative_gain',i) for i in orderneg[mean[orderneg]<-EPS][:10]]
 rows=[]
 for reason,i in selected:
  a=annotation.loc[genes[i]]
  rows.append({**info,'selection':reason,'Feature_GID':genes[i],
   'mean_gain':float(mean[i]),'minimum_library_gain':float(delta[:,i].min()),
   'positive_all3libraries':bool((delta[:,i]>EPS).all()),
   'baseline_gene_mse':float(baseline_loss[:,i].mean()),'centroid_gene_mse':float(centroid_loss[:,i].mean()),
   **{k:a[k] for k in ('EntrezID','Gene_Symbol','Product','GOcovered','BPcovered','ambiguous_entrez','categories')}})
 return rows

def shared_analysis(cache,annotation,sets,names):
 threshold=1.;summary=[];gene_tables=[];function_rows=[];examples=[]
 keys={(x['genotype'],x['comparator']):x for x in cache if x['threshold_fpkm']==threshold}
 ids=np.array(sorted(set(keys[('Car',PRIMARY_COMPARATORS[0])]['genes'])&set(keys[('Cleo',PRIMARY_COMPARATORS[0])]['genes'])))
 masks,eligible,covered,bp=group_masks(ids,annotation,sets,names)
 cores={};ranks={};data={};displacements={}
 for comparator in PRIMARY_COMPARATORS:
  ca=keys[('Car',comparator)];cl=keys[('Cleo',comparator)]
  ia=pd.Index(ca['genes']).get_indexer(ids);il=pd.Index(cl['genes']).get_indexer(ids)
  assert (ia>=0).all() and (il>=0).all()
  da=ca['delta'][:,ia];dl=cl['delta'][:,il];ga=da.mean(axis=0);gl=dl.mean(axis=0)
  ta=ca['target_displacement'][ia];tl=cl['target_displacement'][il]
  displacements={'Car':ta,'Cleo':tl}
  pa=ga>EPS;pl=gl>EPS;core=(da>EPS).all(axis=0)&(dl>EPS).all(axis=0)
  cores[comparator]=core;data[comparator]=(da,dl,ga,gl)
  ra=pd.Series(ga).rank(method='average',pct=True).to_numpy();rl=pd.Series(gl).rank(method='average',pct=True).to_numpy()
  ranks['Car_'+comparator]=ra;ranks['Cleo_'+comparator]=rl
  sa=np.where(ta>EPS,1,np.where(ta < -EPS,-1,0));sl=np.where(tl>EPS,1,np.where(tl < -EPS,-1,0))
  same=(sa==sl)&(sa!=0);opposite=(sa*sl==-1);zero=(sa==0)|(sl==0)
  summary.append({'threshold_fpkm':1.,'primary_threshold':True,'comparator':comparator,
   'n_original_Car_features':len(ca['genes']),'n_original_Cleo_features':len(cl['genes']),
   'n_common_features':len(ids),'pearson_gain':correlation(ga,gl),'spearman_gain':correlation(ga,gl,rank=True),
   'n_positive_both':int((pa&pl).sum()),'n_positive_Car_only':int((pa&~pl).sum()),
   'n_positive_Cleo_only':int((~pa&pl).sum()),'n_nonpositive_both':int((~pa&~pl).sum()),
   'n_with_zero_tolerance_either_gain':int(((np.abs(ga)<=EPS)|(np.abs(gl)<=EPS)).sum()),
   'n_all3_positive_both_core':int(core.sum()),'core_fraction':float(core.mean()),
   'Car_common_net_gain':float(ga.sum()),'Cleo_common_net_gain':float(gl.sum()),
   'Car_core_net_gain':float(ga[core].sum()),'Cleo_core_net_gain':float(gl[core].sum()),
   'Car_core_share_net':ratio(ga[core].sum(),ga.sum()),'Cleo_core_share_net':ratio(gl[core].sum(),gl.sum()),
   'Car_core_share_positive':ratio(ga[core].sum(),ga[ga>0].sum()),'Cleo_core_share_positive':ratio(gl[core].sum(),gl[gl>0].sum()),
   'common_same_nonzero_response_direction':int(same.sum()),'common_opposite_response_direction':int(opposite.sum()),
   'common_zero_response_either':int(zero.sum()),'core_same_nonzero_response_direction':int((core&same).sum()),
   'core_opposite_response_direction':int((core&opposite).sum()),'core_zero_response_either':int((core&zero).sum())})
  frame=pd.DataFrame({'threshold_fpkm':1.,'comparator':comparator,'Feature_GID':ids,
   'Car_mean_gain':ga,'Cleo_mean_gain':gl,'Car_min_library_gain':da.min(axis=0),'Cleo_min_library_gain':dl.min(axis=0),
   'shared_core_all3_each':core,'Car_gain_percentile':ra,'Cleo_gain_percentile':rl,
   'Car_target_minus_control':ta,'Cleo_target_minus_control':tl,'same_nonzero_response_direction':same,'opposite_response_direction':opposite,'zero_response_either':zero})
  gene_tables.append(frame)
 strict=np.logical_and.reduce(list(cores.values()))
 minimum_rank=np.stack(list(ranks.values())).min(axis=0)
 strict_table=pd.DataFrame({'Feature_GID':ids,'strict_core_all3_comparators':strict,
  'minimum_all6_gain_percentile':minimum_rank,**{k+'_percentile':v for k,v in ranks.items()},
  'Car_target_minus_control':displacements['Car'],'Cleo_target_minus_control':displacements['Cleo']})
 for row in strict_table[strict].sort_values(['minimum_all6_gain_percentile','Feature_GID'],ascending=[False,True]).head(20).to_dict('records'):
  a=annotation.loc[row['Feature_GID']]
  examples.append({**row,**{k:a[k] for k in ('EntrezID','Gene_Symbol','Product','GOcovered','BPcovered','ambiguous_entrez','categories')}})
 strict_info={'threshold_fpkm':1.,'n_common_features':len(ids),'n_strict_core_all3_comparators':int(strict.sum()),
  'strict_core_fraction':float(strict.mean()),'n_strict_core_BP_eligible':int((strict&eligible&bp).sum()),
  'n_strict_core_functionally_ineligible':int((strict&~eligible).sum())}
 for comparator in (*PRIMARY_COMPARATORS,'ALL_THREE_COMPARATORS'):
  core=cores[comparator] if comparator in cores else strict
  bg=eligible&bp;bgprev=float(core[bg].mean()) if bg.sum() else None
  allprev=float(core[eligible].mean()) if eligible.sum() else None
  for category,name,role,mask in masks:
   n=int(mask.sum());ncore=int((core&mask).sum());prev=ncore/n if n else None
   function_rows.append({'threshold_fpkm':1.,'comparator':comparator,'category_id':category,
    'category_name':name,'group_role':role,'n_common_eligible_group':n,'n_shared_core':ncore,
    'core_prevalence':prev,'BP_background_core_prevalence':bgprev,
    'relative_core_prevalence_vs_BP_background':ratio(prev,bgprev) if prev is not None and bgprev is not None else None,
    'all_eligible_background_core_prevalence':allprev,
    'relative_core_prevalence_vs_all_eligible':ratio(prev,allprev) if prev is not None and allprev is not None else None,
    'small_group_n_below20':n<20,'n_BP_eligible_background':int(bg.sum())})
 for frame in gene_tables:frame['strict_core_all3_comparators']=strict
 return summary,pd.concat(gene_tables,ignore_index=True),strict_info,function_rows,examples,strict_table

def run(args):
 start=time.perf_counter();fresh(args.output)
 protocol=json.loads(args.protocol.read_text())
 assert protocol['protocol_id']=='CODE9_CITRUS_POSTRESULT_GENE_ATTRIBUTION_V1'
 assert protocol['comparators']==list(COMPARATORS) and protocol['primary_comparators']==list(PRIMARY_COMPARATORS)
 assert protocol['epsilon']==EPS and protocol['concentration']['fractions']==list(FRACTIONS)
 manifest,source=check_parent(args.code8_root,protocol)
 annotation,sets,names,annex=read_annotation(args.annotation_dir,args.annotation_annex)
 write_json(args.output/'RUN_STARTED.json',{'started_utc':now(),'protocol_sha256':sha(args.protocol),
  'annotation_annex_sha256':sha(args.annotation_annex),'annotation_output_hashes':annex['output_hashes'],
  'code_sha256':sha(__file__),'parent_pins':protocol['parent_pins'],
  'scope':'All CODE8 targets already known. Post-result attribution, no new validation.'})
 targets={}
 for g in ('Car','Cleo'):
  item=next(x for x in manifest['scenarios'] if x['genotype']==g)
  targets[g]=read_target(source,item['excluded_target_columns'])
 expected=pd.read_csv(args.code8_root/'results/comparator_summary.csv')
 cache=[];summaries=[];removals=[];critical=[];curves=[];genes_outputs=[];oracles=[]
 functions=[];qc=[];coverage=[];examples=[];identities=[];inventory=[]
 for item in manifest['scenarios']:
  info={'scenario':item['key'],'genotype':item['genotype'],'threshold_fpkm':item['threshold_fpkm'],
        'primary_threshold':item['primary_threshold'],'genes':item['genes']}
  with np.load(args.code8_root/'frozen'/item['assets'],allow_pickle=False) as a:
   genes=a['genes'];pred=a['predictions'];models=a['model_names'].tolist()
   source_ids,source_y=targets[item['genotype']];indices=pd.Index(source_ids).get_indexer(genes)
   assert (indices>=0).all() and len(set(genes))==len(genes) and len(genes)==info['genes']
   y=source_y[:,indices];c=pred[models.index('control')];d=pred[models.index('full_pair_centroid')]
   assert (annotation.loc[genes,'Feature_GID'].to_numpy()==genes).all()
   anns=annotation.loc[genes];amb=anns.ambiguous_entrez.to_numpy();valid=anns.valid_entrez.to_numpy()
   masks,eligible,go,bp=group_masks(genes,annotation,sets,names)
   in_six=np.array([bool(x) for x in anns.categories]);group_counts=anns.categories.str.count(';').to_numpy()+in_six.astype(int)
   coverage.append({**info,'n_valid_entrez':int(valid.sum()),'n_ambiguous_entrez_excluded':int(amb.sum()),
    'n_functionally_eligible':int(eligible.sum()),'n_any_GO_eligible':int((eligible&go).sum()),
    'n_BP_eligible':int((eligible&bp).sum()),'n_any_GO_without_BP_eligible':int((eligible&go&~bp).sum()),
    'n_GO_unannotated_eligible':int((eligible&~go).sum()),'n_six_group_union':int((eligible&in_six).sum()),
    'n_multiple_group_memberships':int((eligible&(group_counts>1)).sum()),
    'any_GO_coverage_of_original_panel':float(go.mean()),'BP_coverage_of_original_panel':float(bp.mean())})
   inventory.append({**info,**{f'target_column_{j+1}':str(a['target_columns'][j]) for j in range(3)}})
   for model in ORACLE_MODELS:
    row=direction_oracle(pred[models.index(model)],c,y)
    row.update({**info,'model':model,'original_centroid_mean_mse':float(((y-d)**2).mean()),
                'rescaled_model_mse_minus_original_centroid_mse':row['target_fitted_mean_mse']-float(((y-d)**2).mean())})
    oracles.append(row)
   for comparator in COMPARATORS:
    b=pred[models.index(comparator)];delta,bl,dl=gain_arrays(b,d,y);g=delta.mean(axis=0)
    base={**info,'comparator':comparator,'comparator_role':'primary' if comparator in PRIMARY_COMPARATORS else 'secondary'}
    stat=sign_summary(delta);summaries.append({**base,**stat,'baseline_mean_mse':float(bl.mean()),'centroid_mean_mse':float(dl.mean())})
    want=expected[(expected.scenario==item['key'])&(expected.comparator==comparator)].iloc[0].mean_comparator_mse_minus_centroid_mse
    discrepancy=abs(float(g.mean())-float(want));assert discrepancy<2e-10
    identities.append({**base,'gene_mean_gain':float(g.mean()),'CODE8_mean_gain':float(want),
       'absolute_difference':discrepancy,'factorized_identity_max_error':float(np.max(np.abs(g-(b-d)*(b+d-2*y.mean(axis=0))))),
       'mass_identity_error':abs(stat['positive_mass']-stat['negative_mass']-stat['net_mass'])})
    rr,curve,cc=removal_analysis(genes,delta,bl,dl)
    removals.extend([{**base,**x} for x in rr]);critical.append({**base,**cc})
    for k,v in base.items():curve[k]=v
    curves.append(curve)
    gf=pd.DataFrame({'Feature_GID':genes,'EntrezID':anns.EntrezID.to_numpy(),'mean_gain':g,
     'gain_library1':delta[0],'gain_library2':delta[1],'gain_library3':delta[2],
     'baseline_prediction':b,'centroid_prediction':d,'control_mean':c,'target_mean':y.mean(axis=0),
     'target_library1':y[0],'target_library2':y[1],'target_library3':y[2],
     'baseline_gene_mse':bl.mean(axis=0),'centroid_gene_mse':dl.mean(axis=0),
     'positive_all3libraries':(delta>EPS).all(axis=0)})
    for k,v in base.items():gf[k]=v
    genes_outputs.append(gf)
    keep=~amb;qr=sign_summary(delta[:,keep]);perlib=delta[:,keep].mean(axis=1)
    qc.append({**base,'n_removed_ambiguous_Entrez_features':int(amb.sum()),**qr,
      'baseline_mse_after_duplicate_only_exclusion':float(bl[:,keep].mean()),
      'centroid_mse_after_duplicate_only_exclusion':float(dl[:,keep].mean()),
      **{f'library_{j+1}_retained_gain':float(perlib[j]) for j in range(3)}})
    functions.extend(functional_rows(base,genes,delta,annotation,sets,names))
    if info['primary_threshold'] and comparator in PRIMARY_COMPARATORS:
     examples.extend(exemplar_rows(base,genes,delta,bl,dl,annotation))
    cache.append({**base,'genes':genes.copy(),'delta':delta.copy(),'target_displacement':y.mean(axis=0)-c})
  print('COMPUTED '+item['key'],flush=True)
 cross,shared,strict,corefunctions,coreexamples,strictgenes=shared_analysis(cache,annotation,sets,names)
 frames={'gene_summary':pd.DataFrame(summaries),'concentration_removal':pd.DataFrame(removals),
  'critical_removal':pd.DataFrame(critical),'direction_oracles':pd.DataFrame(oracles),
  'functional_summary':pd.DataFrame(functions),'annotation_coverage':pd.DataFrame(coverage),
  'duplicate_Entrez_score_QC':pd.DataFrame(qc),'gene_examples':pd.DataFrame(examples),
  'cross_genotype_summary':pd.DataFrame(cross),'strict_core_summary':pd.DataFrame([strict]),
  'functional_core_summary':pd.DataFrame(corefunctions),'strict_core_examples':pd.DataFrame(coreexamples),
  'identity_checks':pd.DataFrame(identities),'scenario_inventory':pd.DataFrame(inventory)}
 for name,frame in frames.items():frame.to_csv(args.output/(name+'.csv'),index=False)
 bigframes={'gene_contributions':pd.concat(genes_outputs,ignore_index=True),
            'concentration_curves':pd.concat(curves,ignore_index=True),
            'shared_gene_gains':shared,'strict_core_all_genes':strictgenes}
 for name,frame in bigframes.items():frame.to_csv(args.output/(name+'.csv.gz'),index=False,compression={'method':'gzip','compresslevel':6,'mtime':0})
 check_parent(args.code8_root,protocol)
 for filename,digest in annex['output_hashes'].items():
  path=args.annotation_dir/filename
  if path.is_file():assert sha(path)==digest
 result={'status':'CODE9_POSTRESULT_ATTRIBUTION_COMPLETE','completed_utc':now(),'runtime_seconds':time.perf_counter()-start,
  'protocol_sha256':sha(args.protocol),'annotation_annex_sha256':sha(args.annotation_annex),'code_sha256':sha(__file__),
  'parent_unchanged':True,'new_heldout_validation':False,'formal_significance_test':False,
  'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,
  'rows':{k:len(v) for k,v in {**frames,**bigframes}.items()},
  'strict_core':strict,'primary_concentration':frames['critical_removal'][frames['critical_removal'].primary_threshold].to_dict('records'),
  'output_hashes':{p.name:sha(p) for p in sorted(args.output.iterdir()) if p.is_file()}}
 write_json(args.output/'RUN_SUMMARY.json',result);print(json.dumps(result,indent=2))

def selftest():
 import itertools
 rng=np.random.default_rng(9012026);b=rng.normal(size=90);d=rng.normal(size=90);y=rng.normal(size=(3,90))
 delta,bl,dl=gain_arrays(b,d,y);assert np.isclose(delta.mean(),bl.mean()-dl.mean())
 centered=y-y.mean(axis=0);altered=y.mean(axis=0)+3*centered
 changed,_,_=gain_arrays(b,d,altered);assert np.allclose(delta.mean(axis=0),changed.mean(axis=0),atol=2e-14)
 vals=np.array([8.,4.,2.,1.,-3.,-1.]);genes=np.array(['b','a','c','d','e','f']);dt=np.tile(vals,(3,1))
 rr,curve,crit=removal_analysis(genes,dt,np.maximum(dt,0),np.maximum(-dt,0))
 assert crit['minimum_top_gain_features_to_erase_net_advantage']==2
 assert rr[0]['retained_net_mass']==3 and rr[0]['retained_mean_gain']==.6
 for k in (1,2,3):
  minimum=min(sum(vals[i] for i in range(6) if i not in removed) for removed in itertools.combinations(range(6),k))
  assert np.isclose(curve.iloc[k-1].remaining_net_gain,minimum)
 equal=np.ones((3,3));_,cv,_=removal_analysis(np.array(['z','a','m']),equal,equal,np.zeros_like(equal))
 assert cv.Feature_GID.tolist()==['a','m','z']
 ss=sign_summary(dt);assert ss['positive_mass']==15 and ss['negative_mass']==4 and ss['net_mass']==11
 tiny=np.full((3,4),.6e-12);_,_,cc=removal_analysis(np.array(['a','b','c','d']),tiny,tiny,np.zeros_like(tiny))
 assert cc['minimum_top_gain_features_to_erase_net_advantage'] is None
 negative=-equal;_,_,cc=removal_analysis(np.array(['a','b','c']),negative,np.zeros_like(negative),-negative)
 assert cc['minimum_top_gain_features_to_erase_net_advantage']==0
 control=rng.normal(size=90);pred=control+2.5*(y.mean(axis=0)-control)
 oo=direction_oracle(pred,control,y);assert np.isclose(oo['alpha_target_fitted'],.4)
 zz=direction_oracle(control,control,y);assert zz['zero_direction_norm'] and zz['rescaling_recoverable_mse']==0
 # Exact feature identity joins; BP-only background; duplicate exclusion.
 ann=pd.DataFrame({'Feature_GID':['a','b','c','d'],'EntrezID':['1','2','3','3'],
  'Gene_Symbol':['a','b','c','d'],'Product':['p']*4,'ambiguous_entrez':[False,False,True,True],
  'valid_entrez':[True]*4,'GOcovered':[True,True,False,False],'BPcovered':[True,False,False,False],
  'categories':['x','','','']}).set_index('Feature_GID',drop=False)
 masks,eligible,go,bp=group_masks(np.array(['a','b','c','d']),ann,{'x':{'a'}},{'x':'test'})
 mm={name:mask for name,_,_,mask in masks}
 assert mm['BP_annotated_background'].tolist()==[True,False,False,False]
 assert mm['GO_without_BP'].tolist()==[False,True,False,False]
 cache=[]
 for genotype in ('Car','Cleo'):
  for j,comparator in enumerate(PRIMARY_COMPARATORS):
   v=np.array([4.,3.,-2.,1.])
   if genotype=='Cleo':v=np.array([4.,-3.,2.,1.])
   if j==2:v[3]=-1
   cache.append({'genotype':genotype,'comparator':comparator,'threshold_fpkm':1.,'genes':np.array(['a','b','c','d']),
    'delta':np.tile(v,(3,1)),'target_displacement':np.array([1.,-1.,1.,0.])})
 _,shared,strict,_,examples,_=shared_analysis(cache,ann,{'x':{'a'}},{'x':'test'})
 assert strict['n_strict_core_all3_comparators']==1 and examples[0]['Feature_GID']=='a'
 assert len([k for k in examples[0] if k.endswith('_percentile')])==7
 print('SELFTEST_PASS: loss identity; target variance cancellation; exhaustive fixed-count deletion; lexical ties; positive/negative mass; prefix edge cases; direction oracle; BP background; duplicate exclusion; strict six-rank core')

def main():
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True);sub.add_parser('selftest')
 q=sub.add_parser('run')
 for arg in ('code8-root','protocol','annotation-dir','annotation-annex','output'):q.add_argument('--'+arg,type=Path,required=True)
 args=p.parse_args()
 if args.command=='selftest':selftest()
 else:run(args)
if __name__=='__main__':main()
