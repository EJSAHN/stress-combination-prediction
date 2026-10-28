"""Replay audited numerical analyses and check untouched historical results.

Scientific engines retain their exact original bytes under provenance so their
recorded hashes remain meaningful. Public stage names describe scientific tasks.
This is a reproduction of known outcomes, never a new blinded validation.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os,platform,subprocess,sys,time
import numpy as np
import pandas as pd

JSON_VECTOR_COLUMNS = {
 'tables/noise/latent_mean_diagnostics.csv': {
  'latent_single_singular_values', 'latent_basis_singular_values'
 }
}

def compare_json_vector(left, right):
 """Compare declared numeric vectors as numbers, never as serialized text."""
 x,y=json.loads(str(left)),json.loads(str(right))
 for value in (x,y):
  if not isinstance(value,list) or not value or any(type(v) not in (int,float) for v in value):
   raise ValueError('Singular values must be a nonempty, one-dimensional numeric JSON array')
 xv,yv=np.asarray(x,dtype=float),np.asarray(y,dtype=float)
 if xv.shape!=yv.shape or not np.isfinite(xv).all() or not np.isfinite(yv).all():
  raise ValueError('Singular-value arrays must have equal shapes and finite entries')
 assert np.allclose(xv,yv,atol=5e-10,rtol=5e-10),('numeric vector mismatch',x,y)
 return len(xv),float(np.max(np.abs(xv-yv)))

def now():return datetime.now(timezone.utc).isoformat()
def digest(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def write_json(path,obj):
 Path(path).parent.mkdir(parents=True,exist_ok=True)
 Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def verify_assets(repo):
 m=json.loads((repo/'config/analysis_assets.json').read_text(encoding='utf-8'))
 for item in m['files']:
  p=repo/item['relative_path']
  if not p.is_file():raise FileNotFoundError('Missing packaged input: '+item['relative_path']+'. For the code-only bundle, run scripts/fetch_inputs.py explicitly or use the complete local bundle.')
  if p.stat().st_size!=item['bytes'] or digest(p)!=item['sha256']:raise ValueError('Input integrity mismatch: '+item['relative_path'])
 return len(m['files'])

def compare_arrays(new,old,stage):
 rows=[]
 with np.load(new,allow_pickle=False) as a,np.load(old,allow_pickle=False) as b:
  assert set(a.files)==set(b.files),(stage,'NPZ keys')
  for name in a.files:
   x,y=a[name],b[name];assert x.shape==y.shape,(stage,name,'shape')
   if x.dtype.kind in 'USb':
    assert np.array_equal(x,y),(stage,name,'exact array mismatch');error=0.;method='exact'
   elif 'basis' in name and x.ndim==2:
    # SVD bases may differ by signs/orthogonal rotations across BLAS builds.
    assert np.allclose(x@x.T,np.eye(len(x)),atol=5e-9),(stage,name,'fresh basis not orthonormal')
    assert np.allclose(y@y.T,np.eye(len(y)),atol=5e-9),(stage,name,'reference basis not orthonormal')
    error=max(float(np.max(np.abs(x-(x@y.T)@y))),float(np.max(np.abs(y-(y@x.T)@x)))) if x.size else 0.
    assert error<5e-9,(stage,name,error);method='same row space'
   else:
    assert np.allclose(x,y,atol=5e-10,rtol=5e-10,equal_nan=True),(stage,name,'numeric mismatch')
    error=float(np.nanmax(np.abs(x-y))) if x.size else 0.;method='numeric allclose'
   rows.append({'stage':stage,'array':name,'shape':str(x.shape),'comparison':method,'max_absolute_error':error,'status':'PASS'})
 return rows

def compare_tables(repo,output):
 refs=json.loads((repo/'config/table_reference_map.json').read_text(encoding='utf-8'))['tables']
 records=[]
 for item in refs:
  expected=repo/item['reference'];actual=output/item['output']
  assert digest(expected)==item['sha256'],'Reference changed: '+item['reference']
  a=pd.read_csv(actual,keep_default_na=False);b=pd.read_csv(expected,keep_default_na=False)
  assert list(a.columns)==list(b.columns),(item['output'],'columns')
  assert a.shape==b.shape,(item['output'],'shape',a.shape,b.shape)
  numeric_cells=0;text_cells=0;max_error=0.
  for col in a:
   x=a[col];y=b[col]
   if col in JSON_VECTOR_COLUMNS.get(item['output'],set()):
    for left,right in zip(x,y):
     n,error=compare_json_vector(left,right)
     numeric_cells+=n;max_error=max(max_error,error)
    continue
   # Empty cells are preserved separately; do not coerce labels into zero.
   xe=x.astype(str).eq('');ye=y.astype(str).eq('');assert xe.equals(ye),(item['output'],col,'missingness')
   xx=pd.to_numeric(x[~xe],errors='coerce');yy=pd.to_numeric(y[~ye],errors='coerce')
   if len(xx) and xx.notna().all() and yy.notna().all():
    xv=xx.to_numpy(float);yv=yy.to_numpy(float)
    assert np.allclose(xv,yv,atol=5e-10,rtol=5e-10,equal_nan=True),(item['output'],col,'numeric mismatch',float(np.max(abs(xv-yv))))
    numeric_cells+=len(xv)
    if len(xv):max_error=max(max_error,float(np.max(abs(xv-yv))))
   else:
    assert x.astype(str).equals(y.astype(str)),(item['output'],col,'text mismatch')
    text_cells+=len(x)
  records.append({'stage':item['stage'],'table':item['output'],'rows':len(a),'numeric_cells':numeric_cells,'text_cells':text_cells,'max_absolute_error':max_error,'status':'PASS'})
 return records

def run(repo,output):
 started=time.perf_counter();repo=repo.resolve();output=output.resolve()
 if output.exists() and any(output.iterdir()):raise ValueError('Use a fresh output directory. Previous results are never overwritten.')
 output.mkdir(parents=True,exist_ok=True)
 tables=output/'tables';work=output/'working';logs=output/'logs';validation=output/'validation'
 for d in (tables,work,logs,validation,work/'temp'):d.mkdir(parents=True,exist_ok=True)
 h=repo/'provenance/historical'
 def historical(n):return h/f'stress_combinations_code{n}'
 env=os.environ.copy()
 env.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')})
 env.update({'TMP':str(work/'temp'),'TEMP':str(work/'temp'),'TMPDIR':str(work/'temp'),'PYTHONDONTWRITEBYTECODE':'1'})
 commands=[];arrays=[]
 def command(label,script,*args):
  argv=[sys.executable,str(script),*map(str,args)]
  print('START '+label,flush=True);t=time.perf_counter()
  with (logs/(label+'.log')).open('w',encoding='utf-8') as f:
   completed=subprocess.run(argv,cwd=repo,env=env,stdout=f,stderr=subprocess.STDOUT,text=True)
  record={'stage':label,'argv':argv,'exit_code':completed.returncode,'seconds':time.perf_counter()-t,'log':str((logs/(label+'.log')).relative_to(output))}
  commands.append(record);write_json(logs/'execution_commands.json',commands)
  if completed.returncode:
   tail=(logs/(label+'.log')).read_text(encoding='utf-8',errors='replace')[-6000:]
   raise RuntimeError(label+' failed. '+tail)
  print('PASS '+label,flush=True)
 try:
  count=verify_assets(repo)
  write_json(validation/'run_started.json',{'started_utc':now(),'known_outcomes_reproduction':True,'new_external_validation':False,'verified_historical_assets':count,'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__})
  # Fresh fits precede evaluation in every prediction stage. Repeating this
  # separation does not make historically opened targets unseen again.
  d=historical(3);source=d/'reproduction/data/GSE147962_normalized_counts.txt.gz';design=d/'reproduction/data/GSE147962_sample_design.csv'
  command('discovery_selftest',d/'code/freeze_predictions.py','--self-test')
  df=work/'discovery_frozen'
  command('discovery_fit',d/'code/freeze_predictions.py','--data',source,'--design',design,'--gene-panel',d/'reproduction/data/gene_panel.csv','--code2-module',d/'code/reference_code2_models.py','--protocol',d/'CODE3_PROTOCOL.json','--output',df)
  arrays+=compare_arrays(df/'frozen_predictions.npz',d/'frozen/models/frozen_predictions.npz','discovery_prediction')
  command('discovery_prediction',d/'code/evaluate_frozen_predictions.py','--data',source,'--design',design,'--protocol',d/'CODE3_PROTOCOL.json','--freeze-dir',df,'--output',tables/'discovery_prediction','--authorize-target-read')
  d=historical(5);df=work/'selection_frozen';script=d/'code/donor_selection.py'
  command('selection_selftest',script,'self-test')
  common=['--data',source,'--design',design,'--protocol',d/'CODE5_PROTOCOL.json']
  command('discovery_selection_fit',script,'freeze',*common,'--output',df)
  arrays+=compare_arrays(df/'predictions.npz',d/'frozen/predictions.npz','discovery_selection')
  command('discovery_selection',script,'evaluate',*common,'--frozen',df,'--output',tables/'discovery_selection')
  d=historical(6);script=d/'code/external_three_factor.py';ef=work/'external_frozen'
  common=['--wt',d/'input/GSE330879_Col0_fpkm.txt.gz','--mutant',d/'input/GSE330879_RAP2.3_fpkm.txt.gz','--design',d/'fallback_metadata/sample_inventory.csv','--protocol',d/'CODE6_PROTOCOL.json','--amendment',d/'INPUT_AMENDMENT_01.json']
  command('external_selftest',script,'selftest')
  command('external_fit',script,'freeze',*common,'--output',ef)
  for f in sorted(ef.glob('*.npz')):arrays+=compare_arrays(f,d/'frozen_verified'/f.name,'external_prediction/'+f.stem)
  command('external_prediction',script,'evaluate',*common,'--frozen',ef,'--output',tables/'external_prediction')
  c=historical(8);script=c/'code/stress_citrus_validation.py';cf=work/'citrus_frozen';poison=work/'citrus_target_exclusion'
  common=['--source',c/'input/GSE203331_Expression_Profile.GCF_000317415.1_Csi_valencia_1.0.gene.txt.gz','--design',c/'input/GSE203331_design_verified.csv','--protocol',c/'CODE8_PROTOCOL.json','--amendment',c/'INPUT_AMENDMENT01.json']
  command('citrus_selftest',script,'selftest')
  command('citrus_fit',script,'freeze',*common,'--output',cf)
  for f in sorted(cf.glob('*.npz')):arrays+=compare_arrays(f,c/'frozen'/f.name,'citrus_prediction/'+f.stem)
  command('citrus_target_exclusion',script,'poison-check',*common,'--frozen',cf,'--output',poison)
  command('citrus_prediction',script,'evaluate',*common,'--frozen',cf,'--poison-report',poison/'POISON_CHECK.json','--output',tables/'citrus_prediction')
  # Diagnostics preserve the historically fixed predictions/feature panels.
  # The preceding source-level fits independently reproduced those arrays.
  d=historical(7)
  for name,module in [('external_span','empirical_span_expansion.py'),('noise','noisy_span_simulation.py')]:
   command(name+'_selftest',d/'code'/module,'selftest')
   extra=['--batch-size','16'] if name=='noise' else []
   command(name,d/'code'/module,'run','--code6-root',historical(6),'--protocol',d/'CODE7_PROTOCOL.json','--output',tables/name,*extra)
   if name=='noise':
    from .noise_diagnostics import canonicalize_noise_diagnostics
    print('START noise_diagnostic_identity_check',flush=True)
    diagnostic_correction=canonicalize_noise_diagnostics(repo,tables/name)
    write_json(validation/'noise_diagnostic_identity_correction.json',diagnostic_correction)
    print('PASS noise_diagnostic_identity_check',flush=True)
  d=historical(9);script=d/'code/gene_gain_attribution.py'
  command('attribution_selftest',script,'selftest')
  command('citrus_attribution',script,'run','--code8-root',historical(8),'--protocol',d/'CODE9_PROTOCOL.json','--annotation-dir',d/'annotation','--annotation-annex',d/'annotation/ANNOTATION_ANNEX.json','--output',tables/'citrus_attribution')
  d=historical(10);script=d/'code/attribution_reliability.py'
  command('repeatability_selftest',script,'selftest')
  command('attribution_repeatability',script,'run','--protocol',d/'CODE10_PROTOCOL.json','--output',tables/'attribution_repeatability')
  refs=compare_tables(repo,output)
  assert verify_assets(repo)==count
  pd.DataFrame(arrays).to_csv(validation/'frozen_array_checks.csv',index=False)
  pd.DataFrame(refs).to_csv(validation/'historical_table_checks.csv',index=False)
  result={'status':'PASS','completed_utc':now(),'runtime_seconds':time.perf_counter()-started,'historical_asset_files_verified':count,'fresh_array_checks':len(arrays),'historical_tables_checked':len(refs),'numeric_cells_checked':sum(x['numeric_cells'] for x in refs),'max_table_absolute_error':max(x['max_absolute_error'] for x in refs),'maximum_array_error':max(x['max_absolute_error'] for x in arrays),'all_scientific_stages_recomputed':True,'simulation_scenarios':36,'simulation_outcomes':36000,'disjoint_refits':1024,'new_external_validation':False,'windows_native_execution':platform.system()=='Windows','platform':platform.platform(),'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__}
  write_json(validation/'analysis_validation.json',result)
  write_json(output/'analysis_manifest.json',{'status':'PASS','tables':[{ 'relative_path':str(p.relative_to(output)),'sha256':digest(p),'bytes':p.stat().st_size} for p in sorted(tables.rglob('*')) if p.is_file()]})
  print(json.dumps(result,indent=2));return result
 except Exception as error:
  write_json(validation/'analysis_validation.json',{'status':'FAIL','failed_utc':now(),'error':str(error),'runtime_seconds':time.perf_counter()-started,'completed_commands':commands})
  raise
