#!/usr/bin/env python3
"""Explicitly fetch missing public GEO inputs for the code-only distribution.

The Windows launcher never invokes this script. The complete local bundle
already includes all four processed matrices and requires no network access.
"""
from pathlib import Path
import argparse,hashlib,json,os,urllib.request
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--check-only',action='store_true')
 args=parser.parse_args()
 sources=json.loads((ROOT/'config/data_sources.json').read_text(encoding='utf-8'))['sources']
 for item in sources:
  p=ROOT/item['relative_path']
  if p.is_file():
   if sha(p)!=item['sha256']:raise ValueError('Existing input hash mismatch; file preserved: '+str(p))
   print('VERIFIED '+p.name);continue
  if args.check_only:raise FileNotFoundError('Missing '+str(p))
  p.parent.mkdir(parents=True,exist_ok=True)
  temporary=p.with_name(p.name+'.download')
  if temporary.exists():raise FileExistsError('Previous partial response preserved: '+str(temporary))
  request=urllib.request.Request(item['source_url'],headers={'User-Agent':'PlantStressReproducibility/1.0'})
  print('DOWNLOAD '+item['source_url'])
  with urllib.request.urlopen(request,timeout=120) as response,temporary.open('xb') as output:
   for block in iter(lambda:response.read(1048576),b''):output.write(block)
  if temporary.stat().st_size!=item['bytes'] or sha(temporary)!=item['sha256']:
   raise ValueError('Downloaded response failed byte/hash verification; preserved for diagnosis: '+str(temporary))
  os.replace(temporary,p);print('VERIFIED '+p.name)
 print('ALL_INPUTS_VERIFIED')
if __name__=='__main__':main()
