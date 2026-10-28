#!/usr/bin/env python3
"""Independent table replay audit; does not import analysis or figure code.

The optional archive root checks references against the untouched research
records in addition to the bundled reference hashes. This is a release audit,
not an independent biological validation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for part in iter(lambda:f.read(1024*1024),b''):
            h.update(part)
    return h.hexdigest()


def numeric_strings(values):
    result=[]
    for v in values:
        if v in ('','NaN','nan','NA'):
            result.append(np.nan)
        else:
            try:result.append(float(v))
            except ValueError:return None
    return np.asarray(result,float)


def compare_csv(actual,expected,atol=1e-10,rtol=1e-9):
    a=pd.read_csv(actual,dtype=str,keep_default_na=False)
    b=pd.read_csv(expected,dtype=str,keep_default_na=False)
    assert list(a.columns)==list(b.columns),(str(actual),'column order/schema')
    assert a.shape==b.shape,(str(actual),'table shape',a.shape,b.shape)
    numeric=0;text=0;maxdiff=0.;fields=[]
    for name in a.columns:
        av=a[name].tolist();bv=b[name].tolist()
        if actual.name=='latent_mean_diagnostics.csv' and name in {
                'latent_single_singular_values','latent_basis_singular_values'}:
            errors=[]
            for left,right in zip(av,bv):
                x,y=json.loads(left),json.loads(right)
                assert isinstance(x,list) and isinstance(y,list) and len(x)==len(y)>0
                assert all(type(v) in (int,float) for v in x+y)
                x,y=np.asarray(x,float),np.asarray(y,float)
                assert np.isfinite(x).all() and np.isfinite(y).all()
                assert np.allclose(x,y,atol=atol,rtol=rtol),(str(actual),name,'numeric vector mismatch')
                numeric+=len(x);errors.append(float(np.max(np.abs(x-y))))
            d=max(errors,default=0.);maxdiff=max(maxdiff,d)
            fields.append({'column':name,'kind':'numeric_json_array','maximum_absolute_difference':d})
            continue
        an=numeric_strings(av);bn=numeric_strings(bv)
        if an is not None and bn is not None:
            equal=np.isclose(an,bn,atol=atol,rtol=rtol,equal_nan=True)
            assert bool(equal.all()),(str(actual),name,'numeric mismatch',np.flatnonzero(~equal)[:8].tolist())
            finite=np.isfinite(an)&np.isfinite(bn)
            d=float(np.max(np.abs(an[finite]-bn[finite]))) if finite.any() else 0.
            maxdiff=max(maxdiff,d);numeric+=len(an)
            fields.append({'column':name,'kind':'numeric','maximum_absolute_difference':d})
        else:
            assert av==bv,(str(actual),name,'text mismatch')
            text+=len(av)
            fields.append({'column':name,'kind':'exact_text'})
    return {'rows':len(a),'columns':len(a.columns),'numeric_cells':numeric,
            'exact_text_cells':text,'maximum_absolute_difference':maxdiff,
            'column_checks':fields}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repository',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--output',type=Path,required=True,help='Freshly executed pipeline output directory')
    p.add_argument('--archive-root',type=Path,help='Optional original research-record parent')
    p.add_argument('--report',type=Path,required=True)
    args=p.parse_args();repo=args.repository.resolve();out=args.output.resolve()
    mapping=repo/'config/table_reference_map.json'
    records=json.loads(mapping.read_text())['tables'];checks=[]
    assert records and len({r['output'] for r in records})==len(records)
    for rec in records:
        reference=repo/rec['reference'];actual=out/rec['output']
        assert digest(reference)==rec['sha256'],('bundled reference changed',rec['reference'])
        expected=reference
        amended=rec.get('correction')=='single_span_identity_diagnostics_v1'
        if amended:
            original=repo/rec['historical_reference']
            assert digest(original)==rec['historical_sha256'],('historical reference changed',str(original))
        if args.archive_root is not None:
            original=args.archive_root/rec['original_relative_path']
            original_hash=rec['historical_sha256'] if amended else rec['sha256']
            assert digest(original)==original_hash,('original reference changed',rec['original_relative_path'])
            if not amended: expected=original
        checked=compare_csv(actual,expected)
        checks.append({'stage':rec['stage'],'output':rec['output'],
            'reference_sha256':rec['sha256'],'output_sha256':digest(actual),
            'documented_diagnostic_correction':amended,**checked})
    report={'status':'PASS_INDEPENDENT_RELEASE_TABLE_AUDIT',
        'completed_utc':datetime.now(timezone.utc).isoformat(),
        'analysis_or_figure_code_imported':False,
        'untouched_original_archive_checked':args.archive_root is not None,
        'comparison_scope':'All cells of all explicitly mapped released CSV tables; not raw-read reprocessing or new independent experimental validation.',
        'amended_diagnostic_tables':sum(x['documented_diagnostic_correction'] for x in checks),
        'absolute_tolerance':1e-10,'relative_tolerance':1e-9,
        'table_count':len(checks),'numeric_cells':sum(v['numeric_cells'] for v in checks),
        'exact_text_cells':sum(v['exact_text_cells'] for v in checks),
        'maximum_absolute_difference':max(v['maximum_absolute_difference'] for v in checks),
        'mapping_sha256':digest(mapping),'audit_code_sha256':digest(__file__),'tables':checks}
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='tables'},indent=2))


if __name__=='__main__':main()
