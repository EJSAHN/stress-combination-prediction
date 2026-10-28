#!/usr/bin/env python3
"""Recompute every released scientific analysis using an existing environment."""
import argparse,os,sys
from pathlib import Path
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
 os.environ[name]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from stress_prediction.replay import run

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--output',type=Path,required=True,help='Fresh output directory; previous results are never overwritten.')
 args=p.parse_args()
 run(ROOT,args.output)
