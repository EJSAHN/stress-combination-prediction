#!/usr/bin/env python3
"""Render publication figures from executed analysis tables, without refitting."""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.ticker import PercentFormatter

WIDTH=180/25.4
COMPS=('second_order_factorial','pair_trained_factor_OLS','within_V3')
SHORT={'second_order_factorial':'Second-order\nmodel','pair_trained_factor_OLS':'Fitted\nsingles','within_V3':'Single-response\nspan'}
COLORS={'second_order_factorial':'#BB5566','pair_trained_factor_OLS':'#228833','within_V3':'#4477AA'}
GROUPS=('WT','RAP2.3','Car','Cleo')
GN={'WT':'Col-0','RAP2.3':r'$\it{rap2.3}$','Car':'Carrizo','Cleo':'Cleopatra'}
GC={'WT':'#4477AA','RAP2.3':'#66CCEE','Car':'#228833','Cleo':'#CCBB44'}
MODEL_LABELS={'control':'Control','full_pair_centroid':'Pair\ncentroid','within_only':'Single\nspan','within_V3':'Single\nspan','single_additive':'Additive\nsingles','second_order_factorial':'Second\norder','single_mean':'Mean\nsingles','pair_trained_scalar':'Scaled\nsingles','pair_trained_factor_OLS':'Fitted\nsingles'}

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

class Figures:
 def __init__(self,results,output):
  self.root=results/'tables';self.output=output;self.cache={};self.inputs={};self.records={};self.outputs=[]
  output.mkdir(parents=True,exist_ok=True);(output/'source_data').mkdir(exist_ok=True)
 def read(self,rel):
  if rel not in self.cache:
   path=self.root/rel
   d=pd.read_csv(path);d['_row']=np.arange(len(d));self.cache[rel]=d
   self.inputs[rel]={'sha256':sha(path),'rows':len(d)}
  return self.cache[rel]
 def rec(self,fig,panel,series,metric,x,y,rel,rows,fields,operation='identity',**extra):
  rowlist=[int(v) for v in (rows if isinstance(rows,(list,tuple,np.ndarray)) else [rows])]
  v={'figure':fig,'panel':panel,'series':series,'metric':metric,'x':x,'y':y,
   'source_file':rel,'source_rows':';'.join(map(str,rowlist)),'source_fields':fields,
   'operation':operation,'value_axis':'y','value':y,**extra}
  self.records.setdefault(fig,[]).append(v)
 def save(self,fig,number,title):
  name='Figure_'+str(number);fig.canvas.draw()
  for ax in fig.axes:
   for collection in ax.collections:collection.set_rasterized(False)
  pdf=self.output/(name+'.pdf');png=self.output/(name+'.png')
  pdf_tmp=pdf.with_suffix('.pdf.tmp');png_tmp=png.with_suffix('.png.tmp')
  fig.savefig(pdf_tmp,format='pdf',metadata={'Title':title,'Author':'','CreationDate':None,'ModDate':None},facecolor='white')
  fig.savefig(png_tmp,format='png',dpi=300,facecolor='white',metadata={'Title':title})
  from PIL import Image
  with Image.open(png_tmp) as check:check.load()
  pdf_tmp.replace(pdf);png_tmp.replace(png)
  src=self.output/'source_data'/(name+'_source.csv')
  data=pd.DataFrame(self.records.get(name,[]));data.to_csv(src,index=False,float_format='%.17g')
  assert len(data)>0
  from PIL import Image
  with Image.open(png) as im:im.load();resolution=im.size;dpi=im.info.get('dpi')
  assert resolution[0]>=1080 and dpi and min(dpi)>=299.9
  raw=pdf.read_bytes();assert raw.startswith(b'%PDF') and b'/Subtype /Image' not in raw
  self.outputs.append({'figure':name,'title':title,'pdf':pdf.name,'pdf_sha256':sha(pdf),'png':png.name,
   'png_sha256':sha(png),'source_csv':str(src.relative_to(self.output)),'source_sha256':sha(src),
   'source_rows':len(data),'png_pixels':resolution,'png_dpi':dpi,'width_mm':180,'vector_pdf_no_image_objects':True})
  plt.close(fig)

def one(d,**kwargs):
 m=np.ones(len(d),bool)
 for k,v in kwargs.items():m &= d[k].eq(v).to_numpy()
 rows=d.loc[m];assert len(rows)==1,(kwargs,len(rows))
 return rows.iloc[0]

def style(ax,title=None,ylabel=None):
 ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=8,width=.7,length=3)
 ax.grid(axis='y',alpha=.16,lw=.5)
 if title:ax.set_title(title,loc='left',fontsize=9,fontweight='bold',pad=8)
 if ylabel:ax.set_ylabel(ylabel,fontsize=8.5)

def configure():
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.labelsize':8.5,
  'axes.titlesize':9,'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,
  'axes.linewidth':.7,'lines.linewidth':1.1,'pdf.fonttype':42,'ps.fonttype':42,
  'mathtext.fontset':'dejavusans','savefig.transparent':False})

def primary_models(f):
 e=f.read('external_prediction/model_summary.csv');c=f.read('citrus_prediction/model_summary.csv')
 return e[e.threshold_fpkm.eq(1.0)],c[c.threshold_fpkm.eq(1.0)]

def figure1(f):
 name='Figure_1';e,c=primary_models(f)
 fig=plt.figure(figsize=(WIDTH,6.25));gs=fig.add_gridspec(2,2,height_ratios=[1.05,1],hspace=.55,wspace=.65)
 ax=fig.add_subplot(gs[0,:]);ax.set_xlim(0,12);ax.set_ylim(0,4.6);ax.axis('off')
 ax.set_title('A  Lower-order observations predict the triple condition',loc='left',fontsize=9,fontweight='bold')
 states=[('Control',(0,0,0)),('A',(1,0,0)),('B',(0,1,0)),('C',(0,0,1)),('A+B',(1,1,0)),('A+C',(1,0,1)),('B+C',(0,1,1)),('A+B+C',(1,1,1))]
 for j,l in enumerate(['A','B','C']):ax.text(1.55+j*.48,4.25,l,ha='center',fontsize=8)
 for i,(label,bits) in enumerate(states):
  y=3.86-i*.43;ax.text(.02,y+.13,label,va='center',fontsize=8)
  for j,b in enumerate(bits):
   ax.add_patch(patches.Rectangle((1.37+j*.48,y),.36,.29,fc=('#EE7733' if i==7 else '#4477AA') if b else '#EEEEEE',ec='white',lw=.5))
 ax.annotate('',xy=(4.0,2.6),xytext=(3.05,2.6),arrowprops={'arrowstyle':'->','lw':1,'color':'#555555'})
 ax.text(4.25,3.83,'Pair-centroid prediction',fontsize=9,fontweight='bold')
 ax.text(4.25,3.12,r'$\hat{y}_{ABC}=(y_{AB}+y_{AC}+y_{BC})/3$',fontsize=13)
 ax.text(4.25,2.55,'Means on the log2(1 + FPKM) scale',fontsize=8.5)
 ax.text(4.25,1.83,r'$\it{Arabidopsis\ thaliana}$: Col-0 and $\it{rap2.3}$',fontsize=8.5)
 ax.text(4.25,1.34,'Citrus: Carrizo and Cleopatra',fontsize=8.5)
 ax.text(4.25,.66,'Three RNA pools per condition and genotype',fontsize=8.5)
 ax.text(4.25,.18,'Triple values excluded from predictor construction',fontsize=8.5)
 for g in GROUPS:
  d,rel,col=(e,'external_prediction/model_summary.csv','mse') if g in GROUPS[:2] else (c,'citrus_prediction/model_summary.csv','mean_library_mse')
  r=one(d,genotype=g,model='control')
  f.rec(name,'A',g,'filtered_features',0,0,rel,r._row,'genes',n=3,n_unit='target RNA pools per genotype',value_axis='annotation',value=int(r.genes))
 ax=fig.add_subplot(gs[1,0]);models=['full_pair_centroid','within_V3','single_additive','second_order_factorial'];values=np.empty((4,4))
 for i,g in enumerate(GROUPS):
  d,rel,col=(e,'external_prediction/model_summary.csv','mse') if i<2 else (c,'citrus_prediction/model_summary.csv','mean_library_mse')
  base=one(d,genotype=g,model='control')
  for j,m in enumerate(models):
   key='within_only' if i<2 and m=='within_V3' else m;r=one(d,genotype=g,model=key);v=r[col]/base[col];values[i,j]=v
   f.rec(name,'B',g,m,j,i,rel,[r._row,base._row],col,'numerator_row / control_row',n=3,n_unit='target RNA pools per genotype',denominator=float(base[col]),value_axis='color',value=v)
 im=ax.pcolormesh(np.arange(5)-.5,np.arange(5)-.5,values,cmap='Blues',vmin=0,vmax=max(1,float(values.max())),shading='flat',rasterized=False,edgecolors='white',linewidth=.3)
 ax.set_ylim(3.5,-.5);ax.set_xlim(-.5,3.5)
 for i in range(4):
  for j in range(4):ax.text(j,i,f'{values[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if values[i,j]>im.norm.vmax*.62 else '#222222')
 ax.set_yticks(range(4),[GN[g] for g in GROUPS]);ax.set_xticks(range(4),[MODEL_LABELS[m] for m in models]);ax.tick_params(length=0,pad=5)
 ax.set_title('B  Primary model comparison',loc='left',fontsize=9,fontweight='bold',pad=8)
 cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.08,pad=.27,aspect=22);cb.set_label('MSE / control MSE',fontsize=8);cb.ax.tick_params(labelsize=8)
 ax=fig.add_subplot(gs[1,1]);extra=['full_pair_centroid','single_mean','pair_trained_scalar','pair_trained_factor_OLS','second_order_factorial']
 for g,offset,marker in [('Car',-.12,'o'),('Cleo',.12,'s')]:
  base=one(c,genotype=g,model='control')
  for j,m in enumerate(extra):
   r=one(c,genotype=g,model=m);v=r.mean_library_mse/base.mean_library_mse
   ax.plot(v,j+offset,marker=marker,color=GC[g],ms=4,label=GN[g] if j==0 else None)
   f.rec(name,'C',g,m,v,j+offset,'citrus_prediction/model_summary.csv',[r._row,base._row],'mean_library_mse','numerator_row / control_row',n=3,n_unit='target RNA pools per genotype',denominator=float(base.mean_library_mse),value_axis='x',value=v)
 ax.set_yticks(range(5),['Pair centroid','Mean singles','Scaled singles','Fitted singles','Second order']);ax.invert_yaxis();ax.set_xlim(0,1.03)
 ax.axvline(1,color='#888888',ls=':',lw=.8);ax.set_xlabel('MSE / control MSE')
 style(ax,'C  Additional citrus baselines');ax.grid(False);ax.grid(axis='x',alpha=.16,lw=.5);ax.legend(frameon=False,loc='upper right',bbox_to_anchor=(1.05,1.3),ncol=2,handletextpad=.3,columnspacing=.7)
 fig.subplots_adjust(left=.09,right=.97,top=.94,bottom=.12)
 f.save(fig,1,'Lower-order observations and held-triple prediction')

def figure2(f):
 name='Figure_2';fig,axs=plt.subplots(2,2,figsize=(WIDTH,6.2))
 for col,(rel,genos,msecol,title) in enumerate([
  ('external_span/model_summary.csv',GROUPS[:2],'mse',r'$\it{A.\ thaliana}$'),
  ('citrus_prediction/model_summary.csv',GROUPS[2:],'mean_library_mse','Citrus')]):
  d=f.read(rel);d=d[d.threshold_fpkm.eq(1)];ax=axs[0,col]
  for g,marker in zip(genos,['o','s']):
   base=one(d,genotype=g,model='full_pair_centroid');vals=[]
   for j,m in enumerate(['within_V3','within_V9','within_V11']):
    r=one(d,genotype=g,model=m);v=float(r[msecol]-base[msecol]);vals.append(v)
    f.rec(name,'AB'[col],g,'outside_span_gain',j,v,rel,[r._row,base._row],msecol,'within_span_MSE - centroid_MSE',n=3,n_unit='target RNA pools per genotype',span=m.removeprefix('within_'))
   ax.plot(range(3),vals,marker=marker,ms=4,color=GC[g],label=GN[g])
   ax.annotate(f'{100*vals[-1]/vals[0]:.1f}% retained',(2,vals[-1]),xytext=(-4,8 if g==genos[0] else -13),textcoords='offset points',ha='right',fontsize=8,color=GC[g])
   first=one(d,genotype=g,model='within_V3');last=one(d,genotype=g,model='within_V11')
   f.rec(name,'AB'[col],g,'retained_outside_span_fraction',2,vals[-1],rel,[last._row,base._row,first._row],msecol,'(V11_MSE - centroid_MSE)/(V3_MSE - centroid_MSE)',n=3,n_unit='target RNA pools per genotype',value_axis='annotation',value=vals[-1]/vals[0],display='percent, one decimal')
  ax.set_xticks(range(3),['Single\nmeans','Single\npools','Single +\ncontrol pools']);ax.set_xlim(-.4,2.4);ax.set_ylim(bottom=0)
  style(ax,'AB'[col]+'  '+title,'MSE gain outside observed span');ax.legend(frameon=False,loc='best',handlelength=1.5)
 sims=f.read('noise/simulation_summary.csv');draws=f.read('noise/simulation_draws.csv.gz')
 for j,g in enumerate(GROUPS[:2]):
  ax=axs[1,j];panel='CD'[j]
  r=one(sims,genotype=g,latent_null='joint_rank3',**{'lambda':0.,'noise_sd_scale':1.,'score':'refitted_basis_estimated_control'})
  subset=draws[draws.genotype.eq(g)&draws.latent_null.eq('joint_rank3')&draws['lambda'].eq(0)&draws.noise_sd_scale.eq(1)]
  vals=subset.refitted_basis_estimated_control.to_numpy();assert len(vals)==1000
  edges=np.linspace(float(vals.min()),float(vals.max()),31);counts,_=np.histogram(vals,bins=edges)
  ax.stairs(counts,edges,fill=True,color='#BBBBBB',edgecolor='#888888',lw=.5)
  for k,count in enumerate(counts):
   mask=(vals>=edges[k])&((vals<edges[k+1]) if k<len(counts)-1 else (vals<=edges[k+1]))
   assert mask.sum()==count
   f.rec(name,panel,g,'conditional_draw_count',(edges[k]+edges[k+1])/2,int(count),'noise/simulation_draws.csv.gz',subset.loc[mask,'_row'].to_numpy(),'refitted_basis_estimated_control','histogram count; final bin includes upper edge',lower=edges[k],upper=edges[k+1],n=1000,n_unit='conditional simulation draws')
  ax.axvline(r.observed_statistic,color='#BB5566',ls='--',lw=1.1,label='Observed')
  ax.axvline(r.q95,color='#333333',ls=':',lw=1,label='Conditional 95th percentile')
  ax.text(.04,.94,f'{int(r.n_at_least_observed)} / 1,000\n≥ observed',transform=ax.transAxes,va='top',fontsize=8)
  for field in ['observed_statistic','q95','n_at_least_observed']:
   f.rec(name,panel,g,field,float(r[field]),0,'noise/simulation_summary.csv',r._row,field,n=1000,n_unit='conditional simulation draws',value_axis='annotation' if field=='n_at_least_observed' else 'x',value=float(r[field]))
  ax.set_ylim(0,float(counts.max())*1.25);ax.set_xlim(min(0,edges[0]),max(edges[-1],r.q95,r.observed_statistic)*1.08)
  ax.set_xlabel('Outside-span MSE gain');style(ax,panel+'  Conditional noise: '+GN[g],'Number of draws')
  if j==0:ax.legend(frameon=False,loc='upper right',bbox_to_anchor=(1,1.33),fontsize=8,handlelength=1.5)
 fig.subplots_adjust(left=.11,right=.98,top=.95,bottom=.10,wspace=.42,hspace=.67)
 f.save(fig,2,'Observed-span sensitivity and conditional noise uncertainty')

def figure3(f):
 name='Figure_3';fig,axs=plt.subplots(2,2,figsize=(WIDTH,6.3));rem=f.read('citrus_attribution/concentration_removal.csv');genes=f.read('citrus_attribution/gene_summary.csv')
 for j,g in enumerate(GROUPS[2:]):
  ax=axs[0,j]
  for comp in COMPS:
   rr=rem[rem.genotype.eq(g)&rem.threshold_fpkm.eq(1)&rem.comparator.eq(comp)].sort_values('requested_fraction')
   r0=one(genes,genotype=g,threshold_fpkm=1.,comparator=comp)
   x=np.r_[0,rr.requested_fraction.to_numpy()];y=np.r_[1,rr.retained_over_original_net_mass.to_numpy()]
   ax.plot(x,y,marker='o',ms=3,color=COLORS[comp],label=SHORT[comp].replace('\n',' '))
   f.rec(name,'AB'[j],comp,'net_gain_retained',0,1,'citrus_attribution/gene_summary.csv',r0._row,'net_mass','net_mass / itself',genotype=g,n=int(r0.n_features),n_unit='features in original panel')
   for _,r in rr.iterrows():f.rec(name,'AB'[j],comp,'net_gain_retained',r.requested_fraction,r.retained_over_original_net_mass,'citrus_attribution/concentration_removal.csv',r._row,'requested_fraction;retained_over_original_net_mass',genotype=g,n=int(r.remaining_count),n_unit='retained features')
  ax.axhline(0,color='#555555',ls='--',lw=.7);ax.set_xlim(-.003,.103);ax.set_ylim(-.48,1.08)
  ax.set_xticks([0,.01,.05,.10]);ax.xaxis.set_major_formatter(PercentFormatter(1,decimals=0));ax.yaxis.set_major_formatter(PercentFormatter(1,decimals=0))
  ax.set_xlabel('Highest-gain features removed');style(ax,'AB'[j]+'  '+GN[g],'Original net gain retained')
  if j==0:
   handles,labels=ax.get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.52,.995),ncol=3,frameon=False,fontsize=8,handlelength=1.3,columnspacing=.7)
 ax=axs[1,0]
 for g,offset,marker in [('Car',-.09,'o'),('Cleo',.09,'s')]:
  vals=[]
  for j,comp in enumerate(COMPS):
   r=one(genes,genotype=g,threshold_fpkm=1.,comparator=comp);vals.append(r.fraction_mean_positive)
   f.rec(name,'C',g,'fraction_mean_positive',j+offset,r.fraction_mean_positive,'citrus_attribution/gene_summary.csv',r._row,'fraction_mean_positive',comparator=comp,n=int(r.n_features),n_unit='features')
  ax.plot(np.arange(3)+offset,vals,ls='none',marker=marker,ms=4,color=GC[g],label=GN[g])
 ax.set_xticks(range(3),[SHORT[k] for k in COMPS]);ax.set_xlim(-.4,2.4);ax.set_ylim(0,1)
 ax.axhline(.5,color='#999999',lw=.7,ls=':');ax.yaxis.set_major_formatter(PercentFormatter(1,decimals=0))
 style(ax,'C  Distribution of improvement','Features with positive mean gain');ax.legend(frameon=False,loc='upper left',ncol=2,handletextpad=.3,columnspacing=.8)
 ax=axs[1,1];ora=f.read('citrus_attribution/direction_oracles.csv');labels=[]
 for i,(g,comp) in enumerate([(g,c) for g in GROUPS[2:] for c in COMPS[:2]]):
  r=one(ora,genotype=g,threshold_fpkm=1.,model=comp);vals=[]
  for fld in ['original_mean_mse','target_fitted_mean_mse']:
   v=r[fld]/r.original_centroid_mean_mse;vals.append(v)
   f.rec(name,'D',g,fld,v,i,'citrus_attribution/direction_oracles.csv',r._row,fld+';original_centroid_mean_mse','model_MSE / original_centroid_MSE',comparator=comp,n=3,n_unit='target RNA pools per genotype',value_axis='x',value=v)
  ax.plot(vals,[i,i],color='#BBBBBB',lw=.9);ax.plot(vals[0],i,'o',ms=4,color='#777777',label='Original' if i==0 else None);ax.plot(vals[1],i,'D',ms=4,color='#4477AA',label='Target-fitted' if i==0 else None)
  labels.append(GN[g]+'\n'+('Second order' if comp==COMPS[0] else 'Fitted singles'))
 ax.axvline(1,color='#BB5566',ls='--',lw=.8);ax.set_yticks(range(4),labels);ax.invert_yaxis();ax.set_xlim(left=0);ax.set_xlabel('MSE / original pair-centroid MSE')
 style(ax,'D  Outcome-fitted rescaling');ax.grid(False);ax.grid(axis='x',alpha=.16,lw=.5);ax.legend(frameon=False,loc='lower left',fontsize=8,handletextpad=.4,labelspacing=.3)
 fig.subplots_adjust(left=.10,right=.98,top=.89,bottom=.12,wspace=.66,hspace=.64)
 f.save(fig,3,'Gene-wise prediction gains and concentration')

def figure4(f):
 name='Figure_4';fig,axs=plt.subplots(2,2,figsize=(WIDTH,6.1));mom=f.read('attribution_repeatability/part_A_moments.csv');cross=f.read('attribution_repeatability/part_A_cross_genotype.csv');summ=f.read('attribution_repeatability/part_B_contrast_summary.csv')
 ax=axs[0,0]
 for g,offset,marker in [('Car',-.08,'o'),('Cleo',.08,'s')]:
  vals=[]
  for j,comp in enumerate(COMPS):
   r=one(mom,genotype=g,comparator=comp);vals.append(r.conditional_reliability)
   f.rec(name,'A',g,'conditional_reliability',j+offset,r.conditional_reliability,'attribution_repeatability/part_A_moments.csv',r._row,'conditional_reliability',n=3,n_unit='target RNA pools per genotype',comparator=comp)
  ax.plot(np.arange(3)+offset,vals,ls='none',marker=marker,ms=4,color=GC[g],label=GN[g])
 ax.set_ylim(0,1.05);ax.axhline(1,color='#999999',ls=':',lw=.7);style(ax,'A  Fixed-predictor repeatability','Conditional stable-moment fraction');ax.legend(frameon=False,loc='lower left',handletextpad=.4)
 ax=axs[0,1]
 for fld,off,col,marker,label in [('raw_Pearson',-.07,'#777777','o','Observed'),('adjusted_Pearson',.07,'#4477AA','D','Moment-adjusted')]:
  vals=[]
  for j,comp in enumerate(COMPS):
   r=one(cross,comparator=comp);vals.append(r[fld]);f.rec(name,'B',fld,comp,j+off,r[fld],'attribution_repeatability/part_A_cross_genotype.csv',r._row,fld,n=3,n_unit='target RNA pools per genotype',n_features=15466)
  ax.plot(np.arange(3)+off,vals,ls='none',marker=marker,ms=4,color=col,label=label)
 for j,comp in enumerate(COMPS):
  r=one(cross,comparator=comp);ax.plot([j-.07,j+.07],[r.raw_Pearson,r.adjusted_Pearson],color='#BBBBBB',lw=.8,zorder=0)
 ax.set_ylim(bottom=0);style(ax,'B  Between-genotype gains','Pearson or adjusted moment ratio');ax.legend(frameon=False,loc='upper left',handletextpad=.4)
 for ax,panel,metric,title,ylabel in [(axs[1,0],'C','Spearman','C  Disjoint refits: ranks','Spearman correlation'),(axs[1,1],'D','top10pct_Jaccard','D  Disjoint refits: top 10%','Top-set Jaccard overlap')]:
  for component,offset,col,marker,label in [('within_mean',-.09,'#228833','o','Within genotype'),('cross_mean',.09,'#BB5566','s','Between genotypes')]:
   values=[];lower=[];upper=[]
   for j,comp in enumerate(COMPS):
    r=one(summ,comparator=comp,metric=metric,component=component);assert r.n_finite==256
    values.append(r['median']);lower.append(r['median']-r.q025);upper.append(r.q975-r['median'])
    f.rec(name,panel,component,metric,j+offset,r['median'],'attribution_repeatability/part_B_contrast_summary.csv',r._row,'median;q025;q975',lower=r.q025,upper=r.q975,n=256,n_unit='overlapping partitions',comparator=comp,range_type='overlapping-partition range; not CI')
   ax.errorbar(np.arange(3)+offset,values,yerr=[lower,upper],fmt=marker,color=col,ms=4,capsize=2,lw=1,label=label)
  ax.set_ylim(bottom=0);style(ax,title,ylabel);ax.legend(frameon=False,loc='upper left',handlelength=1.2,handletextpad=.4,labelspacing=.3)
 for ax in axs.flat:ax.set_xticks(range(3),[SHORT[k] for k in COMPS]);ax.set_xlim(-.4,2.4)
 fig.subplots_adjust(left=.11,right=.98,top=.94,bottom=.12,wspace=.40,hspace=.62)
 f.save(fig,4,'Within- and between-genotype repeatability of gene-gain profiles')

def supplement1(f):
 name='Figure_S1';fig,axs=plt.subplots(1,2,figsize=(WIDTH,4.0));rel='discovery_prediction/target_model_comparison.csv';d=f.read(rel)
 d=d[d.panel.eq('all_filtered')&d.component.eq('expression')];base=one(d,model='control_mean')
 models=['control_mean','single_effect_sum','pairwise','all_triples_mean','quiet_triples_mean','quiet_residual_rank2']
 labels=['Control','Additive singles','Second order','All triple donors','Selected triple mean','Selected residual']
 ax=axs[0]
 for i,m in enumerate(models):
  r=one(d,model=m);v=float(np.log10(r.mse_mean/base.mse_mean));ax.plot(v,i,'o',ms=4,color='#4477AA')
  f.rec(name,'A',m,'log10_MSE_ratio',v,i,rel,[r._row,base._row],'mse_mean','log10(numerator_row / control_row)',n=3,n_unit='target RNA pools per genotype',value_axis='x',value=v)
 ax.axvline(0,color='#888888',ls=':',lw=.8);ax.set_yticks(range(len(models)),labels);ax.invert_yaxis();ax.set_xlabel('log10(MSE / control MSE)');style(ax,'A  Exploratory four-factor target');ax.grid(False);ax.grid(axis='x',alpha=.16,lw=.5)
 rel='discovery_selection/target_model_errors.csv';d=f.read(rel);d=d[d['rank'].eq(2)];targets=d.target.unique().tolist()
 models=['all_donors_mean','common_low2','raw_low2','residual_low2','common_low2_residual_only'];labels=['All\ndonors','Low\ncommon','Low\nraw','Low\nresidual','Low common\nresidual only']
 v=np.empty((len(targets),len(models)))
 for i,t in enumerate(targets):
  base=one(d,target=t,model='control_mean')
  for j,m in enumerate(models):
   r=one(d,target=t,model=m);v[i,j]=np.log10(r.MSE/base.MSE)
   f.rec(name,'B',t,m,j,i,rel,[r._row,base._row],'MSE','log10(numerator_row / control_row)',value_axis='color',value=v[i,j])
 ax=axs[1];lim=max(abs(v.min()),abs(v.max()));im=ax.pcolormesh(np.arange(v.shape[1]+1)-.5,np.arange(v.shape[0]+1)-.5,v,cmap='RdBu_r',vmin=-lim,vmax=lim,shading='flat',rasterized=False,edgecolors='white',linewidth=.3)
 ax.set_ylim(v.shape[0]-.5,-.5);ax.set_xlim(-.5,v.shape[1]-.5)
 for i in range(len(targets)):
  for j in range(len(models)):ax.text(j,i,f'{v[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if abs(v[i,j])>.65*lim else '#222222')
 ax.set_xticks(range(len(models)),labels,rotation=45,ha='right');ax.set_yticks(range(len(targets)),[t.replace('-','+') for t in targets]);ax.tick_params(length=0)
 ax.set_title('B  Donor-selection sensitivity',loc='left',fontsize=9,fontweight='bold',pad=8)
 cb=fig.colorbar(im,ax=ax,orientation='horizontal',pad=.37,fraction=.06,aspect=20);cb.set_label('log10(MSE / control MSE)',fontsize=8)
 fig.subplots_adjust(left=.20,right=.98,top=.88,bottom=.22,wspace=.82)
 f.save(fig,'S1','Exploratory discovery and donor-selection failure')

def supplement2(f):
 name='Figure_S2';fig,axs=plt.subplots(3,2,figsize=(WIDTH,8.2));models=['full_pair_centroid','within_V3','single_additive','second_order_factorial']
 for j,g in enumerate(GROUPS):
  ax=axs.flat[j];rel='external_prediction/library_errors.csv' if j<2 else 'citrus_prediction/library_errors.csv';d=f.read(rel);d=d[d.genotype.eq(g)&d.threshold_fpkm.eq(1)]
  for k,m in enumerate(models):
   key='within_only' if j<2 and m=='within_V3' else m;rr=d[d.model.eq(key)]
   for q,(_,r) in enumerate(rr.iterrows()):
    base=one(d,model='control',target_library=r.target_library);v=r.mse/base.mse;xp=k+(q-1)*.10
    ax.plot(xp,v,'o',ms=3,color=GC[g]);f.rec(name,'ABCD'[j],g,m,xp,v,rel,[r._row,base._row],'mse','library_model_MSE / same_library_control_MSE',n=1,n_unit='target RNA pool',pool=r.target_library)
  ax.axhline(1,color='#888888',ls=':',lw=.7);ax.set_ylim(bottom=0);ax.set_xlim(-.35,3.35);ax.set_xticks(range(4),[MODEL_LABELS[m] for m in models]);style(ax,'ABCD'[j]+'  '+GN[g],'Target-pool MSE / control MSE')
 ax=axs[2,0]
 for g in GROUPS:
  rel='external_prediction/model_summary.csv' if g in GROUPS[:2] else 'citrus_prediction/model_summary.csv';col='mse' if g in GROUPS[:2] else 'mean_library_mse';d=f.read(rel);vals=[]
  for t in [.5,1.,2.]:
   r=one(d,genotype=g,threshold_fpkm=t,model='full_pair_centroid');base=one(d,genotype=g,threshold_fpkm=t,model='control');vals.append(r[col]/base[col]);f.rec(name,'E',g,'centroid_MSE_ratio',t,vals[-1],rel,[r._row,base._row],col,'numerator_row / control_row',n=3,n_unit='target RNA pools per genotype')
  ax.plot([.5,1,2],vals,marker='o',ms=3,color=GC[g],label=GN[g])
 ax.set_ylim(bottom=0);ax.set_xticks([.5,1,2]);ax.set_xlabel('Lower-order mean FPKM filter');style(ax,'E  Expression-filter sensitivity','Pair-centroid MSE / control MSE');ax.legend(frameon=False,fontsize=8,ncol=2,handlelength=1.1,columnspacing=.5)
 ax=axs[2,1]
 for g in GROUPS:
  rel='external_span/model_summary.csv' if g in GROUPS[:2] else 'citrus_prediction/model_summary.csv';col='mse' if g in GROUPS[:2] else 'mean_library_mse';d=f.read(rel);vals=[]
  for t in [.5,1.,2.]:
   a=one(d,genotype=g,threshold_fpkm=t,model='within_V3');b=one(d,genotype=g,threshold_fpkm=t,model='within_V11');base=one(d,genotype=g,threshold_fpkm=t,model='full_pair_centroid');v=(b[col]-base[col])/(a[col]-base[col]);vals.append(v)
   f.rec(name,'F',g,'outside_gain_retention',t,v,rel,[b._row,base._row,a._row],col,'(V11_MSE-centroid_MSE)/(V3_MSE-centroid_MSE)',n=3,n_unit='target RNA pools per genotype')
  ax.plot([.5,1,2],vals,marker='o',ms=3,color=GC[g])
 ax.set_ylim(0,1.05);ax.set_xticks([.5,1,2]);ax.set_xlabel('Lower-order mean FPKM filter');ax.yaxis.set_major_formatter(PercentFormatter(1,decimals=0));style(ax,'F  Expanded-span gain retention','Original outside-span gain retained')
 fig.subplots_adjust(left=.11,right=.98,top=.95,bottom=.08,wspace=.40,hspace=.65)
 f.save(fig,'S2','Target-pool and expression-filter sensitivity')

def supplement3(f):
 name='Figure_S3';rel='noise/simulation_summary.csv';d=f.read(rel);d=d[d.score.eq('refitted_basis_estimated_control')]
 fig,axs=plt.subplots(2,2,figsize=(WIDTH,6.0));low=min(0,float(d.q025.min()));high=float(d.q975.max())*1.1
 palettes={0.:'#4477AA',.5:'#228833',1.:'#AA3377'}
 for i,g in enumerate(GROUPS[:2]):
  for j,null in enumerate(['single_mean_span','joint_rank3']):
   ax=axs[i,j];panel='ABCD'[i*2+j]
   for l,off in [(0.,-.06),(.5,0),(1.,.06)]:
    vals=[];lo=[];hi=[]
    for sd in [.5,1.,2.]:
     r=one(d,genotype=g,latent_null=null,**{'lambda':l,'noise_sd_scale':sd})
     vals.append(r.q50);lo.append(r.q50-r.q025);hi.append(r.q975-r.q50)
     f.rec(name,panel,g,'conditional_gain',sd+off,r.q50,rel,r._row,'q50;q025;q975',lower=r.q025,upper=r.q975,n=1000,n_unit='conditional simulation draws',latent_null=null,diagonal_fraction=l,noise_sd_scale=sd,range_type='conditional simulation range; not CI')
    ax.errorbar(np.array([.5,1,2])+off,vals,yerr=[lo,hi],fmt='o',ms=3,capsize=2,lw=.9,color=palettes[l],label={0.:'0 (empirical)',.5:'0.5',1.:'1 (diagonal)'}[l])
   obs=d[d.genotype.eq(g)].iloc[0];ax.axhline(obs.observed_statistic,color='#BB5566',ls='--',lw=1)
   f.rec(name,panel,g,'observed_statistic',0,obs.observed_statistic,rel,obs._row,'observed_statistic',n=3,n_unit='target RNA pools per genotype')
   ax.set_xlim(0,2.3);ax.set_ylim(low,high);ax.set_xticks([.5,1,2]);ax.set_xlabel('Noise standard-deviation multiplier')
   style(ax,panel+'  '+GN[g]+(': single-mean span' if j==0 else ': joint rank 3'),'Outside-span MSE gain')
 handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.53,.997),ncol=3,frameon=False,title='Diagonal covariance fraction',fontsize=8,title_fontsize=8)
 fig.subplots_adjust(left=.12,right=.98,top=.83,bottom=.10,wspace=.37,hspace=.52)
 f.save(fig,'S3','All conditional noise scenarios')

def supplement4(f):
 name='Figure_S4';rel='citrus_attribution/functional_summary.csv';d=f.read(rel);d=d[d.threshold_fpkm.eq(1)]
 go=[('GO:0006457','Protein folding'),('GO:0015979','Photosynthesis'),('GO:0006412','Translation'),('GO:0098869','Oxidant\ndetoxification*'),('GO:0055085','Membrane transport'),('GO:0071554','Cell wall')]
 chosen=d[d.category_id.isin([x[0] for x in go])&d.comparator.isin(COMPS[:2])]
 xmin=min(0,float(chosen.mean_gain.min()))-.08;xmax=max(0,float(chosen.mean_gain.max()))+.25
 fig,axs=plt.subplots(2,2,figsize=(WIDTH,6.1))
 for i,comp in enumerate(COMPS[:2]):
  for j,g in enumerate(GROUPS[2:]):
   ax=axs[i,j];panel='ABCD'[2*i+j]
   for k,(gid,label) in enumerate(go):
    r=one(d,genotype=g,comparator=comp,category_id=gid);small=bool(r.small_group_n_below20)
    ax.plot(r.mean_gain,k,'o',ms=4,color=COLORS[comp],mfc='white' if small else COLORS[comp]);ax.annotate(f'n={int(r.n_features)}',(r.mean_gain,k),xytext=(5,0),textcoords='offset points',va='center',fontsize=8,bbox={'fc':'white','ec':'none','pad':.4,'alpha':.85})
    f.rec(name,panel,g,gid,r.mean_gain,k,rel,r._row,'mean_gain;n_features;small_group_n_below20',n=int(r.n_features),n_unit='features',comparator=comp,value_axis='x',value=r.mean_gain)
   bg=one(d,genotype=g,comparator=comp,category_id='BP_annotated_background');ax.axvline(bg.mean_gain,color='#777777',ls=':',lw=.8)
   f.rec(name,panel,g,'BP_background',bg.mean_gain,0,rel,bg._row,'mean_gain',n=int(bg.n_features),n_unit='features',comparator=comp,value_axis='x',value=bg.mean_gain)
   ax.axvline(0,color='#555555',ls='--',lw=.7);ax.set_yticks(range(6),[t[1] for t in go] if j==0 else ['']*6);ax.invert_yaxis();ax.set_xlim(xmin,xmax);ax.set_xlabel('Mean error reduction')
   style(ax,panel+'  '+GN[g]+(': second order' if i==0 else ': fitted singles'));ax.grid(False);ax.grid(axis='x',alpha=.16,lw=.5);ax.tick_params(axis='y',length=0)
 cov=f.read('citrus_attribution/annotation_coverage.csv');coverage=[]
 for g in GROUPS[2:]:
  r=one(cov,genotype=g,threshold_fpkm=1.);coverage.append(f'{GN[g]} {100*r.any_GO_coverage_of_original_panel:.1f}%')
  f.rec(name,'annotation',g,'GO_coverage',0,0,'citrus_attribution/annotation_coverage.csv',r._row,'any_GO_coverage_of_original_panel',n=int(r.genes),n_unit='features in original panel',value_axis='annotation',value=r.any_GO_coverage_of_original_panel)
 fig.text(.20,.035,'Any-GO coverage: '+', '.join(coverage)+'.\nOpen markers: fewer than 20 features. Dotted lines: BP-annotated background.',fontsize=8)
 fig.subplots_adjust(left=.20,right=.98,top=.94,bottom=.17,wspace=.25,hspace=.46)
 f.save(fig,'S4','Computational functional annotations of gene-wise gains')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 configure();f=Figures(a.results.resolve(),a.output.resolve())
 for draw in [figure1,figure2,figure3,figure4,supplement1,supplement2,supplement3,supplement4]:draw(f)
 report={'status':'PASS','scope':'Numerical source tracing and export-format checks; visual QA is recorded separately.',
  'figure_count':len(f.outputs),'main_figures':4,'supplementary_figures':4,'input_files':f.inputs,'figures':f.outputs,
  'plot_code_sha256':sha(__file__),'zero_values_and_negative_results_retained':True,
  'axes_and_ranges':'No gene-level uncertainty bars. Simulation and partition intervals are explicitly non-CI ranges.',
  'font':'DejaVu Sans, 8–9pt at180mm page width; panel letters9pt; mathematical formula13pt',
  'visual_review_metadata':'Manual review required after first rendering; current release review is supplied in docs/figure_visual_review.json.'}
 (f.output/'figure_validation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
 print(json.dumps({'status':'PASS','figure_count':len(f.outputs),'output':str(f.output)}))
if __name__=='__main__':main()
