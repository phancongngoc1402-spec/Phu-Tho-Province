from pathlib import Path
import sys
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
import os
P=Path(__file__).resolve().parent;R=Path(os.environ['ESPR_WORK']);sys.path.insert(0,str(P))
from analysis.data import valid_pollutant
from analysis.core import screening_curve
D=pd.read_pickle(R/'work/data.pkl');pols=['COD','TSS','NH4-N'];F=R/'work/figures';F.mkdir(exist_ok=True);sets={(p,y):valid_pollutant(D,p,D.year.eq(y)) for p in pols for y in [2025,2026]}
plt.rcParams.update({'font.size':9,'svg.fonttype':'none'})
def save(f,n):
 for ext in ['png','svg']:f.savefig(F/f'{n}.{ext}',dpi=300,bbox_inches='tight')
 plt.close(f)
f,axs=plt.subplots(3,4,figsize=(11,7.6),sharex=True,sharey=True)
for i,p in enumerate(pols):
 for j,frac in enumerate([.005,.01,.02,.05]):
  a=axs[i,j]
  for y,col,ls in [(2025,'#0072B2','-'),(2026,'#D55E00','--')]:
   c=screening_curve(sets[p,y],load_frac=frac,max_alert=.99);a.plot(c.realized_alert_fraction*100,c.recall*100,color=col,ls=ls,label=str(y))
  a.axhline(90,color='black',ls=':',lw=.7);a.set_title(f'{p}: top {frac*100:g}%');a.set_xlim(0,100);a.set_ylim(0,102);a.grid(alpha=.2)
  if i==2:a.set_xlabel('Flagged intervals (%)')
  if j==0:a.set_ylabel('Recall (%)')
axs[0,0].legend(frameon=False);f.tight_layout();save(f,'Fig1')
m=pd.read_csv(R/'work/new/monthly.csv');f,axs=plt.subplots(2,1,figsize=(9.3,5),sharex=True)
for a,k,lab in [(axs[0],'Recall','Recall (%)'),(axs[1],'Alert','Flagged intervals (%)')]:
 arr=np.array([[m.loc[(m.Channel==p)&(m.Month==f'2026-{mo:02}'),k].iloc[0] for mo in range(1,8)] for p in pols]);im=a.imshow(arr,vmin=0,vmax=100,cmap='cividis',aspect='auto');a.set_yticks(range(3),pols);a.set_title(lab)
 for i in range(3):
  for j in range(7):a.text(j,i,f'{arr[i,j]:.1f}',ha='center',va='center',color='white' if arr[i,j]<50 else 'black')
 f.colorbar(im,ax=a,label=lab,fraction=.03,pad=.02)
axs[1].set_xticks(range(7),['Jan','Feb','Mar','Apr','May','Jun','Jul']);axs[1].set_xlabel('2026');f.tight_layout();save(f,'Fig2')
f,axs=plt.subplots(3,2,figsize=(10,8.5))
for i,p in enumerate(pols):
 po=pd.concat([sets[p,2025],sets[p,2026]]);bounds=np.percentile(po.load_rate,[95,99]);xs=np.geomspace(po.C.min(),po.C.max(),300)
 for j,y in enumerate([2025,2026]):
  a=axs[i,j];dt=sets[p,y];hb=a.hexbin(dt.C,dt.Q,gridsize=45,mincnt=1,xscale='log',yscale='log',norm=LogNorm(vmin=1,vmax=2000),cmap='cividis')
  for b,ls in zip(bounds,['--',':']):a.plot(xs,1000*b/xs,ls=ls,color='black',lw=.8)
  for b,ls in zip(np.percentile(dt.load_rate,[95,99]),['-.','-']):a.plot(xs,1000*b/xs,ls=ls,color='#D55E00',lw=.8)
  a.set_xlim(po.C.min(),po.C.max());a.set_ylim(po.Q.min(),po.Q.max());a.set_title(f'{p}: {y}');a.set_xlabel('C (mg N/L)' if p=='NH4-N' else 'C (mg/L)');a.set_ylabel('Q (m³/h)');f.colorbar(hb,ax=a,label='Count per hexagon',fraction=.04,pad=.02)
f.tight_layout();save(f,'Fig3')
def diagram(n,boxes,arrows):
 f,a=plt.subplots(figsize=(10,5.7));a.set_xlim(0,10);a.set_ylim(0,6);a.axis('off')
 for x,y,w,h,txt,col in boxes:
  a.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.04',facecolor=col,edgecolor='#243746'));a.text(x+w/2,y+h/2,txt,ha='center',va='center',fontsize=9)
 for x,y,u,v in arrows:a.add_patch(FancyArrowPatch((x,y),(u,v),arrowstyle='-|>',mutation_scale=12,color='#243746'))
 f.tight_layout();save(f,n)
diagram('Fig4',[(2.7,4.9,4.6,.8,'Quality-controlled five-minute C and Q\nTimestamps, status and coverage audit','#E5EEF4'),(.2,2.5,2.8,1.4,'CONCENTRATION\nPermit/process threshold\nSeparate response and audit','#E5EEF4'),(3.6,2.5,2.8,1.4,'MASS PRIORITIZATION\nC × Q; interval mass\nLocally approved endpoint','#F9EDDC'),(7,2.5,2.8,1.4,'REVALIDATION\nDistribution + coverage\nRecall, precision, burden','#E7F2EC'),(2,.4,6,1.2,'Human review → process/sensor diagnosis\nDocumented recalibration policy → lock candidate rule\nIndependent future validation; preserve compliance stream','#F1F1F1')],[(5,4.9,1.6,3.9),(5,4.9,5,3.9),(5,4.9,8.4,3.9),(1.6,2.5,5,1.6),(5,2.5,5,1.6),(8.4,2.5,5,1.6)])
diagram('FigS4',[(.2,4.8,4.2,.9,'Supplied analytical workbook\n148,906 unique timestamps','#E5EEF4'),(5.2,4.8,4.6,.9,'Parse, check uniqueness, sort\nConstruct full five-minute calendar','#E5EEF4'),(.2,2.8,4.2,1.1,'Parameter status → positive chemistry\nMutually exclusive exclusion audit\nTables S15–S16','#E5EEF4'),(5.2,2.8,4.6,1.1,'Retain Q=0 for state audit\nQ>0 required for primary ranking\nMissing slots never converted to zero','#F9EDDC'),(.2,.6,4.2,1.1,'Eligible 2025 development\nMass target → rank-screening grid\nAbsolute concentration threshold','#E5EEF4'),(5.2,.6,4.6,1.1,'Lock rule → Jan–Jul 2026\nRelative and fixed historical targets\nDependence-aware resampling','#E7F2EC')],[(4.4,5.25,5.2,5.25),(7.5,4.8,2.3,3.9),(4.4,3.35,5.2,3.35),(7.5,2.8,2.3,1.7),(4.4,1.15,5.2,1.15)])
print('Saved new figures')
