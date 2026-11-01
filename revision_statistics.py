from pathlib import Path
import sys,json,math
import numpy as np,pandas as pd
from urllib.request import Request,urlopen
from concurrent.futures import ThreadPoolExecutor
import os
P=Path(__file__).resolve().parent;R=Path(os.environ['ESPR_WORK']);sys.path.insert(0,str(P))
from analysis.config import *
from analysis.core import *
from analysis.data import valid_pollutant
from analysis.bootstrap import *
from analysis.bootstrap import _bootstrap_weight_context,_weights_from_day_sample,_allocate_top_multiplicity
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=R/'work/new';F=O/'figures';d=pd.read_pickle(R/'work/data.pkl');pols=list(POLLUTANTS);D={(p,y):valid_pollutant(d,p,d.year.eq(y)) for p in pols for y in (2025,2026)}
# Proper pairwise Pearson ACF at exact calendar lags. FFT correlations of sufficient statistics.
def acf_pairwise(x,maxlag=2016):
 x=np.asarray(x,float);ok=np.isfinite(x);x=np.where(ok,x,0);w=ok.astype(float);nfft=1<<(2*len(x)-1).bit_length()
 fw=np.fft.rfft(w,nfft);fx=np.fft.rfft(x,nfft);fx2=np.fft.rfft(x*x,nfft)
 def cc(a,b):return np.fft.irfft(a*np.conj(b),nfft)[:maxlag+1]
 nn=cc(fw,fw);xy=cc(fx,fx);sx=cc(fx,fw);sy=cc(fw,fx);sx2=cc(fx2,fw);sy2=cc(fw,fx2)
 num=xy-sx*sy/nn;den=np.sqrt(np.maximum(0,sx2-sx*sx/nn)*np.maximum(0,sy2-sy*sy/nn));return np.divide(num,den,out=np.full_like(num,np.nan),where=(nn>2)&(den>0))
acs={};rows=[]
for p in pols:
 for y in (2025,2026):
  dt=D[p,y];cal=pd.date_range(f'{y}-01-01','2025-12-31 23:55' if y==2025 else '2026-07-31 23:55',freq='5min')
  for var in ['C','Q']:
   ac=acf_pairwise(dt.set_index('time')[var].reindex(cal).to_numpy());assert np.nanmax(abs(ac))<1.00001;acs[p,y,var]=ac;ze=np.flatnonzero(ac[1:]<=0);cut=int(ze[0]+1) if len(ze) else len(ac);de=np.flatnonzero(ac[1:]<=np.exp(-1));tau=max(1,1+2*np.sum(ac[1:cut]));rows.append({'Channel':p,'Period':y,'Variable':var,'n':len(dt),'Lag1':ac[1],'ACF24h':ac[288],'DecorrMin':int(de[0]+1)*5 if len(de) else np.nan,'ZeroMin':cut*5 if len(ze) else np.nan,'Tau':tau,'ESSproxy':len(dt)/tau,'Truncated':not bool(len(ze))})
pd.DataFrame(rows).to_csv(O/'acf.csv',index=False)
fig,axs=plt.subplots(2,3,figsize=(10,5),sharex=True,sharey=True)
for col,p in enumerate(pols):
 for row,y in enumerate((2025,2026)):
  ax=axs[row,col]
  for var,clr,ls in [('C','#0072B2','-'),('Q','#D55E00','--')]:ax.plot(np.arange(2017)/12,acs[p,y,var],label=var,color=clr,ls=ls,lw=1)
  ax.axvline(24,color='black',ls=':',lw=.7);ax.axhline(0,color='gray',lw=.6);ax.set_title(f'{p}, {y}');ax.set_ylim(-.3,1.05);ax.grid(alpha=.2)
  if row==1:ax.set_xlabel('Calendar lag (h)')
  if col==0:ax.set_ylabel('Pairwise Pearson ACF')
axs[0,0].legend(frameon=False);fig.tight_layout()
for ext in ['png','svg']:fig.savefig(F/f'FigS5_ACF.{ext}',dpi=350,bbox_inches='tight')
plt.close(fig)
# Primary CIs rederived in a shared stream for internal consistency across all four metrics.
# Base published streams remain archived. Separate frozen-label sensitivity uses same day draws.
rows=[];fixed=[];rng=np.random.default_rng(20261026)
for p in pols:
 dt=D[p,2026];th=alert_requirement(D[p,2025])['threshold'];bd=p99_mass_boundary(D[p,2025]);dc,nd,c,ma,co,order,ge=_bootstrap_weight_context(dt);fl=c>=th;obs=np.zeros(len(dt),bool);obs[order[:math.ceil(.01*len(dt))]]=True;fm=ma>=bd;v=[];f=[];s=[]
 for b in range(1000):
  w=_weights_from_day_sample(rng,dc,nd);n=int(w.sum());m=math.ceil(.01*n);tm=_allocate_top_multiplicity(order,w,m);al=int(w[fl].sum());tp=int(tm[fl].sum());nt=int(w[obs].sum());v.append([al/n,tp/m,tp/al if al else np.nan,np.dot(tm[fl],ma[fl])/np.dot(tm,ma),w[obs&fl].sum()/nt if nt else np.nan]);nf=int(w[fm].sum());tf=int(w[fm&fl].sum());f.append([nf,tf/nf if nf else np.nan,tf/al if al else np.nan,np.dot(w[fm&fl],ma[fm&fl])/np.dot(w[fm],ma[fm]) if nf else np.nan]);s.append(np.count_nonzero((w>0)&fm))
 v=np.array(v);f=np.array(f);s=np.array(s);np.savez_compressed(O/f'bootstrap_samples_{p}.npz',primary=v,fixed=f,distinct_target=s)
 q=np.nanpercentile(v,[2.5,50,95,97.5],axis=0);r={'Channel':p,'OneSidedUpper':q[2,1]*100}
 for j,key in enumerate(['Alert','Recall','Precision','Capture','Frozen']):r[key+'Low']=q[0,j]*100;r[key+'Median']=q[1,j]*100;r[key+'High']=q[3,j]*100
 rows.append(r)
 qf=np.nanpercentile(f,[2.5,50,97.5],axis=0);rr={'Channel':p,'Empty':int(np.isnan(f[:,1]).sum()),'UniqueMedian':float(np.median(s)),'UniqueLow':np.percentile(s,2.5),'UniqueHigh':np.percentile(s,97.5)}
 for j,key in enumerate(['Target','Recall','Precision','Capture']):
  factor=1 if key=='Target' else 100
  rr[key+'Low']=qf[0,j]*factor;rr[key+'Median']=qf[1,j]*factor;rr[key+'High']=qf[2,j]*factor
 fixed.append(rr)
pd.DataFrame(rows).to_csv(O/'unified_primary_ci.csv',index=False);pd.DataFrame(fixed).to_csv(O/'unified_fixed_ci.csv',index=False)
# Existing fixed-target primary Tables retain shared original seed; primary table uses new unified stream, explicit disclosed.
print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(fixed).to_string(index=False))
