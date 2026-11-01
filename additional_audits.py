from pathlib import Path
import sys,math,json,hashlib
import numpy as np,pandas as pd
from scipy.spatial.distance import jensenshannon
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
P=Path(__file__).resolve().parent;R=Path(os.environ['ESPR_WORK']);sys.path.insert(0,str(P))
from analysis.config import *
from analysis.data import load_dataset,valid_pollutant,period_masks
from analysis.core import *
from analysis.bootstrap import bootstrap_year_metrics,_bootstrap_weight_context,_weights_from_day_sample,_allocate_top_multiplicity
O=R/'work/new';O.mkdir(parents=True,exist_ok=True);F=O/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'svg.fonttype':'none'})
d=load_dataset(P/'Dataset.xlsx');d.to_pickle(R/'work/data.pkl');pols=list(POLLUTANTS);D={(p,y):valid_pollutant(d,p,d.year.eq(y)) for p in pols for y in (2025,2026)};B=1000;seed=20261026
rows={k:[] for k in ['audit','zero','sample','ties','primary','fixed','monthly','tail','conditional','burden','support','js','alignment','calibration','matched']}
def savecsv():
 for k,v in rows.items():
  if v:pd.DataFrame(v).to_csv(O/(k+'.csv'),index=False)
def figsave(fig,n):
 fig.savefig(F/(n+'.png'),dpi=300,bbox_inches='tight');fig.savefig(F/(n+'.svg'),bbox_inches='tight');plt.close(fig)
for pi,p in enumerate(pols):
 sp=POLLUTANTS[p];cc=sp['value_col'];sc=sp['status_col']
 for y in (2025,2026):
  z=d.loc[d.year.eq(y)];c=pd.to_numeric(z[cc],errors='coerce');q=z[QCOL];s=z[sc];a=s.isin(ACCEPTED_STATUSES);gc=a&c.gt(0)&np.isfinite(c);dt=D[p,y];req=alert_requirement(dt);N=len(dt);k=math.ceil(req['nominal_alert_fraction']*N)
  r={'Channel':p,'Period':y,'Observed':len(z),'Calibration':int(s.eq(CALIBRATION_STATUS).sum()),'Error':int(s.eq(DEVICE_ERROR_STATUS).sum()),'Blank':int(s.isna().sum()),'Other':int((~a&~s.eq(CALIBRATION_STATUS)&~s.eq(DEVICE_ERROR_STATUS)&s.notna()).sum()),'MissingC':int((a&~np.isfinite(c)).sum()),'NonpositiveC':int((a&np.isfinite(c)&c.le(0)).sum()),'ZeroQ':int((gc&q.eq(0)).sum()),'InvalidQ':int((gc&(~np.isfinite(q)|q.lt(0))).sum()),'Primary':N,'AllZeroQ':int(q.eq(0).sum())};assert sum(r[k] for k in ['Calibration','Error','Blank','Other','MissingC','NonpositiveC','ZeroQ','InvalidQ','Primary'])==len(z);rows['audit'].append(r)
  rows['sample'].append({'Channel':p,'Period':y,'n':N,'Cmedian':dt.C.median(),'Cq1':dt.C.quantile(.25),'Cq3':dt.C.quantile(.75),'Qmedian':dt.Q.median(),'Qq1':dt.Q.quantile(.25),'Qq3':dt.Q.quantile(.75),'AllZeroQ':r['AllZeroQ']})
  rows['ties'].append({'Channel':p,'Period':y,'n':N,'Target':req['target_n'],'Nominal':req['nominal_alert_fraction']*100,'k':k,'Threshold':req['threshold'],'CutoffTies':int(dt.C.eq(req['threshold']).sum()),'Extra':req['flag_n']-k,'Alerts':req['flag_n'],'Realized':req['realized_alert_fraction']*100,'ExcessPP':(req['realized_alert_fraction']-req['nominal_alert_fraction'])*100,'Recall':req['recall']*100})
  p95=dt.C.quantile(.95);p99=dt.C.quantile(.99)
  for st,mm in [('Q=0',gc&q.eq(0)),('Q>0',gc&q.gt(0))]:
   x=c.loc[mm];rows['zero'].append({'Channel':p,'Period':y,'State':st,'n':len(x),'Median':x.median(),'q1':x.quantile(.25),'q3':x.quantile(.75),'p95':x.quantile(.95),'PosFlowP95':p95,'AboveP95':int(x.ge(p95).sum()),'AboveP99':int(x.ge(p99).sum()),'AboveP95pct':x.ge(p95).mean()*100})
  n=math.ceil(.01*N);tg=np.zeros(N,bool);tg[stable_top_indices(dt.interval_mass.to_numpy(),.01)]=True;hc=dt.C.ge(p95).to_numpy();hq=dt.Q.ge(dt.Q.quantile(.95)).to_numpy();ma=dt.interval_mass.to_numpy()
  for i,j in [(1,1),(1,0),(0,1),(0,0)]:
   mask=tg&(hc==i)&(hq==j);rows['tail'].append({'Channel':p,'Period':y,'Cclass':'>=P95' if i else '<P95','Qclass':'>=P95' if j else '<P95','n':int(mask.sum()),'Target':n,'Share':mask.sum()/n*100,'MassShare':ma[mask].sum()/ma[tg].sum()*100})
 a=D[p,2025];dt=D[p,2026];th=alert_requirement(a)['threshold'];bd=p99_mass_boundary(a);ev=evaluate_threshold(dt,th);fx=evaluate_fixed_mass_target(dt,th,bd);rows['primary'].append({'Channel':p,'Threshold':th,**ev});rows['fixed'].append({'Channel':p,'Threshold':th,'Boundary':bd,'TargetDays':dt.loc[dt.interval_mass.ge(bd),'date'].nunique(),**fx})
 for mo in range(1,8):
  dm=dt.loc[dt.time.dt.month.eq(mo)].reset_index(drop=True);e=evaluate_threshold(dm,th);dc,nd,c,ma,co,order,ge=_bootstrap_weight_context(dm);flag=c>=th;rng=np.random.default_rng(seed+pi*20+mo);v=[]
  for b in range(B):
   w=_weights_from_day_sample(rng,dc,nd);n=int(w.sum());m=math.ceil(.01*n);tm=_allocate_top_multiplicity(order,w,m);v.append([w[flag].sum()/n,tm[flag].sum()/m])
  ci=np.percentile(v,[2.5,97.5],axis=0);rows['monthly'].append({'Channel':p,'Month':f'2026-{mo:02}','n':e['n'],'Target':e['target_n'],'Alerts':e['alerts'],'TP':e['true_positives'],'Alert':e['alert_rate']*100,'AlertLow':ci[0,0]*100,'AlertHigh':ci[1,0]*100,'Recall':e['recall']*100,'RecallLow':ci[0,1]*100,'RecallHigh':ci[1,1]*100})
 dc,nd,c,ma,co,order,ge=_bootstrap_weight_context(dt);fl=c>=th;obs=np.zeros(len(dt),bool);obs[order[:math.ceil(.01*len(dt))]]=True;rng=np.random.default_rng(FIXED_HOLDOUT_SEEDS[p]);v=[]
 for b in range(B):
  w=_weights_from_day_sample(rng,dc,nd);m=math.ceil(.01*w.sum());tm=_allocate_top_multiplicity(order,w,m);nt=w[obs].sum();v.append([tm[fl].sum()/m,w[obs&fl].sum()/nt if nt else np.nan])
 q=np.nanpercentile(v,[2.5,50,95,97.5],axis=0);rows['conditional'].append({'Channel':p,'JointLow':q[0,0]*100,'JointHigh':q[3,0]*100,'OneSidedUpper':q[2,0]*100,'FrozenMedian':q[1,1]*100,'FrozenLow':q[0,1]*100,'FrozenHigh':q[3,1]*100,'FrozenEmpty':int(np.isnan(np.array(v)[:,1]).sum())})
 b1=bootstrap_year_metrics(a,p,'2025',B=B,return_samples=True);b2=bootstrap_year_metrics(dt,p,'2026',B=B,return_samples=True);diff=(b2['_burden_samples']-b1['_burden_samples'])*100;q=np.percentile(diff,[2.5,50,97.5]);rows['burden'].append({'Channel':p,'ObservedPP':(alert_requirement(dt)['nominal_alert_fraction']-alert_requirement(a)['nominal_alert_fraction'])*100,'Low':q[0],'Median':q[1],'High':q[2]})
 print('QC and bootstrap',p,flush=True);savecsv()
rng=np.random.default_rng(REVISION_METRICS_BOOTSTRAP_SEED)
for p in pols:
 dt=D[p,2026];fm=dt.interval_mass.to_numpy()>=p99_mass_boundary(D[p,2025]);dc,nd,*_=_bootstrap_weight_context(dt);v=[]
 for b in range(B):
  w=_weights_from_day_sample(rng,dc,nd);v.append([np.count_nonzero((w>0)&fm),w[fm].sum()])
 q=np.percentile(v,[2.5,50,97.5],axis=0);rows['support'].append({'Channel':p,'Target':fm.sum(),'Days':dt.loc[fm,'date'].nunique(),'UniqueLow':q[0,0],'UniqueMedian':q[1,0],'UniqueHigh':q[2,0],'FractionMedian':q[1,0]/fm.sum()*100,'Empty':int((np.array(v)[:,1]==0).sum())})
# JS reference: bootstrap day histograms within 2025 at original period day counts; common fixed pooled bins.
for pi,p in enumerate(pols):
 a=D[p,2025];dt=D[p,2026];po=pd.concat([a,dt]);ex=np.unique(np.quantile(np.log1p(po.C),np.linspace(0,1,26)));ey=np.unique(np.quantile(np.log1p(po.Q),np.linspace(0,1,26)))
 def hist(x):return np.stack([np.histogram2d(np.log1p(g.C),np.log1p(g.Q),bins=[ex,ey])[0].ravel() for _,g in x.groupby('date')])
 h1=hist(a);h2=hist(dt);ob=jensenshannon(h1.sum(0)+.5,h2.sum(0)+.5,base=2);rng=np.random.default_rng(seed+200+pi);v=[]
 for b in range(B):
  w1=np.bincount(rng.integers(0,len(h1),len(h1)),minlength=len(h1));w2=np.bincount(rng.integers(0,len(h1),len(h2)),minlength=len(h1));v.append(jensenshannon(w1@h1+.5,w2@h1+.5,base=2))
 v=np.array(v);np.save(O/f'js_{p}.npy',v);q=np.percentile(v,[2.5,50,95,97.5,99]);rows['js'].append({'Channel':p,'Observed':ob,'Low':q[0],'Median':q[1],'P95':q[2],'High':q[3],'P99':q[4],'Exceedances':int((v>=ob).sum()),'Days2025':len(h1),'Days2026':len(h2)})
 print('JS',p,flush=True)
fig,axs=plt.subplots(1,3,figsize=(9,3))
for ax,r in zip(axs,rows['js']):
 ax.hist(np.load(O/f"js_{r['Channel']}.npy"),bins=25,color='#0072B2');ax.axvline(r['Observed'],color='#D55E00',label='Observed');ax.axvline(r['P99'],color='black',ls=':',label='Reference P99');ax.set_xlim(0,.6);ax.set_title(r['Channel']);ax.set_xlabel('JS distance (base 2)')
axs[0].set_ylabel('Replicates');axs[0].legend(frameon=False);fig.tight_layout();figsave(fig,'FigS6_JS_reference')
# Analyses present in manuscript but absent from base runner.
qs=d.set_index('time')[QCOL]
for p in pols:
 sp=POLLUTANTS[p];c=pd.to_numeric(d[sp['value_col']],errors='coerce');acc=d[sp['status_col']].isin(ACCEPTED_STATUSES)&c.gt(0)
 for delta in [-30,-15,-10,-5,0,5,10,15,30]:
  q=qs.reindex(pd.DatetimeIndex(d.time+pd.Timedelta(minutes=delta))).to_numpy();mm=acc&np.isfinite(q)&(q>0);x=d.loc[mm,['time','date','year','month','quarter']].copy();x['C']=c.loc[mm].to_numpy();x['Q']=q[mm];x['load_rate']=x.C*x.Q/1000;x['interval_mass']=x.load_rate/12;a=x.loc[x.year.eq(2025)].reset_index(drop=True);b=x.loc[x.year.eq(2026)].reset_index(drop=True);th=alert_requirement(a)['threshold'];e=evaluate_threshold(b,th);rows['alignment'].append({'Channel':p,'Offset':delta,'Threshold':th,'n':e['n'],'Alert':e['alert_rate']*100,'Recall':e['recall']*100,'Capture':e['target_mass_capture']*100})
 ts=d.time.to_numpy(dtype='datetime64[ns]').astype('int64');ct=np.sort(d.loc[d[sp['status_col']].eq(CALIBRATION_STATUS),'time'].to_numpy(dtype='datetime64[ns]').astype('int64'))
 for hr in [0,12,24]:
  keep=np.ones(len(d),bool)
  if hr and len(ct):
   pos=np.searchsorted(ct,ts);dist=np.minimum(abs(ts-ct[np.clip(pos-1,0,len(ct)-1)]),abs(ts-ct[np.clip(pos,0,len(ct)-1)]));keep=dist>hr*3600*1e9
  a=valid_pollutant(d,p,d.year.eq(2025)&keep);b=valid_pollutant(d,p,d.year.eq(2026)&keep);th=alert_requirement(a)['threshold'];e=evaluate_threshold(b,th);rows['calibration'].append({'Channel':p,'Hours':hr,'Threshold':th,'Alert':e['alert_rate']*100,'Recall':e['recall']*100,'Capture':e['target_mass_capture']*100})
 for mo in range(1,8):
  a=D[p,2025].loc[D[p,2025].time.dt.month.eq(mo)].reset_index(drop=True);b=D[p,2026].loc[D[p,2026].time.dt.month.eq(mo)].reset_index(drop=True);th=alert_requirement(a)['threshold'];e=evaluate_threshold(b,th);rows['matched'].append({'Channel':p,'Month':mo,'Threshold':th,'Alert':e['alert_rate']*100,'Recall':e['recall']*100,'Precision':e['precision']*100,'Capture':e['target_mass_capture']*100})
savecsv();(O/'metadata.json').write_text(json.dumps({'seed':seed,'B':B,'sha256':hashlib.sha256((P/'Dataset.xlsx').read_bytes()).hexdigest()},indent=2));print('DONE',flush=True)
