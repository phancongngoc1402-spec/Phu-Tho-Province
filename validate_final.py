from pathlib import Path
import argparse,math,json,hashlib
import numpy as np,pandas as pd
from analysis.calibration import rank_grid_calibration,rank_grid_sorted
from analysis.bootstrap import _derive_threshold_arrays,_bootstrap_weight_context,_weights_from_day_sample
from analysis.data import load_dataset,valid_pollutant
ROOT=Path(__file__).resolve().parent

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=ROOT/'outputs');ap.add_argument('--write-reference',action='store_true');a=ap.parse_args();F=a.output/'final';C=[]
 def ck(n,v):
  C.append({'name':n,'passed':bool(v)})
  if not v:print('FAIL',n)
 def cl(n,x,y):ck(n,np.allclose(x,y,rtol=1e-9,atol=1e-9,equal_nan=True))
 ck('dataset checksum',hashlib.sha256((ROOT/'Dataset.xlsx').read_bytes()).hexdigest()=='253fd20f2e26f9e5b5cae57a4ed1140433f01ef141805e20b7b1f11b226fa283');meta=json.loads((F/'run_metadata.json').read_text());ck('primary B',meta['primary_B']==1000);ck('secondary B',meta['secondary_B']==500)
 t=pd.read_csv(F/'Table3_locked_threshold_holdout.csv').set_index('Channel');expected={'COD':(47.43,40677,407,489,200),'TSS':(35.44,41854,419,2529,251),'NH4-N':(3.31,41870,419,11831,416)}
 for p,(th,n,nt,al,tp) in expected.items():
  r=t.loc[p]
  for k,v in [('Threshold',th),('n',n),('target_n',nt),('alerts',al),('true_positives',tp)]:cl(p+k,r[k],v)
  cl(p+'recall',r.recall,tp/nt);cl(p+'precision',r.precision,tp/al);cl(p+'alert',r.alert_rate,al/n)
 ref=ROOT/'reference_final';ref.mkdir(exist_ok=True)
 if a.write_reference:
  for f in F.glob('*.csv'):
   if f.name!='OUTPUT_MAP.csv':(ref/f.name).write_bytes(f.read_bytes())
 for f in ref.glob('*.csv'):
  x=pd.read_csv(F/f.name);y=pd.read_csv(f);ck(f.name+'shape',x.shape==y.shape and list(x.columns)==list(y.columns))
  for k in y:
   if pd.api.types.is_numeric_dtype(x[k]) and pd.api.types.is_numeric_dtype(y[k]):cl(f.name+k,x[k],y[k])
   else:ck(f.name+k,x[k].fillna('<NA>').astype(str).equals(y[k].fillna('<NA>').astype(str)))
 rng=np.random.default_rng(20261031)
 for i in range(100):
  n=int(rng.integers(20,180));c=np.round(rng.uniform(.1,10,n),1);m=np.round(rng.uniform(.1,20,n),1);w=rng.integers(0,5,n);w[0]=max(1,w[0]);ce=np.repeat(c,w);me=np.repeat(m,w);N=len(ce);nt=math.ceil(.1*N);co=np.argsort(-ce,kind='mergesort');mo=np.argsort(-me,kind='mergesort')[:nt];rank=np.empty(N,int);rank[co]=np.arange(1,N+1);h=math.ceil(.9*nt);k=np.sort(rank[mo])[h-1];g=min(1,max(.005,math.ceil((k/N-1e-12)/.005)*.005));th=ce[co[math.ceil(g*N)-1]];v=rank_grid_calibration(c,m,frac=.1,weights=w);cl('oraclethreshold'+str(i),v['threshold'],th);cl('oraclefraction'+str(i),v['nominal_alert_fraction'],g);cl('delegate'+str(i),_derive_threshold_arrays(ce,me,frac=.1),th)
 d=load_dataset(ROOT/'Dataset.xlsx')
 for p in expected:
  z=valid_pollutant(d,p,d.year.eq(2025));dc,nd,c,m,co,mo,ge=_bootstrap_weight_context(z);rng=np.random.default_rng(20261031)
  for i in range(50):
   w=_weights_from_day_sample(rng,dc,nd);v=rank_grid_sorted(c,co,mo,w,math.ceil(.01*w.sum()));u=rank_grid_calibration(np.repeat(c,w),np.repeat(m,w));cl(p+'bootstrap equivalence'+str(i),[v['threshold'],v['nominal_alert_fraction']],[u['threshold'],u['nominal_alert_fraction']])
 ac=pd.read_csv(F/'TableS18_acf.csv');ck('Pearson bounds',np.all(abs(ac[['Lag1','ACF24h']].to_numpy())<=1+1e-10));s=pd.read_csv(F/'TableS7c_target_support.csv').set_index('Channel');ck('COD empty6',s.loc['COD','Empty']==6);ck('COD support26/5',s.loc['COD','Target']==26 and s.loc['COD','Days']==5);f=pd.read_csv(F/'TableS7b_fixed_P99_bootstrap.csv').set_index('Channel');cl('COD100-100',[f.loc['COD','Recall CI low'],f.loc['COD','Recall CI high']],[1,1]);ex=pd.read_csv(F/'TableS15_exclusions.csv')
 for _,r in ex.iterrows():cl(str(r.Channel)+str(r.Period)+'partition',sum(r[k] for k in ['Calibration','Error','Blank','Other','MissingC','NonpositiveC','ZeroQ','InvalidQ','Primary']),r.Observed)
 fail=sum(not x['passed'] for x in C);out={'checks':len(C),'passed':len(C)-fail,'failed':fail,'status':'PASS' if not fail else 'FAIL','details':C};(a.output/'validation_report.json').write_text(json.dumps(out,indent=2));print('FINAL VALIDATION:',out['passed'],'/',len(C),'passed;',fail,'failed');return int(fail>0)
if __name__=='__main__':raise SystemExit(main())
