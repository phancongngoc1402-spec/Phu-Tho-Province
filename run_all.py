from pathlib import Path
import argparse,sys,subprocess,tempfile,shutil,json,hashlib,os,platform
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('MPLBACKEND','Agg')
ROOT=Path(__file__).resolve().parent

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,default=ROOT/'Dataset.xlsx');ap.add_argument('--output',type=Path,default=ROOT/'outputs');ap.add_argument('--base-only',action='store_true',help='Internal resumable stage; not full manuscript output');ap.add_argument('--resume',action='store_true',help='Reuse base only when input, code and dependency versions match its signature');ap.add_argument('--skip-figures',action='store_true');a=ap.parse_args();O=a.output.resolve();F=O/'final';F.mkdir(parents=True,exist_ok=True);work=O/'work';(work/'work/new').mkdir(parents=True,exist_ok=True);env=os.environ.copy();env['ESPR_WORK']=str(work)
 if a.input.resolve()!=(ROOT/'Dataset.xlsx').resolve():raise ValueError('Bundled dataset is required; no external replacement is inferred')
 digest=hashlib.sha256(a.input.read_bytes()).hexdigest()
 if digest!='253fd20f2e26f9e5b5cae57a4ed1140433f01ef141805e20b7b1f11b226fa283':raise ValueError('Checksum mismatch')
 if not a.base_only:
  if not (a.resume and (O/'base_complete.txt').exists() and (O/'base_complete.txt').read_text().strip()==stage_signature(a.input)):run_base(a.input,O)
 else:run_base(a.input,O);return
 subprocess.run([sys.executable,str(ROOT/'additional_audits.py')],cwd=ROOT,env=env,check=True)
 subprocess.run([sys.executable,str(ROOT/'revision_statistics.py')],cwd=ROOT,env=env,check=True)
 if not a.skip_figures:subprocess.run([sys.executable,str(ROOT/'make_figures.py')],cwd=ROOT,env=env,check=True)
 finish(O,work,digest,not a.skip_figures)

def stage_signature(inp):
 import numpy,pandas,scipy,matplotlib,openpyxl
 h=hashlib.sha256(Path(inp).read_bytes())
 for p in sorted((ROOT/'analysis').glob('*.py')):h.update(p.name.encode());h.update(p.read_bytes())
 for m in [numpy,pandas,scipy,matplotlib,openpyxl]:h.update((m.__name__+'='+m.__version__).encode())
 return h.hexdigest()

def run_base(inp,O):
 cmd=[sys.executable,'-m','analysis.base_pipeline','--input',str(inp),'--output',str(O/'base'),'--bootstrap','1000','--secondary-bootstrap','500','--skip-figures'];subprocess.run(cmd,cwd=ROOT,check=True);(O/'base_complete.txt').write_text(stage_signature(inp)+'\n')

def finish(O,W,digest,figures):
 import numpy as np,pandas as pd
 from analysis.data import load_dataset,valid_pollutant
 from analysis.core import alert_requirement,evaluate_threshold
 from analysis.robustness import high_coverage_audit
 from analysis.config import POLLUTANTS
 F=O/'final';B=O/'base';N=W/'work/new';save=lambda d,p:d.to_csv(p,index=False,encoding='utf-8-sig')
 for dr in ['main','supplement']:
  for p in (B/dr).glob('*.csv'):shutil.copy2(p,F/p.name)
 for old,new in [('TableS4_expanding_origin.csv','TableS4b_expanding_origin.csv'),('TableS4b_quarter_fixed_threshold.csv','TableS4c_quarter.csv')]:shutil.copy2(B/'supplement'/old,F/new)
 for n in ['Table1_primary_sample.csv','Table3_locked_threshold_holdout.csv','TableS4_expanding_origin.csv','TableS4b_quarter_fixed_threshold.csv','TableS6_influential_day_high_coverage.csv','Month_level_locked_threshold.csv']:(F/n).unlink(missing_ok=True)
 mapping={'sample':'Table1_sample','monthly':'TableS4a_monthly','support':'TableS7c_target_support','alignment':'TableS12_alignment','calibration':'TableS13_calibration','matched':'TableS14_matched','audit':'TableS15_exclusions','zero':'TableS16_zero_flow','ties':'TableS17_ties','acf':'TableS18_acf','burden':'TableS19_burden','unified_primary_ci':'TableS20_conditioning','js':'TableS21_js','tail':'TableS22_tail'}
 for old,new in mapping.items():shutil.copy2(N/(old+'.csv'),F/(new+'.csv'))
 save(pd.read_csv(N/'primary.csv').merge(pd.read_csv(N/'unified_primary_ci.csv'),on='Channel'),F/'Table3_locked_threshold_holdout.csv');save(pd.read_csv(N/'fixed.csv').merge(pd.read_csv(F/'TableS7b_fixed_P99_bootstrap.csv'),on='Channel'),F/'Table5_fixed_historical.csv')
 d=pd.read_pickle(W/'work/data.pkl');inf=[];months=[]
 for p in POLLUTANTS:
  de=valid_pollutant(d,p,d.year.eq(2025));ho=valid_pollutant(d,p,d.year.eq(2026));th=alert_requirement(de)['threshold'];inf.append({'Channel':p,'Threshold':th,**evaluate_threshold(ho[ho.date.ne(pd.Timestamp('2026-07-13'))].reset_index(drop=True),th)});target=np.zeros(len(ho),bool);target[np.argsort(-ho.interval_mass.to_numpy(),kind='mergesort')[:int(np.ceil(.01*len(ho)))]]=True
  for mo in range(1,8):
   mm=ho.time.dt.month.eq(mo).to_numpy();nt=int((mm&target).sum());tp=int((mm&target&ho.C.ge(th).to_numpy()).sum());months.append({'Channel':p,'Month':f'2026-{mo:02}','PooledTarget':nt,'TP':tp,'Recall':100*tp/nt if nt else np.nan})
 save(pd.DataFrame(inf),F/'TableS6a_influential_day_locked.csv');save(high_coverage_audit(d),F/'TableS6b_high_coverage.csv');save(pd.DataFrame(months),F/'Pooled_target_by_month.csv')
 if figures:
  G=F/'figures';G.mkdir(exist_ok=True)
  for dr in [N/'figures',W/'work/figures']:
   for p in dr.glob('*'):shutil.copy2(p,G/p.name)
  from analysis.figures import figure_s1_coverage_status,figure_s2_hourly_error
  figure_s1_coverage_status(d,G/'FigS1.png');figure_s2_hourly_error(pd.read_csv(F/'TableS9_hourly_aggregation.csv'),G/'FigS2.png')
  import matplotlib.pyplot as plt
  t=pd.read_csv(F/'TableS11_event_bootstrap.csv');y=t['Relative change (%)'].to_numpy();fig,ax=plt.subplots(figsize=(6.5,4));ax.bar(t.Channel,y,color='#0072B2');ax.errorbar(np.arange(3),y,yerr=np.vstack([y-t['Relative change CI low'],t['Relative change CI high']-y]),fmt='none',ecolor='black',capsize=4);ax.set_ylabel('Relative change in P90 event mass (%)');ax.set_title('Zero-flow bridging: 30 min versus 0 min');ax.grid(axis='y',alpha=.2);fig.tight_layout()
  for e in ['png','svg']:fig.savefig(G/('FigS3.'+e),dpi=300,bbox_inches='tight')
  plt.close(fig)
 import scipy,matplotlib,openpyxl
 meta={'release':'2026.10.31-canonical-rank-grid','sha256':digest,'rows':len(d),'primary_B':1000,'secondary_B':500,'python':platform.python_version(),'packages':{m.__name__:m.__version__ for m in [np,pd,scipy,matplotlib,openpyxl]},'seeds':{'primary_joint':20261026,'fixed_historical_support':20260828,'monthly':'20261026+20*pollutant_index+month','burden':'20260831+10*pollutant_index+year_offset;2025=0;2026=100','nested':'725911+pollutant_index','moving':'725911+100*block_length+pollutant_index','JS':'20261026+200+pollutant_index'},'calibration':'analysis/calibration.py; same stable rank-grid estimator for observed and all bootstrap calibrations','ACF':'revision_statistics.py:acf_pairwise'};(F/'run_metadata.json').write_text(json.dumps(meta,indent=2))
 rows=[]
 for p in sorted(F.glob('*.csv')):
  if p.name=='OUTPUT_MAP.csv':continue
  loc=p.name.split('_')[0].replace('Table','Table ');seed='See run_metadata.json and analysis/config.py'
  if p.name.startswith(('Table3_','TableS20_')):seed='20261026 shared primary metric stream'
  if p.name.startswith(('Table5_','TableS7b_','TableS7c_')):seed='20260828 shared fixed-target support stream'
  rows.append({'status':'FINAL','file':'outputs/final/'+p.name,'manuscript_location':loc,'seed_note':seed})
 for g in sorted((F/'figures').glob('*.png')):
  rows.append({'status':'FINAL','file':'outputs/final/figures/'+g.name,'manuscript_location':g.stem.replace('Fig','Figure '),'seed_note':'Programmatic final figure; see script and numerical table map'})
 save(pd.DataFrame(rows),ROOT/'OUTPUT_MAP.csv');save(pd.DataFrame(rows),F/'OUTPUT_MAP.csv');print('FINAL OUTPUTS COMPLETE',F,flush=True)
if __name__=='__main__':main()
