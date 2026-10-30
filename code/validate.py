"""Independent checks for core estimates, aggregation, integrity and old outputs."""
from acquire import ROOT,AUDIT,RESULT,SOURCE
from depth import aggregate
from analyse import REF,phase
import pandas as pd,numpy as np,json,hashlib,zipfile,os,sys
from pathlib import Path

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def main():
    checks={}
    # Independently reconstruct design, OLS and the Bartlett-weighted HAC sandwich.
    y=pd.read_csv(RESULT/'hourly_raw.csv',index_col='time',parse_dates=['time'])
    t=y.loc[(y.index>=REF-pd.Timedelta(hours=336))&(y.index<REF+pd.Timedelta(hours=24))]
    hours=(t.index-REF).total_seconds()/3600
    indicators=[((hours>=a)&(hours<b)).astype(float) for a,b in [(-72,-48),(-48,-24),(-24,0),(0,24)]]
    X=np.column_stack([np.ones(360),hours/24,*indicators,*[(t.index.hour==h).astype(float) for h in range(1,24)]])
    effects=pd.read_csv(RESULT/'main_effects.csv');tested=0
    for endpoint in effects.outcome.unique():
        v=t[endpoint].to_numpy();beta=np.linalg.lstsq(X,v,rcond=None)[0];res=v-X@beta
        scores=X*res[:,None];S=scores.T@scores
        for lag in range(1,25):
            cross=scores[lag:].T@scores[:-lag];S+=(1-lag/25)*(cross+cross.T)
        bread=np.linalg.inv(X.T@X);cov=bread@S@bread*len(X)/(len(X)-X.shape[1])
        for k,p in enumerate(['pre48_24','pre24','post24'],start=3):
            r=effects[(effects.outcome==endpoint)&(effects.phase==p)].iloc[0]
            assert np.isclose(beta[k],r.adjusted_difference,atol=1e-9)
            assert np.isclose(np.sqrt(cov[k,k]),r.standard_error,atol=1e-9)
            tested+=1
    checks['independent_ols_and_hac_comparisons']=tested
    # Holm family adjustment checked without statsmodels.
    p=effects.p_value.to_numpy();o=np.argsort(p);adj=np.minimum(1,np.maximum.accumulate(p[o]*(len(p)-np.arange(len(p)))))
    assert np.allclose(adj,effects.p_holm.to_numpy()[o]);checks['holm_family_size']=len(p)
    # Synthetic linear dive has an exact one-metre/minute vertical rate.
    index=pd.date_range('2015-01-01',periods=240,freq='15s',tz='UTC')
    synthetic=pd.DataFrame({'Depth':np.arange(240)/4,'Temperature':20},index=index)
    a=aggregate(synthetic);assert a.n.iloc[0]==240 and a.median_depth_m.iloc[0]==29.875
    assert np.isclose(a.vertical_rate_m_min.iloc[0],1) and np.isclose(a.vertical_rate_5min_m_min.iloc[0],1)
    broken=synthetic.copy();broken.loc[index[:60],'Depth']=np.nan
    b=aggregate(broken);assert pd.isna(b.vertical_rate_m_min.iloc[0]),'Gap must reduce valid movement coverage'
    checks['synthetic_aggregation_and_gap_checks']=True
    # Re-open raw event data and check a full hour with explicit numpy quantiles/differences.
    raw=pd.read_csv(RESULT/'event_depth_subset.csv.gz');raw['time']=pd.to_datetime(raw.time,utc=True)
    for hour in pd.date_range('2015-06-06 12:00',periods=72,freq='h',tz='UTC'):
        v=raw.loc[(raw.time>=hour)&(raw.time<hour+pd.Timedelta(hours=1)),'Depth'].to_numpy()
        assert len(v)==240 and np.isclose(np.median(v),y.loc[hour,'median_depth_m'])
        assert np.isclose(np.subtract(*np.quantile(v,[.75,.25])),y.loc[hour,'depth_iqr_m'])
        assert np.isclose(np.abs(np.diff(np.median(v.reshape(60,4),axis=1))).mean(),y.loc[hour,'vertical_rate_m_min'])
    checks['event_hour_independent_aggregation_checks']=72
    # Check environmental return units, complete event time coverage and threshold outcome.
    e=pd.read_csv(RESULT/'environment_matched_hourly.csv');e['time']=pd.to_datetime(e.time,utc=True)
    e=e[(e.time>='2015-06-05')&(e.time<'2015-06-11')]
    assert len(e)==144 and not e[['wind_speed_10m','pressure_msl','Thgt']].isna().any().any()
    assert (e.position_lag_hours<=6).all() and e.wind_speed_10m.max()<17.5 and e.Thgt.max()<4
    checks['complete_environment_event_hours']=len(e)
    srcs=[]
    for p in SOURCE.iterdir():
        if p.is_file() and not p.name.endswith('.retrieval.json'):srcs.append({'file':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)})
    for p in SOURCE.glob('*.retrieval.json'):
        j=json.loads(p.read_text())
        if j['status']==200:assert sha(SOURCE/j['file'])==j['sha256']
    (ROOT/'SOURCE_MANIFEST.json').write_text(json.dumps(srcs,indent=2))
    archive=AUDIT/'sources/tags/14_07-132980.zip'
    stored=next(r for r in json.loads((AUDIT/'validation/tag_archive_validation.json').read_text()) if r['fileName']==archive.name)
    assert sha(archive)==stored['sha256'];checks['original_shark_zip_unchanged']=True
    old=json.loads((AUDIT/'validation/protected_before.json').read_text());found=set()
    for folder in ['outfolder','All_7_Species_Research','Shark_PreStorm_2026-09-16']:
        for base,dirs,files in os.walk(ROOT.parent/folder,followlinks=False):
            dirs[:]=[d for d in dirs if d not in ['node_modules','__pycache__']]
            found.update((Path(base)/f).relative_to(ROOT.parent).as_posix() for f in files)
    assert found=={r['file'] for r in old}
    for r in old:
        p=ROOT.parent/r['file'];s=p.stat()
        assert sha(p)==r['sha256'] and s.st_size==r['size'] and s.st_mtime_ns==r['mtime_ns'],r['file']
    checks['protected_prior_files_unchanged']=len(old)
    (RESULT/'validation.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
