"""Quality audit and deterministic time-series aggregation."""
from acquire import ROOT,AUDIT,RESULT
import pandas as pd,numpy as np,zipfile,json,hashlib,sys

def aggregate(x,clip=False):
    v=x.Depth.clip(lower=0) if clip else x.Depth
    g=v.resample('1h')
    y=pd.DataFrame({'n':g.count(),'median_depth_m':g.median(),'depth_iqr_m':g.quantile(.75)-g.quantile(.25),'mean_depth_m':g.mean(),'surface_fraction':v.le(5).resample('1h').mean(),'temperature_c':x.Temperature.resample('1h').median()})
    for minutes,need,pairs in [(1,3,45),(5,15,9)]:
        m=v.resample(f'{minutes}min').agg(['median','count'])
        m.loc[m['count']<need,'median']=np.nan
        rate=m['median'].diff().abs()/minutes
        # Both points must be in the same clock hour.
        same=m.index.floor('h')==pd.DatetimeIndex(m.index-pd.Timedelta(minutes=minutes)).floor('h')
        rate=rate.where(same)
        rg=rate.resample('h');col='vertical_rate_m_min' if minutes==1 else 'vertical_rate_5min_m_min'
        y[col]=rg.mean().where(rg.count()>=pairs);y[col+'_pairs']=rg.count()
    y.loc[y.n<180,['median_depth_m','depth_iqr_m','mean_depth_m','surface_fraction']]=np.nan
    return y

def main():
    p=AUDIT/'sources/tags/14_07-132980.zip'
    with zipfile.ZipFile(p) as z:
        with z.open('out-Archive.csv') as f:x=pd.read_csv(f,usecols=['Time','Depth','Temperature'],low_memory=False)
    x['time']=pd.to_datetime(x.Time,format='%H:%M:%S %d-%b-%Y',utc=True,errors='coerce')
    bad=int(x.time.isna().sum());assert bad==0
    x.Depth=pd.to_numeric(x.Depth,errors='coerce');x.Temperature=pd.to_numeric(x.Temperature,errors='coerce')
    dup=int(x.time.duplicated().sum())
    q={'raw_rows':len(x),'invalid_timestamps':bad,'raw_duplicate_timestamp_rows':dup,'raw_zip_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
    q['nonfinite_depth_rows_excluded']=int((~np.isfinite(x.Depth)).sum())
    x=x[np.isfinite(x.Depth)].copy()
    q['identical_finite_rows_removed']=int(x.duplicated(['time','Depth','Temperature']).sum())
    x=x.drop_duplicates(['time','Depth','Temperature'])
    assert not x.time.duplicated().any(),'Conflicting/repeated finite depth rows require review'
    q['finite_depth_duplicate_timestamps']=int(x.time.duplicated().sum())
    x=x[x.time.ge('2014-11-07')&x.time.lt('2015-08-03')].set_index('time').sort_index()
    assert np.isfinite(x.Depth).all()
    dt=x.index.to_series().diff().dt.total_seconds()
    q.update(interior_records=len(x),first_utc=str(x.index.min()),last_utc=str(x.index.max()),median_interval_seconds=float(dt.median()),maximum_interval_seconds=float(dt.max()),negative_depth_fraction=float(x.Depth.lt(0).mean()),minimum_depth_m=float(x.Depth.min()),maximum_depth_m=float(x.Depth.max()),minimum_temperature_c=float(x.Temperature.min()),maximum_temperature_c=float(x.Temperature.max()))
    event=x.loc['2015-05-25':'2015-06-09']
    q['event_negative_depth_fraction']=float(event.Depth.lt(0).mean())
    q['event_minimum_depth_m']=float(event.Depth.min())
    q['event_maximum_vertical_raw_step_m']=float(event.Depth.diff().abs().max())
    for clip,name in [(False,'raw'),(True,'clipped')]:
        y=aggregate(x,clip);y.to_csv(RESULT/f'hourly_{name}.csv',index_label='time')
    daily=x.Depth.resample('D').agg(['count','min','median','max'])
    daily['q01']=x.Depth.resample('D').quantile(.01);daily['negative_fraction']=x.Depth.lt(0).resample('D').mean()
    daily.to_csv(RESULT/'daily_depth_quality.csv',index_label='time')
    # Retain event raw data only as a compressed derived subset, preserving source fields.
    event[['Depth','Temperature']].to_csv(RESULT/'event_depth_subset.csv.gz',index_label='time')
    (RESULT/'depth_quality.json').write_text(json.dumps(q,indent=2))
    print(json.dumps(q,indent=2),flush=True)

if __name__=='__main__':main()
