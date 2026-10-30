"""Match model weather to model shark positions with transparent spatial offsets."""
from acquire import ROOT,SOURCE,RESULT
import pandas as pd,numpy as np,json

def hav(a,b,c,d):
    a,b,c,d=map(np.radians,[a,b,c,d]);v=np.sin((c-a)/2)**2+np.cos(a)*np.cos(c)*np.sin((d-b)/2)**2
    return 6371.0088*2*np.arcsin(np.sqrt(np.clip(v,0,1)))

def onset(frame,wind=17.5,wave=4):
    x=frame.set_index('time').reindex(pd.date_range('2015-06-05','2015-06-11',freq='h',tz='UTC',inclusive='left'))
    yes=(x.wind_speed_10m>=wind)|(x.Thgt>=wave)
    starts=yes&yes.shift(-1,fill_value=False)&yes.shift(-2,fill_value=False)
    return str(starts.index[starts][0]) if starts.any() else None

def main():
    sites=pd.read_csv(RESULT/'weather_requested_sites.csv');frames=[];coords=[]
    for r in sites.to_dict('records'):
        p=SOURCE/(r['site_id']+'.json')
        if not p.exists():continue
        j=json.loads(p.read_text());assert j['hourly_units']['wind_speed_10m']=='m/s' and j['hourly_units']['pressure_msl']=='hPa'
        d=pd.DataFrame(j['hourly']);d['time']=pd.to_datetime(d.time,utc=True);d['site_id']=r['site_id'];frames.append(d)
        coords.append(dict(**r,grid_lat=j['latitude'],grid_lon=j['longitude'],grid_elevation_m=j.get('elevation'),timezone=j['timezone']))
    era=pd.concat(frames,ignore_index=True);coords=pd.DataFrame(coords);coords.to_csv(RESULT/'weather_returned_sites.csv',index=False)
    era.to_csv(RESULT/'era5_hourly_sites.csv.gz',index=False)
    wave=pd.read_csv(SOURCE/'ww3_region.csv',skiprows=[1]);wave['time']=pd.to_datetime(wave.time,utc=True)
    assert not wave.duplicated(['time','latitude','longitude']).any()
    # Grid selection is fixed using missingness alone, independent of wave magnitude.
    wet=wave.groupby(['latitude','longitude']).Thgt.agg(['count','size']).reset_index()
    wet=wet[wet['count']/wet['size']>=.9].copy();assert len(wet)>0
    wet.to_csv(RESULT/'wave_wet_grid_cells.csv',index=False)
    p=pd.read_csv(RESULT/'shark_modelled_positions.csv');p['position_time']=pd.to_datetime(p.time,utc=True);p=p.drop(columns='time')
    hours=pd.DataFrame({'time':pd.date_range('2015-05-01','2015-06-21',freq='h',tz='UTC',inclusive='left')})
    traj=pd.merge_asof(hours,p[['position_time','lat','lon']].sort_values('position_time'),left_on='time',right_on='position_time',direction='nearest',tolerance=pd.Timedelta(hours=6))
    traj['position_lag_hours']=(traj.time-traj.position_time).dt.total_seconds().abs()/3600
    scenarios=[('trajectory',traj)]
    for r in sites[sites.role.eq('spatial_scenario')].to_dict('records'):
        s=hours.copy();s['lat']=r['request_lat'];s['lon']=r['request_lon'];s['position_lag_hours']=np.nan
        scenarios.append((r['site_id'],s))
    output=[]
    for label,t in scenarios:
        matched=[]
        pool=coords[coords.role.eq('trajectory')] if label=='trajectory' else coords[coords.site_id.eq(label)]
        for r in t.to_dict('records'):
            if not np.isfinite(r['lat']):continue
            ds=hav(r['lat'],r['lon'],pool.grid_lat.to_numpy(),pool.grid_lon.to_numpy());s=pool.iloc[np.argmin(ds)]
            dw=hav(r['lat'],r['lon'],wet.latitude.to_numpy(),wet.longitude.to_numpy());w=wet.iloc[np.argmin(dw)]
            matched.append(dict(**r,scenario=label,site_id=s.site_id,era5_grid_lat=s.grid_lat,era5_grid_lon=s.grid_lon,era5_grid_distance_km=float(ds.min()),wave_grid_lat=w.latitude,wave_grid_lon=w.longitude,wave_grid_distance_km=float(dw.min())))
        m=pd.DataFrame(matched).merge(era,on=['time','site_id'],how='left',validate='many_to_one').merge(wave,left_on=['time','wave_grid_lat','wave_grid_lon'],right_on=['time','latitude','longitude'],how='left',validate='many_to_one')
        m.loc[m.era5_grid_distance_km>100,['wind_speed_10m','pressure_msl','temperature_2m']]=np.nan
        m.loc[m.wave_grid_distance_km>100,['Thgt','Tper','shgt','whgt']]=np.nan
        output.append(m)
    allm=pd.concat(output,ignore_index=True);allm.to_csv(RESULT/'environment_all_scenarios.csv.gz',index=False)
    main=allm[allm.scenario.eq('trajectory')];main.to_csv(RESULT/'environment_matched_hourly.csv',index=False)
    drops=main[(main.Thgt.diff()<=-1)&main.Thgt.le(.1)].copy()
    drops[['time','wave_grid_lat','wave_grid_lon','Thgt']].to_csv(RESULT/'wave_source_discontinuities.csv',index=False)
    summaries=[]
    for label,m in allm.groupby('scenario'):
        e=m[m.time.ge('2015-06-05')&m.time.lt('2015-06-11')]
        summaries.append(dict(scenario=label,requested_lat=float(e.lat.iloc[0]) if label!='trajectory' else None,requested_lon=float(e.lon.iloc[0]) if label!='trajectory' else None,event_hours=len(e),wind_hours=int(e.wind_speed_10m.notna().sum()),wave_hours=int(e.Thgt.notna().sum()),max_wind_m_s=float(e.wind_speed_10m.max()),max_wave_m=float(e.Thgt.max()),minimum_pressure_hpa=float(e.pressure_msl.min()),peak_wind_utc=str(e.loc[e.wind_speed_10m.idxmax(),'time']),peak_wave_utc=str(e.loc[e.Thgt.idxmax(),'time']),max_era5_offset_km=float(e.era5_grid_distance_km.max()),max_wave_offset_km=float(e.wave_grid_distance_km.max()),screen_onset_utc=onset(m),wave3m_sensitivity_onset_utc=onset(m,wave=3),wave5m_sensitivity_onset_utc=onset(m,wave=5)))
    d=pd.DataFrame(summaries);d.to_csv(RESULT/'environment_summary.csv',index=False)
    (RESULT/'environment_gate.json').write_text(json.dumps({'model_based_onset_utc':onset(main),'observational_exposure_confirmed':False,'wind_threshold_m_s':17.5,'wave_threshold_m':4,'consecutive_hours':3,'warning':'Models and point matching do not establish measured local exposure. Geographic scenarios are not confidence bounds.'},indent=2))
    print(d.to_string(index=False),flush=True)

if __name__=='__main__':main()
