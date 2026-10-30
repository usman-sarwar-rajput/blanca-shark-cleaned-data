"""Five chronological stages based on storm calendar dates, with quiet-day screening."""
from pathlib import Path
import pandas as pd,numpy as np,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'five_stage_comparison';OUT.mkdir(exist_ok=True)
AUDIT=ROOT.parent/'Shark_Storm_Overlap_Audit_2026-09-16'
STAGES=[('Earlier quiet*','2015-05-25','2015-05-28','#9aaebc'),('Just before','2015-05-28','2015-05-31','#d9ac40'),('Blanca active dates','2015-05-31','2015-06-10','#c56e61'),('Just after','2015-06-10','2015-06-13','#65a598'),('Later quiet*','2015-06-13','2015-06-16','#8b9eb9')]

def distance(a,b,c,d):
    a,b,c,d=map(np.radians,[a,b,c,d]);v=np.sin((c-a)/2)**2+np.cos(a)*np.cos(c)*np.sin((d-b)/2)**2
    return 6371.0088*2*np.arcsin(np.sqrt(np.clip(v,0,1)))

def main():
    y=pd.read_csv(ROOT/'results/hourly_raw.csv',parse_dates=['time']).set_index('time')
    e=pd.read_csv(ROOT/'results/environment_matched_hourly.csv',parse_dates=['time']).set_index('time')
    storm=pd.read_csv(AUDIT/'results/storm_tracks_2001_2020.csv',parse_dates=['storm_time'])
    storm=storm[(storm.SEASON==2015)&storm.USA_WIND.ge(34)].copy()
    positions=pd.read_csv(ROOT/'results/shark_modelled_positions.csv',parse_dates=['time'])[['time','lat','lon']].sort_values('time')
    rows=[];daily=[];screen=[]
    for name,a,b,color in STAGES:
        start=pd.Timestamp(a,tz='UTC');end=pd.Timestamp(b,tz='UTC')
        t=y[(y.index>=start)&(y.index<end)];w=e[(e.index>=start)&(e.index<end)]
        hours=int((end-start).total_seconds()/3600)
        assert len(t)==hours and len(w)==hours and t.n.eq(240).all()
        assert t[['median_depth_m','vertical_rate_m_min']].notna().all().all()
        assert w[['wind_speed_10m','Thgt','lat','lon']].notna().all().all()
        # Search all contemporaneous >=34-knot eastern Pacific cyclone points.
        candidates=storm[(storm.storm_time>=start)&(storm.storm_time<end)]
        minimum=None;nearest=None;nearby=False
        for sid,g in candidates.groupby('SID'):
            j=pd.merge_asof(g.sort_values('storm_time'),positions,left_on='storm_time',right_on='time',direction='nearest',tolerance=pd.Timedelta(minutes=90)).dropna(subset=['lat','lon'])
            if j.empty:continue
            km=distance(j.lat.to_numpy(),j.lon.to_numpy(),j.LAT.to_numpy(),j.LON.to_numpy())
            if minimum is None or km.min()<minimum:minimum=float(km.min());nearest=j.NAME.iloc[0]
            nearby|=bool((km<=500).any())
            screen.append(dict(stage=name,storm=j.NAME.iloc[0],sid=sid,minimum_screened_distance_km=float(km.min()),matched_track_points=len(j),within_500km=bool((km<=500).any())))
        quiet=bool((w.wind_speed_10m<10).all() and (w.Thgt<2).all() and not nearby)
        rows.append(dict(stage=name,start_utc=start,end_utc_exclusive=end,days=hours//24,hours=hours,depth_m=t.median_depth_m.mean(),vertical_rate_m_min=t.vertical_rate_m_min.mean(),depth_iqr_m=t.depth_iqr_m.mean(),near_surface_percent=100*t.surface_fraction.mean(),max_model_wind_m_s=w.wind_speed_10m.max(),max_model_wave_m=w.Thgt.max(),nearest_cyclone=nearest,nearest_cyclone_distance_km=minimum,quiet_proxy_pass=quiet,colour=color))
        for date,g in t.groupby(t.index.floor('D')):
            assert len(g)==24
            daily.append(dict(date=date,stage=name,depth_m=g.median_depth_m.mean(),vertical_rate_m_min=g.vertical_rate_m_min.mean(),depth_iqr_m=g.depth_iqr_m.mean(),hours=len(g)))
    d=pd.DataFrame(rows);dd=pd.DataFrame(daily)
    assert d.iloc[0].quiet_proxy_pass and d.iloc[-1].quiet_proxy_pass
    assert len(dd)==22 and sum(d.hours)==528
    # Period means must equal means of complete daily summaries.
    for _,r in d.iterrows():
        assert np.isclose(r.depth_m,dd[dd.stage.eq(r.stage)].depth_m.mean())
        assert np.isclose(r.vertical_rate_m_min,dd[dd.stage.eq(r.stage)].vertical_rate_m_min.mean())
    d.to_csv(OUT/'five_stage_summary.csv',index=False);dd.to_csv(OUT/'daily_values.csv',index=False);pd.DataFrame(screen).to_csv(OUT/'cyclone_screen.csv',index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.spines.left':False,'text.color':'#243b48','axes.labelcolor':'#243b48'})
    fig,axs=plt.subplots(2,1,figsize=(12,8),layout='constrained',sharex=True)
    fig.suptitle('Shark behaviour across five stages of the Blanca timeline',fontsize=18,fontweight='bold')
    for ax,col,title,label,fmt in [(axs[0],'depth_m','DEPTH: taller bar = deeper water','Mean hourly median depth (m)','{:.2f}'),(axs[1],'vertical_rate_m_min','UP-AND-DOWN MOVEMENT: taller bar = more vertical movement','Vertical movement proxy (m/min)','{:.3f}')]:
        ax.bar(np.arange(5),d[col],color=d.colour,width=.65,zorder=3)
        ax.set_ylim(0,d[col].max()*1.24);ax.set_ylabel(label,fontsize=11);ax.set_title(title,loc='left',fontsize=12,pad=12,fontweight='bold');ax.grid(axis='y',alpha=.15,zorder=0)
        for i,v in enumerate(d[col]):ax.text(i,v+d[col].max()*.025,fmt.format(v),ha='center',fontsize=12,fontweight='bold')
    labels=['1. Earlier quiet*\n25–27 May\n3 days','2. Just before\n28–30 May\n3 days','3. Blanca active dates\n31 May–9 June\n10 days','4. Just after\n10–12 June\n3 days','5. Later quiet*\n13–15 June\n3 days']
    axs[-1].set_xticks(np.arange(5),labels,fontsize=10)
    axs[-1].set_xlabel('2015, UTC calendar days. Bars show averages, not totals; the middle period is longer.\n*Quiet = lower weather activity in the checked models and no screened cyclone within 500 km.\nBlanca active dates do not establish when severe conditions reached the shark.',fontsize=10,labelpad=16)
    for ext in ['png','pdf']:fig.savefig(OUT/f'01_five_stage_comparison.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)
    # Daily view preserves variation within the ten-day storm calendar period.
    fig,axs=plt.subplots(2,1,figsize=(13,8),layout='constrained',sharex=True)
    fig.suptitle('Daily behaviour within the five stages',fontsize=18,fontweight='bold')
    for ax,col,label in [(axs[0],'depth_m','Mean hourly median depth (m)'),(axs[1],'vertical_rate_m_min','Vertical movement proxy (m/min)')]:
        ax.plot(dd.date+pd.Timedelta(hours=12),dd[col],'-o',color='#243b48',lw=1.7,ms=4,zorder=4)
        for name,a,b,color in STAGES:ax.axvspan(pd.Timestamp(a,tz='UTC'),pd.Timestamp(b,tz='UTC'),color=color,alpha=.23)
        ax.set_ylabel(label,fontsize=11);ax.grid(axis='y',alpha=.15);ax.set_ylim(bottom=0)
    axs[0].set_title('Earlier quiet*  →  Just before  →  Blanca active dates  →  Just after  →  Later quiet*',loc='left',fontsize=11)
    for name,a,b,color in STAGES:
        mid=pd.Timestamp(a,tz='UTC')+(pd.Timestamp(b)-pd.Timestamp(a))/2
        short={'Earlier quiet*':'1','Just before':'2','Blanca active dates':'3','Just after':'4','Later quiet*':'5'}[name]
        axs[0].text(mid,.94,short,transform=axs[0].get_xaxis_transform(),ha='center',fontsize=14,fontweight='bold')
    axs[-1].xaxis.set_major_locator(mdates.DayLocator(interval=2));axs[-1].xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
    axs[-1].set_xlabel('2015, UTC. Each point is one complete day. Stage numbers match the summary chart.\n*Candidate quiet days selected using environmental data; no claim that every possible disturbance was absent.',fontsize=10,labelpad=12)
    for ext in ['png','pdf']:fig.savefig(OUT/f'02_daily_five_stage_timeline.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)
    lines=['# Five-stage Blanca comparison','', 'This additional descriptive view follows the user-requested sequence: earlier quiet days, immediately before the cyclone calendar period, during that period, immediately after, and later quiet days. No new significance tests were run.','', '## Calendar definitions','', '| Stage | UTC dates, inclusive | Duration |','|---|---|---|']
    for name,a,b,c in STAGES:lines.append(f'| {name} | {a} to {(pd.Timestamp(b)-pd.Timedelta(days=1)).date()} | {(pd.Timestamp(b)-pd.Timestamp(a)).days} days |')
    lines+=['','Blanca active calendar dates are 31 May–9 June 2015, as documented in the [NHC report](https://www.nhc.noaa.gov/data/tcr/EP022015_Blanca.pdf). Full UTC calendar days are used, not exact genesis/dissipation times. This red period is a regional cyclone calendar, not verified severe local exposure.','', '**Change from the preceding charts:** earlier charts used the June 8 landfall/proximity reference and compared the three days around it. This chart uses the entire cyclone calendar. June 5–7 therefore belong to the red active-period category here, although they appeared yellow before the June 8 reference in the earlier charts. Neither chart establishes local onset.','', '## Meaning of quiet','', 'The first and last periods pass an explicit exploratory proxy: every checked hour has modelled 10-m wind <10 m/s, significant wave height <2 m, complete values, and no contemporaneous >=34-knot eastern Pacific cyclone track point within 500 km of a modelled shark position. The cyclone screen uses the prior downloaded IBTrACS MAIN tracks, matching to hourly positions within 30 minutes. Thresholds describe these candidate control periods and are not a validated ecological definition of normality. Dates were specified from the requested chronology and weather summaries, before computing their behavioural means.','', 'This cannot verify that nothing happened: tides, prey, fishing, other disturbances, weak/faraway systems and model errors are not excluded. GPE3 shark positions and weather are model-derived. The later period also follows the previously identified June 9 wave-source discontinuity, whose cause remains unresolved; lower later model waves therefore require caution. These are candidate quiet periods, not validated storm-free controls.','', 'Andres and Carlos occur elsewhere in the same eastern Pacific season; their contemporaneous tracks are explicitly included in cyclone_screen.csv. Distance >500 km does not eliminate swell or regional influence. The storm-period mean combines ten days with changing conditions and can hide a short response. Use the daily timeline alongside the stage-average bars.','', '## Calculations','', 'Depth is the average of hourly median depths; vertical movement is the average hourly absolute change between adjacent one-minute median depths (m/min). It is a pressure-tag movement proxy. Each three-day period has 72 complete hours and 17,280 sensor samples; the ten-day calendar period has 240 complete hours and 57,600 samples. Periods are disjoint. Summary means were checked against the means of their complete daily values. No observed shift is automatically attributed to Blanca.','', '## Files','', '- 01_five_stage_comparison.png / .pdf: five-stage averages.','- 02_daily_five_stage_timeline.png / .pdf: daily variation with stage shading.','- five_stage_summary.csv: stage means, exact time bounds and environmental checks.','- daily_values.csv: one row per day.','- cyclone_screen.csv: named concurrent cyclone proximity checks.','', 'Inputs are the existing hourly_raw.csv, environment_matched_hourly.csv and the prior archived cyclone tracks. Existing statistical reports and files were not changed.']
    notes='\n'.join(lines).replace('matching to hourly positions within 30 minutes','matching directly to original modelled positions within 90 minutes, without hourly position carry-forward')
    (OUT/'README.md').write_text(notes,encoding='utf-8')
    (OUT/'validation.json').write_text(json.dumps({'complete_days':len(dd),'total_hourly_records':int(sum(d.hours)),'first_and_last_quiet_proxy_pass':True,'period_and_daily_means_agree':True,'new_significance_tests':False},indent=2))
    print(d.drop(columns=['colour']).to_string(index=False));print(pd.DataFrame(screen).to_string(index=False))

if __name__=='__main__':main()
