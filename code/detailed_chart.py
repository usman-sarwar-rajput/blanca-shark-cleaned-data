"""Detailed five-stage figure, with explicit observational/model distinctions."""
from pathlib import Path
import json,hashlib
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'detailed_five_stage_chart';OUT.mkdir(exist_ok=True)
START=pd.Timestamp('2015-05-25',tz='UTC');END=pd.Timestamp('2015-06-16',tz='UTC')
REF=pd.Timestamp('2015-06-08 12:00',tz='UTC');FLAG=pd.Timestamp('2015-06-09 21:00',tz='UTC')

def main():
    sources=[ROOT/'results/hourly_raw.csv',ROOT/'results/environment_matched_hourly.csv',ROOT/'five_stage_comparison/five_stage_summary.csv']
    y=pd.read_csv(sources[0],parse_dates=['time']).set_index('time').loc[START:END-pd.Timedelta(hours=1)]
    e=pd.read_csv(sources[1],parse_dates=['time']).set_index('time').loc[START:END-pd.Timedelta(hours=1)]
    stages=pd.read_csv(sources[2],parse_dates=['start_utc','end_utc_exclusive'])
    assert len(y)==528 and len(e)==528 and y.index.equals(e.index)
    assert y.n.eq(240).all() and not y[['median_depth_m','vertical_rate_m_min']].isna().any().any()
    assert not e[['wind_speed_10m','pressure_msl','Thgt']].isna().any().any()
    y.to_csv(OUT/'plotted_hourly_shark_data.csv',index_label='time');e.to_csv(OUT/'plotted_hourly_weather.csv',index_label='time')
    daily=pd.DataFrame(index=y.resample('D').size().index)
    for col in ['median_depth_m','vertical_rate_m_min']:
        g=y[col].resample('D');daily[col+'_mean']=g.mean();daily[col+'_q25']=g.quantile(.25);daily[col+'_q75']=g.quantile(.75)
    daily.to_csv(OUT/'plotted_daily_summaries.csv',index_label='date_utc')
    for _,r in stages.iterrows():
        v=y[(y.index>=r.start_utc)&(y.index<r.end_utc_exclusive)]
        assert len(v)==r.hours and np.isclose(v.median_depth_m.mean(),r.depth_m) and np.isclose(v.vertical_rate_m_min.mean(),r.vertical_rate_m_min)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'text.color':'#213a49','axes.labelcolor':'#213a49','svg.fonttype':'none'})
    fig=plt.figure(figsize=(13,18),layout='constrained')
    grid=fig.add_gridspec(9,1,height_ratios=[.88,1.65,1.65,1.25,1.35,1.25,.33,1.65,1.04])
    fig.suptitle('Juvenile white shark 14_07 around Hurricane Blanca (2015)\nFive-stage descriptive comparison of recorded diving and modelled weather',fontsize=18,fontweight='bold')
    key=fig.add_subplot(grid[0]);key.axis('off')
    labels=['1  EARLIER QUIET*\n25–27 May | 72 h','2  JUST BEFORE\n28–30 May | 72 h','3  BLANCA ACTIVE DATES\n31 May–9 Jun | 240 h','4  JUST AFTER\n10–12 Jun | 72 h','5  LATER QUIET*\n13–15 Jun | 72 h']
    for i,(_,r) in enumerate(stages.iterrows()):
        key.add_patch(Rectangle((i*.2+.004,.1),.192,.75,facecolor=r.colour,alpha=.25,edgecolor='none',transform=key.transAxes))
        key.text(i*.2+.1,.47,labels[i],transform=key.transAxes,ha='center',va='center',fontsize=10,fontweight='bold')
    axes=[]
    for row in range(1,6):
        ax=fig.add_subplot(grid[row],sharex=axes[0] if axes else None);axes.append(ax)
        for _,r in stages.iterrows():ax.axvspan(r.start_utc,r.end_utc_exclusive,color=r.colour,alpha=.12,zorder=0)
        ax.axvline(REF,color='#97443c',ls='--',lw=1.25,zorder=3)
        ax.set_xlim(START,END);ax.grid(axis='y',alpha=.18)
        if row<5:ax.tick_params(labelbottom=False)
    for ax,col,title,unit in [(axes[0],'median_depth_m','A  RECORDED DEPTH','Depth (m)\nHigher = deeper'),(axes[1],'vertical_rate_m_min','B  RECORDED VERTICAL-MOVEMENT PROXY','Vertical movement\n(m/min)')]:
        t=daily.index+pd.Timedelta(hours=12)
        ax.plot(y.index+pd.Timedelta(minutes=30),y[col],color='#718a98',alpha=.35,lw=.65,label='Hourly summary',zorder=1)
        ax.fill_between(t,daily[col+'_q25'],daily[col+'_q75'],color='#2d617d',alpha=.2,label='Middle 50% of hourly values within each day',zorder=2)
        ax.plot(t,daily[col+'_mean'],color='#174963',lw=2,marker='o',ms=3.8,label='Daily average of hourly summaries',zorder=4)
        ax.set_ylabel(unit);ax.set_title(title,loc='left',fontsize=12,fontweight='bold');ax.set_ylim(bottom=0)
    axes[0].legend(loc='upper left',ncol=1,fontsize=8.5,framealpha=.9)
    axes[1].text(.99,.93,'Shading shows variation, not a confidence interval',transform=axes[1].transAxes,ha='right',fontsize=9)
    axes[2].plot(e.index,e.wind_speed_10m,color='#286b73',lw=1.5)
    axes[2].axhline(10,color='#66717a',ls=':',lw=1);axes[2].axhline(17.5,color='#ad5144',ls=':',lw=1)
    axes[2].text(.01,.88,'17.5 m/s: severe-condition screening threshold',transform=axes[2].transAxes,color='#ad5144',fontsize=9)
    axes[2].set_ylim(0,20);axes[2].set_ylabel('Wind at 10 m\n(m/s)');axes[2].set_title('C  MODELLED WIND — ERA5',loc='left',fontsize=12,fontweight='bold')
    axes[3].plot(e.index[e.index<FLAG],e.loc[e.index<FLAG,'Thgt'],color='#286b73',lw=1.5)
    axes[3].plot(e.index[e.index>=FLAG],e.loc[e.index>=FLAG,'Thgt'],color='#777',ls='--',lw=1.3)
    axes[3].axvspan(FLAG,END,facecolor='none',edgecolor='#aaaaaa',hatch='///',lw=0,alpha=.22)
    axes[3].axhline(2,color='#66717a',ls=':',lw=1);axes[3].axhline(4,color='#ad5144',ls=':',lw=1)
    axes[3].text(.01,.88,'4 m: severe-condition screening threshold',transform=axes[3].transAxes,color='#ad5144',fontsize=9)
    axes[3].annotate('Source discontinuity\n9 Jun, 21:00 UTC',xy=(FLAG,0),xytext=(pd.Timestamp('2015-06-06',tz='UTC'),2.8),fontsize=9,color='#97443c',arrowprops={'arrowstyle':'->','color':'#97443c'})
    axes[3].set_ylim(0,4.7);axes[3].set_ylabel('Significant wave\nheight (m)');axes[3].set_title('D  MODELLED WAVES — WW3; HATCHED TAIL REQUIRES CAUTION',loc='left',fontsize=12,fontweight='bold')
    axes[4].plot(e.index,e.pressure_msl,color='#286b73',lw=1.5);axes[4].set_ylabel('Sea-level pressure\n(hPa)');axes[4].set_title('E  MODELLED PRESSURE — ERA5',loc='left',fontsize=12,fontweight='bold')
    axes[4].xaxis.set_major_locator(mdates.DayLocator(interval=2));axes[4].xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
    axes[4].set_xlabel('Date in 2015 (UTC). Vertical dashed line: 8 June, 12:00 UTC — earlier audit reference, not local storm onset.',fontsize=10)
    spacer=fig.add_subplot(grid[6]);spacer.axis('off')
    ta=fig.add_subplot(grid[7]);ta.axis('off');ta.set_title('F  STAGE AVERAGES AND ENVIRONMENTAL CHECKS',loc='left',fontsize=12,fontweight='bold')
    cells=[]
    for i,(_,r) in enumerate(stages.iterrows()):cells.append([str(i+1),str(int(r.hours)),f'{r.depth_m:.2f}',f'{r.vertical_rate_m_min:.3f}',f'{r.near_surface_percent:.1f}',f'{r.max_model_wind_m_s:.2f}',f'{r.max_model_wave_m:.2f}'])
    tab=ta.table(cellText=cells,colLabels=['Stage','Hours','Depth¹\n(m)','Movement²\n(m/min)','At ≤5 m³\n(%)','Maximum wind\n(m/s)','Maximum waves\n(m)'],cellLoc='center',colLoc='center',bbox=[0,0,1,.91],colWidths=[.075,.075,.13,.16,.14,.21,.21])
    tab.auto_set_font_size(False);tab.set_fontsize(10)
    for (r,c),cell in tab.get_celld().items():
        cell.set_edgecolor('#d7e0e4');cell.set_linewidth(.5)
        if r==0:cell.set_facecolor('#e8eef0');cell.get_text().set_fontweight('bold')
        else:cell.set_facecolor(matplotlib.colors.to_rgba(stages.iloc[r-1].colour,.17))
    notes=fig.add_subplot(grid[8]);notes.axis('off')
    notes.text(0,1,'READING AND LIMITATIONS',fontweight='bold',fontsize=11,va='top')
    note=('¹ Average of hourly median depths. ² Mean absolute change between adjacent one-minute median depths.\n'
          '³ Fraction of 15-second readings at or above 5 m depth (depth ≤5 m). One shark; 528 complete hours; 126,720 readings.\n'
          '*Quiet means candidate lower-weather periods: model wind <10 m/s, waves <2 m, and no screened ≥34-knot\n'
          'cyclone within 500 km. Later waves follow an unresolved source discontinuity; neither quiet period is certified disturbance-free.\n'
          'The red stage covers Blanca’s calendar dates, not confirmed local severe exposure. Differences are descriptive, not causal.\n'
          'Sources: shark DOI 10.24431/rw1k6c3; NOAA IBTrACS/NHC; ERA5 via Open-Meteo; PacIOOS WW3 via NOAA AOML.')
    notes.text(0,.82,note,fontsize=9.5,va='top',linespacing=1.5)
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'DETAILED_FIVE_STAGE_CHART.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)
    provenance=[{'file':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources]
    (OUT/'validation_and_sources.json').write_text(json.dumps({'hours':528,'raw_measurements':126720,'days':22,'stages':5,'hourly_coverage_complete':True,'stage_means_verified':True,'daily_band':'25th–75th percentiles of hourly values; not confidence interval','sources':provenance},indent=2))
    explanation='''# Detailed five-stage chart: interpretation and justification

This figure shows one juvenile white shark across 22 complete UTC calendar days. It is a descriptive comparison, not a new confirmatory analysis.

## Why these five stages?

1. Earlier candidate quiet period: 25–27 May, 72 hours.
2. Immediately before Blanca's calendar period: 28–30 May, 72 hours.
3. Blanca's active calendar dates: 31 May–9 June, 240 hours.
4. Immediately after: 10–12 June, 72 hours.
5. Later candidate quiet period: 13–15 June, 72 hours.

The calendar boundaries follow the [official NHC report](https://www.nhc.noaa.gov/data/tcr/EP022015_Blanca.pdf), using whole days rather than exact genesis/dissipation times. Each flank is three days to match the user's comparison request. The active period is ten days; averages rather than totals prevent its longer duration mechanically inflating the measurements, but precision and within-period variability are not equal. These boundaries describe the regional cyclone; local shark exposure and local storm onset remain unverified. June 5–7 are during the active calendar period, even though they were before the June 8 reference in previous figures.

## Measurements and display

Panel A: grey hourly median pressure-derived depth, blue daily mean of those hourly medians. Higher on this axis means deeper water. Panel B: hourly mean absolute change between consecutive one-minute median depths, in metres per minute. It measures recorded vertical displacement, not swimming effort, horizontal speed or number of dives. The blue daily bands show the middle 50% of hourly summaries within a day. They are variability bands, not confidence intervals, and contain neither 50% of raw samples nor a range of all possible true effects.

Panels C–E show hourly modelled wind, significant wave height and sea-level pressure. They are not measurements taken alongside the shark. All times are UTC. Vertical dashed lines retain the previous audit reference (8 June 12:00) only as orientation; they do not mark independently confirmed local onset. Grey dotted lines in wind and wave panels identify the quiet-screen limits (10 m/s and 2 m); red dotted lines identify the earlier investigator-defined severe-condition limits (17.5 m/s and 4 m). These are screening choices, not demonstrated shark-response thresholds.

Panel F reports averages, sample-hour counts, near-surface percentage and environmental maxima for each stage. Near-surface time uses raw depth <=5 m, including small negative pressure offsets, and is descriptive. Negative-depth clipping did not alter the original primary median-depth result; calibration uncertainty remains.

## What does quiet mean?

Both candidate quiet periods have model wind <10 m/s and significant waves <2 m at every checked hour, with no screened >=34-knot eastern Pacific cyclone point within 500 km. This is an environmental screening rule, not proof that nothing happened or that the animal was in a normal physiological state. Tides, prey, habitat, human activity, other hazards and distant swell are not excluded. Andres and Carlos are included in the preceding cyclone screen and were distant in the selected flank periods. The original GPE3 positions were paired to cyclone points within 90 minutes; see ../five_stage_comparison/cyclone_screen.csv.

The ERA5 and WW3 grids are coarse. During June 5–10 the maximum grid offsets were approximately 37 km and 44 km, respectively; the shark positions themselves are model-derived. The WW3 series drops abruptly at 9 June 21:00. Its cause remains unresolved. The subsequent wave trace is dashed and hatched to prevent treating it as independently validated calm water. This affects confidence in later candidate quiet days. No source values were silently corrected or replaced.

## What can be concluded?

The daily view shows short-term variation that the ten-day average can hide. Movement also increases during the later candidate quiet period, so the observations are not a unique pre-Blanca signal. Stage averages do not establish causal effects, anticipation or recovery. Three-day selection differs from the original 11-day baseline and therefore gives different averages; neither should be substituted for the other without explanation.

The original model-dependent statistical findings and multiple-comparison limitations remain. No new p-values or significance stars were added. One shark is one biological subject: 126,720 measurements do not constitute 126,720 independent replicates. Replicated events, credible local exposure, matched habitat controls and validated sensor/geolocation uncertainty would be needed for stronger inference.

## Reproducibility and sources

The chart is generated by ../code/detailed_chart.py from existing hourly results. All 528 hourly records were present, with 240 raw depth samples each. The five phase averages were independently reconciled to those hourly values. Plotted hourly and daily data accompany the figure, and input SHA-256 hashes are saved in validation_and_sources.json. PNG is for viewing; PDF and SVG are vector exports for resizing. Existing analysis files were not changed.

- [Original shark archive](https://doi.org/10.24431/rw1k6c3)
- [Dataset descriptor](https://pmc.ncbi.nlm.nih.gov/articles/PMC8976051/)
- [ERA5 access and documentation](https://open-meteo.com/en/docs/historical-weather-api)
- [PacIOOS WW3 2015, NOAA AOML](https://erddap.aoml.noaa.gov/hdb/erddap/griddap/WaveWatch_2015.html)
- [NOAA IBTrACS](https://www.ncei.noaa.gov/products/international-best-track-archive)

The complete analysis report and original retrieval metadata remain in the parent folder. This figure is suitable for discussing the exploratory evidence and its limitations; it does not convert the case into proof of disaster detection.
'''
    (OUT/'CHART_EXPLANATION.md').write_text(explanation,encoding='utf-8')
    print(json.dumps({'chart':str(OUT/'DETAILED_FIVE_STAGE_CHART.png'),'hours_checked':len(y),'period_means_match':True}))

if __name__=='__main__':main()
