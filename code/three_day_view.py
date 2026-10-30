"""Additional descriptive view; preserves the original inferential analysis."""
from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'three_day_comparison';OUT.mkdir(exist_ok=True)
REF=pd.Timestamp('2015-06-08 12:00',tz='UTC')
YELLOW='#d9ae43';GREEN='#4d9585';GREY='#9babb5';INK='#233846'

def main():
    y=pd.read_csv(ROOT/'results/hourly_raw.csv',parse_dates=['time']).set_index('time')
    e=pd.read_csv(ROOT/'results/environment_matched_hourly.csv',parse_dates=['time']).set_index('time')
    periods=[('Baseline (11 days)',-336,-72)]+[(label,a,a+24) for label,a in [('3rd day before',-72),('2nd day before',-48),('Last day before',-24),('1st day after',0),('2nd day after',24),('3rd day after',48)]]
    rows=[]
    for label,a,b in periods:
        start=REF+pd.Timedelta(hours=a);end=REF+pd.Timedelta(hours=b)
        t=y[(y.index>=start)&(y.index<end)]
        assert len(t)==b-a and t.n.eq(240).all()
        assert t[['median_depth_m','vertical_rate_m_min','depth_iqr_m']].notna().all().all()
        rows.append(dict(period=label,start_utc=start,end_utc_exclusive=end,hours=len(t),depth_m=t.median_depth_m.mean(),vertical_rate_m_min=t.vertical_rate_m_min.mean(),depth_iqr_m=t.depth_iqr_m.mean(),near_surface_percent=100*t.surface_fraction.mean()))
    d=pd.DataFrame(rows);d.to_csv(OUT/'daily_comparison.csv',index=False)
    combined=[]
    for label,a,b in [('11-day baseline',-336,-72),('72 hours before',-72,0),('72 hours after',0,72)]:
        t=y[(y.index>=REF+pd.Timedelta(hours=a))&(y.index<REF+pd.Timedelta(hours=b))]
        combined.append(dict(period=label,hours=len(t),depth_m=t.median_depth_m.mean(),vertical_rate_m_min=t.vertical_rate_m_min.mean()))
    pd.DataFrame(combined).to_csv(OUT/'combined_comparison.csv',index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.spines.top':False,'axes.spines.right':False,'axes.spines.left':False,'text.color':INK,'axes.labelcolor':INK})
    fig,axs=plt.subplots(2,1,figsize=(12,8.2),layout='constrained',sharex=True)
    fig.suptitle('Shark 14_07: baseline versus 3 days before and after',fontsize=19,fontweight='bold')
    x=np.array([-1.6,0,1,2,3,4,5]);colors=[GREY]+[YELLOW]*3+[GREEN]*3
    for ax,column,title,unit,fmt in [(axs[0],'depth_m','DEPTH: taller bar = deeper water','Average hourly median depth (m)','{:.2f}'),(axs[1],'vertical_rate_m_min','UP-AND-DOWN MOVEMENT: taller bar = more vertical movement','Vertical movement proxy (m/min)','{:.3f}')]:
        ax.bar(x,d[column],color=colors,width=.66,zorder=3)
        ax.axhline(d[column].iloc[0],color='#667782',ls='--',lw=1.3,zorder=2)
        ax.axvline(2.5,color='#ab493b',ls=':',lw=1.7)
        ax.set_ylim(0,d[column].max()*1.38);ax.set_ylabel(unit,fontsize=11)
        ax.set_title(title,loc='left',fontsize=12,pad=12,fontweight='bold');ax.grid(axis='y',alpha=.15,zorder=0)
        for xi,v in zip(x,d[column]):ax.text(xi,v+d[column].max()*.025,fmt.format(v),ha='center',va='bottom',fontsize=12,fontweight='bold')
    axs[0].text(.48,.94,'YELLOW = BEFORE',transform=axs[0].transAxes,ha='center',color='#956e0d',fontsize=11)
    axs[0].text(.83,.94,'GREEN = AFTER',transform=axs[0].transAxes,ha='center',color='#286d60',fontsize=11)
    axs[0].text(.07,.94,'Dashed horizontal line = baseline',transform=axs[0].transAxes,fontsize=10)
    axs[1].set_xticks(x,['Baseline*\n11 days','3rd day\nbefore','2nd day\nbefore','Last day\nbefore','1st day\nafter','2nd day\nafter','3rd day\nafter'])
    axs[1].set_xlabel('Each coloured bar = 24 hours. Red dotted divider = 8 June 2015, 12:00 UTC.\n*Baseline is a comparison period, not verified storm-free days. Reference is not verified local storm onset.',fontsize=10,labelpad=16)
    fig.savefig(OUT/'01_simple_daily_comparison.png',dpi=170,bbox_inches='tight');fig.savefig(OUT/'01_simple_daily_comparison.pdf',bbox_inches='tight');plt.close(fig)
    # Complementary weather context for exactly the same six days.
    w=e[(e.index>=REF-pd.Timedelta(hours=72))&(e.index<REF+pd.Timedelta(hours=72))]
    fig,axs=plt.subplots(3,1,figsize=(12,8),layout='constrained',sharex=True)
    fig.suptitle('Modelled weather during the same 3 days before and after',fontsize=17,fontweight='bold')
    for ax,col,label in [(axs[0],'wind_speed_10m','Wind at 10 m (m/s)'),(axs[1],'Thgt','Significant wave height (m)'),(axs[2],'pressure_msl','Sea-level pressure (hPa)')]:
        ax.plot(w.index,w[col],color='#286581',lw=1.8)
        ax.axvspan(REF-pd.Timedelta(hours=72),REF,color=YELLOW,alpha=.16)
        ax.axvspan(REF,REF+pd.Timedelta(hours=72),color=GREEN,alpha=.12)
        ax.axvline(REF,color='#ab493b',ls='--',lw=1.5);ax.set_ylabel(label,fontsize=11);ax.grid(axis='y',alpha=.15)
    axs[1].annotate('Flagged source discontinuity\n(not an observed sea-state collapse)',xy=(pd.Timestamp('2015-06-09 21:00',tz='UTC'),0),xytext=(pd.Timestamp('2015-06-08 18:00',tz='UTC'),.7),fontsize=9,arrowprops={'arrowstyle':'->','color':'#ab493b'},color='#ab493b')
    axs[0].set_title('Yellow: preceding 72 h | Green: following 72 h | Dashed line: reference time',loc='left',fontsize=11)
    axs[-1].xaxis.set_major_locator(mdates.DayLocator());axs[-1].xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
    axs[-1].set_xlabel('2015, UTC. Coarse weather models; local severe storm onset was not established.',fontsize=10)
    fig.savefig(OUT/'02_six_day_weather.png',dpi=170,bbox_inches='tight');fig.savefig(OUT/'02_six_day_weather.pdf',bbox_inches='tight');plt.close(fig)
    lines=['# Three-day descriptive comparison','', 'Additional exploratory visualization. No new significance tests were run. The original analysis remains unchanged.','', 'Reference: 8 June 2015 at 12:00 UTC. Each daily interval starts and ends at 12:00 UTC, rather than midnight. Baseline: 25 May at 12:00 through 5 June at 12:00. These comparison days are not certified storm-free.','', '| Period | Depth (m) | Vertical movement proxy (m/min) |','|---|---:|---:|']
    for r in rows:lines.append(f"| {r['period']} | {r['depth_m']:.2f} | {r['vertical_rate_m_min']:.3f} |")
    lines+=['','Depth is the average of hourly median depths; the movement proxy is the hourly average of absolute changes between adjacent one-minute median depths. Taller depth bars mean deeper water, unlike the inverted depth axis in the earlier time-series chart. Dashed horizontal lines mark baseline averages.','', 'The movement increase is concentrated in the final day before the reference. A wider view describes persistence and recovery but does not determine whether the storm caused the change or whether it began before local environmental cues. Severe local exposure remains unverified; the existing sensitivity and multiple-testing limitations still apply.','', 'Weather panel: the abrupt wave-model change at 9 June 21:00 is a previously flagged source discontinuity. It falls inside this extended after-period, so interpreting the extended wave comparison requires particular caution. Source data have not been corrected or replaced.','', 'Files: daily_comparison.csv gives exact interval boundaries; combined_comparison.csv gives the 72-hour summaries. Both charts are exported as PNG and vector PDF. All hourly data and models are inherited from the completed analysis.']
    (OUT/'README.md').write_text('\n'.join(lines),encoding='utf-8')
    print(d[['period','depth_m','vertical_rate_m_min']].to_string(index=False));print(pd.DataFrame(combined).to_string(index=False))

if __name__=='__main__':main()
