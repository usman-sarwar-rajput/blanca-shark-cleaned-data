"""Generate transparent report and exportable scientific figures from computed tables."""
from acquire import ROOT,RESULT
from analyse import REF
import pandas as pd,numpy as np,json,html,sys
sys.path.insert(0,str(ROOT/'dependencies'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','savefig.facecolor':'white'})

def save(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight');plt.close(fig)

def table(d,columns,formats=None):
    formats=formats or {};lines=['| '+' | '.join(columns.values())+' |','| '+' | '.join(['---']*len(columns))+' |']
    for r in d.to_dict('records'):
        v=[]
        for k in columns:
            value=r[k]
            if k in formats:value=formats[k].format(value)
            v.append(str(value))
        lines.append('| '+' | '.join(v)+' |')
    return '\n'.join(lines)

def main():
    y=pd.read_csv(RESULT/'hourly_raw.csv',parse_dates=['time']).set_index('time')
    e=pd.read_csv(RESULT/'environment_matched_hourly.csv',parse_dates=['time']).set_index('time')
    es=pd.read_csv(RESULT/'environment_summary.csv');m=pd.read_csv(RESULT/'main_effects.csv');s=pd.read_csv(RESULT/'sensitivity_effects.csv')
    ranks=pd.read_csv(RESULT/'contextual_ranks.csv');q=json.loads((RESULT/'depth_quality.json').read_text());checks=json.loads((RESULT/'validation.json').read_text())
    plot=y.loc['2015-05-25':'2015-06-10'];weather=e.loc[plot.index.min():plot.index.max()]
    fig,axs=plt.subplots(6,1,figsize=(11,13),sharex=True,layout='constrained')
    axs[0].plot(weather.index,weather.wind_speed_10m,color='#256582');axs[0].axhline(17.5,color='#a64b3c',ls=':',lw=1);axs[0].set_ylabel('Wind at 10 m\n(m/s)');axs[0].set_ylim(0,19)
    axs[1].plot(weather.index,weather.Thgt,color='#256582');axs[1].axhline(4,color='#a64b3c',ls=':',lw=1);axs[1].set_ylabel('Significant wave\nheight (m)');axs[1].set_ylim(0,4.5)
    axs[1].annotate('Source discontinuity',xy=(pd.Timestamp('2015-06-09 21:00',tz='UTC'),0),xytext=(pd.Timestamp('2015-06-06 00:00',tz='UTC'),3.1),fontsize=9,color='#a64b3c',arrowprops={'arrowstyle':'->','color':'#a64b3c'})
    axs[2].plot(weather.index,weather.pressure_msl,color='#256582');axs[2].set_ylabel('Mean sea-level\npressure (hPa)')
    axs[3].plot(plot.index,plot.median_depth_m,color='#1d3834',lw=1);axs[3].set_ylabel('Hourly median\ndepth (m)');axs[3].invert_yaxis()
    axs[4].plot(plot.index,plot.depth_iqr_m,color='#1d3834',lw=1);axs[4].set_ylabel('Hourly depth\nIQR (m)')
    axs[5].plot(plot.index,plot.vertical_rate_m_min,color='#1d3834',lw=1);axs[5].set_ylabel('Vertical-rate proxy\n(m/min)')
    for ax in axs:
        ax.axvspan(REF-pd.Timedelta(hours=48),REF,color='#d6a846',alpha=.15)
        ax.axvspan(REF,REF+pd.Timedelta(hours=24),color='#78a59e',alpha=.16)
        ax.axvline(REF,color='#a64b3c',ls='--',lw=1);ax.grid(axis='y',alpha=.15)
    axs[0].set_title('Blanca: modelled environment and recorded shark diving',loc='left',fontsize=14,pad=32)
    axs[0].text(0,1.02,'Shark 14_07 | Dashed vertical line: audit reference, not verified storm onset',transform=axs[0].transAxes,fontsize=10)
    axs[-1].xaxis.set_major_locator(mdates.DayLocator(interval=2));axs[-1].xaxis.set_major_formatter(mdates.DateFormatter('%d %b'));axs[-1].set_xlabel('2015, UTC | Yellow: preceding 48 h; green: following 24 h | Dotted: screening thresholds')
    save(fig,'01_environment_and_diving')
    fig,axs=plt.subplots(1,3,figsize=(12,4.6),layout='constrained')
    names={'median_depth_m':'Hourly median depth (m)','depth_iqr_m':'Depth IQR (m)','vertical_rate_m_min':'Vertical-rate proxy (m/min)'}
    for ax,(outcome,title) in zip(axs,names.items()):
        for offset,data,color,label in [(-.10,m,'#256582','Linear trend'),(.10,s[s.variant.eq('no_trend')],'#b76440','No trend')]:
            rows=data[data.outcome.eq(outcome)].set_index('phase').loc[['pre48_24','pre24','post24']]
            positions=np.arange(3)+offset
            ax.errorbar(rows.adjusted_difference,positions,xerr=[rows.adjusted_difference-rows.ci_low,rows.ci_high-rows.adjusted_difference],fmt='o',ms=5,capsize=3,color=color,label=label)
        ax.axvline(0,color='#666',lw=1);ax.set_yticks(np.arange(3),['48–24 h before','24 h before','24 h after']);ax.invert_yaxis();ax.set_title(title,fontsize=11);ax.set_xlabel('Difference and nominal 95% HAC interval')
    handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=2,fontsize=9)
    fig.suptitle('Sensitivity to the assumed baseline trend',fontsize=14)
    save(fig,'02_effects_and_trend_sensitivity')
    p=pd.read_csv(RESULT/'contextual_pseudo_references.csv')
    fig,axs=plt.subplots(1,3,figsize=(12,3.8),layout='constrained')
    for ax,phase,label in zip(axs,['pre48_24','pre24','post24'],['48–24 h before','24 h before','24 h after']):
        values=p[p.phase.eq(phase)].matched_depth_difference_m
        v=float(ranks[ranks.phase.eq(phase)].event_matched_difference_m.iloc[0])
        ax.hist(values,bins=16,color='#a9bcb9',edgecolor='white');ax.axvline(v,color='#a64b3c',lw=2,label='Blanca reference');ax.set_title(label);ax.set_xlabel('UTC-hour-matched depth difference (m)');ax.set_ylabel('Contextual windows');ax.legend(fontsize=8)
    fig.suptitle('Ordinary-period context: overlapping windows, not independent controls',fontsize=13)
    save(fig,'03_contextual_variation')
    # Export a compact descriptive table including explicitly secondary measurements.
    periods={'baseline':(-336,-72),'pre48_24':(-48,-24),'pre24':(-24,0),'post24':(0,24)};desc=[]
    for name,(a,b) in periods.items():
        t=y[(y.index>=REF+pd.Timedelta(hours=a))&(y.index<REF+pd.Timedelta(hours=b))]
        d={'period':name,'hours':len(t),'median_depth_hourly_mean_m':t.median_depth_m.mean(),'depth_iqr_hourly_mean_m':t.depth_iqr_m.mean(),'vertical_rate_mean_m_min':t.vertical_rate_m_min.mean(),'near_surface_percent':100*t.surface_fraction.mean(),'median_in_situ_temperature_c':t.temperature_c.median()};desc.append(d)
    desc=pd.DataFrame(desc);desc.to_csv(RESULT/'descriptive_periods.csv',index=False)
    environment=es[es.scenario.eq('trajectory')].iloc[0]
    significant=m[m.p_holm<.05]
    main_depth=m[m.outcome.eq('median_depth_m')].iloc[0];no=s[(s.variant.eq('no_trend'))&(s.outcome.eq('median_depth_m'))&(s.phase.eq('pre48_24'))].iloc[0]
    text=f'''# Blanca and juvenile white shark diving: exploratory statistical report

Completed 16 September 2026. Animal 14_07, PAT 14_07-132980. This is a reproducible single-animal case analysis, not a submission-ready manuscript.

## Finding

**Diving varied around Blanca, but the results do not provide robust evidence that this shark anticipated severe local storm conditions.** Modelled local wind and waves did not cross the environmental screening thresholds. One depth comparison passed the main model's multiple-testing adjustment, but the finding depended on assuming a continuing linear depth trend and did not retain significance in the prespecified no-trend sensitivity analysis.

The work now includes actual behavioural calculations, local environmental matching, quality checks, uncertainty inspection, sensitivity analyses and exportable figures. This goes beyond the previous availability audit.

## Data and recording quality

The unchanged original ZIP is `../Shark_Storm_Overlap_Audit_2026-09-16/sources/tags/14_07-132980.zip`, member `out-Archive.csv`. [The original shark archive](https://doi.org/10.24431/rw1k6c3) and its [peer-reviewed descriptor](https://pmc.ncbi.nlm.nih.gov/articles/PMC8976051/) supply provenance. These are recovered pressure-derived depth and temperature observations; locations are GPE3 estimates rather than GPS.

The source has {q['raw_rows']:,} rows. We excluded {q['nonfinite_depth_rows_excluded']:,} rows without finite depth and removed one identical finite duplicate at the August detachment time. There were no remaining duplicated measurement timestamps. Restricting to the conservative deployment interior, 7 November 2014 through 2 August 2015, leaves **{q['interior_records']:,} measurements at uninterrupted 15-second intervals**. The 360-hour modelling interval contains 86,400 depth observations. Raw readings after the deployment interior were excluded.

Depth ranges from {q['minimum_depth_m']:.1f} to {q['maximum_depth_m']:.1f} m over the retained deployment. {q['negative_depth_fraction']*100:.2f}% of retained readings are negative; in the exported May 25–June 9 subset the fraction is {q['event_negative_depth_fraction']*100:.3f}%. Surface offset remains a calibration limitation. Clipping negative depths to zero did not change the principal median-depth result. The supplied clock-correction rows have identical tag and real times; no additional clock adjustment was inferred. The release date is months after Blanca. Original metadata also contain a release-location comment inconsistent with the first modelled coordinate; this should be clarified with the data custodians before publication, although it does not interrupt the June sensor record.

## Local environmental assessment

[ERA5 via Open-Meteo](https://open-meteo.com/en/docs/historical-weather-api) supplies hourly 10-m wind and sea-level pressure, with the ERA5 model explicitly selected and sea-grid preference. [NOAA AOML's preserved 2015 PacIOOS WW3 archive](https://erddap.aoml.noaa.gov/hdb/erddap/griddap/WaveWatch_2015.html) supplies significant wave height. WW3 is a model product forced by GFS winds; it is not an observed buoy series. All request URLs, retrieval dates, response bytes and hashes are stored with the raw files.

Hourly weather was paired with the nearest recorded GPE3 location within six hours, without interpolating a shark track. Wave cells were selected from persistently wet cells using availability, independently of wave magnitude. All 144 hours from June 5 through June 10 have matched wind, pressure and wave data. The event's maximum grid offsets are **{environment.max_era5_offset_km:.1f} km for ERA5 and {environment.max_wave_offset_km:.1f} km for WW3**. These offsets and coarse coastal grids limit local interpretation.

| Quantity, June 5–10 | Matched model estimate |
|---|---|
| Maximum wind at 10 m | {environment.max_wind_m_s:.2f} m/s |
| Wind maximum time | {environment.peak_wind_utc} |
| Maximum significant wave height | {environment.max_wave_m:.2f} m |
| Wave maximum time | {environment.peak_wave_utc} |
| Minimum sea-level pressure | {environment.minimum_pressure_hpa:.1f} hPa |

The protocol's operational screen required wind at least 17.5 m/s or significant wave height at least 4 m for three consecutive hours. **No model-based onset was identified.** The wave-height sensitivities of 3 and 5 m also produced no onset. Fixed +/-0.5-degree scenarios had maximum winds of {es.max_wind_m_s.min():.2f}–{es.max_wind_m_s.max():.2f} m/s and maximum waves of {es.max_wave_m.min():.2f}–{es.max_wave_m.max():.2f} m; none met the main rule. These thresholds are screening choices, not validated biological thresholds. Failure to cross them does not establish calm conditions, absence of a response to weaker cues, or absence of unmodelled extremes.

Visual quality review identified an abrupt wave-source drop from 1.30 m to 0.00 m at June 9, 21:00 UTC in the same grid cell, followed by low values. It may represent a model/archive discontinuity; its cause has not been verified. The source values are retained and flagged, not interpreted as an observed collapse of the sea state. This occurs nine hours after the analysed post-reference window ends, so it does not affect those behavioural comparisons or the reported earlier wave maximum. It further limits confidence in the extended weather series; nonmissing model values are not automatically valid observations.

The [official NHC Blanca report](https://www.nhc.noaa.gov/data/tcr/EP022015_Blanca.pdf) confirms regional storm conditions and weakening around Baja California landfall on June 8. Its land-station observations are regional context; no collocated observational series was obtained for this shark. Severe local exposure remains unverified.

## Location uncertainty

The supplied GPE3 NetCDF contains twelve-hour location likelihood grids. After normalizing the nonnegative weights, the highest-weight 0.25-degree cell contains at least 95% of model mass at all June 5–10 twelve-hour steps. At the audit reference, the cell centred at 26.875 N, 113.125 W contains approximately 97.4%. This concentration is **conditional on the original model and its assumptions**, not externally validated positional accuracy; within-cell location and model misspecification remain unresolved. GPE3 uses depth, temperature and bathymetry, so reconstructed location is not fully independent of the biologging measurements. [Manufacturer description](https://wildlifecomputers.com/data/technologies/geolocation/).

## Statistical design

The reference remains **8 June 2015 at 12:00 UTC**, fixed by the preceding proximity audit. It is not local onset. The baseline is May 25 at 12:00 through June 5 at 12:00 (264 hours). A 24-hour guard interval follows. We compare June 6 at 12:00–June 7 at 12:00, June 7 at 12:00–June 8 at 12:00, and June 8 at 12:00–June 9 at 12:00.

The primary endpoint is hourly median depth. Secondary endpoints are hourly depth IQR and a vertical-rate proxy: mean absolute change between adjacent one-minute median depths, in m/min. This measures vertical displacement recorded by the tag, not horizontal swimming speed, energy use, awareness or formally defined dive frequency. Temperature and time at <=5 m are descriptive outputs only. Tag temperature varies with the occupied depth and is not independent sea-surface temperature.

OLS models use UTC-hour effects, a linear time trend and period indicators, including a nuisance guard indicator. All 360 hours are retained in sequence. Newey–West HAC errors use lag 24 with a finite-sample correction; nominal intervals use the normal approximation. [HAC method implementation](https://www.statsmodels.org/stable/generated/statsmodels.stats.sandwich_covariance.cov_hac.html). We apply [Holm adjustment](https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html) across all nine endpoint-by-period comparisons. Hourly records from one animal are temporally dependent and cannot be counted as independent animals or replicated storms. Intervals are conditional on this fitted time-series model and may be optimistic for a single short response window.

## Observed behaviour

Values below are averages of the hourly summaries, not the median of all pooled readings.

{table(desc,{'period':'Period','hours':'Hours','median_depth_hourly_mean_m':'Hourly median depth, mean (m)','depth_iqr_hourly_mean_m':'Depth IQR, mean (m)','vertical_rate_mean_m_min':'Vertical rate (m/min)'},{'median_depth_hourly_mean_m':'{:.2f}','depth_iqr_hourly_mean_m':'{:.2f}','vertical_rate_mean_m_min':'{:.3f}'})}

The raw period summaries show a shallower 48–24-hour period and a higher vertical-rate proxy in the last 24 hours before the reference. These are descriptions of the recordings; a storm-related explanation requires stronger evidence.

## Adjusted comparisons

Negative depth differences mean shallower than the model's continuation of the baseline pattern. Nominal 95% intervals are not simultaneous intervals; use the Holm column when assessing the nine-test family.

{table(m,{'outcome':'Endpoint','phase':'Period','adjusted_difference':'Adjusted difference','ci_low':'95% lower','ci_high':'95% upper','p_holm':'Holm p'},{'adjusted_difference':'{:.3f}','ci_low':'{:.3f}','ci_high':'{:.3f}','p_holm':'{:.4f}'})}

Only the 48–24-hour median-depth contrast passes Holm adjustment in the main model: **{main_depth.adjusted_difference:.2f} m**, nominal 95% interval [{main_depth.ci_low:.2f}, {main_depth.ci_high:.2f}], Holm p={main_depth.p_holm:.6f}. The assumed increasing baseline trend makes the fitted counterfactual deeper than the observed baseline average. It therefore explains why the adjusted estimate differs substantially from the raw period difference.

## Robustness and ordinary variation

Removing the linear trend changes that same estimate to **{no.adjusted_difference:.2f} m**, nominal interval [{no.ci_low:.2f}, {no.ci_high:.2f}], Holm p={no.p_holm_within_variant:.3f}. It no longer passes the family adjustment. The immediate pre-reference depth contrast also changes sign when the trend is removed. HAC lags 12 and 48 and clipping negative readings do not resolve this model dependence. The full sensitivity file includes five-minute movement aggregation; tests are adjusted within each sensitivity variant and are not independent confirmations.

There are 64 eligible contextual pseudo-reference dates. Their windows overlap, and other weather and habitat changes were not excluded, so they are **not certified non-storm controls or independent replication**. For the 48–24-hour UTC-hour-matched depth shift, {int(ranks.iloc[0].context_abs_ge_event)} of 64 contextual windows have an equally large or larger absolute shift. The corresponding counts for the immediate preceding and following day are {int(ranks.iloc[1].context_abs_ge_event)} and {int(ranks.iloc[2].context_abs_ge_event)}. These are descriptive ranks, not p-values; habitat changes can make distant periods poorly comparable.

## Scientific conclusion and publication implications

This selected shark provides a useful exploratory time series, but **it does not currently support a claim of reliable pre-disaster detection**. Local severe conditions were not established, the main depth finding depends on baseline-trend assumptions, and the other prespecified main comparisons do not survive multiple-testing adjustment. A response to weaker contemporaneous wind, pressure or waves remains possible, but this analysis does not identify its cause or timing. A before-peak change is not necessarily advance sensing: environmental cues can already be present before the peak.

These findings could form an exploratory case-study component or feasibility/limitations section, subject to expert ecological and statistical review. They do not establish a repeatable result across sharks or storms. Stronger inference needs independent exposed events, credible local observations or validated higher-resolution environmental reconstruction, better matched habitat controls, and explicit sensor-calibration documentation. A confirmatory design should be fixed before examining new event outcomes. No synthetic observations, ML classifier or forecasting accuracy was used.

## Figures

![Environment and diving](figures/01_environment_and_diving.png)

Figure 1. Recorded hourly depth summaries alongside modelled environment. The vertical line is the audit reference; depth increases downward. Dotted horizontal lines show the operational environmental thresholds. Grid offsets and location uncertainty apply throughout.

![Trend sensitivity](figures/02_effects_and_trend_sensitivity.png)

Figure 2. Model-based period differences with and without a linear trend. Error bars are nominal 95% HAC intervals; family-adjusted p-values are in the tables.

![Contextual variation](figures/03_contextual_variation.png)

Figure 3. Descriptive contextual-window distributions. Overlapping windows and changing habitat prevent interpreting their frequency as independent-event evidence.

## Reproduction, files and validation

Run the Python scripts in this order: `acquire.py`, `retry_weather.py`, `depth.py`, `uncertainty.py`, `environment.py`, `analyse.py`, `validate.py`, `report.py`. Execute from any directory using each script's path. The scripts resolve this folder automatically; the original shark ZIP and prior audit tables are expected in the sibling audit folder. For offline replication, preserve the downloaded snapshots. The validator's previous-folder check is workspace-specific. Packages: numpy, pandas, scipy, statsmodels, matplotlib, h5py and Markdown (the last two installed only in this analysis folder). Exact versions are pinned in `requirements.txt` and analysis versions are also in `results/software_versions.json`.

Main reading files: `results/main_effects.csv`, `results/descriptive_periods.csv`, `results/sensitivity_effects.csv`, `results/environment_summary.csv`, `results/depth_quality.json`, `results/location_uncertainty.csv`, and the figures exported as both PNG and vector PDF. CSV tables open in Excel. The modelled-position matching details and all requested/returned weather coordinates are retained. `SOURCE_MANIFEST.json` records environmental source checksums; the raw ZIP hash links to the prior publisher-verified archive. Retain source author attribution and applicable licences when sharing. Failed acquisition attempts are recorded; the initial combined weather request failed with HTTP 502, after which all 19 individual site requests succeeded.

Validation passed: {checks['independent_ols_and_hac_comparisons']} independent OLS/HAC coefficient and standard-error checks; nine-test Holm verification; synthetic aggregation and missing-data checks; {checks['event_hour_independent_aggregation_checks']} independently aggregated event hours; 144 complete environmental event hours; source-hash checks; and preservation of all {checks['protected_prior_files_unchanged']} previously protected files. These validate data handling and calculations, not the causal interpretation. The SmartCodex advisory checker was unavailable because its yaml dependency was missing; no automated complexity or credit figure was invented.
'''
    (ROOT/'RESEARCH_REPORT.md').write_text(text,encoding='utf-8')
    # A lightweight self-contained reading page; local images remain beside the report.
    try:
        import markdown
        body=markdown.markdown(text,extensions=['tables'])
    except ImportError:
        body='<pre style="white-space:pre-wrap">'+html.escape(text)+'</pre>'
    css='body{font-family:Georgia,serif;max-width:1000px;margin:45px auto;padding:0 28px;color:#172d34;line-height:1.55}h1{font-size:30px}h2{margin-top:2em;font-size:22px}table{border-collapse:collapse;width:100%;font:13px Arial,sans-serif}th,td{padding:8px;border-bottom:1px solid #ccd7d8;text-align:left}th{background:#edf3f3}img{max-width:100%}code{overflow-wrap:anywhere}@media print{body{margin:0;font-size:10pt}h2{break-after:avoid}tr,img{break-inside:avoid}}'
    (ROOT/'RESEARCH_REPORT.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Blanca shark exploratory analysis</title><style>'+css+'</style><body>'+body+'</body></html>',encoding='utf-8')
    print('Report and three PNG/PDF figures generated',flush=True)

if __name__=='__main__':main()
