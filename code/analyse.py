"""Prespecified exploratory comparisons with time dependence and sensitivity checks."""
from acquire import ROOT,RESULT
import pandas as pd,numpy as np,json,sys,platform
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy.stats import norm

REF=pd.Timestamp('2015-06-08 12:00',tz='UTC')
OUTCOMES=['median_depth_m','depth_iqr_m','vertical_rate_m_min']
PHASES=['pre48_24','pre24','post24']

def phase(h):
    return np.select([h<-72,h<-48,h<-24,h<0],['baseline','guard','pre48_24','pre24'],default='post24')

def fit(y,ref,lag=24,trend=True,five=False,label='main'):
    t=y[(y.index>=ref-pd.Timedelta(hours=336))&(y.index<ref+pd.Timedelta(hours=24))].copy()
    assert len(t)==360 and (t.index.to_series().diff().dropna()==pd.Timedelta(hours=1)).all(),'HAC needs complete regular hours'
    if five:t['vertical_rate_m_min']=t['vertical_rate_5min_m_min']
    assert not t[OUTCOMES].isna().any().any(),'Missing endpoint hours: model needs review'
    t['phase']=phase((t.index-ref).total_seconds()/3600)
    X=pd.DataFrame({'const':np.ones(len(t))},index=t.index)
    if trend:X['trend_days']=(t.index-ref).total_seconds()/86400
    for p in ['guard']+PHASES:X[p]=(t.phase==p).astype(float)
    for h in range(1,24):X[f'utc_hour_{h}']=(t.index.hour==h).astype(float)
    assert np.linalg.matrix_rank(X)==X.shape[1]
    rows=[];pred=[];diagnostics=[]
    for outcome in OUTCOMES:
        model=sm.OLS(t[outcome],X).fit(cov_type='HAC',cov_kwds={'maxlags':lag,'use_correction':True},use_t=False)
        # Counterfactual includes trend and hour pattern but removes period indicators.
        xc=X.copy();xc[['guard']+PHASES]=0
        diagnostics.append(dict(reference_utc=str(ref),variant=label,outcome=outcome,n_hours=len(t),parameters=X.shape[1],residual_lag1=model.resid.autocorr(1),residual_lag24=model.resid.autocorr(24),condition_number=float(model.condition_number),r_squared=float(model.rsquared)))
        for p in PHASES:
            i=t.phase.eq(p);ci=model.conf_int().loc[p]
            rows.append(dict(reference_utc=str(ref),variant=label,outcome=outcome,phase=p,n_hours=int(i.sum()),baseline_hours=int(t.phase.eq('baseline').sum()),period_mean=float(t.loc[i,outcome].mean()),baseline_mean=float(t.loc[t.phase.eq('baseline'),outcome].mean()),adjusted_difference=float(model.params[p]),ci_low=float(ci.iloc[0]),ci_high=float(ci.iloc[1]),standard_error=float(model.bse[p]),p_value=float(model.pvalues[p]),hac_lag_hours=lag,trend=trend))
        pred.append(pd.DataFrame({'time':t.index,'outcome':outcome,'observed':t[outcome].to_numpy(),'fitted':model.predict(X).to_numpy(),'baseline_counterfactual':model.predict(xc).to_numpy(),'phase':t.phase.to_numpy(),'reference_utc':str(ref)}))
    return pd.DataFrame(rows),pd.concat(pred),pd.DataFrame(diagnostics)

def matched_difference(y,ref,p):
    starts={'pre48_24':(-48,-24),'pre24':(-24,0),'post24':(0,24)}
    baseline=y[(y.index>=ref-pd.Timedelta(hours=336))&(y.index<ref-pd.Timedelta(hours=72))]['median_depth_m']
    a,b=starts[p];event=y[(y.index>=ref+pd.Timedelta(hours=a))&(y.index<ref+pd.Timedelta(hours=b))]['median_depth_m']
    if len(baseline)!=264 or len(event)!=24 or baseline.isna().any() or event.isna().any():return np.nan
    hour=baseline.groupby(baseline.index.hour).mean()
    return float(np.mean(event.to_numpy()-hour.reindex(event.index.hour).to_numpy()))

def main():
    y=pd.read_csv(RESULT/'hourly_raw.csv',index_col='time',parse_dates=['time'])
    clipped=pd.read_csv(RESULT/'hourly_clipped.csv',index_col='time',parse_dates=['time'])
    gate=json.loads((RESULT/'environment_gate.json').read_text());refs=[REF]
    if gate['model_based_onset_utc'] and pd.Timestamp(gate['model_based_onset_utc'])!=REF:refs.append(pd.Timestamp(gate['model_based_onset_utc']))
    mainrows=[];sens=[];pred=[];diags=[]
    for ref in refs:
        m,p,d=fit(y,ref);mainrows.append(m);pred.append(p);diags.append(d)
        for name,lag,trend,data,five in [('hac12',12,True,y,False),('hac48',48,True,y,False),('no_trend',24,False,y,False),('clip_negative',24,True,clipped,False),('five_minute_proxy',24,True,y,True)]:
            m,_,d=fit(data,ref,lag,trend,five,name);sens.append(m);diags.append(d)
    m=pd.concat(mainrows,ignore_index=True);m['p_holm']=multipletests(m.p_value,method='holm')[1];m.to_csv(RESULT/'main_effects.csv',index=False)
    s=pd.concat(sens,ignore_index=True);s['p_holm_within_variant']=s.groupby('variant').p_value.transform(lambda x:multipletests(x,method='holm')[1]);s.to_csv(RESULT/'sensitivity_effects.csv',index=False)
    pd.concat(pred).to_csv(RESULT/'model_predictions.csv',index=False);pd.concat(diags).to_csv(RESULT/'model_diagnostics.csv',index=False)
    controls=[]
    for ref in pd.date_range('2015-04-16 12:00','2015-07-15 12:00',freq='D',tz='UTC'):
        if ref+pd.Timedelta(days=1)>pd.Timestamp('2015-06-01',tz='UTC') and ref-pd.Timedelta(days=14)<pd.Timestamp('2015-06-13',tz='UTC'):continue
        for p in PHASES:
            value=matched_difference(y,ref,p)
            if np.isfinite(value):controls.append(dict(reference_utc=ref,phase=p,matched_depth_difference_m=value))
    c=pd.DataFrame(controls);c.to_csv(RESULT/'contextual_pseudo_references.csv',index=False)
    ranks=[]
    for p in PHASES:
        v=matched_difference(y,REF,p);control=c[c.phase.eq(p)].matched_depth_difference_m
        ranks.append(dict(phase=p,event_matched_difference_m=v,contextual_windows=len(control),context_abs_ge_event=int((control.abs()>=abs(v)).sum()),context_abs_ge_event_fraction=float((control.abs()>=abs(v)).mean()),context_min_difference=float(control.min()),context_max_difference=float(control.max()),meaning='Descriptive rank among overlapping contextual periods; not a p-value or non-storm control'))
    ranks=pd.DataFrame(ranks);ranks.to_csv(RESULT/'contextual_ranks.csv',index=False)
    version={'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'statsmodels':sm.__version__}
    (RESULT/'software_versions.json').write_text(json.dumps(version,indent=2))
    print(m[['outcome','phase','period_mean','baseline_mean','adjusted_difference','ci_low','ci_high','p_holm']].to_string(index=False),flush=True)
    print(ranks.to_string(index=False),flush=True)

if __name__=='__main__':main()
