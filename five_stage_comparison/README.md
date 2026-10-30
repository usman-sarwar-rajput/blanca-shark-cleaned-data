# Five-stage Blanca comparison

This additional descriptive view follows the user-requested sequence: earlier quiet days, immediately before the cyclone calendar period, during that period, immediately after, and later quiet days. No new significance tests were run.

## Calendar definitions

| Stage | UTC dates, inclusive | Duration |
|---|---|---|
| Earlier quiet* | 2015-05-25 to 2015-05-27 | 3 days |
| Just before | 2015-05-28 to 2015-05-30 | 3 days |
| Blanca active dates | 2015-05-31 to 2015-06-09 | 10 days |
| Just after | 2015-06-10 to 2015-06-12 | 3 days |
| Later quiet* | 2015-06-13 to 2015-06-15 | 3 days |

Blanca active calendar dates are 31 May–9 June 2015, as documented in the [NHC report](https://www.nhc.noaa.gov/data/tcr/EP022015_Blanca.pdf). Full UTC calendar days are used, not exact genesis/dissipation times. This red period is a regional cyclone calendar, not verified severe local exposure.

**Change from the preceding charts:** earlier charts used the June 8 landfall/proximity reference and compared the three days around it. This chart uses the entire cyclone calendar. June 5–7 therefore belong to the red active-period category here, although they appeared yellow before the June 8 reference in the earlier charts. Neither chart establishes local onset.

## Meaning of quiet

The first and last periods pass an explicit exploratory proxy: every checked hour has modelled 10-m wind <10 m/s, significant wave height <2 m, complete values, and no contemporaneous >=34-knot eastern Pacific cyclone track point within 500 km of a modelled shark position. The cyclone screen uses the prior downloaded IBTrACS MAIN tracks, matching directly to original modelled positions within 90 minutes, without hourly position carry-forward. Thresholds describe these candidate control periods and are not a validated ecological definition of normality. Dates were specified from the requested chronology and weather summaries, before computing their behavioural means.

This cannot verify that nothing happened: tides, prey, fishing, other disturbances, weak/faraway systems and model errors are not excluded. GPE3 shark positions and weather are model-derived. The later period also follows the previously identified June 9 wave-source discontinuity, whose cause remains unresolved; lower later model waves therefore require caution. These are candidate quiet periods, not validated storm-free controls.

Andres and Carlos occur elsewhere in the same eastern Pacific season; their contemporaneous tracks are explicitly included in cyclone_screen.csv. Distance >500 km does not eliminate swell or regional influence. The storm-period mean combines ten days with changing conditions and can hide a short response. Use the daily timeline alongside the stage-average bars.

## Calculations

Depth is the average of hourly median depths; vertical movement is the average hourly absolute change between adjacent one-minute median depths (m/min). It is a pressure-tag movement proxy. Each three-day period has 72 complete hours and 17,280 sensor samples; the ten-day calendar period has 240 complete hours and 57,600 samples. Periods are disjoint. Summary means were checked against the means of their complete daily values. No observed shift is automatically attributed to Blanca.

## Files

- 01_five_stage_comparison.png / .pdf: five-stage averages.
- 02_daily_five_stage_timeline.png / .pdf: daily variation with stage shading.
- five_stage_summary.csv: stage means, exact time bounds and environmental checks.
- daily_values.csv: one row per day.
- cyclone_screen.csv: named concurrent cyclone proximity checks.

Inputs are the existing hourly_raw.csv, environment_matched_hourly.csv and the prior archived cyclone tracks. Existing statistical reports and files were not changed.