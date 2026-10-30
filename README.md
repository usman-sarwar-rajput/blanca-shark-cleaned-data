# Hurricane Blanca juvenile white shark analysis data

This repository contains cleaned and derived analysis tables for an exploratory comparison of juvenile white shark `14_07` around the 2015 Hurricane Blanca calendar period. The original telemetry was collected by the Monterey Bay Aquarium collaboration; this repository's contribution is data preparation, period definitions, environmental matching, and descriptive summaries.

## Original data and attribution

- Original telemetry: https://doi.org/10.24431/rw1k6c3 (CC BY; retain attribution to the original data creators).
- Data description: https://doi.org/10.1038/s41597-022-01235-3
- Hurricane Blanca report: https://www.nhc.noaa.gov/data/tcr/EP022015_Blanca.pdf

Please cite the original dataset and data-description article when using these derived tables. The original tag records are available at the source DOI and are not claimed as newly collected here.

## Contents

- `results/`: cleaned hourly tables, environmental matches, model outputs, diagnostics, and validation files.
- `five_stage_comparison/`: daily values, cyclone screen, and the five-period summary used for the paper's descriptive table.
- `three_day_comparison/`: an alternate shorter comparison.
- `code/`: processing and analysis scripts.
- `ANALYSIS_PROTOCOL.md`, `SOURCE_MANIFEST.json`, `requirements.txt`: methods, source checksums, and software requirements.

The five-stage table compares 25–27 May, 28–30 May, 31 May–9 June, 10–12 June, and 13–15 June 2015. The active period is a cyclone calendar interval, not verified local storm exposure at the animal. This single-animal descriptive analysis does not establish anticipation, causation, or predictive performance.

Some results are model-derived and should be read with the methods and validation files. The original archived telemetry and third-party weather source files are not bundled here; see the source DOI and manifest. The data files and code have different potential rights, so no new blanket licence is asserted over third-party material.
