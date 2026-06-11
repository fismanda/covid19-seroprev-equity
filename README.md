# covid19-seroprev-equity

Analysis code and data accompanying:

**Hassan A, Nassrallah M, Fisman DN.** *Seroprevalence convergence masks
persistent socioeconomic disparities in SARS-CoV-2 infection risk in Canada.*
*Epidemics*, 2026 (revision under review).

medRxiv preprint: [10.64898/2026.01.02.26343363](https://doi.org/10.64898/2026.01.02.26343363)

Zenodo archive: [10.5281/zenodo.18064548](https://doi.org/10.5281/zenodo.18064548)

## Overview

We analyzed serial cross-sectional Canadian Blood Services anti-nucleocapsid
seroprevalence data (April 2021 to April 2023) stratified by area-level material
deprivation quintile (Canadian Index of Multiple Deprivation). We fit a
two-compartment seroNegative–seroPositive (NP) dynamic model with seroreversion
to estimate quintile-specific forces of infection before and after Omicron
emergence (January 1, 2022).

Key findings:

- Pre-Omicron, the most deprived quintile (Q5) experienced 71% higher force of
  infection than the least deprived (Q1; IRR 1.71, 95% CI 1.60 to 1.83).
- Following Omicron emergence, the relative gradient compressed (Q5 vs Q1 IRR
  1.12, 95% CI 1.11 to 1.14) but the absolute force of infection remained higher
  in more-deprived quintiles throughout the period.
- The compression was driven by differential amplification (Q1 Omicron multiplier
  48.5×; Q5 Omicron multiplier 31.8×; non-overlapping CIs).
- Apparent seroprevalence convergence by late 2022 masked these persistent
  disparities in the underlying hazard of infection.

## Repository structure

- **Data/** — Aggregate analytic dataset and documentation
  - `cbs_seroprev_by_quintile.csv` (deposited analytic dataset)
  - `README.md` (column definitions and source)
- **Analysis/revision_2026/** — Revised analysis code (June 2026)
  - `01_fit_saturated_model.R`
  - `02_fit_parsimonious_model.R`
  - `03_compare_models.R`
  - `04_calculate_irr_with_ci.R`
  - `05_calculate_omicron_multipliers_and_rho.R`
  - `06_calculate_parsimonious_rho_ci.R`
  - `07_generate_sample_size_table_S1B.R`
  - `README.md`
- **Figures/revision_2026/** — Revised figure code (June 2026)
  - `figure1_two_panel.R`
  - `figure2_two_panel.R`
  - `figure3_label_fix.R`
  - `figureS1_with_CI.R`
  - `README.md`

Original scripts from the initial submission are retained in the `Analysis/` and
`Figures/` parent directories for transparency.

## Reproducibility

1. **Data.** Load `Data/cbs_seroprev_by_quintile.csv`. This is the aggregate
   quintile-stratified dataset that drives the analysis (410 observations: 82
   sampling dates × 5 quintiles; cumulative 849,452 donor-tests across April 21,
   2021 to April 26, 2023). See `Data/README.md` for column definitions.
2. **Model fitting.** Run scripts in `Analysis/revision_2026/` in numerical order
   (01 through 07). Scripts 01 and 02 produce model fit objects that 03–07
   consume.
3. **Figures.** Run scripts in `Figures/revision_2026/` after the model fits
   exist. Each script saves PNG (300 dpi) and PDF versions of its figure.

R packages required: `bbmle`, `ggplot2`, `dplyr`, `lubridate`, `patchwork`,
`MASS`, `tidyr`.

## Citation

If you use this code or data, please cite the manuscript (DOI to be added upon
publication) and the Zenodo archive linked at the top of this README.

## Data source acknowledgment

Seroprevalence data are derived from the publicly released Canadian Blood Services
serosurveillance program data aggregated by the COVID-19 Immunity Task Force
(CITF). CITF and CBS serosurveillance methodology described in O'Brien et al.
(*Microbiology Spectrum*, 2023) and Murphy et al. (*CMAJ*, 2023).

## Contact

David Fisman, Dalla Lana School of Public Health, University of Toronto.
[david.fisman@utoronto.ca](mailto:david.fisman@utoronto.ca)
