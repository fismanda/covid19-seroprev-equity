# Data

This directory contains the aggregate quintile-stratified anti-N seroprevalence data
used to fit the seroNegative–seroPositive (NP) compartmental model reported in:

Hassan A, Nassrallah M, Fisman DN. Seroprevalence convergence masks persistent
socioeconomic disparities in SARS-CoV-2 infection risk in Canada. *Epidemics*
(under revision, 2026).

## File

`cbs_seroprev_by_quintile.csv` — Canadian Blood Services anti-nucleocapsid
seroprevalence estimates, stratified by area-level material deprivation quintile
(Q1 = least deprived, Q5 = most deprived), April 21, 2021 to April 26, 2023.

## Columns

| Column | Description |
|---|---|
| `date` | Sampling date (YYYY-MM-DD) |
| `quintile` | Material deprivation quintile (Q1 – Q5) |
| `n_donors` | Number of CBS donors tested in the serosurvey |
| `seroprev_est` | Crude observed seroprevalence (proportion) |
| `seroprev_lo95` | Lower 95% binomial confidence bound |
| `seroprev_hi95` | Upper 95% binomial confidence bound |

## Source

Aggregated from publicly available COVID-19 Immunity Task Force (CITF) serosurvey
data: https://www.covid19immunitytaskforce.ca/

## Usage in code

The R scripts in `/Analysis/` and `/Figures/` load this dataset from
`Data/cbs_seroprev_by_quintile.csv` relative to the repository root.

## Citation

If you use this dataset, please cite the original CITF data release and the
manuscript above. A Zenodo archive with persistent DOI is associated with this
repository.
