# Figures — Revision 2026

R scripts that generate the figures for the manuscript submitted to *Epidemics* in
June 2026.

## Contents

- `figure1_two_panel.R` — **Figure 1.** Two-panel display of infection-induced
  seroprevalence trajectories by material deprivation quintile.
  - Panel A: observed (points) vs model-predicted (lines) seroprevalence, with
    parametric bootstrap 95% CI bands (1,000 draws from the multivariate normal
    approximation to the parameter sampling distribution at the MLE).
  - Panel B: Q5 minus Q1 absolute difference in seroprevalence over time, with
    bootstrap 95% confidence bounds shown as dashed lines.

- `figure2_two_panel.R` — **Figure 2.** Forest-plot-style display of socioeconomic
  gradients in force of infection, with quintiles on the y-axis (Q1 top, Q5 bottom).
  - Panel A: incidence rate ratios (relative force of infection) by quintile and
    period, with delta-method 95% CIs.
  - Panel B: absolute weekly force of infection by quintile and period, on a log
    scale, with delta-method 95% CIs.

- `figure3_label_fix.R` — **Figure 3.** Bar chart of Omicron multipliers
  (fold-change in force of infection during Omicron relative to pre-Omicron) by
  material deprivation quintile, with 95% CI whiskers. Multiplier text labels are
  positioned above the upper whiskers to eliminate overlap.

- `figureS1_with_CI.R` — **Supplementary Figure S1.** Parsimonious-model fit
  (same display style as Figure 1 Panel A, but using the model with a single
  shared Omicron multiplier across all quintiles).

## Dependencies

Scripts read intermediate model objects from the working directory. Before running
any figure script, ensure the corresponding analysis script has produced its
output:

| Figure script | Required input |
|---|---|
| `figure1_two_panel.R` | `saturated_model_estimate_rho.rds` |
| `figure2_two_panel.R` | `irr_estimates_with_ci.rds`, `saturated_model_estimate_rho.rds` |
| `figure3_label_fix.R` | `omicron_multipliers_with_ci.rds` |
| `figureS1_with_CI.R` | `parsimonious_model_estimate_rho.rds` |

R packages required: `ggplot2`, `dplyr`, `lubridate`, `patchwork`, `MASS`,
`tidyr`. The bootstrap-based scripts (Figure 1 and S1) take 30 to 60 seconds to
run due to the 1,000-draw bootstrap. Other figures render in seconds.

## Outputs

Each script saves both PNG (300 dpi) and PDF (vector) versions of its figure to the
working directory.
