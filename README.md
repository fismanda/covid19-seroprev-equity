# Seroprevalence Convergence Does Not Reflect Transmission Equity

This repository contains code for the analysis presented in:

**Hassan A, Nasrallah EI, Fisman DN (2025).** *Seroprevalence Convergence Does Not Reflect Transmission Equity: Persistent Socioeconomic Disparities in COVID-19 Force of Infection in Canada.* [Journal TBD]

## Overview

This analysis demonstrates that apparent convergence in SARS-CoV-2 seroprevalence across socioeconomic strata masked persistent inequalities in force of infection. Using dynamic compartmental modeling of serial seroprevalence data from Canadian blood donors, we show that materially deprived populations maintained elevated transmission risk throughout the pandemic, despite reaching similar cumulative infection levels by 2023.

### Key Findings

- **Pre-Omicron period (2020-2021):** The most deprived quintile experienced 71% higher force of infection than the least deprived quintile (IRR: 1.71; 95% CI: 1.60-1.83)
- **Omicron period (2022-2023):** Force of infection increased dramatically in all quintiles but by differing magnitudes
  - Least deprived quintile: 48.5-fold increase
  - Most deprived quintile: 31.8-fold increase
  - Resulting gradient compression: Q5 vs Q1 IRR = 1.12 (95% CI: 1.11-1.14)
- **Interpretation:** Seroprevalence convergence arose from differential amplification during Omicron rather than elimination of underlying socioeconomic disparities

## Repository Structure

```
├── README.md                           # This file
├── analysis/                           # Main analysis scripts
│   ├── 01_fit_saturated_model.R       # Fit saturated model (quintile-specific λ pre/post Omicron)
│   ├── 02_fit_parsimonious_model.R    # Fit parsimonious model (shared Omicron multiplier)
│   ├── 03_compare_models.R            # Model comparison (AIC, BIC, likelihood ratio test)
│   ├── 04_calculate_irr_with_ci.R     # Calculate incidence rate ratios with 95% CIs
│   └── 05_calculate_omicron_multipliers_and_rho.R  # Calculate Omicron multipliers
└── figures/                            # Figure generation scripts
    ├── figure1_model_fit_FINAL.R      # Figure 1: Model fit to seroprevalence data
    ├── figure2_irr_forest_plot.R      # Figure 2: Socioeconomic gradient (forest plot)
    ├── figure3_omicron_multipliers.R  # Figure 3: Differential amplification
    └── figureS1_parsimonious_fit.R    # Figure S1: Parsimonious model fit
```

## Data Source

Seroprevalence data are publicly available from the **COVID-19 Immunity Task Force (CITF) Data Portal**:  
https://portal.citf.mcgill.ca

The analysis uses Canadian Blood Services (CBS) serosurvey data measuring infection-induced antibodies (anti-nucleocapsid) stratified by material deprivation quintile, collected between April 2021 and April 2023.

## Methods

### Mathematical Model

We developed a two-compartment susceptible-infected (SI) model with seroreversion to estimate force of infection from observed seroprevalence trajectories:

```
S(t+1) = S(t) × exp(-λ(t)) + I(t) × (1 - exp(-ρ))
I(t+1) = I(t) × exp(-ρ) + S(t) × (1 - exp(-λ(t)))
```

Where:
- `S(t)` = proportion seronegative at time t
- `I(t)` = proportion seropositive at time t  
- `λ(t)` = time-varying force of infection
- `ρ` = seroreversion rate

### Model Specifications

**Saturated Model:**
- Quintile-specific forces of infection for both pre-Omicron and Omicron periods
- 11 parameters: 10 λ values (5 quintiles × 2 periods) + ρ

**Parsimonious Model:**
- Quintile-specific baseline forces of infection
- Shared Omicron multiplier applied to all quintiles
- 7 parameters: 5 baseline λ values + 1 multiplier + ρ

The saturated model was strongly preferred (ΔAIC = -194.17, p < 0.001).

## Requirements

### R Version
R ≥ 4.0.0

### Required Packages
```r
install.packages(c("bbmle", "dplyr", "tidyr", "ggplot2", "readr"))
```

## Usage

### Running the Full Analysis

1. **Download seroprevalence data** from CITF Data Portal
2. **Fit models:**
   ```r
   source("analysis/01_fit_saturated_model.R")
   source("analysis/02_fit_parsimonious_model.R")
   ```
3. **Compare models:**
   ```r
   source("analysis/03_compare_models.R")
   ```
4. **Calculate estimates:**
   ```r
   source("analysis/04_calculate_irr_with_ci.R")
   source("analysis/05_calculate_omicron_multipliers_and_rho.R")
   ```
5. **Generate figures:**
   ```r
   source("figures/figure1_model_fit_FINAL.R")
   source("figures/figure2_irr_forest_plot.R")
   source("figures/figure3_omicron_multipliers.R")
   source("figures/figureS1_parsimonious_fit.R")
   ```

### Key Outputs

- `saturated_model_estimate_rho.rds` - Fitted saturated model
- `parsimonious_model_estimate_rho.rds` - Fitted parsimonious model
- `irr_estimates_with_ci.csv` - Incidence rate ratios with 95% CIs
- `omicron_multipliers.csv` - Quintile-specific Omicron amplification factors
- Figure files (PNG and PDF formats)

## Citation

If you use this code, please cite:

```
Hassan A, Nasrallah EI, Fisman DN (2025). Seroprevalence Convergence Does Not 
Reflect Transmission Equity: Persistent Socioeconomic Disparities in COVID-19 
Force of Infection in Canada. [Journal TBD].
```

## Authors

- **Abdullah Hassan** - Dalla Lana School of Public Health, University of Toronto
- **Emmanuel Issa Nasrallah** - Dalla Lana School of Public Health, University of Toronto
- **David N. Fisman** - Dalla Lana School of Public Health, University of Toronto

*Both first authors contributed equally to this work.

## License

This project is licensed under the MIT License - see below for details.

```
MIT License

Copyright (c) 2025 Hassan, Nasrallah, Fisman

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Funding

This research was supported by grants to DNF from the Canadian Institutes for Health Research (2019 COVID-19 rapid research funding OV4-170360 and operating grant #518192).

## Contact

For questions about the analysis or code, please contact:
- David Fisman: david.fisman@utoronto.ca

## Acknowledgments

We acknowledge the COVID-19 Immunity Task Force and Canadian Blood Services for making seroprevalence data publicly available.
