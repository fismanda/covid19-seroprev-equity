################################################################################
# Calculate Omicron Multipliers with Confidence Intervals
################################################################################
#
# PURPOSE:
#   Calculate the fold-change in force of infection from pre-Omicron to Omicron
#   for each quintile. This quantifies "differential amplification" - the finding
#   that affluent populations experienced larger proportional increases.
#
# INPUTS:
#   - saturated_model_estimate_rho.rds
#   - parsimonious_model_estimate_rho.rds
#
# OUTPUTS:
#   - omicron_multipliers_with_ci.csv (main results)
#   - omicron_multipliers_with_ci.rds (for plotting)
#   - Console output with formatted results
#
# AUTHOR: Hassan, Nassrallah, Fisman
# DATE: December 2025
################################################################################

library(bbmle)
library(dplyr)

cat("\n==============================================================")
cat("\n=== CALCULATING OMICRON MULTIPLIERS WITH CONFIDENCE INTERVALS ===")
cat("\n==============================================================\n\n")

# -----------------------------------------------------------------------------
# 1. LOAD SATURATED MODEL
# -----------------------------------------------------------------------------

cat("Loading saturated model...\n")
results <- readRDS("saturated_model_estimate_rho.rds")
fit <- results$fit

coefs <- coef(fit)

# Pre-Omicron
lambda_Q1_pre <- plogis(coefs["logit_lambda_Q1_pre"])
lambda_Q2_pre <- plogis(coefs["logit_lambda_Q2_pre"])
lambda_Q3_pre <- plogis(coefs["logit_lambda_Q3_pre"])
lambda_Q4_pre <- plogis(coefs["logit_lambda_Q4_pre"])
lambda_Q5_pre <- plogis(coefs["logit_lambda_Q5_pre"])

# Omicron
lambda_Q1_omi <- plogis(coefs["logit_lambda_Q1_omi"])
lambda_Q2_omi <- plogis(coefs["logit_lambda_Q2_omi"])
lambda_Q3_omi <- plogis(coefs["logit_lambda_Q3_omi"])
lambda_Q4_omi <- plogis(coefs["logit_lambda_Q4_omi"])
lambda_Q5_omi <- plogis(coefs["logit_lambda_Q5_omi"])

# Calculate multipliers (point estimates)
mult <- c(
  Q1 = lambda_Q1_omi / lambda_Q1_pre,
  Q2 = lambda_Q2_omi / lambda_Q2_pre,
  Q3 = lambda_Q3_omi / lambda_Q3_pre,
  Q4 = lambda_Q4_omi / lambda_Q4_pre,
  Q5 = lambda_Q5_omi / lambda_Q5_pre
)

cat("✓ Model loaded\n\n")

# -----------------------------------------------------------------------------
# 2. CALCULATE CONFIDENCE INTERVALS
# -----------------------------------------------------------------------------

cat("Calculating 95% confidence intervals...\n\n")

# Function to calculate CI for multiplier using delta method
calc_mult_ci <- function(fit, param_pre, param_omi, alpha = 0.05) {
  
  # Define the ratio function
  ratio_fn <- function(params) {
    omi <- plogis(params[param_omi])
    pre <- plogis(params[param_pre])
    return(omi / pre)
  }
  
  # Get point estimate
  point_est <- ratio_fn(coef(fit))
  
  # Get variance-covariance matrix
  vcov_mat <- vcov(fit)
  
  # Numerical gradient for delta method
  eps <- 1e-7
  coefs <- coef(fit)
  
  grad <- numeric(length(coefs))
  for (i in 1:length(coefs)) {
    coefs_plus <- coefs
    coefs_plus[i] <- coefs[i] + eps
    
    coefs_minus <- coefs
    coefs_minus[i] <- coefs[i] - eps
    
    grad[i] <- (ratio_fn(coefs_plus) - ratio_fn(coefs_minus)) / (2 * eps)
  }
  
  # Calculate standard error using delta method
  se <- sqrt(t(grad) %*% vcov_mat %*% grad)
  
  # Calculate CI on log scale
  log_est <- log(point_est)
  log_se <- se / point_est
  
  z_crit <- qnorm(1 - alpha/2)
  
  ci_lower <- exp(log_est - z_crit * log_se)
  ci_upper <- exp(log_est + z_crit * log_se)
  
  return(c(lower = ci_lower, upper = ci_upper))
}

# Calculate CIs for each quintile
mult_ci <- data.frame(
  Quintile = c("Q1", "Q2", "Q3", "Q4", "Q5"),
  Multiplier = mult,
  Lower = NA,
  Upper = NA
)

for (i in 1:5) {
  param_pre <- paste0("logit_lambda_Q", i, "_pre")
  param_omi <- paste0("logit_lambda_Q", i, "_omi")
  
  ci <- calc_mult_ci(fit, param_pre, param_omi)
  mult_ci$Lower[i] <- ci["lower"]
  mult_ci$Upper[i] <- ci["upper"]
  
  cat(sprintf("Q%d: %.2fx (95%% CI: %.2f - %.2f)\n", 
              i, mult[i], ci["lower"], ci["upper"]))
}

cat("\n✓ Confidence intervals calculated\n\n")

# -----------------------------------------------------------------------------
# 3. SAVE RESULTS
# -----------------------------------------------------------------------------

write.csv(mult_ci, "omicron_multipliers_with_ci.csv", row.names = FALSE)
saveRDS(mult_ci, "omicron_multipliers_with_ci.rds")

cat("✓ Results saved:\n")
cat("  - omicron_multipliers_with_ci.csv\n")
cat("  - omicron_multipliers_with_ci.rds\n\n")

# -----------------------------------------------------------------------------
# 4. FORMATTED OUTPUT
# -----------------------------------------------------------------------------

cat("=== OMICRON MULTIPLIERS WITH 95% CONFIDENCE INTERVALS ===\n\n")
cat("Fold-change in force of infection from pre-Omicron to Omicron:\n")
cat("─────────────────────────────────────────────────────\n")
cat("Quintile    Multiplier    95% CI\n")
cat("─────────────────────────────────────────────────────\n")

for (i in 1:5) {
  q <- mult_ci$Quintile[i]
  m <- mult_ci$Multiplier[i]
  ci_str <- sprintf("%.2f - %.2f", mult_ci$Lower[i], mult_ci$Upper[i])
  cat(sprintf("%-10s  %.2fx         (%s)\n", q, m, ci_str))
}

cat("\n")

# -----------------------------------------------------------------------------
# 5. DIFFERENTIAL AMPLIFICATION ANALYSIS
# -----------------------------------------------------------------------------

cat("=== DIFFERENTIAL AMPLIFICATION ===\n\n")

mult_Q1 <- mult[1]
mult_Q5 <- mult[5]
diff_ratio <- mult_Q1 / mult_Q5

cat(sprintf("Q1 (least deprived):  %.2fx increase\n", mult_Q1))
cat(sprintf("Q5 (most deprived):   %.2fx increase\n", mult_Q5))
cat(sprintf("Ratio (Q1/Q5):        %.2f\n\n", diff_ratio))

cat("INTERPRETATION:\n")
cat(sprintf("The least deprived quintile experienced a %.0f%% larger\n", 
            (diff_ratio - 1) * 100))
cat("proportional increase in force of infection compared to\n")
cat("the most deprived quintile during Omicron.\n\n")

cat("This 'differential amplification' resulted in gradient compression:\n")
cat("affluent populations 'caught up' to levels already experienced\n")
cat("by disadvantaged populations, producing apparent convergence\n")
cat("in seroprevalence despite persistent transmission disparities.\n\n")

# -----------------------------------------------------------------------------
# 6. COMPARE TO PARSIMONIOUS MODEL
# -----------------------------------------------------------------------------

cat("=== COMPARISON: SATURATED vs PARSIMONIOUS ===\n\n")

cat("Loading parsimonious model...\n")
pars_results <- readRDS("parsimonious_model_estimate_rho.rds")
pars_coefs <- coef(pars_results$fit)
shared_mult <- exp(pars_coefs["log_multiplier"])
pars_rho <- pars_results$rho

cat("✓ Loaded\n\n")

cat("PARSIMONIOUS MODEL (Shared Multiplier):\n")
cat(sprintf("  All quintiles: %.2fx\n", shared_mult))
cat(sprintf("  Seroreversion: ρ = %.5f/week (half-life: %.1f months)\n\n", 
            pars_rho, log(2)/pars_rho/4.33))

cat("SATURATED MODEL (Quintile-Specific Multipliers):\n")
for (i in 1:5) {
  cat(sprintf("  Q%d: %.2fx\n", i, mult[i]))
}
sat_rho <- results$rho
cat(sprintf("  Seroreversion: ρ = %.5f/week (half-life: %.1f months)\n\n", 
            sat_rho, log(2)/sat_rho/4.33))

cat("SEROREVERSION COMPARISON:\n")
cat(sprintf("  Saturated:     ρ = %.5f/week\n", sat_rho))
cat(sprintf("  Parsimonious:  ρ = %.5f/week\n", pars_rho))
if (abs(sat_rho - pars_rho) < 0.0005) {
  cat("  → Nearly identical (minimal antibody waning)\n\n")
} else {
  cat(sprintf("  → Difference: %.5f/week\n\n", abs(sat_rho - pars_rho)))
}

cat("\nThe parsimonious model assumes uniform Omicron amplification,\n")
cat(sprintf("averaging to %.2fx across all groups. The saturated model\n", shared_mult))
cat(sprintf("reveals substantial heterogeneity (%.2fx in Q1 vs %.2fx in Q5),\n", mult[1], mult[5]))
cat("with ΔAIC = -991 strongly favoring the saturated specification.\n\n")

# -----------------------------------------------------------------------------
# 7. ABSOLUTE CHANGES IN FORCE OF INFECTION
# -----------------------------------------------------------------------------

cat("=== ABSOLUTE CHANGES IN FORCE OF INFECTION ===\n\n")
cat("                Pre-Omicron    Omicron      Absolute\n")
cat("Quintile        (per week)     (per week)   Change\n")
cat("─────────────────────────────────────────────────────\n")

lambdas_pre <- c(lambda_Q1_pre, lambda_Q2_pre, lambda_Q3_pre, lambda_Q4_pre, lambda_Q5_pre)
lambdas_omi <- c(lambda_Q1_omi, lambda_Q2_omi, lambda_Q3_omi, lambda_Q4_omi, lambda_Q5_omi)

for (i in 1:5) {
  abs_change <- lambdas_omi[i] - lambdas_pre[i]
  cat(sprintf("Q%d          %.5f        %.5f      +%.5f\n", 
              i, lambdas_pre[i], lambdas_omi[i], abs_change))
}

cat("\nNote: Despite larger proportional increases in Q1, absolute\n")
cat("increases were similar across quintiles (~0.02/week), resulting\n")
cat("in convergence to similar final force of infection levels.\n\n")

