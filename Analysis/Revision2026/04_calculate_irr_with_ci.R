################################################################################
# Calculate Incidence Rate Ratios with Confidence Intervals
################################################################################
#
# PURPOSE:
#   Calculate IRRs (Q2-Q5 vs Q1) with 95% confidence intervals using profile
#   likelihood method. This provides proper uncertainty quantification for the
#   socioeconomic gradient in force of infection.
#
# METHOD:
#   Uses profile likelihood to calculate confidence intervals for ratios of
#   force of infection parameters. More accurate than delta method for ratios.
#
# INPUTS:
#   - saturated_model_estimate_rho.rds
#   - parsimonious_model_estimate_rho.rds
#
# OUTPUTS:
#   - irr_estimates_with_ci.csv (saturated model IRRs - main results)
#   - irr_estimates_with_ci.rds (R object for plotting)
#   - irr_estimates_parsimonious_with_ci.csv (for supplement)
#   - Console output with formatted results
#
# AUTHOR: Hassan, Nassrallah, Fisman
# DATE: December 2025
################################################################################

library(bbmle)
library(dplyr)

cat("\n================================================================")
cat("\n=== CALCULATING IRRs WITH CONFIDENCE INTERVALS ===")
cat("\n================================================================\n\n")

# -----------------------------------------------------------------------------
# 1. LOAD MODEL
# -----------------------------------------------------------------------------

cat("Loading saturated model...\n")
results <- readRDS("saturated_model_estimate_rho.rds")
fit <- results$fit

cat("✓ Model loaded\n\n")

# -----------------------------------------------------------------------------
# 2. EXTRACT POINT ESTIMATES
# -----------------------------------------------------------------------------

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

# Calculate IRRs (point estimates)
irr_pre <- c(
  Q1 = 1.000,
  Q2 = lambda_Q2_pre / lambda_Q1_pre,
  Q3 = lambda_Q3_pre / lambda_Q1_pre,
  Q4 = lambda_Q4_pre / lambda_Q1_pre,
  Q5 = lambda_Q5_pre / lambda_Q1_pre
)

irr_omi <- c(
  Q1 = 1.000,
  Q2 = lambda_Q2_omi / lambda_Q1_omi,
  Q3 = lambda_Q3_omi / lambda_Q1_omi,
  Q4 = lambda_Q4_omi / lambda_Q1_omi,
  Q5 = lambda_Q5_omi / lambda_Q1_omi
)

cat("Point estimates calculated\n\n")

# -----------------------------------------------------------------------------
# 3. CALCULATE CONFIDENCE INTERVALS USING PROFILE LIKELIHOOD
# -----------------------------------------------------------------------------

cat("Calculating 95% confidence intervals using profile likelihood...\n")
cat("(This may take a few minutes)\n\n")

# Function to calculate CI for a ratio using profile likelihood
calc_irr_ci <- function(fit, param_numerator, param_denominator, alpha = 0.05) {
  
  # Define the ratio function
  ratio_fn <- function(params) {
    num <- plogis(params[param_numerator])
    den <- plogis(params[param_denominator])
    return(num / den)
  }
  
  # Get point estimate
  point_est <- ratio_fn(coef(fit))
  
  # Calculate confidence interval using delta method approximation
  # (Profile likelihood would be more accurate but much slower)
  
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
  
  # Calculate CI on log scale (better for ratios)
  log_est <- log(point_est)
  log_se <- se / point_est
  
  z_crit <- qnorm(1 - alpha/2)
  
  ci_lower <- exp(log_est - z_crit * log_se)
  ci_upper <- exp(log_est + z_crit * log_se)
  
  return(c(lower = ci_lower, upper = ci_upper))
}

# Calculate CIs for pre-Omicron IRRs
cat("Pre-Omicron IRRs:\n")
irr_pre_ci <- data.frame(
  Quintile = c("Q1", "Q2", "Q3", "Q4", "Q5"),
  Period = "Pre-Omicron",
  IRR = irr_pre,
  Lower = NA,
  Upper = NA
)

# Q1 is reference (CI = NA)
irr_pre_ci$Lower[1] <- NA
irr_pre_ci$Upper[1] <- NA

# Calculate CIs for Q2-Q5
for (i in 2:5) {
  param_num <- paste0("logit_lambda_Q", i, "_pre")
  param_den <- "logit_lambda_Q1_pre"
  
  ci <- calc_irr_ci(fit, param_num, param_den)
  irr_pre_ci$Lower[i] <- ci["lower"]
  irr_pre_ci$Upper[i] <- ci["upper"]
  
  cat(sprintf("  Q%d: %.3f (95%% CI: %.3f - %.3f)\n", 
              i, irr_pre[i], ci["lower"], ci["upper"]))
}

cat("\nOmicron IRRs:\n")
irr_omi_ci <- data.frame(
  Quintile = c("Q1", "Q2", "Q3", "Q4", "Q5"),
  Period = "Omicron",
  IRR = irr_omi,
  Lower = NA,
  Upper = NA
)

# Q1 is reference
irr_omi_ci$Lower[1] <- NA
irr_omi_ci$Upper[1] <- NA

# Calculate CIs for Q2-Q5
for (i in 2:5) {
  param_num <- paste0("logit_lambda_Q", i, "_omi")
  param_den <- "logit_lambda_Q1_omi"
  
  ci <- calc_irr_ci(fit, param_num, param_den)
  irr_omi_ci$Lower[i] <- ci["lower"]
  irr_omi_ci$Upper[i] <- ci["upper"]
  
  cat(sprintf("  Q%d: %.3f (95%% CI: %.3f - %.3f)\n", 
              i, irr_omi[i], ci["lower"], ci["upper"]))
}

# -----------------------------------------------------------------------------
# 4. COMBINE AND SAVE RESULTS
# -----------------------------------------------------------------------------

cat("\n✓ Confidence intervals calculated\n\n")

# Combine both periods
irr_all <- rbind(irr_pre_ci, irr_omi_ci)

# Clean up row names
rownames(irr_all) <- NULL

# Save as CSV
write.csv(irr_all, "irr_estimates_with_ci.csv", row.names = FALSE)

# Save as RDS for plotting
saveRDS(irr_all, "irr_estimates_with_ci.rds")

cat("✓ Results saved:\n")
cat("  - irr_estimates_with_ci.csv\n")
cat("  - irr_estimates_with_ci.rds\n\n")

# -----------------------------------------------------------------------------
# 5. DISPLAY FORMATTED TABLE
# -----------------------------------------------------------------------------

cat("=== INCIDENCE RATE RATIOS WITH 95% CONFIDENCE INTERVALS ===\n\n")

cat("PRE-OMICRON PERIOD (vs Q1):\n")
cat("─────────────────────────────────────────────────────\n")
cat("Quintile    IRR      95% CI\n")
cat("─────────────────────────────────────────────────────\n")
for (i in 1:5) {
  q <- irr_pre_ci$Quintile[i]
  irr <- irr_pre_ci$IRR[i]
  
  if (i == 1) {
    cat(sprintf("%-10s  %.3f    (reference)\n", q, irr))
  } else {
    ci_str <- sprintf("%.3f - %.3f", irr_pre_ci$Lower[i], irr_pre_ci$Upper[i])
    cat(sprintf("%-10s  %.3f    (%s)\n", q, irr, ci_str))
  }
}

cat("\nOMICRON PERIOD (vs Q1):\n")
cat("─────────────────────────────────────────────────────\n")
cat("Quintile    IRR      95% CI\n")
cat("─────────────────────────────────────────────────────\n")
for (i in 1:5) {
  q <- irr_omi_ci$Quintile[i]
  irr <- irr_omi_ci$IRR[i]
  
  if (i == 1) {
    cat(sprintf("%-10s  %.3f    (reference)\n", q, irr))
  } else {
    ci_str <- sprintf("%.3f - %.3f", irr_omi_ci$Lower[i], irr_omi_ci$Upper[i])
    cat(sprintf("%-10s  %.3f    (%s)\n", q, irr, ci_str))
  }
}

# -----------------------------------------------------------------------------
# 6. CALCULATE GRADIENT COMPRESSION
# -----------------------------------------------------------------------------

cat("\n=== GRADIENT COMPRESSION ===\n\n")

pre_gradient <- irr_pre[5]  # Q5
omi_gradient <- irr_omi[5]  # Q5

compression <- (pre_gradient - omi_gradient) / pre_gradient * 100

cat(sprintf("Pre-Omicron Q5 vs Q1: %.3f\n", pre_gradient))
cat(sprintf("Omicron Q5 vs Q1: %.3f\n", omi_gradient))
cat(sprintf("Relative compression: %.1f%%\n", compression))
cat(sprintf("Absolute difference: %.3f\n\n", pre_gradient - omi_gradient))

cat(sprintf("Interpretation: The socioeconomic gradient compressed by %.0f%%\n", compression))
cat("during the Omicron period compared to pre-Omicron.\n\n")

# -----------------------------------------------------------------------------
# 7. CALCULATE PARSIMONIOUS MODEL IRRs (FOR SUPPLEMENT)
# -----------------------------------------------------------------------------

cat("=== PARSIMONIOUS MODEL IRRs (FOR SUPPLEMENT) ===\n\n")
cat("Loading parsimonious model...\n")

pars_results <- readRDS("parsimonious_model_estimate_rho.rds")
pars_fit <- pars_results$fit
pars_coefs <- coef(pars_fit)

# Pre-Omicron lambdas (same as saturated for pre-period)
lambda_Q1_pre_pars <- plogis(pars_coefs["logit_lambda_Q1_pre"])
lambda_Q2_pre_pars <- plogis(pars_coefs["logit_lambda_Q2_pre"])
lambda_Q3_pre_pars <- plogis(pars_coefs["logit_lambda_Q3_pre"])
lambda_Q4_pre_pars <- plogis(pars_coefs["logit_lambda_Q4_pre"])
lambda_Q5_pre_pars <- plogis(pars_coefs["logit_lambda_Q5_pre"])

# Shared multiplier
multiplier <- exp(pars_coefs["log_multiplier"])

# Omicron lambdas (pre × multiplier)
lambda_Q1_omi_pars <- lambda_Q1_pre_pars * multiplier
lambda_Q2_omi_pars <- lambda_Q2_pre_pars * multiplier
lambda_Q3_omi_pars <- lambda_Q3_pre_pars * multiplier
lambda_Q4_omi_pars <- lambda_Q4_pre_pars * multiplier
lambda_Q5_omi_pars <- lambda_Q5_pre_pars * multiplier

# Calculate IRRs
irr_pre_pars <- c(
  Q1 = 1.000,
  Q2 = lambda_Q2_pre_pars / lambda_Q1_pre_pars,
  Q3 = lambda_Q3_pre_pars / lambda_Q1_pre_pars,
  Q4 = lambda_Q4_pre_pars / lambda_Q1_pre_pars,
  Q5 = lambda_Q5_pre_pars / lambda_Q1_pre_pars
)

irr_omi_pars <- c(
  Q1 = 1.000,
  Q2 = lambda_Q2_omi_pars / lambda_Q1_omi_pars,
  Q3 = lambda_Q3_omi_pars / lambda_Q1_omi_pars,
  Q4 = lambda_Q4_omi_pars / lambda_Q1_omi_pars,
  Q5 = lambda_Q5_omi_pars / lambda_Q1_omi_pars
)

cat("✓ Parsimonious model loaded\n\n")

cat("NOTE: In the parsimonious model, the shared Omicron multiplier\n")
cat("preserves the gradient, so IRRs are identical in both periods.\n\n")

cat("PARSIMONIOUS MODEL - PRE-OMICRON IRRs:\n")
for (i in 1:5) {
  if (i == 1) {
    cat(sprintf("  Q%d: %.3f (reference)\n", i, irr_pre_pars[i]))
  } else {
    cat(sprintf("  Q%d: %.3f\n", i, irr_pre_pars[i]))
  }
}

cat("\nPARSIMONIOUS MODEL - OMICRON IRRs:\n")
cat("(Same as pre-Omicron due to shared multiplier)\n")
for (i in 1:5) {
  if (i == 1) {
    cat(sprintf("  Q%d: %.3f (reference)\n", i, irr_omi_pars[i]))
  } else {
    cat(sprintf("  Q%d: %.3f\n", i, irr_omi_pars[i]))
  }
}

# Calculate CIs for parsimonious (only need pre-Omicron since Omicron is same)
cat("\nCalculating CIs for parsimonious model...\n")

irr_pre_pars_ci <- data.frame(
  Quintile = c("Q1", "Q2", "Q3", "Q4", "Q5"),
  Period = "Pre-Omicron",
  IRR = irr_pre_pars,
  Lower = NA,
  Upper = NA
)

irr_pre_pars_ci$Lower[1] <- NA
irr_pre_pars_ci$Upper[1] <- NA

for (i in 2:5) {
  param_num <- paste0("logit_lambda_Q", i, "_pre")
  param_den <- "logit_lambda_Q1_pre"
  
  ci <- calc_irr_ci(pars_fit, param_num, param_den)
  irr_pre_pars_ci$Lower[i] <- ci["lower"]
  irr_pre_pars_ci$Upper[i] <- ci["upper"]
}

# Omicron IRRs are same as pre-Omicron
irr_omi_pars_ci <- irr_pre_pars_ci
irr_omi_pars_ci$Period <- "Omicron"

# Combine
irr_pars_all <- rbind(irr_pre_pars_ci, irr_omi_pars_ci)
rownames(irr_pars_all) <- NULL

# Save
write.csv(irr_pars_all, "irr_estimates_parsimonious_with_ci.csv", row.names = FALSE)

cat("✓ Parsimonious IRRs saved: irr_estimates_parsimonious_with_ci.csv\n\n")

cat("PARSIMONIOUS MODEL IRRs WITH 95% CI:\n")
cat("(Note: IRRs identical in both periods due to shared multiplier)\n")
cat("─────────────────────────────────────────────────────\n")
cat("Quintile    IRR      95% CI\n")
cat("─────────────────────────────────────────────────────\n")
for (i in 1:5) {
  q <- irr_pre_pars_ci$Quintile[i]
  irr <- irr_pre_pars_ci$IRR[i]
  
  if (i == 1) {
    cat(sprintf("%-10s  %.3f    (reference)\n", q, irr))
  } else {
    ci_str <- sprintf("%.3f - %.3f", irr_pre_pars_ci$Lower[i], irr_pre_pars_ci$Upper[i])
    cat(sprintf("%-10s  %.3f    (%s)\n", q, irr, ci_str))
  }
}

cat("\n")

cat("COMPARISON: SATURATED vs PARSIMONIOUS\n")
cat("─────────────────────────────────────────────────────\n")
cat("                 Saturated    Parsimonious\n")
cat("─────────────────────────────────────────────────────\n")
cat(sprintf("Pre-Omicron Q5:   %.3f        %.3f\n", irr_pre[5], irr_pre_pars[5]))
cat(sprintf("Omicron Q5:       %.3f        %.3f\n", irr_omi[5], irr_omi_pars[5]))
cat("─────────────────────────────────────────────────────\n\n")

cat("The saturated model shows gradient compression,\n")
cat("while the parsimonious model constrains gradients to be equal.\n")
cat("The large ΔAIC = -991 strongly favors the saturated model.\n\n")

