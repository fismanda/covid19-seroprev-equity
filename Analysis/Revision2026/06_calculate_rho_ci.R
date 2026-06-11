################################################################################
# Calculate 95% Confidence Interval for Seroreversion Rate (rho)
################################################################################

library(bbmle)

cat("\n=== CALCULATING RHO CONFIDENCE INTERVAL ===\n\n")

# Load saturated model
results <- readRDS("saturated_model_estimate_rho.rds")
fit <- results$fit

# Get coefficient estimates
coefs <- coef(fit)

# Get rho (seroreversion rate)
# It's stored on logit scale and scaled to (0, 0.01)
logit_rho_scaled <- coefs["logit_rho_scaled"]

# Transform back to rho
rho <- plogis(logit_rho_scaled) * 0.01

cat(sprintf("Point estimate: ρ = %.5f per week\n", rho))

# Get 95% CI using delta method
vcov_mat <- vcov(fit)
se_logit_rho <- sqrt(vcov_mat["logit_rho_scaled", "logit_rho_scaled"])

# CI on logit scale
z_crit <- qnorm(0.975)
logit_lower <- logit_rho_scaled - z_crit * se_logit_rho
logit_upper <- logit_rho_scaled + z_crit * se_logit_rho

# Transform back to rho scale
rho_lower <- plogis(logit_lower) * 0.01
rho_upper <- plogis(logit_upper) * 0.01

cat(sprintf("95%% CI: (%.5f, %.5f)\n\n", rho_lower, rho_upper))

# Calculate half-life for each
halflife_point <- log(2) / rho
halflife_lower <- log(2) / rho_upper  # inverted because higher rho = shorter half-life
halflife_upper <- log(2) / rho_lower

cat("HALF-LIFE ESTIMATES:\n")
cat(sprintf("  Weeks:  %.1f (95%% CI: %.1f, %.1f)\n", 
            halflife_point, halflife_lower, halflife_upper))
cat(sprintf("  Months: %.1f (95%% CI: %.1f, %.1f)\n", 
            halflife_point/4.33, halflife_lower/4.33, halflife_upper/4.33))
cat(sprintf("  Years:  %.1f (95%% CI: %.1f, %.1f)\n\n", 
            halflife_point/52, halflife_lower/52, halflife_upper/52))

# For manuscript text
cat("FOR MANUSCRIPT:\n")
cat(sprintf("ρ = %.5f per week (95%% CI: %.5f, %.5f)\n", rho, rho_lower, rho_upper))
cat(sprintf("Half-life: %.0f months (%.1f years; 95%% CI: %.0f-%.0f months)\n\n", 
            halflife_point/4.33, halflife_point/52, 
            halflife_lower/4.33, halflife_upper/4.33))


# Get rho
rho <- results$rho
cat(sprintf("Point estimate: ρ = %.5f per week\n", rho))

# Get 95% CI using delta method
vcov_mat <- vcov(fit)
se_logit_rho <- sqrt(vcov_mat["logit_rho", "logit_rho"])

# Get the logit_rho value
logit_rho <- coef(fit)["logit_rho"]

# CI on logit scale
z_crit <- qnorm(0.975)
logit_lower <- logit_rho - z_crit * se_logit_rho
logit_upper <- logit_rho + z_crit * se_logit_rho

# Transform back to rho scale (assuming it was scaled between 0 and 0.01)
rho_lower <- plogis(logit_lower) * 0.01
rho_upper <- plogis(logit_upper) * 0.01

cat(sprintf("95%% CI: (%.5f, %.5f)\n\n", rho_lower, rho_upper))

# Calculate half-life
halflife_point <- log(2) / rho
halflife_lower <- log(2) / rho_upper
halflife_upper <- log(2) / rho_lower

cat("HALF-LIFE ESTIMATES:\n")
cat(sprintf("  Weeks:  %.1f (95%% CI: %.1f, %.1f)\n", 
            halflife_point, halflife_lower, halflife_upper))
cat(sprintf("  Months: %.1f (95%% CI: %.1f, %.1f)\n", 
            halflife_point/4.33, halflife_lower/4.33, halflife_upper/4.33))
cat(sprintf("  Years:  %.1f (95%% CI: %.1f, %.1f)\n\n", 
            halflife_point/52, halflife_lower/52, halflife_upper/52))

# For manuscript
cat("FOR MANUSCRIPT:\n")
cat(sprintf("ρ = %.5f per week (95%% CI: %.5f, %.5f)\n", rho, rho_lower, rho_upper))

