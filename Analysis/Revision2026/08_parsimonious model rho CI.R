# calculate_parsimonious_rho_ci.R

library(bbmle)

results <- readRDS("parsimonious_model_estimate_rho.rds")
fit <- results$fit
coefs <- coef(fit)
vcov_mat <- vcov(fit)

# Parsimonious model uses 'logit_rho' (not 'logit_rho_scaled')
logit_rho <- coefs["logit_rho"]
se_logit_rho <- sqrt(vcov_mat["logit_rho", "logit_rho"])

rho <- plogis(logit_rho) * 0.01

z_crit <- qnorm(0.975)
rho_lower <- plogis(logit_rho - z_crit * se_logit_rho) * 0.01
rho_upper <- plogis(logit_rho + z_crit * se_logit_rho) * 0.01

cat(sprintf("\nParsimonious model ρ:\n"))
cat(sprintf("  Point estimate: %.5f per week\n", rho))
cat(sprintf("  95%% CI: (%.5f, %.5f)\n\n", rho_lower, rho_upper))

cat(sprintf("FOR MANUSCRIPT (line 105):\n"))
cat(sprintf("  ρ = %.5f per week, 95%% CI: %.5f, %.5f\n", rho, rho_lower, rho_upper))

