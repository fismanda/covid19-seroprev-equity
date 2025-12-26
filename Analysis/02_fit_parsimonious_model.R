################################################################################
# Fit Parsimonious Model: Shared Omicron Multiplier
################################################################################
#
# PURPOSE:
#   Estimate force of infection assuming Omicron affected all quintiles equally.
#   This "parsimonious" model has fewer parameters than the saturated model and
#   serves as a null hypothesis: did Omicron amplify transmission uniformly?
#
# MODEL STRUCTURE:
#   - 7 parameters total:
#     * 5 pre-Omicron λ (Q1-Q5)
#     * 1 shared Omicron multiplier m
#     * 1 seroreversion rate ρ
#   - SI compartmental model with seroreversion
#   - Model starts March 1, 2020 (S=1, I=0)
#   - Fits to data from April 21, 2021 onwards
#
# INPUTS:
#   - ~/Dropbox/Family Room/CITF data/dataverse_files/seroprevalence_by_social_determinant.csv
#
# OUTPUTS:
#   - parsimonious_model_estimate_rho.rds
#
# AUTHOR: Hassan, Nassrallah, Fisman
# DATE: December 2025
################################################################################

library(bbmle)
library(dplyr)
library(lubridate)

cat("\n================================================================")
cat("\n=== FITTING PARSIMONIOUS MODEL: SHARED OMICRON MULTIPLIER ===")
cat("\n================================================================\n\n")

# -----------------------------------------------------------------------------
# 1. LOAD AND PREPARE DATA (same as saturated model)
# -----------------------------------------------------------------------------

ses_data <- read.csv("~/Dropbox/Family Room/CITF data/dataverse_files/seroprevalence_by_social_determinant.csv")
ses_data$date <- mdy(ses_data$samplingdate)

cbs_data <- ses_data %>%
  filter(ab_target == "N",
         grepl("^Q[1-5]", population),
         project == "CBS") %>%
  arrange(date)

cbs_data$quintile_clean <- case_when(
  cbs_data$population == "Q1 (least deprived)" ~ "Q1",
  cbs_data$population == "Q2" ~ "Q2",
  cbs_data$population == "Q3" ~ "Q3",
  cbs_data$population == "Q4" ~ "Q4",
  cbs_data$population == "Q5 (most deprived)" ~ "Q5",
  TRUE ~ NA_character_
)

fitting_start_date <- as.Date("2021-04-21")
fitting_data <- cbs_data %>%
  filter(date >= fitting_start_date,
         !is.na(quintile_clean))

cat("Data loaded and filtered\n")
cat(sprintf("Observations for fitting: %d\n\n", nrow(fitting_data)))

# Define key dates
model_start_date <- as.Date("2020-03-01")
omicron_date <- as.Date("2022-01-01")

# Prepare data
quintiles <- c("Q1", "Q2", "Q3", "Q4", "Q5")
all_observed <- list()
all_dates <- list()
all_sample_sizes <- list()

for (i in 1:5) {
  qdata <- fitting_data %>%
    filter(quintile_clean == quintiles[i]) %>%
    arrange(date)
  all_observed[[i]] <- qdata$seroprev_est
  all_dates[[i]] <- qdata$date
  all_sample_sizes[[i]] <- qdata$sero_denom
}

# -----------------------------------------------------------------------------
# 2. DEFINE SI MODEL (with shared multiplier)
# -----------------------------------------------------------------------------

simulate_si_parsimonious <- function(lambda_pre, multiplier, rho,
                                      target_dates, model_start_date, omicron_date) {
  
  all_dates_seq <- seq(model_start_date, max(target_dates), by = "7 days")
  num_steps <- length(all_dates_seq)
  
  susceptible <- numeric(num_steps)
  infected <- numeric(num_steps)
  
  susceptible[1] <- 1
  infected[1] <- 0
  
  for (t in 2:num_steps) {
    # Apply shared multiplier after Omicron
    if (all_dates_seq[t] < omicron_date) {
      r1 <- lambda_pre
    } else {
      r1 <- lambda_pre * multiplier  # Shared multiplier
    }
    r2 <- rho
    
    new_susceptible <- susceptible[t - 1] * exp(-r1) + infected[t - 1] * (1 - exp(-r2))
    new_infected <- infected[t - 1] * exp(-r2) + susceptible[t - 1] * (1 - exp(-r1))
    
    susceptible[t] <- new_susceptible
    infected[t] <- new_infected
  }
  
  predictions <- numeric(length(target_dates))
  for (i in 1:length(target_dates)) {
    idx <- which(all_dates_seq == target_dates[i])
    if (length(idx) > 0) {
      predictions[i] <- infected[idx]
    } else {
      closest_idx <- which.min(abs(all_dates_seq - target_dates[i]))
      predictions[i] <- infected[closest_idx]
    }
  }
  
  return(predictions)
}

# -----------------------------------------------------------------------------
# 3. DEFINE LIKELIHOOD FUNCTION
# -----------------------------------------------------------------------------

joint_loglik_parsimonious <- function(lambda_Q1_pre, lambda_Q2_pre, lambda_Q3_pre,
                                       lambda_Q4_pre, lambda_Q5_pre,
                                       multiplier, rho) {
  
  lambdas_pre <- c(lambda_Q1_pre, lambda_Q2_pre, lambda_Q3_pre, lambda_Q4_pre, lambda_Q5_pre)
  
  total_ll <- 0
  
  for (i in 1:5) {
    pred <- simulate_si_parsimonious(
      lambdas_pre[i], multiplier, rho,
      all_dates[[i]], model_start_date, omicron_date
    )
    
    obs_prop <- all_observed[[i]] / 100
    n <- all_sample_sizes[[i]]
    
    # Binomial likelihood with actual sample sizes
    obs_positive <- round(n * obs_prop)
    ll <- sum(dbinom(x = obs_positive,
                     size = n,
                     prob = pred,
                     log = TRUE))
    total_ll <- total_ll + ll
  }
  
  return(total_ll)
}

# Objective function
f_parsimonious <- function(logit_lambda_Q1_pre, logit_lambda_Q2_pre, logit_lambda_Q3_pre,
                            logit_lambda_Q4_pre, logit_lambda_Q5_pre,
                            log_multiplier, logit_rho) {
  
  lambda_Q1_pre <- plogis(logit_lambda_Q1_pre)
  lambda_Q2_pre <- plogis(logit_lambda_Q2_pre)
  lambda_Q3_pre <- plogis(logit_lambda_Q3_pre)
  lambda_Q4_pre <- plogis(logit_lambda_Q4_pre)
  lambda_Q5_pre <- plogis(logit_lambda_Q5_pre)
  
  multiplier <- exp(log_multiplier)  # Ensure positive
  rho <- plogis(logit_rho) * 0.01
  
  -joint_loglik_parsimonious(lambda_Q1_pre, lambda_Q2_pre, lambda_Q3_pre,
                              lambda_Q4_pre, lambda_Q5_pre,
                              multiplier, rho)
}

# -----------------------------------------------------------------------------
# 4. SET INITIAL GUESSES
# -----------------------------------------------------------------------------

guess_pars <- list(
  logit_lambda_Q1_pre = qlogis(0.0005),
  logit_lambda_Q2_pre = qlogis(0.0006),
  logit_lambda_Q3_pre = qlogis(0.0006),
  logit_lambda_Q4_pre = qlogis(0.0008),
  logit_lambda_Q5_pre = qlogis(0.001),
  log_multiplier = log(35),  # Starts at ~35x multiplier
  logit_rho = qlogis(0.1)
)

# -----------------------------------------------------------------------------
# 5. FIT MODEL
# -----------------------------------------------------------------------------

cat("=== FITTING PARSIMONIOUS MODEL ===\n")
cat("Parameters: 5 λ_pre + 1 shared multiplier + 1 ρ = 7 parameters\n")
cat("Starting maximum likelihood estimation...\n\n")

fit_parsimonious <- mle2(
  function(logit_lambda_Q1_pre, logit_lambda_Q2_pre, logit_lambda_Q3_pre,
           logit_lambda_Q4_pre, logit_lambda_Q5_pre,
           log_multiplier, logit_rho)
    f_parsimonious(logit_lambda_Q1_pre, logit_lambda_Q2_pre, logit_lambda_Q3_pre,
                   logit_lambda_Q4_pre, logit_lambda_Q5_pre,
                   log_multiplier, logit_rho),
  start = guess_pars,
  method = "L-BFGS-B"
)

cat("✓ Fitting complete!\n\n")

# -----------------------------------------------------------------------------
# 6. EXTRACT AND DISPLAY RESULTS
# -----------------------------------------------------------------------------

coefs <- coef(fit_parsimonious)

lambda_Q1_pre <- plogis(coefs["logit_lambda_Q1_pre"])
lambda_Q2_pre <- plogis(coefs["logit_lambda_Q2_pre"])
lambda_Q3_pre <- plogis(coefs["logit_lambda_Q3_pre"])
lambda_Q4_pre <- plogis(coefs["logit_lambda_Q4_pre"])
lambda_Q5_pre <- plogis(coefs["logit_lambda_Q5_pre"])

multiplier <- exp(coefs["log_multiplier"])
rho <- plogis(coefs["logit_rho"]) * 0.01

cat("=== ESTIMATED SEROREVERSION RATE ===\n\n")
cat(sprintf("ρ = %.5f per week (%.3f%% per week)\n", rho, rho * 100))
cat(sprintf("Half-life = %.1f weeks (%.1f months)\n\n", log(2)/rho, log(2)/rho/4.33))

cat("=== SHARED OMICRON MULTIPLIER ===\n\n")
cat(sprintf("Multiplier = %.2fx (applied to all quintiles)\n\n", multiplier))

cat("=== PRE-OMICRON FORCE OF INFECTION ===\n\n")
cat(sprintf("Q1: %.5f per week\n", lambda_Q1_pre))
cat(sprintf("Q2: %.5f per week\n", lambda_Q2_pre))
cat(sprintf("Q3: %.5f per week\n", lambda_Q3_pre))
cat(sprintf("Q4: %.5f per week\n", lambda_Q4_pre))
cat(sprintf("Q5: %.5f per week\n", lambda_Q5_pre))

cat("\n=== OMICRON FORCE OF INFECTION (λ_pre × multiplier) ===\n\n")
cat(sprintf("Q1: %.5f per week\n", lambda_Q1_pre * multiplier))
cat(sprintf("Q2: %.5f per week\n", lambda_Q2_pre * multiplier))
cat(sprintf("Q3: %.5f per week\n", lambda_Q3_pre * multiplier))
cat(sprintf("Q4: %.5f per week\n", lambda_Q4_pre * multiplier))
cat(sprintf("Q5: %.5f per week\n", lambda_Q5_pre * multiplier))

cat("\n=== PRE-OMICRON IRRs (vs Q1) ===\n\n")
cat(sprintf("Q1: 1.000 (reference)\n"))
cat(sprintf("Q2: %.3f\n", lambda_Q2_pre / lambda_Q1_pre))
cat(sprintf("Q3: %.3f\n", lambda_Q3_pre / lambda_Q1_pre))
cat(sprintf("Q4: %.3f\n", lambda_Q4_pre / lambda_Q1_pre))
cat(sprintf("Q5: %.3f\n", lambda_Q5_pre / lambda_Q1_pre))

cat("\n=== OMICRON IRRs (vs Q1) - SAME AS PRE-OMICRON ===\n")
cat("(Because shared multiplier preserves gradient)\n\n")
cat(sprintf("Q1: 1.000 (reference)\n"))
cat(sprintf("Q2: %.3f\n", lambda_Q2_pre / lambda_Q1_pre))
cat(sprintf("Q3: %.3f\n", lambda_Q3_pre / lambda_Q1_pre))
cat(sprintf("Q4: %.3f\n", lambda_Q4_pre / lambda_Q1_pre))
cat(sprintf("Q5: %.3f\n", lambda_Q5_pre / lambda_Q1_pre))

# -----------------------------------------------------------------------------
# 7. SAVE RESULTS
# -----------------------------------------------------------------------------

results <- list(
  fit = fit_parsimonious,
  loglik = logLik(fit_parsimonious),
  k = 7,
  aic = AIC(fit_parsimonious),
  bic = BIC(fit_parsimonious),
  rho = rho,
  multiplier = multiplier,
  model_start_date = model_start_date,
  fitting_start_date = fitting_start_date,
  omicron_date = omicron_date
)

saveRDS(results, "parsimonious_model_estimate_rho.rds")

cat("\n✓ Results saved: parsimonious_model_estimate_rho.rds\n\n")
cat("Model fit statistics:\n")
cat(sprintf("  Log-likelihood: %.2f\n", results$loglik))
cat(sprintf("  AIC: %.2f\n", results$aic))
cat(sprintf("  Parameters: %d\n\n", results$k))

