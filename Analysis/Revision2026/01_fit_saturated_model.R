################################################################################
# Fit Saturated Model: Quintile-Specific Force of Infection (Pre-Omicron & Omicron)
################################################################################
#
# PURPOSE:
#   Estimate force of infection for each material deprivation quintile during
#   pre-Omicron (March 2020 - December 2021) and Omicron (January 2022+) periods.
#   This "saturated" model allows Omicron's effect to vary by socioeconomic status.
#
# MODEL STRUCTURE:
#   - 11 parameters total:
#     * 5 pre-Omicron λ (Q1-Q5)
#     * 5 Omicron λ (Q1-Q5)  
#     * 1 seroreversion rate ρ
#   - SI compartmental model with seroreversion
#   - Model starts March 1, 2020 (S=1, I=0)
#   - Fits to data from April 21, 2021 onwards (when all quintiles available)
#
# INPUTS:
#   - ~/Dropbox/Family Room/CITF data/dataverse_files/seroprevalence_by_social_determinant.csv
#
# OUTPUTS:
#   - saturated_model_estimate_rho.rds (fitted model object with parameters)
#
# AUTHOR: Hassan, Nassrallah, Fisman
# DATE: December 2025
################################################################################

library(bbmle)
library(dplyr)
library(lubridate)

cat("\n====================================================================")
cat("\n=== FITTING SATURATED MODEL: QUINTILE-SPECIFIC OMICRON EFFECTS ===")
cat("\n====================================================================\n\n")

# -----------------------------------------------------------------------------
# 1. LOAD AND PREPARE DATA
# -----------------------------------------------------------------------------

ses_data <- read.csv("~/Dropbox/Family Room/CITF data/dataverse_files/seroprevalence_by_social_determinant.csv")
ses_data$date <- mdy(ses_data$samplingdate)

# Filter to CBS, anti-N, quintile data
cbs_data <- ses_data %>%
  filter(ab_target == "N",
         grepl("^Q[1-5]", population),
         project == "CBS") %>%
  arrange(date)

# Create clean quintile labels
cbs_data$quintile_clean <- case_when(
  cbs_data$population == "Q1 (least deprived)" ~ "Q1",
  cbs_data$population == "Q2" ~ "Q2",
  cbs_data$population == "Q3" ~ "Q3",
  cbs_data$population == "Q4" ~ "Q4",
  cbs_data$population == "Q5 (most deprived)" ~ "Q5",
  TRUE ~ NA_character_
)

cat("Data loaded successfully\n")
cat(sprintf("Total observations: %d\n", nrow(cbs_data)))
cat(sprintf("Date range: %s to %s\n\n", min(cbs_data$date), max(cbs_data$date)))

# Filter to April 2021 onwards (when all quintiles available)
fitting_start_date <- as.Date("2021-04-21")
fitting_data <- cbs_data %>%
  filter(date >= fitting_start_date,
         !is.na(quintile_clean))

cat("=== DATA FOR FITTING ===\n")
cat(sprintf("Fitting start date: %s\n", fitting_start_date))
cat(sprintf("Observations used for fitting: %d\n\n", nrow(fitting_data)))

for (q in c("Q1", "Q2", "Q3", "Q4", "Q5")) {
  qdata <- fitting_data %>% filter(quintile_clean == q)
  cat(sprintf("  %s: %d observations from %s to %s\n", 
              q, nrow(qdata), min(qdata$date), max(qdata$date)))
}
cat("\n")

# -----------------------------------------------------------------------------
# 2. DEFINE KEY DATES
# -----------------------------------------------------------------------------

model_start_date <- as.Date("2020-03-01")  # Model starts here (S=1, I=0)
omicron_date <- as.Date("2022-01-01")      # Omicron emergence

cat(sprintf("Model start date: %s (S=1, I=0)\n", model_start_date))
cat(sprintf("Omicron emergence: %s\n\n", omicron_date))

# -----------------------------------------------------------------------------
# 3. PREPARE DATA FOR LIKELIHOOD CALCULATION
# -----------------------------------------------------------------------------

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
# 4. DEFINE SI MODEL WITH SEROREVERSION
# -----------------------------------------------------------------------------

# SI model simulation from March 2020
simulate_si_saturated <- function(lambda_pre, lambda_omi, rho, 
                                   target_dates, model_start_date, omicron_date) {
  
  # Create complete date sequence from model start to last observation
  all_dates_seq <- seq(model_start_date, max(target_dates), by = "7 days")
  num_steps <- length(all_dates_seq)
  
  susceptible <- numeric(num_steps)
  infected <- numeric(num_steps)
  
  # Initialize at March 2020
  susceptible[1] <- 1
  infected[1] <- 0
  
  # Run model forward from March 2020
  for (t in 2:num_steps) {
    # Switch force of infection at Omicron emergence
    if (all_dates_seq[t] < omicron_date) {
      r1 <- lambda_pre
    } else {
      r1 <- lambda_omi
    }
    r2 <- rho
    
    # SI dynamics with seroreversion
    new_susceptible <- susceptible[t - 1] * exp(-r1) + infected[t - 1] * (1 - exp(-r2))
    new_infected <- infected[t - 1] * exp(-r2) + susceptible[t - 1] * (1 - exp(-r1))
    
    susceptible[t] <- new_susceptible
    infected[t] <- new_infected
  }
  
  # Extract predictions for target dates only
  predictions <- numeric(length(target_dates))
  for (i in 1:length(target_dates)) {
    idx <- which(all_dates_seq == target_dates[i])
    if (length(idx) > 0) {
      predictions[i] <- infected[idx]
    } else {
      # Interpolate if exact date not in sequence
      closest_idx <- which.min(abs(all_dates_seq - target_dates[i]))
      predictions[i] <- infected[closest_idx]
    }
  }
  
  return(predictions)
}

# -----------------------------------------------------------------------------
# 5. DEFINE LIKELIHOOD FUNCTION
# -----------------------------------------------------------------------------

# Joint log-likelihood across all quintiles
joint_loglik_saturated <- function(lambda_Q1_pre, lambda_Q2_pre, lambda_Q3_pre, 
                                     lambda_Q4_pre, lambda_Q5_pre,
                                     lambda_Q1_omi, lambda_Q2_omi, lambda_Q3_omi,
                                     lambda_Q4_omi, lambda_Q5_omi, rho) {
  
  lambdas_pre <- c(lambda_Q1_pre, lambda_Q2_pre, lambda_Q3_pre, lambda_Q4_pre, lambda_Q5_pre)
  lambdas_omi <- c(lambda_Q1_omi, lambda_Q2_omi, lambda_Q3_omi, lambda_Q4_omi, lambda_Q5_omi)
  
  total_ll <- 0
  
  for (i in 1:5) {
    # Simulate model for this quintile
    pred <- simulate_si_saturated(
      lambdas_pre[i], lambdas_omi[i], rho,
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

# Objective function (negative log-likelihood)
f_saturated <- function(logit_lambda_Q1_pre, logit_lambda_Q2_pre, logit_lambda_Q3_pre,
                         logit_lambda_Q4_pre, logit_lambda_Q5_pre,
                         logit_lambda_Q1_omi, logit_lambda_Q2_omi, logit_lambda_Q3_omi,
                         logit_lambda_Q4_omi, logit_lambda_Q5_omi,
                         logit_rho) {
  
  # Transform parameters from logit scale
  lambda_Q1_pre <- plogis(logit_lambda_Q1_pre)
  lambda_Q2_pre <- plogis(logit_lambda_Q2_pre)
  lambda_Q3_pre <- plogis(logit_lambda_Q3_pre)
  lambda_Q4_pre <- plogis(logit_lambda_Q4_pre)
  lambda_Q5_pre <- plogis(logit_lambda_Q5_pre)
  
  lambda_Q1_omi <- plogis(logit_lambda_Q1_omi)
  lambda_Q2_omi <- plogis(logit_lambda_Q2_omi)
  lambda_Q3_omi <- plogis(logit_lambda_Q3_omi)
  lambda_Q4_omi <- plogis(logit_lambda_Q4_omi)
  lambda_Q5_omi <- plogis(logit_lambda_Q5_omi)
  
  # Transform rho to [0, 0.01] range (0-1% seroreversion per week)
  rho <- plogis(logit_rho) * 0.01
  
  # Return negative log-likelihood
  -joint_loglik_saturated(lambda_Q1_pre, lambda_Q2_pre, lambda_Q3_pre, 
                           lambda_Q4_pre, lambda_Q5_pre,
                           lambda_Q1_omi, lambda_Q2_omi, lambda_Q3_omi,
                           lambda_Q4_omi, lambda_Q5_omi, rho)
}

# -----------------------------------------------------------------------------
# 6. SET INITIAL PARAMETER GUESSES
# -----------------------------------------------------------------------------

guess_sat <- list(
  logit_lambda_Q1_pre = qlogis(0.0005),
  logit_lambda_Q2_pre = qlogis(0.0006),
  logit_lambda_Q3_pre = qlogis(0.0006),
  logit_lambda_Q4_pre = qlogis(0.0008),
  logit_lambda_Q5_pre = qlogis(0.001),
  logit_lambda_Q1_omi = qlogis(0.025),
  logit_lambda_Q2_omi = qlogis(0.025),
  logit_lambda_Q3_omi = qlogis(0.025),
  logit_lambda_Q4_omi = qlogis(0.028),
  logit_lambda_Q5_omi = qlogis(0.028),
  logit_rho = qlogis(0.1)  # Starts at rho = 0.001
)

# -----------------------------------------------------------------------------
# 7. FIT MODEL USING MAXIMUM LIKELIHOOD ESTIMATION
# -----------------------------------------------------------------------------

cat("=== FITTING SATURATED MODEL ===\n")
cat("Parameters: 5 λ_pre + 5 λ_omi + 1 ρ = 11 parameters\n")
cat("Starting maximum likelihood estimation...\n\n")

fit_saturated <- mle2(
  function(logit_lambda_Q1_pre, logit_lambda_Q2_pre, logit_lambda_Q3_pre,
           logit_lambda_Q4_pre, logit_lambda_Q5_pre,
           logit_lambda_Q1_omi, logit_lambda_Q2_omi, logit_lambda_Q3_omi,
           logit_lambda_Q4_omi, logit_lambda_Q5_omi, logit_rho)
    f_saturated(logit_lambda_Q1_pre, logit_lambda_Q2_pre, logit_lambda_Q3_pre,
                logit_lambda_Q4_pre, logit_lambda_Q5_pre,
                logit_lambda_Q1_omi, logit_lambda_Q2_omi, logit_lambda_Q3_omi,
                logit_lambda_Q4_omi, logit_lambda_Q5_omi, logit_rho),
  start = guess_sat,
  method = "L-BFGS-B"
)

cat("✓ Fitting complete!\n\n")

# -----------------------------------------------------------------------------
# 8. EXTRACT AND DISPLAY RESULTS
# -----------------------------------------------------------------------------

coefs <- coef(fit_saturated)

# Transform back to natural scale
lambda_Q1_pre <- plogis(coefs["logit_lambda_Q1_pre"])
lambda_Q2_pre <- plogis(coefs["logit_lambda_Q2_pre"])
lambda_Q3_pre <- plogis(coefs["logit_lambda_Q3_pre"])
lambda_Q4_pre <- plogis(coefs["logit_lambda_Q4_pre"])
lambda_Q5_pre <- plogis(coefs["logit_lambda_Q5_pre"])

lambda_Q1_omi <- plogis(coefs["logit_lambda_Q1_omi"])
lambda_Q2_omi <- plogis(coefs["logit_lambda_Q2_omi"])
lambda_Q3_omi <- plogis(coefs["logit_lambda_Q3_omi"])
lambda_Q4_omi <- plogis(coefs["logit_lambda_Q4_omi"])
lambda_Q5_omi <- plogis(coefs["logit_lambda_Q5_omi"])

rho <- plogis(coefs["logit_rho"]) * 0.01

cat("=== ESTIMATED SEROREVERSION RATE ===\n\n")
cat(sprintf("ρ = %.5f per week (%.3f%% per week)\n", rho, rho * 100))
cat(sprintf("Half-life = %.1f weeks (%.1f months)\n\n", log(2)/rho, log(2)/rho/4.33))

cat("=== PRE-OMICRON FORCE OF INFECTION ===\n\n")
cat(sprintf("Q1: %.5f per week\n", lambda_Q1_pre))
cat(sprintf("Q2: %.5f per week\n", lambda_Q2_pre))
cat(sprintf("Q3: %.5f per week\n", lambda_Q3_pre))
cat(sprintf("Q4: %.5f per week\n", lambda_Q4_pre))
cat(sprintf("Q5: %.5f per week\n", lambda_Q5_pre))

cat("\n=== OMICRON FORCE OF INFECTION ===\n\n")
cat(sprintf("Q1: %.5f per week\n", lambda_Q1_omi))
cat(sprintf("Q2: %.5f per week\n", lambda_Q2_omi))
cat(sprintf("Q3: %.5f per week\n", lambda_Q3_omi))
cat(sprintf("Q4: %.5f per week\n", lambda_Q4_omi))
cat(sprintf("Q5: %.5f per week\n", lambda_Q5_omi))

cat("\n=== OMICRON MULTIPLIERS ===\n\n")
cat(sprintf("Q1: %.2fx\n", lambda_Q1_omi / lambda_Q1_pre))
cat(sprintf("Q2: %.2fx\n", lambda_Q2_omi / lambda_Q2_pre))
cat(sprintf("Q3: %.2fx\n", lambda_Q3_omi / lambda_Q3_pre))
cat(sprintf("Q4: %.2fx\n", lambda_Q4_omi / lambda_Q4_pre))
cat(sprintf("Q5: %.2fx\n", lambda_Q5_omi / lambda_Q5_pre))

cat("\n=== PRE-OMICRON IRRs (vs Q1) ===\n\n")
cat(sprintf("Q1: 1.000 (reference)\n"))
cat(sprintf("Q2: %.3f\n", lambda_Q2_pre / lambda_Q1_pre))
cat(sprintf("Q3: %.3f\n", lambda_Q3_pre / lambda_Q1_pre))
cat(sprintf("Q4: %.3f\n", lambda_Q4_pre / lambda_Q1_pre))
cat(sprintf("Q5: %.3f\n", lambda_Q5_pre / lambda_Q1_pre))

cat("\n=== OMICRON IRRs (vs Q1) ===\n\n")
cat(sprintf("Q1: 1.000 (reference)\n"))
cat(sprintf("Q2: %.3f\n", lambda_Q2_omi / lambda_Q1_omi))
cat(sprintf("Q3: %.3f\n", lambda_Q3_omi / lambda_Q1_omi))
cat(sprintf("Q4: %.3f\n", lambda_Q4_omi / lambda_Q1_omi))
cat(sprintf("Q5: %.3f\n", lambda_Q5_omi / lambda_Q1_omi))

# -----------------------------------------------------------------------------
# 9. SAVE RESULTS
# -----------------------------------------------------------------------------

results <- list(
  fit = fit_saturated,
  loglik = logLik(fit_saturated),
  k = 11,
  aic = AIC(fit_saturated),
  bic = BIC(fit_saturated),
  rho = rho,
  model_start_date = model_start_date,
  fitting_start_date = fitting_start_date,
  omicron_date = omicron_date
)

saveRDS(results, "saturated_model_estimate_rho.rds")

cat("\n✓ Results saved: saturated_model_estimate_rho.rds\n\n")
cat("Model fit statistics:\n")
cat(sprintf("  Log-likelihood: %.2f\n", results$loglik))
cat(sprintf("  AIC: %.2f\n", results$aic))
cat(sprintf("  Parameters: %d\n\n", results$k))

