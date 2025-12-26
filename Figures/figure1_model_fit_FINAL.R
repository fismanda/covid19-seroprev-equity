################################################################################
# Figure 1: Model Fit - Observed vs Predicted Seroprevalence
################################################################################
# 
# PURPOSE:
#   Generate Figure 1 showing observed CBS seroprevalence data (points) and
#   saturated model predictions (lines) for all five material deprivation 
#   quintiles from April 2021 through April 2023.
#
# INPUTS:
#   - ~/Dropbox/Family Room/CITF data/dataverse_files/seroprevalence_by_social_determinant.csv
#   - saturated_model_estimate_rho.rds (fitted model object)
#
# OUTPUTS:
#   - Figure1_Model_Fit.png (300 dpi)
#   - Figure1_Model_Fit.pdf (vector)
#
# AUTHOR: Hassan, Nassrallah, Fisman
# DATE: December 2025
################################################################################

library(ggplot2)
library(dplyr)
library(lubridate)

cat("\n=== GENERATING FIGURE 1: MODEL FIT ===\n\n")

# -----------------------------------------------------------------------------
# 1. LOAD AND PREPARE DATA
# -----------------------------------------------------------------------------

# Load CBS seroprevalence data
ses_data <- read.csv("~/Dropbox/Family Room/CITF data/dataverse_files/seroprevalence_by_social_determinant.csv")
ses_data$date <- mdy(ses_data$samplingdate)

# Filter to CBS, anti-N, quintile data
cbs_data <- ses_data %>%
  filter(ab_target == "N",
         grepl("^Q[1-5]", population),
         project == "CBS") %>%
  arrange(date)

# Create clean quintile labels FIRST (before filtering)
cbs_data$quintile_clean <- case_when(
  cbs_data$population == "Q1 (least deprived)" ~ "Q1",
  cbs_data$population == "Q2" ~ "Q2",
  cbs_data$population == "Q3" ~ "Q3",
  cbs_data$population == "Q4" ~ "Q4",
  cbs_data$population == "Q5 (most deprived)" ~ "Q5",
  TRUE ~ NA_character_  # Catch any unexpected values
)

# Check that all quintiles are present
cat("Quintiles in data:\n")
print(table(cbs_data$quintile_clean, useNA = "ifany"))

# Restrict to April 2021+ (when all quintiles have data)
fitting_start_date <- as.Date("2021-04-21")
plot_data_obs <- cbs_data %>%
  filter(date >= fitting_start_date,
         !is.na(quintile_clean)) %>%
  select(date, seroprev_est, quintile_clean)

cat(sprintf("Loaded %d observations from %s to %s\n", 
            nrow(plot_data_obs), min(plot_data_obs$date), max(plot_data_obs$date)))

# -----------------------------------------------------------------------------
# 2. LOAD MODEL AND EXTRACT PARAMETERS
# -----------------------------------------------------------------------------

results <- readRDS("saturated_model_estimate_rho.rds")
coefs <- coef(results$fit)

# Extract force of infection parameters
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

rho <- results$rho
omicron_date <- results$omicron_date
model_start_date <- results$model_start_date

cat("Parameters extracted successfully\n")
cat(sprintf("  ρ = %.5f per week\n", rho))
cat(sprintf("  Model starts: %s\n", model_start_date))
cat(sprintf("  Omicron date: %s\n\n", omicron_date))

# -----------------------------------------------------------------------------
# 3. SIMULATE MODEL PREDICTIONS
# -----------------------------------------------------------------------------

# SI model simulation function
simulate_si_saturated <- function(lambda_pre, lambda_omi, rho, 
                                   target_dates, model_start_date, omicron_date) {
  # Create complete date sequence from model start
  all_dates_seq <- seq(model_start_date, max(target_dates), by = "7 days")
  num_steps <- length(all_dates_seq)
  
  susceptible <- numeric(num_steps)
  infected <- numeric(num_steps)
  
  # Initialize at March 2020
  susceptible[1] <- 1
  infected[1] <- 0
  
  # Run model forward
  for (t in 2:num_steps) {
    if (all_dates_seq[t] < omicron_date) {
      r1 <- lambda_pre
    } else {
      r1 <- lambda_omi
    }
    r2 <- rho
    
    new_susceptible <- susceptible[t - 1] * exp(-r1) + infected[t - 1] * (1 - exp(-r2))
    new_infected <- infected[t - 1] * exp(-r2) + susceptible[t - 1] * (1 - exp(-r1))
    
    susceptible[t] <- new_susceptible
    infected[t] <- new_infected
  }
  
  # Extract predictions for target dates
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
  
  return(predictions * 100)  # Convert to percentage
}

# Generate predictions for each quintile
quintiles <- c("Q1", "Q2", "Q3", "Q4", "Q5")
lambdas_pre <- c(lambda_Q1_pre, lambda_Q2_pre, lambda_Q3_pre, lambda_Q4_pre, lambda_Q5_pre)
lambdas_omi <- c(lambda_Q1_omi, lambda_Q2_omi, lambda_Q3_omi, lambda_Q4_omi, lambda_Q5_omi)

prediction_data <- data.frame()

for (i in 1:5) {
  # Get observed dates for this quintile
  obs_dates <- plot_data_obs %>%
    filter(quintile_clean == quintiles[i]) %>%
    pull(date)
  
  if (length(obs_dates) > 0) {
    # Create dense date sequence for smooth lines
    date_seq <- seq(min(obs_dates), max(obs_dates), by = "7 days")
    
    # Simulate model
    pred <- simulate_si_saturated(lambdas_pre[i], lambdas_omi[i], rho,
                                   date_seq, model_start_date, omicron_date)
    
    pred_df <- data.frame(
      date = date_seq,
      predicted = pred,
      quintile_clean = quintiles[i]
    )
    
    prediction_data <- rbind(prediction_data, pred_df)
    
    cat(sprintf("Generated predictions for %s: %d timepoints\n", quintiles[i], length(date_seq)))
  }
}

# -----------------------------------------------------------------------------
# 4. CREATE LABELS FOR PLOTTING
# -----------------------------------------------------------------------------

plot_data_obs$Quintile <- factor(plot_data_obs$quintile_clean,
                                  levels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                                  labels = c("Q1 (Least Deprived)", "Q2", "Q3", "Q4", "Q5 (Most Deprived)"))

prediction_data$Quintile <- factor(prediction_data$quintile_clean,
                                    levels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                                    labels = c("Q1 (Least Deprived)", "Q2", "Q3", "Q4", "Q5 (Most Deprived)"))

# -----------------------------------------------------------------------------
# 5. CREATE FIGURE
# -----------------------------------------------------------------------------

p <- ggplot() +
  # Model predictions (smooth lines)
  geom_line(data = prediction_data,
            aes(x = date, y = predicted, color = Quintile),
            linewidth = 1.2, alpha = 0.9) +
  # Observed data (points)
  geom_point(data = plot_data_obs,
             aes(x = date, y = seroprev_est, color = Quintile),
             size = 2.5, alpha = 0.7) +
  # Omicron emergence line
  geom_vline(xintercept = as.numeric(omicron_date),
             linetype = "dashed", color = "gray30", linewidth = 0.8) +
  annotate("text", x = omicron_date, y = 85,
           label = "Omicron Emergence\n(Jan 1, 2022)",
           hjust = -0.05, size = 4, color = "gray20") +
  # Color scale: blue (Q1) to red (Q5)
  scale_color_manual(values = c("#2166ac", "#67a9cf", "#d1e5f0", "#fddbc7", "#b2182b")) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y",
               limits = c(as.Date("2021-04-01"), as.Date("2023-05-01"))) +
  scale_y_continuous(limits = c(0, 90), breaks = seq(0, 90, 15)) +
  labs(
    title = "Infection-Induced SARS-CoV-2 Seroprevalence by Material Deprivation Quintile",
    subtitle = "Observed data (points) and model predictions (lines) from Canadian Blood Services, April 2021-April 2023",
    x = "Date",
    y = "Seroprevalence (%)",
    color = "Material Deprivation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30", lineheight = 1.2),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  guides(color = guide_legend(nrow = 1))

# -----------------------------------------------------------------------------
# 6. SAVE FIGURE
# -----------------------------------------------------------------------------

ggsave("Figure1_Model_Fit.png", p,
       width = 12, height = 8, dpi = 300, bg = "white")

ggsave("Figure1_Model_Fit.pdf", p,
       width = 12, height = 8, bg = "white")

cat("\n✓ Figure 1 saved:\n")
cat("  - Figure1_Model_Fit.png (300 dpi)\n")
cat("  - Figure1_Model_Fit.pdf (vector)\n")

# -----------------------------------------------------------------------------
# 7. CALCULATE AND DISPLAY FIT STATISTICS
# -----------------------------------------------------------------------------

cat("\n=== MODEL FIT QUALITY ===\n\n")

for (i in 1:5) {
  obs <- plot_data_obs %>%
    filter(quintile_clean == quintiles[i])
  
  if (nrow(obs) > 0) {
    pred <- simulate_si_saturated(lambdas_pre[i], lambdas_omi[i], rho,
                                   obs$date, model_start_date, omicron_date)
    
    ss_res <- sum((obs$seroprev_est - pred)^2)
    ss_tot <- sum((obs$seroprev_est - mean(obs$seroprev_est))^2)
    r2 <- 1 - (ss_res / ss_tot)
    rmse <- sqrt(mean((obs$seroprev_est - pred)^2))
    
    cat(sprintf("%s: R² = %.4f, RMSE = %.2f%%, n = %d\n",
                quintiles[i], r2, rmse, nrow(obs)))
  }
}

cat("\n✓ Done!\n\n")

