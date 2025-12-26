################################################################################
# Figure S1: Parsimonious Model Fit
################################################################################
#
# PURPOSE:
#   Create supplementary figure showing model fit for the parsimonious model
#   (shared Omicron multiplier). For comparison with Figure 1 (saturated model).
#   Demonstrates that while parsimonious model fits reasonably, it cannot
#   capture gradient compression as well as the saturated model.
#
# INPUTS:
#   - ~/Dropbox/Family Room/CITF data/dataverse_files/seroprevalence_by_social_determinant.csv
#   - parsimonious_model_estimate_rho.rds
#
# OUTPUTS:
#   - FigureS1_Parsimonious_Fit.png (300 dpi)
#   - FigureS1_Parsimonious_Fit.pdf (vector)
#
# AUTHOR: Hassan, Nassrallah, Fisman
# DATE: December 2025
################################################################################

library(dplyr)
library(lubridate)
library(ggplot2)

cat("\n========================================================")
cat("\n=== CREATING FIGURE S1: PARSIMONIOUS MODEL FIT ===")
cat("\n========================================================\n\n")

# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------

ses_data <- read.csv("~/Dropbox/Family Room/CITF data/dataverse_files/seroprevalence_by_social_determinant.csv")
ses_data$date <- mdy(ses_data$samplingdate)

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

fitting_start_date <- as.Date("2021-04-21")
plot_data_obs <- cbs_data %>%
  filter(date >= fitting_start_date,
         !is.na(quintile_clean)) %>%
  select(date, seroprev_est, quintile_clean)

cat("✓ Observed data loaded\n")
cat(sprintf("  Observations: %d\n", nrow(plot_data_obs)))

# -----------------------------------------------------------------------------
# 2. LOAD PARSIMONIOUS MODEL AND GENERATE PREDICTIONS
# -----------------------------------------------------------------------------

results <- readRDS("parsimonious_model_estimate_rho.rds")
coefs <- coef(results$fit)

# Extract parameters
lambda_Q1_pre <- plogis(coefs["logit_lambda_Q1_pre"])
lambda_Q2_pre <- plogis(coefs["logit_lambda_Q2_pre"])
lambda_Q3_pre <- plogis(coefs["logit_lambda_Q3_pre"])
lambda_Q4_pre <- plogis(coefs["logit_lambda_Q4_pre"])
lambda_Q5_pre <- plogis(coefs["logit_lambda_Q5_pre"])

multiplier <- exp(coefs["log_multiplier"])
rho <- plogis(coefs["logit_rho"]) * 0.01

cat("✓ Parsimonious model loaded\n")
cat(sprintf("  Shared multiplier: %.2fx\n", multiplier))
cat(sprintf("  Seroreversion: %.5f/week\n\n", rho))

# Model dates
model_start_date <- results$model_start_date
omicron_date <- results$omicron_date

# SI model function
simulate_si_parsimonious <- function(lambda_pre, multiplier, rho,
                                      model_start_date, omicron_date, end_date) {
  
  all_dates_seq <- seq(model_start_date, end_date, by = "7 days")
  num_steps <- length(all_dates_seq)
  
  susceptible <- numeric(num_steps)
  infected <- numeric(num_steps)
  
  susceptible[1] <- 1
  infected[1] <- 0
  
  for (t in 2:num_steps) {
    if (all_dates_seq[t] < omicron_date) {
      r1 <- lambda_pre
    } else {
      r1 <- lambda_pre * multiplier
    }
    r2 <- rho
    
    new_susceptible <- susceptible[t - 1] * exp(-r1) + infected[t - 1] * (1 - exp(-r2))
    new_infected <- infected[t - 1] * exp(-r2) + susceptible[t - 1] * (1 - exp(-r1))
    
    susceptible[t] <- new_susceptible
    infected[t] <- new_infected
  }
  
  data.frame(
    date = all_dates_seq,
    seroprevalence = infected * 100
  )
}

# Generate predictions for all quintiles
lambdas_pre <- c(lambda_Q1_pre, lambda_Q2_pre, lambda_Q3_pre, lambda_Q4_pre, lambda_Q5_pre)
quintiles <- c("Q1", "Q2", "Q3", "Q4", "Q5")

end_date <- max(plot_data_obs$date)

plot_data_pred <- data.frame()
for (i in 1:5) {
  pred <- simulate_si_parsimonious(lambdas_pre[i], multiplier, rho,
                                    model_start_date, omicron_date, end_date)
  pred$quintile_clean <- quintiles[i]
  plot_data_pred <- rbind(plot_data_pred, pred)
}

cat("✓ Model predictions generated\n\n")

# -----------------------------------------------------------------------------
# 3. CALCULATE R² FOR EACH QUINTILE
# -----------------------------------------------------------------------------

r_squared <- data.frame()
for (q in quintiles) {
  obs <- plot_data_obs %>% 
    filter(quintile_clean == q) %>%
    arrange(date)
  
  pred <- plot_data_pred %>% 
    filter(quintile_clean == q) %>%
    arrange(date)
  
  # Match prediction dates to observation dates
  pred_matched <- pred %>%
    filter(date %in% obs$date) %>%
    arrange(date)
  
  obs_matched <- obs %>%
    filter(date %in% pred_matched$date) %>%
    arrange(date)
  
  # Calculate R²
  ss_res <- sum((obs_matched$seroprev_est - pred_matched$seroprevalence)^2)
  ss_tot <- sum((obs_matched$seroprev_est - mean(obs_matched$seroprev_est))^2)
  r2 <- 1 - (ss_res / ss_tot)
  
  r_squared <- rbind(r_squared, data.frame(quintile = q, r2 = r2))
}

cat("Model fit (R²):\n")
for (i in 1:nrow(r_squared)) {
  cat(sprintf("  %s: %.4f\n", r_squared$quintile[i], r_squared$r2[i]))
}
cat("\n")

# -----------------------------------------------------------------------------
# 4. CREATE FIGURE
# -----------------------------------------------------------------------------

cat("Creating figure...\n")

# Prepare data with labels
plot_data_obs_labeled <- plot_data_obs %>%
  mutate(quintile_label = case_when(
    quintile_clean == "Q1" ~ "Q1 (Least Deprived)",
    quintile_clean == "Q2" ~ "Q2",
    quintile_clean == "Q3" ~ "Q3",
    quintile_clean == "Q4" ~ "Q4",
    quintile_clean == "Q5" ~ "Q5 (Most Deprived)"
  ))

plot_data_pred_labeled <- plot_data_pred %>%
  mutate(quintile_label = case_when(
    quintile_clean == "Q1" ~ "Q1 (Least Deprived)",
    quintile_clean == "Q2" ~ "Q2",
    quintile_clean == "Q3" ~ "Q3",
    quintile_clean == "Q4" ~ "Q4",
    quintile_clean == "Q5" ~ "Q5 (Most Deprived)"
  ))

p <- ggplot() +
  
  # Model predictions (lines)
  geom_line(data = plot_data_pred_labeled,
            aes(x = date, y = seroprevalence, color = quintile_label),
            linewidth = 1) +
  
  # Observed data (points)
  geom_point(data = plot_data_obs_labeled,
             aes(x = date, y = seroprev_est, color = quintile_label),
             size = 2, alpha = 0.7) +
  
  # Omicron emergence line
  geom_vline(xintercept = as.numeric(omicron_date),
             linetype = "dashed", color = "gray40", linewidth = 0.7) +
  
  # Omicron label
  annotate("text", x = omicron_date, y = 90,
           label = "Omicron Emergence\n(Jan 1, 2022)",
           hjust = -0.05, vjust = 1, size = 3.5, color = "gray40") +
  
  # Color scheme (blue to red gradient)
  scale_color_manual(
    values = c(
      "Q1 (Least Deprived)" = "#2E86AB",
      "Q2" = "#5FA5C7",
      "Q3" = "#8BC4E3",
      "Q4" = "#C77DA5",
      "Q5 (Most Deprived)" = "#A23B72"
    ),
    name = "Material Deprivation"
  ) +
  
  # Labels
  labs(
    x = "Date",
    y = "Seroprevalence (%)",
    title = "Parsimonious Model Fit: Shared Omicron Multiplier",
    subtitle = sprintf("Observed data (points) and model predictions (lines) - Shared multiplier: %.2fx", multiplier),
    caption = "Canadian Blood Services data, April 2021 - April 2023"
  ) +
  
  # Scales
  scale_y_continuous(breaks = seq(0, 90, 15), limits = c(0, 95)) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    plot.caption = element_text(size = 9, color = "gray50", hjust = 0),
    plot.margin = margin(10, 10, 10, 10)
  )

# Save as PNG (300 dpi)
ggsave("FigureS1_Parsimonious_Fit.png",
       plot = p,
       width = 10,
       height = 6,
       dpi = 300,
       bg = "white")

# Save as PDF (vector)
ggsave("FigureS1_Parsimonious_Fit.pdf",
       plot = p,
       width = 10,
       height = 6)

cat("✓ Figure saved:\n")
cat("  - FigureS1_Parsimonious_Fit.png (300 dpi)\n")
cat("  - FigureS1_Parsimonious_Fit.pdf (vector)\n\n")

cat("NOTE: Compare to Figure 1 (saturated model, ΔAIC = -991).\n")
cat("The parsimonious model constrains gradient to be preserved,\n")
cat("while saturated model allows differential amplification.\n\n")

