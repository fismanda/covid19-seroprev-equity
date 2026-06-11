################################################################################
# Figure S1: Parsimonious model fit with parametric bootstrap CI bands
################################################################################

library(ggplot2)
library(dplyr)
library(lubridate)
library(MASS)

cat("\n=== CREATING FIGURE S1: PARSIMONIOUS MODEL FIT WITH CI BANDS ===\n\n")

# ---- 1. LOAD DATA --------------------------------------------------------

ses_data <- read.csv("~/Dropbox/Family Room/CITF data/dataverse_files/seroprevalence_by_social_determinant.csv")
ses_data$date <- mdy(ses_data$samplingdate)

cbs_data <- ses_data %>%
  filter(ab_target == "N", grepl("^Q[1-5]", population), project == "CBS") %>%
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
plot_data_obs <- cbs_data %>%
  filter(date >= fitting_start_date, !is.na(quintile_clean)) %>%
  dplyr::select(date, seroprev_est, quintile_clean)

# ---- 2. LOAD PARSIMONIOUS MODEL ------------------------------------------

results <- readRDS("parsimonious_model_estimate_rho.rds")
fit <- results$fit
coefs <- coef(fit)
vcov_mat <- vcov(fit)

model_start_date <- results$model_start_date
omicron_date <- results$omicron_date

# ---- 3. NP MODEL SIMULATION (parsimonious specification) -----------------
# Parsimonious: each quintile has its own lambda_pre, but a SHARED multiplier m
# Lambda_omi = lambda_pre * m

simulate_np_parsim <- function(params, target_dates, model_start_date, omicron_date, quintile_idx) {
  lambda_pre <- plogis(params[paste0("logit_lambda_Q", quintile_idx, "_pre")])
  multiplier <- exp(params["log_multiplier"])
  lambda_omi <- lambda_pre * multiplier
  rho <- plogis(params["logit_rho"]) * 0.01

  all_dates_seq <- seq(model_start_date, max(target_dates), by = "7 days")
  num_steps <- length(all_dates_seq)
  N_t <- numeric(num_steps); P_t <- numeric(num_steps)
  N_t[1] <- 1; P_t[1] <- 0

  for (t in 2:num_steps) {
    lam <- if (all_dates_seq[t] < omicron_date) lambda_pre else lambda_omi
    N_t[t] <- N_t[t-1] * exp(-lam) + P_t[t-1] * (1 - exp(-rho))
    P_t[t] <- P_t[t-1] * exp(-rho) + N_t[t-1] * (1 - exp(-lam))
  }

  out <- numeric(length(target_dates))
  for (i in seq_along(target_dates)) {
    idx <- which.min(abs(all_dates_seq - target_dates[i]))
    out[i] <- P_t[idx]
  }
  out * 100
}

# ---- 4. PARAMETRIC BOOTSTRAP ---------------------------------------------

set.seed(2026)
n_boot <- 1000
cat(sprintf("Running %d-draw bootstrap for CI bands...\n", n_boot))

param_draws <- mvrnorm(n_boot, mu = coefs, Sigma = vcov_mat)
colnames(param_draws) <- names(coefs)

plot_dates <- seq(as.Date("2021-04-01"), as.Date("2023-05-01"), by = "7 days")

quintiles <- c("Q1", "Q2", "Q3", "Q4", "Q5")
prediction_data <- data.frame()

for (i in seq_along(quintiles)) {
  point_pred <- simulate_np_parsim(coefs, plot_dates, model_start_date, omicron_date, i)
  boot_preds <- t(apply(param_draws, 1, function(p) {
    simulate_np_parsim(p, plot_dates, model_start_date, omicron_date, i)
  }))
  lower <- apply(boot_preds, 2, quantile, probs = 0.025, na.rm = TRUE)
  upper <- apply(boot_preds, 2, quantile, probs = 0.975, na.rm = TRUE)

  prediction_data <- rbind(prediction_data, data.frame(
    date = plot_dates,
    predicted = point_pred,
    lower = lower,
    upper = upper,
    quintile_clean = quintiles[i]
  ))
  cat(sprintf("  %s: done\n", quintiles[i]))
}

plot_data_obs$Quintile <- factor(plot_data_obs$quintile_clean,
                                  levels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                                  labels = c("Q1 (Least Deprived)", "Q2", "Q3", "Q4", "Q5 (Most Deprived)"))
prediction_data$Quintile <- factor(prediction_data$quintile_clean,
                                    levels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                                    labels = c("Q1 (Least Deprived)", "Q2", "Q3", "Q4", "Q5 (Most Deprived)"))

quintile_colors <- c(
  "Q1 (Least Deprived)" = "#003f5c",
  "Q2"                  = "#58508d",
  "Q3"                  = "#bc5090",
  "Q4"                  = "#ff6361",
  "Q5 (Most Deprived)"  = "#ffa600"
)

p <- ggplot() +
  geom_ribbon(data = prediction_data,
              aes(x = date, ymin = lower, ymax = upper, fill = Quintile),
              alpha = 0.18, color = NA) +
  geom_line(data = prediction_data,
            aes(x = date, y = predicted, color = Quintile),
            linewidth = 1.3, alpha = 0.95) +
  geom_point(data = plot_data_obs,
             aes(x = date, y = seroprev_est, color = Quintile),
             size = 2.8, alpha = 0.85) +
  geom_vline(xintercept = as.numeric(omicron_date),
             linetype = "dashed", color = "gray30", linewidth = 0.7) +
  annotate("text", x = omicron_date, y = 87,
           label = "Omicron Emergence\n(Jan 1, 2022)",
           hjust = -0.05, size = 3.5, color = "gray20") +
  scale_color_manual(values = quintile_colors) +
  scale_fill_manual(values = quintile_colors, guide = "none") +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y",
               limits = c(as.Date("2021-04-01"), as.Date("2023-05-01"))) +
  scale_y_continuous(limits = c(0, 95), breaks = seq(0, 90, 15)) +
  labs(x = "Date", y = "Seroprevalence (%)",
       color = "Material Deprivation") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
  ) +
  guides(color = guide_legend(nrow = 1))

ggsave("FigureS1_Parsimonious_Fit.png", p,
       width = 11, height = 7, dpi = 300, bg = "white")
ggsave("FigureS1_Parsimonious_Fit.pdf", p,
       width = 11, height = 7, bg = "white")

cat("\n✓ Figure S1 saved with CI bands\n\n")
print(p)
