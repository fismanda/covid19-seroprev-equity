################################################################################
# Figure 1: Two-panel display
#   Panel A: Observed vs predicted seroprevalence with parametric bootstrap CI bands
#   Panel B: Q5 - Q1 absolute difference in observed seroprevalence over time
################################################################################

library(ggplot2)
library(dplyr)
library(lubridate)
library(patchwork)
library(MASS)  # for mvrnorm

cat("\n=== GENERATING FIGURE 1: TWO-PANEL MODEL FIT WITH CI BANDS ===\n\n")

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

# ---- 2. LOAD MODEL --------------------------------------------------------

results <- readRDS("saturated_model_estimate_rho.rds")
fit <- results$fit
coefs <- coef(fit)
vcov_mat <- vcov(fit)

model_start_date <- results$model_start_date
omicron_date <- results$omicron_date

# ---- 3. NP MODEL SIMULATION FUNCTION -------------------------------------
# Takes logit-transformed parameters, returns seropositive proportion at each date

simulate_np <- function(params, target_dates, model_start_date, omicron_date, quintile_idx) {
  # Extract logit-transformed lambdas for this quintile and ρ
  lambda_pre <- plogis(params[paste0("logit_lambda_Q", quintile_idx, "_pre")])
  lambda_omi <- plogis(params[paste0("logit_lambda_Q", quintile_idx, "_omi")])
  # Find the rho parameter name (saturated model uses "logit_rho_scaled" usually,
  # parsimonious uses "logit_rho")
  rho_name <- intersect(c("logit_rho_scaled", "logit_rho"), names(params))[1]
  rho <- plogis(params[rho_name]) * 0.01

  all_dates_seq <- seq(model_start_date, max(target_dates), by = "7 days")
  num_steps <- length(all_dates_seq)
  N_t <- numeric(num_steps); P_t <- numeric(num_steps)
  N_t[1] <- 1; P_t[1] <- 0

  for (t in 2:num_steps) {
    lam <- if (all_dates_seq[t] < omicron_date) lambda_pre else lambda_omi
    N_t[t] <- N_t[t-1] * exp(-lam) + P_t[t-1] * (1 - exp(-rho))
    P_t[t] <- P_t[t-1] * exp(-rho) + N_t[t-1] * (1 - exp(-lam))
  }

  # Match target dates
  out <- numeric(length(target_dates))
  for (i in seq_along(target_dates)) {
    idx <- which.min(abs(all_dates_seq - target_dates[i]))
    out[i] <- P_t[idx]
  }
  out * 100  # percentage
}

# ---- 4. PARAMETRIC BOOTSTRAP FOR CI BANDS --------------------------------

set.seed(2026)
n_boot <- 1000
cat(sprintf("Running %d-draw parametric bootstrap for CI bands...\n", n_boot))

# Sample parameters from MVN(theta_hat, Sigma)
param_draws <- mvrnorm(n_boot, mu = coefs, Sigma = vcov_mat)
colnames(param_draws) <- names(coefs)

# Common date sequence for smooth lines and bands
plot_dates <- seq(as.Date("2021-04-01"), as.Date("2023-05-01"), by = "7 days")

quintiles <- c("Q1", "Q2", "Q3", "Q4", "Q5")
prediction_data <- data.frame()
boot_preds_list <- vector("list", length(quintiles))  # store Q1/Q5 matrices for Panel B

for (i in seq_along(quintiles)) {
  # Point estimate
  point_pred <- simulate_np(coefs, plot_dates, model_start_date, omicron_date, i)

  # Bootstrap draws (matrix: n_boot rows x length(plot_dates) cols)
  boot_preds <- t(apply(param_draws, 1, function(p) {
    simulate_np(p, plot_dates, model_start_date, omicron_date, i)
  }))

  # Save Q1 and Q5 matrices for cross-quintile difference computation
  if (quintiles[i] %in% c("Q1", "Q5")) {
    boot_preds_list[[i]] <- boot_preds
  }

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

# Compute Q5 - Q1 predicted-difference CI from saved bootstrap matrices
boot_q1 <- boot_preds_list[[1]]   # Q1
boot_q5 <- boot_preds_list[[5]]   # Q5
boot_diff <- boot_q5 - boot_q1    # element-wise: n_boot x length(plot_dates)
diff_lower <- apply(boot_diff, 2, quantile, probs = 0.025, na.rm = TRUE)
diff_upper <- apply(boot_diff, 2, quantile, probs = 0.975, na.rm = TRUE)

# Apply factor for plotting
plot_data_obs$Quintile <- factor(plot_data_obs$quintile_clean,
                                  levels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                                  labels = c("Q1 (Least Deprived)", "Q2", "Q3", "Q4", "Q5 (Most Deprived)"))
prediction_data$Quintile <- factor(prediction_data$quintile_clean,
                                    levels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                                    labels = c("Q1 (Least Deprived)", "Q2", "Q3", "Q4", "Q5 (Most Deprived)"))

# Better-discrimination palette (viridis-like)
quintile_colors <- c(
  "Q1 (Least Deprived)" = "#003f5c",
  "Q2"                  = "#58508d",
  "Q3"                  = "#bc5090",
  "Q4"                  = "#ff6361",
  "Q5 (Most Deprived)"  = "#ffa600"
)

# ---- 5. PANEL A: TRAJECTORIES WITH CI BANDS ------------------------------

panel_A <- ggplot() +
  geom_ribbon(data = prediction_data,
              aes(x = date, ymin = lower, ymax = upper, fill = Quintile),
              alpha = 0.18, color = NA) +
  geom_line(data = prediction_data,
            aes(x = date, y = predicted, color = Quintile),
            linewidth = 1.4, alpha = 0.95) +
  geom_point(data = plot_data_obs,
             aes(x = date, y = seroprev_est, color = Quintile),
             size = 3, alpha = 0.85) +
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
  labs(x = NULL, y = "Seroprevalence (%)",
       color = "Material Deprivation",
       tag = "A") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    plot.tag = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
  ) +
  guides(color = guide_legend(nrow = 1))

# ---- 6. PANEL B: Q5 - Q1 ABSOLUTE DIFFERENCE -----------------------------

# Compute observed Q5-Q1 difference at dates where both are observed
# Aggregate any duplicate (date, quintile) rows by simple mean before pivoting
diff_obs <- plot_data_obs %>%
  dplyr::select(date, quintile_clean, seroprev_est) %>%
  filter(quintile_clean %in% c("Q1", "Q5")) %>%
  group_by(date, quintile_clean) %>%
  summarise(seroprev_est = mean(seroprev_est, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = quintile_clean, values_from = seroprev_est) %>%
  filter(!is.na(Q1) & !is.na(Q5)) %>%
  mutate(diff = Q5 - Q1)

# Also compute predicted Q5-Q1 difference from the model, with bootstrap CIs
pred_wide <- prediction_data %>%
  dplyr::select(date, quintile_clean, predicted) %>%
  tidyr::pivot_wider(names_from = quintile_clean, values_from = predicted) %>%
  mutate(diff_pred = Q5 - Q1,
         diff_lo = diff_lower,
         diff_hi = diff_upper)

panel_B <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray60") +
  # Bootstrap CI bounds as dashed lines
  geom_line(data = pred_wide,
            aes(x = date, y = diff_lo),
            color = "#003f5c", linewidth = 0.6, linetype = "dashed") +
  geom_line(data = pred_wide,
            aes(x = date, y = diff_hi),
            color = "#003f5c", linewidth = 0.6, linetype = "dashed") +
  # Model-predicted difference (point estimate)
  geom_line(data = pred_wide,
            aes(x = date, y = diff_pred),
            color = "#003f5c", linewidth = 1.2) +
  # Observed Q5 - Q1 difference
  geom_point(data = diff_obs,
             aes(x = date, y = diff),
             color = "#bc5090", size = 2.5, alpha = 0.85) +
  geom_vline(xintercept = as.numeric(omicron_date),
             linetype = "dashed", color = "gray30", linewidth = 0.7) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y",
               limits = c(as.Date("2021-04-01"), as.Date("2023-05-01"))) +
  labs(x = "Date",
       y = "Q5 - Q1 seroprevalence difference (percentage points)",
       tag = "B") +
  theme_minimal(base_size = 12) +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
  )

# ---- 7. COMPOSE AND SAVE -------------------------------------------------

fig1 <- panel_A / panel_B + plot_layout(heights = c(2, 1))

ggsave("Figure1_Two_Panel.png", fig1,
       width = 11, height = 11, dpi = 300, bg = "white")
ggsave("Figure1_Two_Panel.pdf", fig1,
       width = 11, height = 11, bg = "white")

cat("\n✓ Figure 1 saved:\n")
cat("  - Figure1_Two_Panel.png (300 dpi)\n")
cat("  - Figure1_Two_Panel.pdf (vector)\n\n")
