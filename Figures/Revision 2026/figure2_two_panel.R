################################################################################
# Figure 2: Two-panel display, axes flipped (quintile on y-axis throughout)
#   Panel A: Forest plot of IRRs by quintile and period (existing, restyled)
#   Panel B: Absolute weekly force of infection by quintile and period
################################################################################

library(ggplot2)
library(dplyr)
library(patchwork)

cat("\n=== GENERATING FIGURE 2: TWO-PANEL IRR + ABSOLUTE FOI ===\n\n")

# ---- 1. LOAD DATA --------------------------------------------------------

irr_data <- readRDS("irr_estimates_with_ci.rds")

# Extract saturated fit so we can compute absolute lambda CIs via delta method
results <- readRDS("saturated_model_estimate_rho.rds")
fit <- results$fit
coefs <- coef(fit)
vcov_mat <- vcov(fit)

# ---- 2. ABSOLUTE LAMBDA WITH 95% CIs --------------------------------------
# lambda = plogis(x), where x is the logit-transformed parameter
# 95% CI on lambda by transforming the CI endpoints on x

z_crit <- qnorm(0.975)
lambda_rows <- list()

for (q in 1:5) {
  for (period in c("pre", "omi")) {
    par_name <- paste0("logit_lambda_Q", q, "_", period)
    x_hat <- coefs[par_name]
    se_x  <- sqrt(vcov_mat[par_name, par_name])
    lambda_est <- plogis(x_hat)
    lambda_lo  <- plogis(x_hat - z_crit * se_x)
    lambda_hi  <- plogis(x_hat + z_crit * se_x)
    lambda_rows[[length(lambda_rows) + 1]] <- data.frame(
      Quintile = paste0("Q", q),
      Period = ifelse(period == "pre", "Pre-Omicron", "Omicron"),
      lambda = lambda_est,
      lambda_lo = lambda_lo,
      lambda_hi = lambda_hi
    )
  }
}
lambda_df <- do.call(rbind, lambda_rows)
lambda_df$Period <- factor(lambda_df$Period, levels = c("Pre-Omicron", "Omicron"))

# ---- 3. PREPARE IRR DATA FOR PANEL A --------------------------------------

irr_plot <- irr_data %>%
  mutate(Period = factor(Period, levels = c("Pre-Omicron", "Omicron")))

# Add Q1 reference rows (IRR = 1.00 exactly, no CI)
q1_rows <- data.frame(
  Quintile = c("Q1", "Q1"),
  Period = factor(c("Pre-Omicron", "Omicron"), levels = c("Pre-Omicron", "Omicron")),
  IRR = c(1, 1),
  Lower = c(NA, NA),
  Upper = c(NA, NA)
)
irr_plot <- bind_rows(irr_plot, q1_rows)

period_colors <- c("Pre-Omicron" = "#2E86AB", "Omicron" = "#A23B72")

# Y-axis: Q1 top, Q5 bottom (forest-plot convention)
y_levels <- c("Q5", "Q4", "Q3", "Q2", "Q1")
y_labels <- c("Q5 (most deprived)", "Q4", "Q3", "Q2", "Q1 (least deprived)")

# ---- 4. PANEL A: IRR FOREST PLOT ------------------------------------------

panel_A <- ggplot(irr_plot, aes(x = IRR, y = Quintile, color = Period)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.2, linewidth = 0.8,
                 position = position_dodge(width = 0.5),
                 na.rm = TRUE) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = period_colors, name = "Period") +
  scale_x_log10(breaks = c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8)) +
  scale_y_discrete(limits = y_levels, labels = y_labels) +
  labs(x = "Incidence Rate Ratio (vs Q1)",
       y = NULL,
       tag = "A") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.tag = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11)
  ) +
  coord_cartesian(xlim = c(0.95, 1.85))

# ---- 5. PANEL B: ABSOLUTE WEEKLY FOI --------------------------------------

panel_B <- ggplot(lambda_df,
                  aes(x = lambda, y = Quintile, color = Period)) +
  geom_errorbarh(aes(xmin = lambda_lo, xmax = lambda_hi),
                 height = 0.2, linewidth = 0.8,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = period_colors, name = "Period") +
  scale_x_log10() +
  scale_y_discrete(limits = y_levels, labels = y_labels) +
  labs(x = "Weekly Force of Infection (λ, log scale)",
       y = NULL,
       tag = "B") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.tag = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11)
  )

# ---- 6. COMPOSE AND SAVE --------------------------------------------------

fig2 <- panel_A / panel_B + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("Figure2_Two_Panel.png", fig2,
       width = 9, height = 9, dpi = 300, bg = "white")
ggsave("Figure2_Two_Panel.pdf", fig2,
       width = 9, height = 9, bg = "white")

cat("\n✓ Figure 2 saved:\n")
cat("  - Figure2_Two_Panel.png (300 dpi)\n")
cat("  - Figure2_Two_Panel.pdf (vector)\n\n")

lambda_df %>%
  filter(Period == "Omicron") %>%
  mutate(ci_width = lambda_hi - lambda_lo) %>%
  dplyr::select(Quintile, lambda, lambda_lo, lambda_hi, ci_width)
