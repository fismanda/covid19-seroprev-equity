################################################################################
# Figure 2: Forest Plot of Incidence Rate Ratios
################################################################################
#
# PURPOSE:
#   Create forest plot comparing IRRs (vs Q1) for pre-Omicron and Omicron periods.
#   Visualizes gradient compression: strong gradient pre-Omicron (Q5=1.71) 
#   compresses to weak gradient during Omicron (Q5=1.12).
#
# INPUTS:
#   - irr_estimates_with_ci.rds
#
# OUTPUTS:
#   - Figure2_IRR_Forest_Plot.png (300 dpi)
#   - Figure2_IRR_Forest_Plot.pdf (vector)
#
# AUTHOR: Hassan, Nassrallah, Fisman
# DATE: December 2025
################################################################################

library(ggplot2)
library(dplyr)

cat("\n==============================================")
cat("\n=== CREATING FIGURE 2: IRR FOREST PLOT ===")
cat("\n==============================================\n\n")

# Load IRR data
irr_data <- readRDS("irr_estimates_with_ci.rds")

cat("✓ Data loaded\n")
cat(sprintf("  Total data points: %d\n\n", nrow(irr_data)))

# Prepare data for plotting
plot_data <- irr_data %>%
  filter(Quintile != "Q1") %>%  # Exclude reference group
  mutate(
    Quintile_num = as.numeric(gsub("Q", "", Quintile)),
    Period = factor(Period, levels = c("Pre-Omicron", "Omicron"))
  )

cat("Creating forest plot...\n")

# Create forest plot
p <- ggplot(plot_data, aes(x = IRR, y = Quintile, color = Period)) +
  
  # Reference line at IRR = 1
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  
  # Confidence intervals
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.2, linewidth = 0.8,
                 position = position_dodge(width = 0.5)) +
  
  # Point estimates
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  
  # Color scheme
  scale_color_manual(
    values = c("Pre-Omicron" = "#2E86AB", "Omicron" = "#A23B72"),
    name = "Period"
  ) +
  
  # Axis labels
  labs(
    x = "Incidence Rate Ratio (vs Q1, least deprived)",
    y = "Material Deprivation Quintile",
    title = "Socioeconomic Gradient in SARS-CoV-2 Force of Infection",
    subtitle = "Gradient compression from pre-Omicron to Omicron period"
  ) +
  
  # Log scale for x-axis
  scale_x_log10(
    breaks = c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8),
    labels = c("1.0", "1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8")
  ) +
  
  # Y-axis: reverse order so Q2 at top, Q5 at bottom
  scale_y_discrete(
    limits = rev(c("Q2", "Q3", "Q4", "Q5")),
    labels = rev(c("Q2", "Q3", "Q4", "Q5 (most deprived)"))
  ) +
  
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  
  # Coordinate system
  coord_cartesian(xlim = c(0.95, 1.85))

# Save as PNG (300 dpi)
ggsave("Figure2_IRR_Forest_Plot.png", 
       plot = p, 
       width = 8, 
       height = 5, 
       dpi = 300,
       bg = "white")

# Save as PDF (vector)
ggsave("Figure2_IRR_Forest_Plot.pdf", 
       plot = p, 
       width = 8, 
       height = 5)

cat("✓ Figure saved:\n")
cat("  - Figure2_IRR_Forest_Plot.png (300 dpi)\n")
cat("  - Figure2_IRR_Forest_Plot.pdf (vector)\n\n")

# Print summary statistics
cat("=== GRADIENT COMPRESSION SUMMARY ===\n\n")

pre_q5 <- irr_data %>% filter(Period == "Pre-Omicron", Quintile == "Q5")
omi_q5 <- irr_data %>% filter(Period == "Omicron", Quintile == "Q5")

cat(sprintf("Pre-Omicron Q5 vs Q1: %.3f (95%% CI: %.3f - %.3f)\n",
            pre_q5$IRR, pre_q5$Lower, pre_q5$Upper))
cat(sprintf("Omicron Q5 vs Q1:     %.3f (95%% CI: %.3f - %.3f)\n\n",
            omi_q5$IRR, omi_q5$Lower, omi_q5$Upper))

compression <- (pre_q5$IRR - omi_q5$IRR) / pre_q5$IRR * 100
cat(sprintf("Gradient compression: %.1f%%\n", compression))
cat(sprintf("Absolute reduction:   %.3f\n\n", pre_q5$IRR - omi_q5$IRR))

cat("The socioeconomic gradient in force of infection compressed\n")
cat("substantially during Omicron, as affluent populations experienced\n")
cat("larger proportional increases in transmission.\n\n")

