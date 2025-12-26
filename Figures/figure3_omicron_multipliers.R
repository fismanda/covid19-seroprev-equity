################################################################################
# Figure 3: Omicron Multipliers by Quintile
################################################################################
#
# PURPOSE:
#   Create bar chart showing fold-change in force of infection from pre-Omicron
#   to Omicron by material deprivation quintile. Visualizes "differential
#   amplification" - affluent populations experienced larger proportional increases.
#
# INPUTS:
#   - omicron_multipliers_with_ci.rds
#
# OUTPUTS:
#   - Figure3_Omicron_Multipliers.png (300 dpi)
#   - Figure3_Omicron_Multipliers.pdf (vector)
#
# AUTHOR: Hassan, Nassrallah, Fisman
# DATE: December 2025
################################################################################

library(ggplot2)
library(dplyr)

cat("\n===================================================")
cat("\n=== CREATING FIGURE 3: OMICRON MULTIPLIERS ===")
cat("\n===================================================\n\n")

# Load multiplier data
mult_data <- readRDS("omicron_multipliers_with_ci.rds")

cat("✓ Data loaded\n")
cat(sprintf("  Quintiles: %d\n\n", nrow(mult_data)))

# Prepare data for plotting
plot_data <- mult_data %>%
  mutate(
    Quintile_label = case_when(
      Quintile == "Q1" ~ "Q1\n(least deprived)",
      Quintile == "Q5" ~ "Q5\n(most deprived)",
      TRUE ~ Quintile
    ),
    Quintile_label = factor(Quintile_label, 
                           levels = c("Q1\n(least deprived)", "Q2", "Q3", "Q4", "Q5\n(most deprived)"))
  )

cat("Creating bar chart...\n")

# Create bar chart with error bars
p <- ggplot(plot_data, aes(x = Quintile_label, y = Multiplier)) +
  
  # Bars with gradient fill
  geom_col(aes(fill = Multiplier), width = 0.7, color = "black", linewidth = 0.3) +
  
  # Error bars
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.3, linewidth = 0.6) +
  
  # Add value labels on bars
  geom_text(aes(label = sprintf("%.1fx", Multiplier)),
            vjust = -0.5, size = 3.5, fontface = "bold") +
  
  # Color gradient (blue to red, high to low)
  scale_fill_gradient(
    low = "#A23B72",    # Red/pink for lower multipliers (Q5)
    high = "#2E86AB",   # Blue for higher multipliers (Q1)
    guide = "none"      # Hide legend
  ) +
  
  # Labels
  labs(
    x = "Material Deprivation Quintile",
    y = "Omicron Multiplier\n(Fold-Change in Force of Infection)",
    title = "Differential Amplification During Omicron",
    subtitle = "Affluent populations experienced larger proportional increases in transmission"
  ) +
  
  # Y-axis
  scale_y_continuous(
    limits = c(0, 52),
    breaks = seq(0, 50, 10),
    expand = c(0, 0)
  ) +
  
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    plot.margin = margin(10, 10, 10, 10)
  )

# Save as PNG (300 dpi)
ggsave("Figure3_Omicron_Multipliers.png", 
       plot = p, 
       width = 8, 
       height = 6, 
       dpi = 300,
       bg = "white")

# Save as PDF (vector)
ggsave("Figure3_Omicron_Multipliers.pdf", 
       plot = p, 
       width = 8, 
       height = 6)

cat("✓ Figure saved:\n")
cat("  - Figure3_Omicron_Multipliers.png (300 dpi)\n")
cat("  - Figure3_Omicron_Multipliers.pdf (vector)\n\n")

# Print summary
cat("=== DIFFERENTIAL AMPLIFICATION SUMMARY ===\n\n")

q1_mult <- mult_data %>% filter(Quintile == "Q1")
q5_mult <- mult_data %>% filter(Quintile == "Q5")

cat(sprintf("Q1 (least deprived):  %.2fx (95%% CI: %.2f - %.2f)\n",
            q1_mult$Multiplier, q1_mult$Lower, q1_mult$Upper))
cat(sprintf("Q5 (most deprived):   %.2fx (95%% CI: %.2f - %.2f)\n\n",
            q5_mult$Multiplier, q5_mult$Lower, q5_mult$Upper))

ratio <- q1_mult$Multiplier / q5_mult$Multiplier
pct_larger <- (ratio - 1) * 100

cat(sprintf("Q1/Q5 ratio: %.2f\n", ratio))
cat(sprintf("Q1 experienced %.0f%% larger proportional increase than Q5.\n\n", pct_larger))

cat("This differential amplification produced gradient compression:\n")
cat("affluent populations 'caught up' through policy-driven increases,\n")
cat("resulting in apparent convergence despite persistent disparities.\n\n")

