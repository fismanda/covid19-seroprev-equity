################################################################################
# Figure 3: Omicron Multipliers - label positioning fix
#
# Change from original: multiplier text labels positioned above the Upper CI
# whisker instead of just above the bar top, eliminating overlap.
################################################################################

library(ggplot2)
library(dplyr)

cat("\n=== CREATING FIGURE 3: OMICRON MULTIPLIERS (LABELS REPOSITIONED) ===\n\n")

# Load multiplier data with CIs (same source as original script)
mult_data <- readRDS("omicron_multipliers_with_ci.rds")

# Create labels for x-axis with quintile descriptors
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

# Vertical offset for labels — proportional to plotting range
label_offset <- 1.5  # in y-axis units; sits cleanly above the upper whisker

p <- ggplot(plot_data, aes(x = Quintile_label, y = Multiplier)) +
  geom_col(aes(fill = Multiplier), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.3, linewidth = 0.6) +
  # CHANGED: label positioned at Upper + offset, not just above the bar
  geom_text(aes(y = Upper + label_offset,
                label = sprintf("%.1fx", Multiplier)),
            size = 3.8, fontface = "bold") +
  scale_fill_gradient(low = "#A23B72", high = "#2E86AB", guide = "none") +
  labs(x = "Material Deprivation Quintile",
       y = "Omicron Multiplier\n(Fold-Change in Force of Infection)") +
  scale_y_continuous(limits = c(0, 55),  # extended slightly to accommodate labels
                     breaks = seq(0, 50, 10),
                     expand = c(0, 0)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave("Figure3_Omicron_Multipliers.png", p,
       width = 8, height = 6, dpi = 300, bg = "white")
ggsave("Figure3_Omicron_Multipliers.pdf", p,
       width = 8, height = 6, bg = "white")

cat("✓ Figure 3 saved (labels above upper CI whiskers)\n\n")


