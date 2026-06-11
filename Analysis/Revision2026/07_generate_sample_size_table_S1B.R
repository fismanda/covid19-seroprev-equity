# 06_generate_sample_size_table_S1B.R
#
# Generates Supplementary Table S1B (sample sizes by material deprivation
# quintile and observation period) from the deposited analytic dataset.
#
# Run from the repository root.

library(dplyr)
library(lubridate)

# Load the deposited analytic dataset
ses_data <- read.csv("Data/cbs_seroprev_by_quintile.csv")
ses_data$date <- as.Date(ses_data$date)

ss_table <- ses_data %>%
  mutate(period = case_when(
    date < as.Date("2022-01-01") ~ "Pre-Omicron (Apr-Dec 2021)",
    date < as.Date("2023-01-01") ~ "Omicron 2022",
    TRUE                          ~ "Omicron 2023"
  )) %>%
  group_by(period, quintile) %>%
  summarise(
    n_obs          = n(),
    n_donors_total = sum(n_donors, na.rm = TRUE),
    n_min          = min(n_donors, na.rm = TRUE),
    n_max          = max(n_donors, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(period, quintile)

print(ss_table)
write.csv(ss_table, "Analysis/revision_2026/Table_S1B_sample_sizes.csv",
          row.names = FALSE)