library(readxl)
library(tidyverse)

# Import our main data
data_aa <- read_excel("for_jedediah_AAs.xlsx")
data_ea <- read_excel("for_jedediah_EAs.xlsx")

# Tidy up the AA Data
data_aa_clean <- data_aa %>%
  
  # Make column names lowercase and rename for readability
  rename_with(tolower) %>%
  rename(
    incident_stroke = stroke22,
    baseline_ecb = baseline_impair, 
    baseline_stroke = stroke_sr) %>%
  
  # Rearrange columns so that demographic data is first
  select(
    id_num, age, gender, 
    incident_ecb, incident_stroke,   # Incidents can be covariates in the model.
    baseline_ecb, baseline_stroke,   # Baselines can be used to filter participants.
    everything(), -race) %>%         # Remove unused, non-numeric column.
  
  # Drop columns with missing incident_ecb values
  drop_na(incident_ecb) %>%
  
  # Drop columns where baseline_ecb = 1
  filter(incident_ecb != 1) %>%
  
  # Recode non-numeric demographics
  mutate(
    gender = case_match(
      gender,                        # SELF-REPORTED
      "F" ~ 0,
      "M" ~ 1),
    baseline_stroke = case_match(
      baseline_stroke,
      "N" ~ 0,
      "Y" ~ 1)) %>%
  
  # Recode APOE4 carrier status as binary
  select(-rs429358, -rs7412, -gt) %>%           # Remove excess non-numeric columns
  mutate(
    apoe_4 = case_when(                         # Create new binary column
      is.na(apoe_gwas) ~ NA_integer_,           # Preserve NAs
      grepl("4", as.character(apoe_gwas)) ~ 1L,
      TRUE ~ 0L),
    across(everything(), as.numeric)) %>%      # Make everything numeric
  select(-apoe_gwas)                           # Remove leftover column

# Save for later
write.csv(data_aa_clean, "aa_clean.csv", row.names = FALSE)

# Tidy up the EA Data
data_ea_clean <- data_ea %>%
  
  # Make column names lowercase and rename for readability
  rename_with(tolower) %>%
  rename(
    incident_stroke = stroke22,
    baseline_ecb = baseline_impair, 
    baseline_stroke = stroke_sr) %>%
  
  # Rearrange columns so that demographic data is first
  select(
    id_num, age, gender, 
    incident_ecb, incident_stroke,         # Incidents can be covariates in the model.
    baseline_ecb, baseline_stroke,         # Baselines can be used to filter participants.
    everything(), -race, -population) %>%  # Remove unused, non-numeric columns.
  
  # Drop columns with missing incident_ecb values
  drop_na(incident_ecb) %>%
  
  # Drop columns where baseline_ecb = 1
  filter(incident_ecb != 1) %>%
  
  # Recode non-numeric demographics
  mutate(
    gender = case_match(
      gender,                        # SELF-REPORTED
      "F" ~ 0,
      "M" ~ 1),
    baseline_stroke = case_match(
      baseline_stroke,
      "N" ~ 0,
      "Y" ~ 1)) %>%
  
  # Recode APOE4 carrier status as binary
  select(-rs429358, -rs7412) %>%               # Remove excess non-numeric columns
  mutate(
    apoe_4 = case_when(                         # Create new binary column
      is.na(apoe_gwas) ~ NA_integer_,           # Preserve NAs
      grepl("4", as.character(apoe_gwas)) ~ 1L,
      TRUE ~ 0L),
    across(everything(), as.numeric)) %>%      # Make everything numeric
  select(-apoe_gwas)                           # Remove leftover column

# Save for later
write.csv(data_ea_clean, "ea_clean.csv", row.names = FALSE)

# Calculate summary statistics for AA Data
summary_aa <- data_aa_clean %>%
  mutate(group = ifelse(incident_ecb == 1, "Cases", "Controls")) %>%
  group_by(group) %>%
  reframe(
    Count = n(),
    Percent = n() / nrow(data_aa_clean) * 100,
    Median_Age = median(age, na.rm = TRUE),
    Age_Q1 = quantile(age, 0.25, na.rm = TRUE),
    Age_Q3 = quantile(age, 0.75, na.rm = TRUE),
    Male_n = sum(gender == 1, na.rm = TRUE),
    Male_percent = sum(gender == 1, na.rm = TRUE) / n() * 100,
    Female_n = sum(gender == 0, na.rm = TRUE),
    Female_percent = sum(gender == 0, na.rm = TRUE) / n() * 100
  ) %>%
  mutate(
    Count_fmt = sprintf("%d (%.1f%%)", Count, Percent),
    Age_fmt = sprintf("%d (%.0f, %.0f)", Median_Age, Age_Q1, Age_Q3),
    Male_fmt = sprintf("%d (%.1f%%)", Male_n, Male_percent),
    Female_fmt = sprintf("%d (%.1f%%)", Female_n, Female_percent)
  ) %>%
  select(group, Count_fmt, Age_fmt, Male_fmt, Female_fmt)

# Format summary statistics as table
table_aa <- summary_aa %>%
  pivot_longer(-group, names_to = "Variable", values_to = "Value") %>%
  pivot_wider(names_from = group, values_from = Value) %>%
  mutate(
    Variable = recode(Variable,
                      Count_fmt = "Count (n%)",
                      Age_fmt = "Med. Age (IQR)",
                      Male_fmt = "Male (n%)",
                      Female_fmt = "Female (n%)")
  )

# Print table
print(table_aa)

# Calculate summary statistics for AE Data
summary_ea <- data_ea_clean %>%
  mutate(group = ifelse(incident_ecb == 1, "Cases", "Controls")) %>%
  group_by(group) %>%
  reframe(
    Count = n(),
    Percent = n() / nrow(data_ea_clean) * 100,
    Median_Age = median(age, na.rm = TRUE),
    Age_Q1 = quantile(age, 0.25, na.rm = TRUE),
    Age_Q3 = quantile(age, 0.75, na.rm = TRUE),
    Male_n = sum(gender == 1, na.rm = TRUE),
    Male_percent = sum(gender == 1, na.rm = TRUE) / n() * 100,
    Female_n = sum(gender == 0, na.rm = TRUE),
    Female_percent = sum(gender == 0, na.rm = TRUE) / n() * 100
  ) %>%
  mutate(
    Count_fmt = sprintf("%d (%.1f%%)", Count, Percent),
    Age_fmt = sprintf("%d (%.0f, %.0f)", Median_Age, Age_Q1, Age_Q3),
    Male_fmt = sprintf("%d (%.1f%%)", Male_n, Male_percent),
    Female_fmt = sprintf("%d (%.1f%%)", Female_n, Female_percent)
  ) %>%
  select(group, Count_fmt, Age_fmt, Male_fmt, Female_fmt)

# Format summary statistics as table
table_ea <- summary_ea %>%
  pivot_longer(-group, names_to = "Variable", values_to = "Value") %>%
  pivot_wider(names_from = group, values_from = Value) %>%
  mutate(
    Variable = recode(Variable,
                      Count_fmt = "Count (n%)",
                      Age_fmt = "Med. Age (IQR)",
                      Male_fmt = "Male (n%)",
                      Female_fmt = "Female (n%)")
  )

# Print table
print(table_ea)

# Create a combined dataset
data_combined <- bind_rows(data_aa_clean, data_ea_clean)

# Save for later
write.csv(data_combined, "combo_clean.csv", row.names = FALSE)

# Calculate summary statistics for AA Data
summary_combo <- data_combined %>%
  mutate(group = ifelse(incident_ecb == 1, "Cases", "Controls")) %>%
  group_by(group) %>%
  reframe(
    Count = n(),
    Percent = n() / nrow(data_combined) * 100,
    Median_Age = median(age, na.rm = TRUE),
    Age_Q1 = quantile(age, 0.25, na.rm = TRUE),
    Age_Q3 = quantile(age, 0.75, na.rm = TRUE),
    Male_n = sum(gender == 1, na.rm = TRUE),
    Male_percent = sum(gender == 1, na.rm = TRUE) / n() * 100,
    Female_n = sum(gender == 0, na.rm = TRUE),
    Female_percent = sum(gender == 0, na.rm = TRUE) / n() * 100
  ) %>%
  mutate(
    Count_fmt = sprintf("%d (%.1f%%)", Count, Percent),
    Age_fmt = sprintf("%d (%.0f, %.0f)", Median_Age, Age_Q1, Age_Q3),
    Male_fmt = sprintf("%d (%.1f%%)", Male_n, Male_percent),
    Female_fmt = sprintf("%d (%.1f%%)", Female_n, Female_percent)
  ) %>%
  select(group, Count_fmt, Age_fmt, Male_fmt, Female_fmt)

# Format summary statistics as table
table_combo <- summary_combo %>%
  pivot_longer(-group, names_to = "Variable", values_to = "Value") %>%
  pivot_wider(names_from = group, values_from = Value) %>%
  mutate(
    Variable = recode(Variable,
                      Count_fmt = "Count (n%)",
                      Age_fmt = "Med. Age (IQR)",
                      Male_fmt = "Male (n%)",
                      Female_fmt = "Female (n%)")
  )

# Print table
print(table_combo)
