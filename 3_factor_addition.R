library(readxl)
library(tidyverse)

# Import our main data
data_aa <- read.csv("aa_clean.csv", header = TRUE)
data_ea <- read.csv("ea_clean.csv", header = TRUE)
data_aa_norm <- read.csv("aa_norm.csv", header = TRUE)
data_ea_norm <- read.csv("ea_norm.csv", header = TRUE)

# Import our factors too
factor_eq <- read_excel("for_jedediah_equamax.xlsx")
factor_vr <- read_excel("for_jedediah_varimax.xlsx")

# Rename equamax column headers
factor_eq_clean <- factor_eq %>%
  
  # Make column names lowercase and rename for readability
  rename_with(tolower) %>%
  rename_with(function(x) str_replace(x, "^factor", "equamax"), matches("^factor[0-9]+$")) %>%
  
  # Rearrange columns so that demographic data is first
  select(id_num, everything()) %>%
  
  # Made the id column a number
  mutate(id_num = as.integer(id_num)) %>%
  
  # Remove duplicate rows by id
  distinct(id_num, .keep_all = TRUE)

# Rename varimax column headers
factor_vr_clean <- factor_vr %>%
  
  # Make column names lowercase and rename for readability
  rename_with(tolower) %>%
  rename_with(function(x) str_replace(x, "^factor", "varimax"), matches("^factor[0-9]+$")) %>%
  
  # Rearrange columns so that demographic data is first
  select(id_num, everything()) %>%
  
  # Made the id column a number
  mutate(id_num = as.integer(id_num)) %>%
  
  # Remove duplicate rows by id
  distinct(id_num, .keep_all = TRUE)

# Grab colnames
eq_cols <- colnames(factor_eq_clean)[-1]
vr_cols <- colnames(factor_vr_clean)[-1]

# Merge datasets
data_aa_factors <- data_aa_norm %>%
  distinct(id_num, .keep_all = TRUE) %>%
  left_join(factor_eq_clean, by = "id_num") %>%
  left_join(factor_vr_clean, by = "id_num") %>%
  relocate(all_of(eq_cols), .after = xanthurenic_acid) %>%
  relocate(all_of(vr_cols), .after = equamax13) %>%
  relocate(s1p_, .before = acetoacetic_acid)
data_ea_factors <- data_ea_norm %>%
  distinct(id_num, .keep_all = TRUE) %>%
  left_join(factor_eq_clean, by = "id_num") %>%
  left_join(factor_vr_clean, by = "id_num") %>%
  relocate(all_of(eq_cols), .after = xanthurenic_acid) %>%
  relocate(all_of(vr_cols), .after = equamax13) %>%
  relocate(s1p_, .before = acetoacetic_acid)

# Save for later
write.csv(data_aa_factors, "aa_combo.csv", row.names = FALSE)
write.csv(data_ea_factors, "ea_combo.csv", row.names = FALSE)

# Create a combined dataset
data_combined <- bind_rows(data_aa_factors, data_ea_factors)

# Save for later
write.csv(data_combined, "combo_combo.csv", row.names = FALSE)