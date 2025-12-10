library(tidyverse)
library(BiocManager)

# Import data
data_aa <- read.csv("aa_clean.csv", header = TRUE)
data_ea <- read.csv("ea_clean.csv", header = TRUE)

# Separate metabolites
data_aa_mtb <- data_aa %>%
  select(id_num, s1p_, anandamide:xanthurenic_acid)
data_ea_mtb <- data_ea %>%
  select(id_num, s1p_, anandamide:xanthurenic_acid)

# Separate lipids
data_aa_lpd <- data_aa %>%
  select(id_num, lpc226_4:tg5203rd_84)
data_ea_lpd <- data_ea %>%
  select(id_num, lpc226_4:tg5203rd_84)

# 10% Missingness filter by metabolite
missing_percent <- colMeans(is.na(data_aa_mtb))*100                             # Calculate the % missing
missing_percent[missing_percent > 0]                                            # Sanity check, none removed!
data_aa_mtb <- data_aa_mtb[, missing_percent <= 10]                             # Filter them if needed
missing_percent <- colMeans(is.na(data_ea_mtb))*100
missing_percent[missing_percent > 0]
data_ea_mtb <- data_ea_mtb[, missing_percent <= 10]

# 10% Missingness filter by lipid
missing_percent <- colMeans(is.na(data_aa_lpd))*100
missing_percent[missing_percent > 0]
data_aa_lpd <- data_aa_lpd[, missing_percent <= 10]
missing_percent <- colMeans(is.na(data_ea_lpd))*100
missing_percent[missing_percent > 0]
data_ea_lpd <- data_ea_lpd[, missing_percent <= 10]

# Create function for inverse-normal transformation
inverse_normal_transform <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

# Apply inverse-normal transformation to metabolite data
normalized_matrix <- data_aa_mtb[, -1]                                          # Drop ID column
normalized_matrix <- apply(normalized_matrix, 2, inverse_normal_transform)      # Normalize data
data_aa_mtb <- data.frame(id_num = data_aa_mtb[, 1], normalized_matrix)         # Add ID column
normalized_matrix <- data_ea_mtb[, -1]
normalized_matrix <- apply(normalized_matrix, 2, inverse_normal_transform)
data_ea_mtb <- data.frame(id_num = data_ea_mtb[, 1], normalized_matrix)

# Apply inverse-normal transformation to lipid data
normalized_matrix <- data_aa_lpd[, -1]
normalized_matrix <- apply(normalized_matrix, 2, inverse_normal_transform)
data_aa_lpd <- data.frame(id_num = data_aa_lpd[, 1], normalized_matrix)
normalized_matrix <- data_ea_lpd[, -1]
normalized_matrix <- apply(normalized_matrix, 2, inverse_normal_transform)
data_ea_lpd <- data.frame(id_num = data_ea_lpd[, 1], normalized_matrix)

# Merge data back together
data_aa_norm <- data_aa                                                         # Create new object
data_aa_norm[, names(data_aa_mtb)[-1]] <- data_aa_mtb[, -1]                     # Merge metabolites
data_aa_norm[, names(data_aa_lpd)[-1]] <- data_aa_lpd[, -1]                     # Merge lipids
data_ea_norm <- data_ea
data_ea_norm[, names(data_ea_mtb)[-1]] <- data_ea_mtb[, -1]
data_ea_norm[, names(data_ea_lpd)[-1]] <- data_ea_lpd[, -1]

# Save for later
write.csv(data_aa_norm, "aa_norm.csv", row.names = FALSE)
write.csv(data_ea_norm, "ea_norm.csv", row.names = FALSE)

# Create a combined dataset
data_combined <- bind_rows(data_aa_norm, data_ea_norm)

# Save for later
write.csv(data_combined, "combo_norm.csv", row.names = FALSE)