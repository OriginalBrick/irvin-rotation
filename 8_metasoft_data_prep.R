library(tidyverse)

# Metabolite Data
data_aa_mtb <- read.csv("mwas_results/mwas_interaction_aa_metabolites.csv", header = TRUE) %>%
  select(metabolite, estimate, std_err) %>%
  rename(estimate_aa = estimate, std_err_aa = std_err)

data_ea_mtb <- read.csv("mwas_results/mwas_interaction_ea_metabolites.csv", header = TRUE) %>%
  select(metabolite, estimate, std_err) %>%
  rename(estimate_ea = estimate, std_err_ea = std_err)

data_mtb <- merge(data_aa_mtb, data_ea_mtb, by = "metabolite")
#write.csv(data_mtb, file = "../metabolites.csv", row.names = FALSE)

output_mtb <- apply(data_mtb, 1, function(row) paste(row, collapse = " "))
writeLines(output_mtb, "../metabolites.txt")

# Lipid Data
data_aa_lpd <- read.csv("mwas_results/mwas_interaction_aa_lipids.csv", header = TRUE) %>%
  select(metabolite, estimate, std_err) %>%
  rename(lipid = metabolite, estimate_aa = estimate, std_err_aa = std_err)

data_ea_lpd <- read.csv("mwas_results/mwas_interaction_ea_lipids.csv", header = TRUE) %>%
  select(metabolite, estimate, std_err) %>%
  rename(lipid = metabolite, estimate_ea = estimate, std_err_ea = std_err)

data_lpd <- merge(data_aa_lpd, data_ea_lpd, by = "lipid")
#write.csv(data_lpd, file = "../lipids.csv", row.names = FALSE)

output_lpd <- apply(data_lpd, 1, function(row) paste(row, collapse = " "))
writeLines(output_lpd, "../lipids.txt")