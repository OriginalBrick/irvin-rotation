library(tidyverse)
library(car)

# Import data
data_aa <- read.csv("aa_combo.csv", header = TRUE) %>%
  column_to_rownames("id_num")
data_ea <- read.csv("ea_combo.csv", header = TRUE) %>%
  column_to_rownames("id_num")
# data_combo <- read.csv("combo_combo.csv", header = TRUE) %>%
#  column_to_rownames("id_num")

# Create metabolite list
metabolites <- data_aa %>% select(lpc226_4:xanthurenic_acid)
metabolites <- colnames(metabolites)

# Create factor lists
# equamax <- paste0("equamax", 1:13)
# varimax <- paste0("varimax", 1:12)

# Create curated factor datasets
# data_aa_factors <- data_aa %>%
#   filter(if_all(equamax1:varimax12, ~ !is.na(.)))
# data_ea_factors <- data_ea %>%
#   filter(if_all(equamax1:varimax12, ~ !is.na(.)))

# Create interaction model function
run_mwas_interaction <- function(main_data, metabolite_data, output_file) {
  results <- list()
  metabolites <- metabolite_data
  
  for (met in metabolites) {
    formula <- as.formula(paste("incident_ecb ~", met, "* apoe_4 + age + gender + baseline_stroke"))
    model <- glm(formula, data = main_data, family = binomial)
    
    interaction_term <- paste(met, ":apoe_4", sep = "")
    coef_summary <- summary(model)$coefficients
    
    # Calculate joint p-value
    joint_test <- linearHypothesis(model, c(paste0(met, " = 0"), paste0(met, ":apoe_4 = 0")), test = "Chisq")
    joint_p <- joint_test$"Pr(>Chisq)"[2]
  
    # Extract interaction effect
    if (interaction_term %in% rownames(coef_summary)) {
      stats <- coef_summary[interaction_term, c("Estimate", "Std. Error", "Pr(>|z|)")]
      names(stats) <- c("estimate", "std_err", "wald_p")
    } else {
      stats <- c(estimate = NA, std_error = NA, wald_p = NA)
    }
    
    # Add joint p-value to results
    results[[met]] <- c(stats, joint_p = joint_p)
  }
  
  # Convert to data frame
  df <- do.call(rbind, results)
  df <- as.data.frame(df)
  df$metabolite <- rownames(df)
  df$wald_p <- as.numeric(df$wald_p)
  df$log_wald_p <- -log10(df$wald_p)
  df$log_joint_p <- -log10(df$joint_p)
  df <- df[, c("metabolite", "estimate", "std_err", "wald_p", "log_wald_p", "joint_p", "log_joint_p")]
  
  # Save to file
  write.csv(df, file = output_file, row.names = FALSE)
  
  return(df)
}

# Create a stratum-specific interaction model function
run_mwas_specific_interaction <- function(main_data, metabolite_data, output_file) {
  results <- list()
  metabolites <- metabolite_data
  
  for (met in metabolites) {
    formula <- as.formula(paste("incident_ecb ~", met, "* apoe_4 + age + gender + baseline_stroke"))
    model <- glm(formula, data = data_aa, family = binomial) #main data?
    
    coef_summary <- summary(model)$coefficients
    coefs <- coef(model)
    vc <- vcov(model)  # variance-covariance matrix
    
    main_term <- met
    interaction_term <- paste0(met, ":apoe_4")
    
    # Now safely extract
    beta_m <- coefs[main_term]
    var_m  <- vc[main_term, main_term]
    
    # APOE4 = 0
    est0 <- beta_m
    se0  <- sqrt(var_m)
    
    # APOE4 = 1
    if (interaction_term %in% names(coefs)) {
      beta_i <- coefs[interaction_term]
      var_i  <- vc[interaction_term, interaction_term]
      cov_mi <- vc[main_term, interaction_term]
      
      est1 <- beta_m + beta_i
      se1  <- sqrt(var_m + var_i + 2 * cov_mi)
    } else {
      est1 <- NA_real_
      se1  <- NA_real_
    }
    
    # Store both strata
    results[[met]] <- data.frame(
      metabolite = met,
      apoe4 = c(0, 1),
      estimate = c(est0, est1),
      std_err = c(se0, se1),
      lower_ci = c(est0 - 1.96*se0, if (!is.na(se1)) est1 - 1.96*se1 else NA_real_),
      upper_ci = c(est0 + 1.96*se0, if (!is.na(se1)) est1 + 1.96*se1 else NA_real_)
    )
  }
  
  df <- do.call(rbind, results)
  write.csv(df, file = output_file, row.names = FALSE)
  return(df)
}

# Create factor model function
# Made a separate model without baseline_stroke covariate because this broke the factors
run_factor_interaction <- function(main_data, metabolite_data, output_file) {
  results <- list()
  metabolites <- metabolite_data
  
  for (met in metabolites) {
    formula <- as.formula(paste("incident_ecb ~", met, "* apoe_4 + age + gender"))
    model <- glm(formula, data = main_data, family = binomial)
    
    interaction_term <- paste(met, ":apoe_4", sep = "")
    coef_summary <- summary(model)$coefficients
    
    # Calculate joint p-value
    joint_test <- linearHypothesis(model, c(paste0(met, " = 0"), paste0(met, ":apoe_4 = 0")), test = "Chisq")
    joint_p <- joint_test$"Pr(>Chisq)"[2]
    
    # Extract interaction effect
    if (interaction_term %in% rownames(coef_summary)) {
      stats <- coef_summary[interaction_term, c("Estimate", "Std. Error", "Pr(>|z|)")]
      names(stats) <- c("estimate", "std_err", "wald_p")
    } else {
      stats <- c(estimate = NA, std_error = NA, wald_p = NA)
    }
    
    # Add joint p-value to results
    results[[met]] <- c(stats, joint_p = joint_p)
  }
  
  # Convert to data frame
  df <- do.call(rbind, results)
  df <- as.data.frame(df)
  df$metabolite <- rownames(df)
  df$wald_p <- as.numeric(df$wald_p)
  df$log_wald_p <- -log10(df$wald_p)
  df$log_joint_p <- -log10(df$joint_p)
  df <- df[, c("metabolite", "estimate", "std_err", "wald_p", "log_wald_p", "joint_p", "log_joint_p")]
  
  # Save to file
  write.csv(df, file = output_file, row.names = FALSE)
  
  return(df)
}

# Create main model function
run_mwas_main <- function(main_data, metabolite_data, output_file) {
  results <- list()
  metabolites <- metabolite_data
  
  for (met in metabolites) {
    # Main effects model: no interaction
    formula <- as.formula(paste("incident_ecb ~", met, "+ apoe_4 + age + gender + baseline_stroke"))
    model <- glm(formula, data = main_data, family = binomial)
    
    coef_summary <- summary(model)$coefficients
    
    # Extract the metabolite's main effect
    if (met %in% rownames(coef_summary)) {
      stats <- coef_summary[met, c("Estimate", "Std. Error", "Pr(>|z|)")]
      names(stats) <- c("estimate", "std_err", "wald_p")
      results[[met]] <- stats
    } else {
      results[[met]] <- c(estimate = NA, std_err = NA, wald_p = NA)
    }
  }
  
  # Convert to data frame
  df <- do.call(rbind, results)
  df <- as.data.frame(df)
  df$metabolite <- rownames(df)
  df$wald_p <- as.numeric(df$wald_p)
  df$log_wald_p <- -log10(df$wald_p)
  df <- df[, c("metabolite", "estimate", "std_err", "wald_p", "log_wald_p")]
  
  # Save to file
  write.csv(df, file = output_file, row.names = FALSE)
  
  return(df)
}

# Create main model function
# Made a separate model without baseline_stroke covariate because this broke the factors
run_factor_main <- function(main_data, metabolite_data, output_file) {
  results <- list()
  metabolites <- metabolite_data
  
  for (met in metabolites) {
    # Main effects model: no interaction
    formula <- as.formula(paste("incident_ecb ~", met, "+ apoe_4 + age + gender"))
    model <- glm(formula, data = main_data, family = binomial)
    
    coef_summary <- summary(model)$coefficients
    
    # Extract the metabolite's main effect
    if (met %in% rownames(coef_summary)) {
      stats <- coef_summary[met, c("Estimate", "Std. Error", "Pr(>|z|)")]
      names(stats) <- c("estimate", "std_err", "wald_p")
      results[[met]] <- stats
    } else {
      results[[met]] <- c(estimate = NA, std_err = NA, wald_p = NA)
    }
  }
  
  # Convert to data frame
  df <- do.call(rbind, results)
  df <- as.data.frame(df)
  df$metabolite <- rownames(df)
  df$wald_p <- as.numeric(df$wald_p)
  df$log_wald_p <- -log10(df$wald_p)
  df <- df[, c("metabolite", "estimate", "std_err", "wald_p", "log_wald_p")]
  
  # Save to file
  write.csv(df, file = output_file, row.names = FALSE)
  
  return(df)
}

# Run our metabolite functions!
run_mwas_interaction(data_aa, metabolites, "mwas_results/mwas_interaction_aa.csv")
run_mwas_interaction(data_ea, metabolites, "mwas_results/mwas_interaction_ea.csv")
run_mwas_interaction(data_combo, metabolites, "mwas_results/mwas_interaction_combo.csv")
run_mwas_main(data_aa, metabolites, "mwas_results/mwas_main_aa.csv")
run_mwas_main(data_ea, metabolites, "mwas_results/mwas_main_ea.csv")
run_mwas_main(data_combo, metabolites, "mwas_results/mwas_main_combo.csv")

# Run our factor functions!
# run_factor_interaction(data_aa, equamax, "equamax_results/equamax_interaction_aa.csv")
# run_factor_interaction(data_ea, equamax, "equamax_results/equamax_interaction_ea.csv")
# run_factor_main(data_aa, equamax, "equamax_results/equamax_main_aa.csv")
# run_factor_main(data_ea, equamax, "equamax_results/equamax_main_ea.csv")
# run_factor_interaction(data_aa, varimax, "varimax_results/varimax_interaction_aa.csv")
# run_factor_interaction(data_ea, varimax, "varimax_results/varimax_interaction_ea.csv")
# run_factor_main(data_aa, varimax, "varimax_results/varimax_main_aa.csv")
# run_factor_main(data_ea, varimax, "varimax_results/varimax_main_ea.csv")

# Run the stratum-specific functions
run_mwas_specific_interaction(data_aa, metabolites, "mwas_results/mwas_specific_interaction_aa.csv")
run_mwas_specific_interaction(data_ea, metabolites, "mwas_results/mwas_specific_interaction_ea.csv")
run_mwas_specific_interaction(data_combo, metabolites, "mwas_results/mwas_specific_interaction_combo.csv")