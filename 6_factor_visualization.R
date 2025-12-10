library(tidyverse)
library(ggplot2)
library(gridExtra)

# Import equamax data
equamax_int_aa <- read.csv("equamax_results/equamax_interaction_aa.csv", header = TRUE)
equamax_int_ea <- read.csv("equamax_results/equamax_interaction_ea.csv", header = TRUE)
equamax_main_aa <- read.csv("equamax_results/equamax_main_aa.csv", header = TRUE)
equamax_main_ea <- read.csv("equamax_results/equamax_main_ea.csv", header = TRUE)

# Import varimax data
varimax_int_aa <- read.csv("varimax_results/varimax_interaction_aa.csv", header = TRUE)
varimax_int_ea <- read.csv("varimax_results/varimax_interaction_ea.csv", header = TRUE)
varimax_main_aa <- read.csv("varimax_results/varimax_main_aa.csv", header = TRUE)
varimax_main_ea <- read.csv("varimax_results/varimax_main_ea.csv", header = TRUE)

# Add equamax confidence intervals
equamax_int_aa <- equamax_int_aa %>%
  mutate(
    lower_ci = estimate - 1.96 * std_err,
    upper_ci = estimate + 1.96 * std_err
  )
equamax_int_ea <- equamax_int_ea %>%
  mutate(
    lower_ci = estimate - 1.96 * std_err,
    upper_ci = estimate + 1.96 * std_err  )

# Add varimax confidence intervals
varimax_int_aa <- varimax_int_aa %>%
  mutate(
    lower_ci = estimate - 1.96 * std_err,
    upper_ci = estimate + 1.96 * std_err
  )
varimax_int_ea <- varimax_int_ea %>%
  mutate(
    lower_ci = estimate - 1.96 * std_err,
    upper_ci = estimate + 1.96 * std_err  )

# Equamax Manhattan Plots
equamax_manhattan_int_aa <- ggplot(equamax_int_aa, aes(x = metabolite, y = log_joint_p)) +
  geom_point(color = "steelblue") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +  # p = 0.05
  scale_y_continuous(limits = c(0, 5)) +  # sets y-axis from 0 to 4
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "Equamax Manhattan Plot (AA w/ Metabolite-APOE4 Interaction)",
       x = "Factor",
       y = "-log10(p-value)"
  )
equamax_manhattan_main_aa <- ggplot(equamax_main_aa, aes(x = metabolite, y = log_wald_p)) +
  geom_point(color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +  # p = 0.05
  scale_y_continuous(limits = c(0, 5)) +  # sets y-axis from 0 to 4
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "Equamax Manhattan Plot (AA w/ APOE4 Covariate)",
       x = "Factor",
       y = "-log10(p-value)"
  )
equamax_manhattan_int_ea <- ggplot(equamax_int_ea, aes(x = metabolite, y = log_joint_p)) +
  geom_point(color = "steelblue") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +  # p = 0.05
  scale_y_continuous(limits = c(0, 5)) +  # sets y-axis from 0 to 4
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "Equamax Manhattan Plot (EA w/ Metabolite-APOE4 Interaction)",
       x = "Factor",
       y = "-log10(p-value)"
  )
equamax_manhattan_main_ea <- ggplot(equamax_main_ea, aes(x = metabolite, y = log_wald_p)) +
  geom_point(color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +  # p = 0.05
  scale_y_continuous(limits = c(0, 5)) +  # sets y-axis from 0 to 4
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "Equamax Manhattan Plot (EA w/ APOE4 Covariate)",
       x = "Factor",
       y = "-log10(p-value)"
  )

pdf("equamax_results/equamax_manhattan_plots.pdf", width = 12, height = 6)
grid.arrange(equamax_manhattan_main_aa, equamax_manhattan_main_ea, ncol = 2)
grid.arrange(equamax_manhattan_int_aa, equamax_manhattan_int_ea, ncol = 2)
dev.off()

# Varimax Manhattan Plots
varimax_manhattan_int_aa <- ggplot(varimax_int_aa, aes(x = metabolite, y = log_joint_p)) +
  geom_point(color = "steelblue") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +  # p = 0.05
  scale_y_continuous(limits = c(0, 5)) +  # sets y-axis from 0 to 4
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "Varimax Manhattan Plot (AA w/ Metabolite-APOE4 Interaction)",
       x = "Factor",
       y = "-log10(p-value)"
  )

varimax_manhattan_main_aa <- ggplot(varimax_main_aa, aes(x = metabolite, y = log_wald_p)) +
  geom_point(color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +  # p = 0.05
  scale_y_continuous(limits = c(0, 5)) +  # sets y-axis from 0 to 4
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "Varimax Manhattan Plot (AA w/ APOE4 Covariate)",
       x = "Factor",
       y = "-log10(p-value)"
  )

varimax_manhattan_int_ea <- ggplot(varimax_int_ea, aes(x = metabolite, y = log_joint_p)) +
  geom_point(color = "steelblue") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +  # p = 0.05
  scale_y_continuous(limits = c(0, 5)) +  # sets y-axis from 0 to 4
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "Varimax Manhattan Plot (EA w/ Metabolite-APOE4 Interaction)",
       x = "Factor",
       y = "-log10(p-value)"
  )

varimax_manhattan_main_ea <- ggplot(varimax_main_ea, aes(x = metabolite, y = log_wald_p)) +
  geom_point(color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +  # p = 0.05
  scale_y_continuous(limits = c(0, 5)) +  # sets y-axis from 0 to 4
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "Varimax Manhattan Plot (EA w/ APOE4 Covariate)",
       x = "Factor",
       y = "-log10(p-value)"
  )

pdf("varimax_results/varimax_manhattan_plots.pdf", width = 12, height = 6)
grid.arrange(varimax_manhattan_main_aa, varimax_manhattan_main_ea, ncol = 2)
grid.arrange(varimax_manhattan_int_aa, varimax_manhattan_int_ea, ncol = 2)
dev.off()

# Forest Plots
forest_aa <- ggplot(top_int_aa, aes(x = metabolite, y = estimate)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, color = "gray40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(-2.2, 2.2)) +
  coord_flip() +
  labs(
    title = "MWAS Forest Plot (AA w/ Metabolite-APOE4 Interaction)",
    x = "Metabolite",
    y = "Log Odds Ratio (95% CI)"
  ) +
  theme_minimal()
forest_ea <- ggplot(top_int_ea, aes(x = metabolite, y = estimate)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, color = "gray40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(-2.2, 2.2)) +
  coord_flip() +
  labs(
    title = "MWAS Forest Plot (AA w/ Metabolite-APOE4 Interaction)",
    x = "Metabolite",
    y = "Log Odds Ratio (95% CI)"
  ) +
  theme_minimal()

pdf("mwas_forest_plots.pdf", width = 12, height = 6)
print(forest_aa)
print(forest_ea)
dev.off()

# Volcano Plots
volcano_aa <- ggplot(mwas_int_aa, aes(x = estimate, y = log_joint_p)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text(aes(label = metabolite), hjust = 0, vjust = 1.2, size = 3, na.rm = TRUE) +
  # scale_y_continuous(limits = c(0, 2.5)) +
  # scale_x_continuous(limits = c(-1.5, 1.5)) +
  labs(
    title = "MWAS Volcano Plot (AA w/ Metabolite-APOE4 Interaction)",
    x = "Log Odds Ratio",
    y = "-log10(p-value)"
  ) +
  theme_minimal()
volcano_ea <- ggplot(mwas_int_ea, aes(x = estimate, y = log_joint_p)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text(aes(label = metabolite), hjust = 0, vjust = 1.2, size = 3, na.rm = TRUE) +
  # scale_y_continuous(limits = c(0, 2.2)) +
  # scale_x_continuous(limits = c(-1.2, 1.2)) +
  labs(
    title = "MWAS Volcano Plot (EA w/ Metabolite-APOE4 Interaction)",
    x = "Log Odds Ratio",
    y = "-log10(p-value)"
  ) +
  theme_minimal()

pdf("mwas_volcano_plots.pdf", width = 12, height = 6)
print(volcano_aa)
print(volcano_ea)
dev.off()