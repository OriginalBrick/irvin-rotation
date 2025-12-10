library(tidyverse)
library(ggplot2)
library(gridExtra)
library(writexl)

# Import data
mwas_int_aa <- read.csv("mwas_results/mwas_interaction_aa.csv", header = TRUE)
mwas_int_ea <- read.csv("mwas_results/mwas_interaction_ea.csv", header = TRUE)
# mwas_int_combo <- read.csv("mwas_results/mwas_interaction_combo.csv", header = TRUE)
mwas_main_aa <- read.csv("mwas_results/mwas_main_aa.csv", header = TRUE)
mwas_main_ea <- read.csv("mwas_results/mwas_main_ea.csv", header = TRUE)
# mwas_main_combo <- read.csv("mwas_results/mwas_main_combo.csv", header = TRUE)
mwas_sint_aa <- read.csv("mwas_results/mwas_specific_interaction_aa.csv", header = TRUE)
mwas_sint_ea <- read.csv("mwas_results/mwas_specific_interaction_ea.csv", header = TRUE)
# mwas_sint_combo <- read.csv("mwas_results/mwas_specific_interaction_combo.csv", header = TRUE)

# Grab the top metabolites
top_int_aa <- mwas_int_aa %>%
  # Add FDR correction
  mutate(joint_fdr = p.adjust(joint_p, method = "fdr")) %>%
  # Filter on p-value threshold
  filter(joint_p < 0.05) %>%
  # Add confidence intervals and order metabolites
  mutate(
    lower_ci = estimate - 1.96 * std_err,
    upper_ci = estimate + 1.96 * std_err,
    metabolite = factor(metabolite, levels = metabolite[order(estimate)])
  )
top_int_ea <- mwas_int_ea %>%
  # Add FDR correction
  mutate(joint_fdr = p.adjust(joint_p, method = "fdr")) %>%
  # Filter on p-value threshold
  filter(joint_p < 0.05) %>%
  # Add confidence intervals and order metabolites
  mutate(
    lower_ci = estimate - 1.96 * std_err,
    upper_ci = estimate + 1.96 * std_err,
    metabolite = factor(metabolite, levels = metabolite[order(estimate)])
  )
top_int_combo <- mwas_int_combo %>%
  # Add FDR correction
  mutate(joint_fdr = p.adjust(joint_p, method = "fdr")) %>%
  # Filter on p-value threshold
  filter(joint_p < 0.05) %>%
  # Add confidence intervals and order metabolites
  mutate(
    lower_ci = estimate - 1.96 * std_err,
    upper_ci = estimate + 1.96 * std_err,
    metabolite = factor(metabolite, levels = metabolite[order(estimate)])
  )
top_main_aa <- mwas_int_aa[mwas_main_aa$wald_p < 0.05, ]
top_main_ea <- mwas_int_ea[mwas_main_ea$wald_p < 0.05, ]
top_main_combo <- mwas_int_combo[mwas_main_combo$wald_p < 0.05, ]
top_sint_aa <- mwas_sint_aa %>%
  filter(metabolite %in% top_int_aa$metabolite) %>%
  mutate(apoe4 = factor(apoe4, levels = c(0,1),labels = c("Non-Carrier","Carrier")))
top_sint_ea <- mwas_sint_ea %>%
  filter(metabolite %in% top_int_ea$metabolite) %>%
  mutate(apoe4 = factor(apoe4, levels = c(0,1),labels = c("Non-Carrier","Carrier")))
top_sint_combo <- mwas_sint_combo %>%
  filter(metabolite %in% top_int_combo$metabolite) %>%
  mutate(apoe4 = factor(apoe4, levels = c(0,1),labels = c("Non-Carrier","Carrier")))

# Create comparison tables
mwas_int_aa <- mwas_int_aa %>%
  mutate(joint_fdr = p.adjust(joint_p, method = "fdr"))
mwas_int_ea <- mwas_int_ea %>%
  mutate(joint_fdr = p.adjust(joint_p, method = "fdr"))
mwas_int_combo <- mwas_int_combo %>%
  mutate(joint_fdr = p.adjust(joint_p, method = "fdr"))

aa_table <- mwas_int_aa %>%
  select(metabolite, joint_p, joint_fdr) %>%
  filter(joint_p < 0.05)
write_xlsx(aa_table, "mwas_results/aa_results.xlsx")

ea_table <- mwas_int_ea %>%
  select(metabolite, joint_p, joint_fdr) %>%
  filter(joint_p < 0.05)
write_xlsx(ea_table, "mwas_results/ea_results.xlsx")

combo_table <- mwas_int_combo %>%
  select(metabolite, joint_p, joint_fdr) %>%
  filter(joint_p < 0.05)
write_xlsx(combo_table, "mwas_results/combo_results.xlsx")

combined_table <- mwas_int_aa %>%
  select(metabolite, joint_p, joint_fdr) %>%
  rename(joint_p_aa = joint_p, joint_fdr_aa = joint_fdr) %>%
  full_join(
    mwas_int_ea %>%
      select(metabolite, joint_p, joint_fdr) %>%
      rename(joint_p_ea = joint_p, joint_fdr_ea = joint_fdr),
    by = "metabolite") %>%
  full_join(
    mwas_int_combo %>%
      select(metabolite, joint_p, joint_fdr) %>%
      rename(joint_p_cb = joint_p, joint_fdr_cb = joint_fdr),
    by = "metabolite") %>%
  filter(joint_p_aa < 0.05 | joint_p_ea < 0.05 | joint_p_cb < 0.05)
write_xlsx(combined_table, "mwas_results/fulljoin_results.xlsx")

# Manhattan Plots
manhattan_int_aa <- ggplot(mwas_int_aa, aes(x = metabolite, y = log_joint_p)) +
  geom_point(color = "steelblue") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +  # p = 0.05
  scale_y_continuous(limits = c(0, 5)) +  # sets y-axis from 0 to 4
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "AA w/ Metabolite-APOE4 Interaction",
       x = "Metabolite",
       y = "-log10(p-value)"
  )
manhattan_main_aa <- ggplot(mwas_main_aa, aes(x = metabolite, y = log_wald_p)) +
  geom_point(color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +  # p = 0.05
  scale_y_continuous(limits = c(0, 5)) +  # sets y-axis from 0 to 4
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "AA w/ APOE4 Covariate",
       x = "Metabolite",
       y = "-log10(p-value)"
  )
manhattan_int_ea <- ggplot(mwas_int_ea, aes(x = metabolite, y = log_joint_p)) +
  geom_point(color = "steelblue") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +  # p = 0.05
  scale_y_continuous(limits = c(0, 5)) +  # sets y-axis from 0 to 4
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "EA w/ Metabolite-APOE4 Interaction",
       x = "Metabolite",
       y = "-log10(p-value)"
  )
manhattan_main_ea <- ggplot(mwas_main_ea, aes(x = metabolite, y = log_wald_p)) +
  geom_point(color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +  # p = 0.05
  scale_y_continuous(limits = c(0, 5)) +  # sets y-axis from 0 to 4
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "EA w/ APOE4 Covariate",
       x = "Metabolite",
       y = "-log10(p-value)"
  )
manhattan_int_combo <- ggplot(mwas_int_combo, aes(x = metabolite, y = log_joint_p)) +
  geom_point(color = "steelblue") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +  # p = 0.05
  scale_y_continuous(limits = c(0, 5)) +  # sets y-axis from 0 to 4
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "Combo w/ Metabolite-APOE4 Interaction",
       x = "Metabolite",
       y = "-log10(p-value)"
  )
manhattan_main_combo <- ggplot(mwas_main_combo, aes(x = metabolite, y = log_wald_p)) +
  geom_point(color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +  # p = 0.05
  scale_y_continuous(limits = c(0, 5)) +  # sets y-axis from 0 to 4
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "Combo w/ APOE4 Covariate",
       x = "Metabolite",
       y = "-log10(p-value)"
  )

pdf("mwas_results/mwas_manhattan_plots.pdf", width = 12, height = 6)
grid.arrange(manhattan_main_aa, manhattan_main_ea, manhattan_main_combo, ncol = 3)
grid.arrange(manhattan_int_aa, manhattan_int_ea, manhattan_int_combo, ncol = 3)
dev.off()

# Forest Plots
forest_aa <- ggplot(top_int_aa, aes(x = metabolite, y = estimate)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, color = "gray40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(-2.2, 2.2)) +
  coord_flip() +
  labs(
    title = "AA w/ Metabolite-APOE4 Interaction",
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
    title = "EA w/ Metabolite-APOE4 Interaction",
    x = "Metabolite",
    y = "Log Odds Ratio (95% CI)"
  ) +
  theme_minimal()
forest_combo <- ggplot(top_int_combo, aes(x = metabolite, y = estimate)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, color = "gray40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(-2.2, 2.2)) +
  coord_flip() +
  labs(
    title = "Combo w/ Metabolite-APOE4 Interaction",
    x = "Metabolite",
    y = "Log Odds Ratio (95% CI)"
  ) +
  theme_minimal()

pdf("mwas_results/mwas_forest_plots.pdf", width = 12, height = 6)
print(forest_aa)
print(forest_ea)
print(forest_combo)
dev.off()

# Stratum Specific Forest Plots
ss_forest_aa <- ggplot(top_sint_aa,
       aes(x = metabolite,
           y = estimate,
           color = apoe4)) +
  geom_point(size = 3,
             position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(-2.2, 2.2)) +
  coord_flip() +
  labs(
    title = "AA w/ Metabolite-APOE4 Stratified Interaction",
    x = "Metabolite",
    y = "Log Odds Ratio (95% CI)",
    color = "APOE4 Status"
  ) +
  theme_minimal()
ss_forest_ea <- ggplot(top_sint_ea,
       aes(x = metabolite,
           y = estimate,
           color = apoe4)) +
  geom_point(size = 3,
             position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(-2.2, 2.2)) +
  coord_flip() +
  labs(
    title = "EA w/ Metabolite-APOE4 Stratified Interaction",
    x = "Metabolite",
    y = "Log Odds Ratio (95% CI)",
    color = "APOE4 Status"
  ) +
  theme_minimal()
ss_forest_combo <- ggplot(top_sint_combo,
                       aes(x = metabolite,
                           y = estimate,
                           color = apoe4)) +
  geom_point(size = 3,
             position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(-2.2, 2.2)) +
  coord_flip() +
  labs(
    title = "Combo w/ Metabolite-APOE4 Stratified Interaction",
    x = "Metabolite",
    y = "Log Odds Ratio (95% CI)",
    color = "APOE4 Status"
  ) +
  theme_minimal()

pdf("mwas_results/mwas_specific_forest_plots.pdf", width = 12, height = 6)
print(ss_forest_aa)
print(ss_forest_ea)
print(ss_forest_combo)
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
    title = "AA w/ Metabolite-APOE4 Interaction",
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
    title = "EA w/ Metabolite-APOE4 Interaction",
    x = "Log Odds Ratio",
    y = "-log10(p-value)"
  ) +
  theme_minimal()
volcano_combo <- ggplot(mwas_int_combo, aes(x = estimate, y = log_joint_p)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text(aes(label = metabolite), hjust = 0, vjust = 1.2, size = 3, na.rm = TRUE) +
  # scale_y_continuous(limits = c(0, 2.2)) +
  # scale_x_continuous(limits = c(-1.2, 1.2)) +
  labs(
    title = "Combo w/ Metabolite-APOE4 Interaction",
    x = "Log Odds Ratio",
    y = "-log10(p-value)"
  ) +
  theme_minimal()

pdf("mwas_results/mwas_volcano_plots.pdf", width = 12, height = 6)
print(volcano_aa)
print(volcano_ea)
print(volcano_combo)
dev.off()
