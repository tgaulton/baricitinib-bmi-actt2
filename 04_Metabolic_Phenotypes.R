################################################################################
# 04_Metabolic_Phenotypes.R
#
# Purpose: Test treatment effect heterogeneity across four metabolic
#          phenotype groups defined by BMI and comorbidity status.
#
# Phenotype groups:
#   MHNW: BMI <30 kg/m^2 without diabetes or hypertension
#   MUNW: BMI <30 kg/m^2 with diabetes or hypertension
#   MHO:  BMI >=30 kg/m^2 without diabetes or hypertension
#   MUO:  BMI >=30 kg/m^2 with diabetes or hypertension
#
# Specification:
#   - Fine-Gray competing risk regression
#   - Phenotype x treatment interaction
#   - Adjusted for age, sex, baseline severity (factor)
#   - LRT for interaction (df = 3)
#   - Phenotype-specific HRs and 95% CIs via delta method
#
# Produces Figure 2 of the manuscript.
#
# NOTE: Requires 01_Analytic_Sample.R to have been run first.
################################################################################

library(tidyverse)
library(cmprsk)
library(ggplot2)

cat("================================================================================\n")
cat("METABOLIC PHENOTYPE ANALYSIS\n")
cat("================================================================================\n\n")

################################################################################
# LOAD DATA AND PREPARE
################################################################################

data_final <- readRDS("02_Data/processed/analytic_dataset.rds")

# Complete cases for diabetes/hypertension (1.8% missing)
data_pheno <- data_final %>%
  filter(!is.na(Diabetes), !is.na(Hypertension), !is.na(Phenotype_factor))

cat("Metabolic phenotype cohort: n =", nrow(data_pheno), "\n\n")

cat("Phenotype distribution:\n")
print(table(data_pheno$Phenotype_factor))
cat("\n")

################################################################################
# FIT MODELS (full and reduced)
################################################################################

# Build model matrices (drop the intercept since crr does not include one)
mm_full <- model.matrix(
  ~ TRTP_bin * Phenotype_factor + AGE + SEX_M + sev5 + sev6 + sev7,
  data = data_pheno
)[, -1]

mm_reduced <- model.matrix(
  ~ TRTP_bin + Phenotype_factor + AGE + SEX_M + sev5 + sev6 + sev7,
  data = data_pheno
)[, -1]

fg_full <- crr(
  ftime    = data_pheno$ftime,
  fstatus  = data_pheno$fstatus,
  cov1     = mm_full,
  failcode = 1,
  cencode  = 0
)

fg_reduced <- crr(
  ftime    = data_pheno$ftime,
  fstatus  = data_pheno$fstatus,
  cov1     = mm_reduced,
  failcode = 1,
  cencode  = 0
)

# LRT (df = 3 for the three interaction terms)
lrt_chi <- 2 * (fg_full$loglik - fg_reduced$loglik)
lrt_p   <- pchisq(lrt_chi, df = 3, lower.tail = FALSE)

cat("Phenotype x Treatment Interaction (Likelihood Ratio Test):\n")
cat("  chi-square:", round(lrt_chi, 3), "\n")
cat("  df: 3\n")
cat("  p-value:", signif(lrt_p, 3), "\n\n")

################################################################################
# COMPUTE PHENOTYPE-SPECIFIC HRs (delta method)
################################################################################

coef_full  <- fg_full$coef
vcov_full  <- fg_full$var
coef_names <- colnames(mm_full)

trtp_idx <- which(coef_names == "TRTP_bin")
int_idx  <- grep("TRTP_bin:Phenotype_factor", coef_names)

phenotypes <- c("MHNW", "MUNW", "MHO", "MUO")
phenotype_labels <- c(
  "BMI <30 kg/m\u00B2\nWithout DM/HTN",
  "BMI <30 kg/m\u00B2\nWith DM/HTN",
  "BMI \u226530 kg/m\u00B2\nWithout DM/HTN",
  "BMI \u226530 kg/m\u00B2\nWith DM/HTN"
)
n_pheno <- as.numeric(table(data_pheno$Phenotype_factor))

# Phenotype-specific HRs
results <- tibble(
  phenotype = phenotypes,
  label = phenotype_labels,
  n = n_pheno,
  log_HR = NA_real_,
  SE = NA_real_
)

# Reference category (MHNW)
results$log_HR[1] <- coef_full[trtp_idx]
results$SE[1]     <- sqrt(vcov_full[trtp_idx, trtp_idx])

# Other phenotypes (MUNW, MHO, MUO)
for (i in 2:4) {
  log_hr <- coef_full[trtp_idx] + coef_full[int_idx[i - 1]]
  var_total <- vcov_full[trtp_idx, trtp_idx] +
               vcov_full[int_idx[i - 1], int_idx[i - 1]] +
               2 * vcov_full[trtp_idx, int_idx[i - 1]]
  
  results$log_HR[i] <- log_hr
  results$SE[i]     <- sqrt(var_total)
}

results <- results %>%
  mutate(
    HR = exp(log_HR),
    HR_lower = exp(log_HR - 1.96 * SE),
    HR_upper = exp(log_HR + 1.96 * SE),
    HR_CI = sprintf("%.2f (%.2f-%.2f)", HR, HR_lower, HR_upper)
  )

cat("Phenotype-specific HRs:\n")
print(results %>% select(phenotype, n, HR, HR_lower, HR_upper))
cat("\n")

################################################################################
# FIGURE 2: METABOLIC SUBGROUP FOREST PLOT (BLACK AND WHITE)
################################################################################

results_plot <- results %>%
  mutate(
    label_with_n = sprintf("%s\n(n=%d)", label, n),
    label_with_n = fct_rev(fct_inorder(label_with_n))
  )

p_figure2 <- ggplot(results_plot, aes(y = label_with_n, x = HR)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "gray40", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper),
                 height = 0.25, linewidth = 0.9, color = "black") +
  geom_point(size = 3, shape = 21, fill = "white",
             color = "black", stroke = 1.2) +
  geom_text(aes(label = HR_CI),
            x = 1.65, size = 4, hjust = 0, color = "black") +
  scale_x_continuous(
    breaks = seq(0.5, 1.55, length.out = 4),
    limits = c(0.5, 2.0),
    expand = c(0, 0)
  ) +
  labs(x = "Hazard Ratio", y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y = element_text(hjust = 0, size = 10, color = "black"),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold",
                                 margin = margin(t = 10)),
    plot.margin = margin(15, 10, 10, 10)
  )

################################################################################
# SAVE OUTPUTS
################################################################################

dir.create("04_Results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Results/tables",  recursive = TRUE, showWarnings = FALSE)

write_csv(results %>% select(phenotype, n, HR, HR_lower, HR_upper, HR_CI),
          "04_Results/tables/figure2_metabolic_subgroup_hrs.csv")

ggsave("04_Results/figures/figure2_metabolic_subgroups.tiff",
       p_figure2, width = 10, height = 5.5, dpi = 500,
       device = "tiff", compression = "lzw")

cat("Figure 2 saved to 04_Results/figures/figure2_metabolic_subgroups.tiff\n")
cat("Phenotype HRs saved to 04_Results/tables/figure2_metabolic_subgroup_hrs.csv\n")
