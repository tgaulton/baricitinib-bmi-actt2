################################################################################
# 06_Model_Adequacy.R
#
# Purpose: Model adequacy diagnostics for the primary Fine-Gray model.
#
# Three diagnostics:
#   1. AIC comparison of linear vs spline BMI
#   2. Likelihood ratio test for non-linearity in BMI
#   3. Proportional subdistribution hazards test for treatment
#   4. Influence diagnostics using standardized DFBETAs
#
# NOTE: Requires 01_Analytic_Sample.R to have been run first.
################################################################################

library(tidyverse)
library(cmprsk)
library(splines)
library(survival)

cat("================================================================================\n")
cat("MODEL ADEQUACY DIAGNOSTICS\n")
cat("================================================================================\n\n")

################################################################################
# LOAD DATA
################################################################################

data_final <- readRDS("02_Data/processed/analytic_dataset.rds")

bmi_mean <- mean(data_final$BMI_kgm2)

data_model <- data_final %>%
  mutate(
    BMI_c = BMI_kgm2 - bmi_mean,
    int_linear = TRTP_bin * BMI_c
  )

################################################################################
# DIAGNOSTIC 1 AND 2: AIC AND LRT FOR LINEAR VS SPLINE BMI
################################################################################

cat("Diagnostics 1-2: Linear vs spline BMI\n\n")

# Linear specification (primary)
mm_linear <- as.matrix(data_model[, c("TRTP_bin", "BMI_c", "int_linear",
                                       "AGE", "SEX_M",
                                       "sev5", "sev6", "sev7")])

fg_linear <- crr(
  ftime    = data_model$ftime,
  fstatus  = data_model$fstatus,
  cov1     = mm_linear,
  failcode = 1, cencode = 0
)

# Spline specification (alternative)
bmi_spline <- ns(data_model$BMI_c, df = 3)
spline_df <- as.data.frame(bmi_spline)
colnames(spline_df) <- paste0("spl", 1:3)

data_model_spl <- bind_cols(data_model, spline_df) %>%
  mutate(
    int_spl1 = TRTP_bin * spl1,
    int_spl2 = TRTP_bin * spl2,
    int_spl3 = TRTP_bin * spl3
  )

mm_spline <- as.matrix(data_model_spl[, c("TRTP_bin",
                                            "spl1", "spl2", "spl3",
                                            "int_spl1", "int_spl2", "int_spl3",
                                            "AGE", "SEX_M",
                                            "sev5", "sev6", "sev7")])

fg_spline <- crr(
  ftime    = data_model$ftime,
  fstatus  = data_model$fstatus,
  cov1     = mm_spline,
  failcode = 1, cencode = 0
)

# AIC computation
AIC_crr <- function(fit) {
  -2 * fit$loglik + 2 * length(fit$coef)
}

aic_linear <- AIC_crr(fg_linear)
aic_spline <- AIC_crr(fg_spline)
delta_aic <- aic_linear - aic_spline

cat("AIC linear:", round(aic_linear, 2), "\n")
cat("AIC spline:", round(aic_spline, 2), "\n")
cat("Delta AIC (linear - spline):", round(delta_aic, 2), "\n")
cat("Lower AIC favored:", ifelse(delta_aic < 0, "linear", "spline"), "\n\n")

# LRT for non-linearity (4 df: 3 spline terms + 3 spline interactions vs 1 BMI + 1 int_linear)
lrt_nonlin <- 2 * (fg_spline$loglik - fg_linear$loglik)
df_nonlin  <- 4
p_nonlin   <- pchisq(lrt_nonlin, df = df_nonlin, lower.tail = FALSE)

cat("LRT for non-linearity:\n")
cat("  chi-square:", round(lrt_nonlin, 2), "\n")
cat("  df:", df_nonlin, "\n")
cat("  p-value:", signif(p_nonlin, 3), "\n\n")

################################################################################
# DIAGNOSTIC 3: PROPORTIONAL SUBDISTRIBUTION HAZARDS FOR TREATMENT
################################################################################

cat("Diagnostic 3: Proportional subdistribution hazards test for treatment\n\n")

# Create expanded Fine-Gray dataset (Austin & Fine 2017)
fg_data <- data_model %>%
  select(ftime, fstatus, TRTP_bin, BMI_c, int_linear,
         AGE, SEX_M, sev5, sev6, sev7) %>%
  mutate(fstatus_f = factor(fstatus, levels = c(0, 1, 2),
                            labels = c("censored", "recovery", "death")))

fg_expanded <- finegray(
  Surv(ftime, fstatus_f) ~ .,
  data = fg_data,
  etype = "recovery"
)

# Fit base model (no time-varying term)
cox_base <- coxph(
  Surv(fgstart, fgstop, fgstatus) ~ TRTP_bin + BMI_c + int_linear +
    AGE + SEX_M + sev5 + sev6 + sev7,
  data = fg_expanded,
  weights = fgwt
)

# Fit model with time x treatment interaction
cox_tt <- coxph(
  Surv(fgstart, fgstop, fgstatus) ~ TRTP_bin + BMI_c + int_linear +
    AGE + SEX_M + sev5 + sev6 + sev7 + tt(TRTP_bin),
  data = fg_expanded,
  weights = fgwt,
  tt = function(x, t, ...) x * log(t)
)

lrt_ph <- 2 * (cox_tt$loglik[2] - cox_base$loglik[2])
p_ph   <- pchisq(lrt_ph, df = 1, lower.tail = FALSE)

cat("Proportional hazards test (time x treatment interaction):\n")
cat("  chi-square:", round(lrt_ph, 2), "\n")
cat("  df: 1\n")
cat("  p-value:", signif(p_ph, 3), "\n\n")

################################################################################
# DIAGNOSTIC 4: INFLUENTIAL OBSERVATIONS (DFBETAs)
################################################################################

cat("Diagnostic 4: Influence diagnostics (standardized DFBETAs)\n\n")

# Standardized DFBETAs from the expanded Cox model
dfbeta_full <- residuals(cox_base, type = "dfbeta")
se_coefs <- summary(cox_base)$coefficients[, "se(coef)"]
dfbeta_std <- sweep(dfbeta_full, 2, se_coefs, FUN = "/")

# Examine treatment, BMI, and interaction coefficients
key_indices <- 1:3
names(key_indices) <- c("TRTP_bin", "BMI_c", "int_linear")

threshold <- 2 / sqrt(nrow(fg_expanded))

cat("Threshold for influence: |DFBETA| >", round(threshold, 3), "\n\n")

for (idx in key_indices) {
  coef_name <- names(key_indices)[which(key_indices == idx)]
  max_dfbeta <- max(abs(dfbeta_std[, idx]))
  n_above <- sum(abs(dfbeta_std[, idx]) > threshold)
  
  cat(sprintf("%-12s  max |DFBETA| = %.3f   obs above threshold: %d\n",
              coef_name, max_dfbeta, n_above))
}

cat("\n")

################################################################################
# SAVE DIAGNOSTICS SUMMARY
################################################################################

dir.create("04_Results/tables", recursive = TRUE, showWarnings = FALSE)

diagnostics <- tibble(
  Diagnostic = c("AIC linear", "AIC spline", "Delta AIC",
                 "LRT non-linearity chi-square",
                 "LRT non-linearity df",
                 "LRT non-linearity p",
                 "PH test chi-square",
                 "PH test df",
                 "PH test p"),
  Value = c(round(aic_linear, 2),
            round(aic_spline, 2),
            round(delta_aic, 2),
            round(lrt_nonlin, 2),
            df_nonlin,
            signif(p_nonlin, 3),
            round(lrt_ph, 2),
            1,
            signif(p_ph, 3))
)

write_csv(diagnostics, "04_Results/tables/model_adequacy_diagnostics.csv")
cat("Diagnostics saved to 04_Results/tables/model_adequacy_diagnostics.csv\n")
