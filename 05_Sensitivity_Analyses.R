################################################################################
# 05_Sensitivity_Analyses.R
#
# Purpose: Pre-specified sensitivity analyses for the primary outcome.
#
# Three sensitivities:
#   1. Race and ethnicity adjustment (n=694)
#   2. Inclusion of participants with BMI >=50 kg/m^2 (n=976)
#   3. Natural cubic spline BMI (n=942)
#
# Produces Table S.5 of the supplement.
#
# NOTE: Requires 01_Analytic_Sample.R to have been run first.
################################################################################

library(tidyverse)
library(cmprsk)
library(splines)

cat("================================================================================\n")
cat("SENSITIVITY ANALYSES (PRIMARY OUTCOME)\n")
cat("================================================================================\n\n")

################################################################################
# LOAD DATA AND HELPERS
################################################################################

data_final <- readRDS("02_Data/processed/analytic_dataset.rds")

bmi_mean <- mean(data_final$BMI_kgm2)

data_model <- data_final %>%
  mutate(
    BMI_c = BMI_kgm2 - bmi_mean,
    int_linear = TRTP_bin * BMI_c
  )

# Helper: LRT between nested crr models
lrt_crr <- function(full, reduced, df) {
  chi <- 2 * (full$loglik - reduced$loglik)
  p <- pchisq(chi, df = df, lower.tail = FALSE)
  list(chi = chi, df = df, p = p)
}

################################################################################
# SENSITIVITY 1: RACE AND ETHNICITY ADJUSTMENT
################################################################################

cat("Sensitivity 1: Race and ethnicity adjustment\n")

# Need to re-derive from raw data since race/ethnicity not in analytic dataset
actt2_raw <- read_csv("02_Data/raw/ACTT2.csv", show_col_types = FALSE)

data_race <- actt2_raw %>%
  filter(!is.na(BMI), BMI >= 18.5, BMI < 50,
         !is.na(RACE), !is.na(ETHNIC),
         RACE != "UNKNOWN",
         ETHNIC != "UNKNOWN", ETHNIC != "NOT REPORTED") %>%
  mutate(
    TRTP_bin = ifelse(TRTP == "Baricitinib + RDV", 1, 0),
    SEX_M = ifelse(SEX == "M", 1L, 0L),
    ORDSCRG_num = as.integer(str_extract(ORDSCRG, "\\d+")),
    BMI_c = BMI - mean(BMI),
    int_linear = TRTP_bin * BMI_c,
    sev5 = as.integer(ORDSCRG_num == 5),
    sev6 = as.integer(ORDSCRG_num == 6),
    sev7 = as.integer(ORDSCRG_num == 7),
    RACE_factor = factor(RACE),
    ETHNIC_factor = factor(ETHNIC),
    ftime = case_when(
      RECCNSR == 0 & DTHCNSR == 0 ~ pmin(TTRECOV, TTDEATH, na.rm = TRUE),
      RECCNSR == 0 & DTHCNSR == 1 ~ TTRECOV,
      RECCNSR == 1 & DTHCNSR == 0 ~ TTDEATH,
      RECCNSR == 1 & DTHCNSR == 1 ~ ifelse(!is.na(TTRECOV), TTRECOV, 29),
      TRUE ~ NA_real_
    ),
    fstatus = case_when(
      DTHCNSR == 0 & (RECCNSR == 1 | TTDEATH <= TTRECOV) ~ 2L,
      RECCNSR == 0 & (DTHCNSR == 1 | TTRECOV  < TTDEATH) ~ 1L,
      RECCNSR == 1 & DTHCNSR == 1 ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(ftime > 0 | is.na(ftime))

mm_race_full <- model.matrix(
  ~ TRTP_bin + BMI_c + int_linear + AGE + SEX_M +
    sev5 + sev6 + sev7 + RACE_factor + ETHNIC_factor,
  data = data_race
)[, -1]

mm_race_reduced <- model.matrix(
  ~ TRTP_bin + BMI_c + AGE + SEX_M +
    sev5 + sev6 + sev7 + RACE_factor + ETHNIC_factor,
  data = data_race
)[, -1]

fg_race_full <- crr(ftime = data_race$ftime, fstatus = data_race$fstatus,
                    cov1 = mm_race_full, failcode = 1, cencode = 0)
fg_race_reduced <- crr(ftime = data_race$ftime, fstatus = data_race$fstatus,
                       cov1 = mm_race_reduced, failcode = 1, cencode = 0)

lrt_race <- lrt_crr(fg_race_full, fg_race_reduced, df = 1)

cat("  n:", nrow(data_race),
    " chi-square:", round(lrt_race$chi, 2),
    " p:", signif(lrt_race$p, 3), "\n\n")

################################################################################
# SENSITIVITY 2: INCLUSION OF BMI >=50 kg/m^2
################################################################################

cat("Sensitivity 2: Inclusion of BMI >=50 kg/m^2\n")

data_expanded <- actt2_raw %>%
  filter(!is.na(BMI), BMI >= 18.5) %>%   # No upper BMI exclusion
  mutate(
    TRTP_bin = ifelse(TRTP == "Baricitinib + RDV", 1, 0),
    SEX_M = ifelse(SEX == "M", 1L, 0L),
    ORDSCRG_num = as.integer(str_extract(ORDSCRG, "\\d+")),
    BMI_c = BMI - mean(BMI),
    int_linear = TRTP_bin * BMI_c,
    sev5 = as.integer(ORDSCRG_num == 5),
    sev6 = as.integer(ORDSCRG_num == 6),
    sev7 = as.integer(ORDSCRG_num == 7),
    ftime = case_when(
      RECCNSR == 0 & DTHCNSR == 0 ~ pmin(TTRECOV, TTDEATH, na.rm = TRUE),
      RECCNSR == 0 & DTHCNSR == 1 ~ TTRECOV,
      RECCNSR == 1 & DTHCNSR == 0 ~ TTDEATH,
      RECCNSR == 1 & DTHCNSR == 1 ~ ifelse(!is.na(TTRECOV), TTRECOV, 29),
      TRUE ~ NA_real_
    ),
    fstatus = case_when(
      DTHCNSR == 0 & (RECCNSR == 1 | TTDEATH <= TTRECOV) ~ 2L,
      RECCNSR == 0 & (DTHCNSR == 1 | TTRECOV  < TTDEATH) ~ 1L,
      RECCNSR == 1 & DTHCNSR == 1 ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(ftime > 0 | is.na(ftime))

mm_exp_full <- as.matrix(data_expanded[, c("TRTP_bin", "BMI_c", "int_linear",
                                            "AGE", "SEX_M",
                                            "sev5", "sev6", "sev7")])
mm_exp_reduced <- as.matrix(data_expanded[, c("TRTP_bin", "BMI_c",
                                                "AGE", "SEX_M",
                                                "sev5", "sev6", "sev7")])

fg_exp_full <- crr(ftime = data_expanded$ftime, fstatus = data_expanded$fstatus,
                   cov1 = mm_exp_full, failcode = 1, cencode = 0)
fg_exp_reduced <- crr(ftime = data_expanded$ftime, fstatus = data_expanded$fstatus,
                      cov1 = mm_exp_reduced, failcode = 1, cencode = 0)

lrt_exp <- lrt_crr(fg_exp_full, fg_exp_reduced, df = 1)

cat("  n:", nrow(data_expanded),
    " chi-square:", round(lrt_exp$chi, 2),
    " p:", signif(lrt_exp$p, 3), "\n\n")

################################################################################
# SENSITIVITY 3: NATURAL CUBIC SPLINE BMI
################################################################################

cat("Sensitivity 3: Natural cubic spline BMI\n")

# 3-df natural cubic spline (default knot placement)
bmi_spline <- ns(data_model$BMI_c, df = 3)
spline_df <- as.data.frame(bmi_spline)
colnames(spline_df) <- paste0("spl", 1:3)

data_model_spl <- bind_cols(data_model, spline_df) %>%
  mutate(
    int_spl1 = TRTP_bin * spl1,
    int_spl2 = TRTP_bin * spl2,
    int_spl3 = TRTP_bin * spl3
  )

mm_spl_full <- as.matrix(data_model_spl[, c("TRTP_bin",
                                              "spl1", "spl2", "spl3",
                                              "int_spl1", "int_spl2", "int_spl3",
                                              "AGE", "SEX_M",
                                              "sev5", "sev6", "sev7")])
mm_spl_reduced <- as.matrix(data_model_spl[, c("TRTP_bin",
                                                 "spl1", "spl2", "spl3",
                                                 "AGE", "SEX_M",
                                                 "sev5", "sev6", "sev7")])

fg_spl_full <- crr(ftime = data_model$ftime, fstatus = data_model$fstatus,
                   cov1 = mm_spl_full, failcode = 1, cencode = 0)
fg_spl_reduced <- crr(ftime = data_model$ftime, fstatus = data_model$fstatus,
                      cov1 = mm_spl_reduced, failcode = 1, cencode = 0)

lrt_spl <- lrt_crr(fg_spl_full, fg_spl_reduced, df = 3)

cat("  n:", nrow(data_model),
    " chi-square:", round(lrt_spl$chi, 2),
    " p:", signif(lrt_spl$p, 3), "\n\n")

################################################################################
# COMPILE TABLE S.5
################################################################################

primary <- readRDS("04_Results/primary_model.rds")

table_s5 <- tibble(
  Sensitivity_Analysis = c(
    "Primary model (linear BMI, factor severity)",
    "Race and ethnicity adjustment",
    "Inclusion of BMI >=50 kg/m^2",
    "Natural cubic spline BMI"
  ),
  n = c(942, nrow(data_race), nrow(data_expanded), 942),
  chi_square = c(round(primary$lrt_chi, 2),
                 round(lrt_race$chi, 2),
                 round(lrt_exp$chi, 2),
                 round(lrt_spl$chi, 2)),
  df = c(1, 1, 1, 3),
  P_value = c(signif(primary$lrt_p, 3),
              signif(lrt_race$p, 3),
              signif(lrt_exp$p, 3),
              signif(lrt_spl$p, 3)),
  Notes = c(
    "Reference for comparison",
    "Additional adjustment for self-reported race and ethnicity",
    "Expanded cohort",
    "3-df natural cubic spline"
  )
)

cat("Table S.5: Sensitivity Analyses\n")
print(table_s5)

write_csv(table_s5, "04_Results/tables/table_s5_sensitivity_analyses.csv")
cat("\nTable S.5 saved to 04_Results/tables/table_s5_sensitivity_analyses.csv\n")
