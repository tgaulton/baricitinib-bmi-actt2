################################################################################
# 01_Analytic_Sample.R
#
# Purpose: Create analytic dataset from ACTT-2 raw data
# Final N: 942
#
# Exclusions (applied sequentially):
#   - Missing BMI: 43
#   - BMI <18.5 kg/m^2: 8
#   - BMI >=50 kg/m^2: 35
#   - Invalid follow-up time (ftime = 0): 5
#   Total excluded: 91 of 1,033 randomized
#
# NOTE: Set working directory to the project root before running.
################################################################################

library(tidyverse)

cat("================================================================================\n")
cat("CREATING ACTT-2 ANALYTIC DATASET\n")
cat("================================================================================\n\n")

################################################################################
# CONFIGURATION
################################################################################

BMI_LOWER <- 18.5
BMI_UPPER <- 50

################################################################################
# LOAD RAW DATA
################################################################################

if (!file.exists("02_Data/raw/ACTT2.csv")) {
  stop("ACTT2.csv not found in 02_Data/raw/")
}

actt2_raw <- read_csv("02_Data/raw/ACTT2.csv", show_col_types = FALSE)
cat("Raw data loaded:", nrow(actt2_raw), "participants\n\n")

################################################################################
# APPLY EXCLUSIONS
################################################################################

exclusions <- tibble(
  step = 0L,
  reason = "Randomized in ACTT-2",
  n_excluded = 0L,
  n_remaining = nrow(actt2_raw)
)

add_exclusion <- function(df, data_before, data_after, reason) {
  n_before <- nrow(data_before)
  n_after <- nrow(data_after)
  n_excl <- n_before - n_after
  
  cat(sprintf("  %s: excluded %d (%.1f%%), remaining %d\n",
              reason, n_excl, 100 * n_excl / 1033, n_after))
  
  bind_rows(df, tibble(
    step = nrow(df),
    reason = reason,
    n_excluded = n_excl,
    n_remaining = n_after
  ))
}

cat("Applying exclusion criteria:\n")

data_step1 <- actt2_raw %>% filter(!is.na(BMI))
exclusions <- add_exclusion(exclusions, actt2_raw, data_step1, "Missing BMI")

data_step2 <- data_step1 %>% filter(BMI >= BMI_LOWER)
exclusions <- add_exclusion(exclusions, data_step1, data_step2,
                            sprintf("BMI < %.1f kg/m2", BMI_LOWER))

data_step3 <- data_step2 %>% filter(BMI < BMI_UPPER)
exclusions <- add_exclusion(exclusions, data_step2, data_step3,
                            sprintf("BMI >= %.0f kg/m2", BMI_UPPER))

################################################################################
# CREATE DERIVED VARIABLES
################################################################################

data_step4 <- data_step3 %>%
  mutate(
    # Treatment
    TRTP_bin = case_when(
      TRTP == "Baricitinib + RDV" ~ 1,
      TRTP == "Placebo + RDV"     ~ 0,
      TRUE ~ NA_real_
    ),
    TRTP_label = factor(TRTP_bin, levels = c(0, 1),
                        labels = c("Placebo + RDV", "Baricitinib + RDV")),
    
    # Demographics
    SEX_M = case_when(
      SEX == "M" ~ 1L,
      SEX == "F" ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    # BMI variables
    BMI_kgm2 = BMI,
    BMI_cat = cut(
      BMI,
      breaks = c(BMI_LOWER, 25, 30, 35, 40, BMI_UPPER),
      labels = c("Normal Weight (18.5-24.9)",
                 "Overweight (25-29.9)",
                 "Obesity Class I (30-34.9)",
                 "Obesity Class II (35-39.9)",
                 "Obesity Class III (40-49.9)"),
      right = FALSE,
      include.lowest = TRUE
    ),
    Obese = BMI >= 30,
    
    # Baseline severity (ordinal score 4-7)
    ORDSCRG_num = as.integer(str_extract(ORDSCRG, "\\d+")),
    
    # Severity dummies for factor specification (primary)
    sev5 = as.integer(ORDSCRG_num == 5),
    sev6 = as.integer(ORDSCRG_num == 6),
    sev7 = as.integer(ORDSCRG_num == 7),
    
    # Comorbidities
    Diabetes = case_when(
      is.na(DIAB2FL) ~ NA,
      DIAB2FL == "Y" ~ TRUE,
      TRUE ~ FALSE
    ),
    Hypertension = case_when(
      is.na(HYPFL) ~ NA,
      HYPFL == "Y" ~ TRUE,
      TRUE ~ FALSE
    ),
    Cardiometabolic = case_when(
      is.na(Diabetes) | is.na(Hypertension) ~ NA,
      Diabetes | Hypertension ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # Metabolic phenotype groups
    Metabolic_Phenotype = case_when(
      is.na(Cardiometabolic) ~ NA_character_,
      !Obese & !Cardiometabolic ~ "MHNW",
      !Obese &  Cardiometabolic ~ "MUNW",
       Obese & !Cardiometabolic ~ "MHO",
       Obese &  Cardiometabolic ~ "MUO"
    ),
    Phenotype_factor = factor(
      Metabolic_Phenotype,
      levels = c("MHNW", "MUNW", "MHO", "MUO")
    ),
    
    # Competing risk outcome variables
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
    ),
    fstatus_label = factor(fstatus, levels = c(0, 1, 2),
                           labels = c("Censored", "Recovered", "Died"))
  )

################################################################################
# FINAL EXCLUSION: Invalid follow-up time
################################################################################

data_final <- data_step4 %>% filter(ftime > 0 | is.na(ftime))
exclusions <- add_exclusion(exclusions, data_step4, data_final,
                            "Invalid follow-up time")

cat("\nFinal analytic cohort: n =", nrow(data_final), "\n\n")

if (nrow(data_final) != 942) {
  warning("Expected n=942, got n=", nrow(data_final))
}

################################################################################
# QUALITY CHECKS
################################################################################

cat("Treatment distribution:\n")
print(table(data_final$TRTP_label))
cat("\nOutcome distribution:\n")
print(table(data_final$fstatus_label))
cat("\nBMI category distribution:\n")
print(table(data_final$BMI_cat))
cat("\nMetabolic phenotype distribution:\n")
print(table(data_final$Phenotype_factor, useNA = "ifany"))

################################################################################
# SAVE
################################################################################

dir.create("02_Data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Results/tables", recursive = TRUE, showWarnings = FALSE)

saveRDS(data_final, "02_Data/processed/analytic_dataset.rds")
write_csv(data_final, "02_Data/processed/analytic_dataset.csv")
write_csv(exclusions, "04_Results/tables/exclusion_flowchart.csv")

cat("\nAnalytic dataset saved to 02_Data/processed/analytic_dataset.rds\n")
cat("Exclusion flowchart saved to 04_Results/tables/exclusion_flowchart.csv\n")
