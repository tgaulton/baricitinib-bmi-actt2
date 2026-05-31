################################################################################
# 02_Primary_Analysis.R
#
# Purpose: Primary analysis of treatment heterogeneity by BMI.
#
# Specification (per revised manuscript):
#   - Fine-Gray competing risk regression
#   - Time to recovery, with death as a competing event
#   - Linear BMI with BMI x treatment interaction
#   - Adjusted for age (continuous), sex, baseline severity (four-level factor)
#   - Heterogeneity tested by likelihood ratio test
#   - HRs at BMI 20, 25, 30, 35, 40, 45 kg/m^2 via delta method
#
# Produces Figure 1 and Table 2 of the manuscript.
#
# NOTE: Requires 01_Analytic_Sample.R to have been run first.
################################################################################

library(tidyverse)
library(cmprsk)
library(ggplot2)

cat("================================================================================\n")
cat("PRIMARY ANALYSIS: BMI x TREATMENT INTERACTION\n")
cat("================================================================================\n\n")

################################################################################
# LOAD DATA AND PREPARE
################################################################################

data_final <- readRDS("02_Data/processed/analytic_dataset.rds")

bmi_mean <- mean(data_final$BMI_kgm2)

data_model <- data_final %>%
  mutate(
    BMI_c = BMI_kgm2 - bmi_mean,
    int_linear = TRTP_bin * BMI_c
  )

cat("Analytic cohort: n =", nrow(data_model), "\n")
cat("Mean BMI:", round(bmi_mean, 2), "kg/m^2\n\n")

################################################################################
# FIT PRIMARY MODEL (Fine-Gray with linear BMI x treatment interaction)
################################################################################

# Full model with interaction
mm_full <- as.matrix(data_model[, c("TRTP_bin", "BMI_c", "int_linear",
                                     "AGE", "SEX_M",
                                     "sev5", "sev6", "sev7")])

fg_primary <- crr(
  ftime    = data_model$ftime,
  fstatus  = data_model$fstatus,
  cov1     = mm_full,
  failcode = 1,
  cencode  = 0
)

# Reduced model without interaction
mm_reduced <- as.matrix(data_model[, c("TRTP_bin", "BMI_c",
                                        "AGE", "SEX_M",
                                        "sev5", "sev6", "sev7")])

fg_reduced <- crr(
  ftime    = data_model$ftime,
  fstatus  = data_model$fstatus,
  cov1     = mm_reduced,
  failcode = 1,
  cencode  = 0
)

# LRT for interaction
lrt_chi <- 2 * (fg_primary$loglik - fg_reduced$loglik)
lrt_p   <- pchisq(lrt_chi, df = 1, lower.tail = FALSE)

cat("BMI x Treatment Interaction (Likelihood Ratio Test):\n")
cat("  chi-square:", round(lrt_chi, 3), "\n")
cat("  df: 1\n")
cat("  p-value:", signif(lrt_p, 3), "\n\n")

################################################################################
# COMPUTE HR AT EACH BMI USING DELTA METHOD
################################################################################

compute_hr_at_bmi <- function(bmi_value, fit, bmi_mean) {
  bmi_c_val <- bmi_value - bmi_mean
  coef_fit  <- as.numeric(fit$coef)
  vcov_fit  <- fit$var
  
  # Positions in model matrix: TRTP_bin = 1, BMI_c = 2, int_linear = 3
  trtp_idx <- 1
  int_idx  <- 3
  
  log_hr <- coef_fit[trtp_idx] + coef_fit[int_idx] * bmi_c_val
  
  var_log_hr <- vcov_fit[trtp_idx, trtp_idx] +
                bmi_c_val^2 * vcov_fit[int_idx, int_idx] +
                2 * bmi_c_val * vcov_fit[trtp_idx, int_idx]
  
  se_log_hr <- sqrt(var_log_hr)
  
  c(exp(log_hr),
    exp(log_hr - 1.96 * se_log_hr),
    exp(log_hr + 1.96 * se_log_hr))
}

# Continuous curve for figure
bmi_grid <- seq(18.5, 50, by = 0.5)
hr_curve <- tibble(bmi = bmi_grid, hr = NA_real_, lower = NA_real_, upper = NA_real_)

for (i in seq_along(bmi_grid)) {
  result <- compute_hr_at_bmi(bmi_grid[i], fg_primary, bmi_mean)
  hr_curve$hr[i]    <- result[1]
  hr_curve$lower[i] <- result[2]
  hr_curve$upper[i] <- result[3]
}

# Point estimates at clinically relevant BMI values
bmi_points <- c(20, 25, 30, 35, 40, 45)
hr_points <- tibble(bmi = bmi_points, hr = NA_real_, lower = NA_real_, upper = NA_real_)

for (i in seq_along(bmi_points)) {
  result <- compute_hr_at_bmi(bmi_points[i], fg_primary, bmi_mean)
  hr_points$hr[i]    <- result[1]
  hr_points$lower[i] <- result[2]
  hr_points$upper[i] <- result[3]
}

cat("HR at clinically relevant BMI values (Table 2):\n")
print(hr_points %>% mutate(across(c(hr, lower, upper), ~ round(.x, 2))))
cat("\n")

################################################################################
# FIGURE 1: HR ACROSS BMI (BLACK AND WHITE)
################################################################################

p_figure1 <- ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "gray40", linewidth = 0.7) +
  geom_ribbon(data = hr_curve,
              aes(x = bmi, ymin = lower, ymax = upper),
              fill = "gray80", alpha = 0.5) +
  geom_line(data = hr_curve,
            aes(x = bmi, y = hr),
            color = "black", linewidth = 1.0) +
  geom_errorbar(data = hr_points,
                aes(x = bmi, ymin = lower, ymax = upper),
                width = 0.5, linewidth = 0.8, color = "black") +
  geom_point(data = hr_points,
             aes(x = bmi, y = hr),
             size = 3, shape = 21, fill = "white",
             color = "black", stroke = 1.2) +
  scale_x_continuous(breaks = c(20, 25, 30, 35, 40, 45), limits = c(18, 47)) +
  scale_y_continuous(breaks = seq(0.5, 2.0, 0.25), limits = c(0.5, 2.0)) +
  labs(x = "Body Mass Index (kg/m\u00B2)",
       y = "Hazard Ratio") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold",
                                 margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

################################################################################
# SAVE OUTPUTS
################################################################################

dir.create("04_Results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Results/tables",  recursive = TRUE, showWarnings = FALSE)

# Table 2
write_csv(hr_points, "04_Results/tables/table2_hr_by_bmi.csv")

# Figure 1
ggsave("04_Results/figures/figure1_hr_by_bmi.tiff",
       p_figure1, width = 7, height = 5, dpi = 500,
       device = "tiff", compression = "lzw")

cat("Figure 1 saved to 04_Results/figures/figure1_hr_by_bmi.tiff\n")
cat("Table 2 saved to 04_Results/tables/table2_hr_by_bmi.csv\n")

# Save model object for downstream scripts
saveRDS(list(fit = fg_primary, bmi_mean = bmi_mean, lrt_chi = lrt_chi, lrt_p = lrt_p),
        "04_Results/primary_model.rds")
