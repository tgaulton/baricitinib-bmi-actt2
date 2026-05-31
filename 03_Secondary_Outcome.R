################################################################################
# 03_Secondary_Outcome.R
#
# Purpose: Secondary outcome analysis for progression to mechanical ventilation,
#          ECMO, or death by day 15.
#
# Specification:
#   - Logistic regression
#   - Cohort: participants not on invasive MV or ECMO at baseline (n=842)
#   - Linear BMI with BMI x treatment interaction
#   - Adjusted for age, sex, baseline severity (factor)
#   - ORs at BMI 20, 25, 30, 35, 40, 45 via delta method
#
# Produces Figure 3 of the manuscript.
#
# NOTE: Requires 01_Analytic_Sample.R to have been run first.
################################################################################

library(tidyverse)
library(ggplot2)

cat("================================================================================\n")
cat("SECONDARY OUTCOME ANALYSIS: BMI x TREATMENT INTERACTION\n")
cat("================================================================================\n\n")

################################################################################
# LOAD DATA AND PREPARE
################################################################################

data_final <- readRDS("02_Data/processed/analytic_dataset.rds")

# Restrict to participants not requiring invasive MV or ECMO at baseline
data_secondary <- data_final %>%
  filter(ORDSCRG_num <= 6) %>%
  mutate(
    ORDSCR15_num = as.numeric(ORDSCR15),
    composite_outcome = as.integer(ORDSCR15_num >= 7),
    BMI_c = BMI_kgm2 - mean(BMI_kgm2),
    int_linear = TRTP_bin * BMI_c
  )

bmi_mean_secondary <- mean(data_secondary$BMI_kgm2)

cat("Secondary outcome cohort: n =", nrow(data_secondary), "\n")
cat("Composite events:", sum(data_secondary$composite_outcome), "\n\n")

# Event rates by arm
events_by_arm <- table(data_secondary$TRTP_bin, data_secondary$composite_outcome)
cat("Event rates by treatment arm:\n")
cat("  Placebo:", events_by_arm[1, 2],
    "of", sum(events_by_arm[1, ]),
    sprintf("(%.1f%%)\n", 100 * events_by_arm[1, 2] / sum(events_by_arm[1, ])))
cat("  Baricitinib:", events_by_arm[2, 2],
    "of", sum(events_by_arm[2, ]),
    sprintf("(%.1f%%)\n\n", 100 * events_by_arm[2, 2] / sum(events_by_arm[2, ])))

################################################################################
# FIT MODELS (full and reduced)
################################################################################

# In the secondary cohort, only sev5 and sev6 are present (ORDSCRG_num <= 6)
logit_full <- glm(
  composite_outcome ~ TRTP_bin + BMI_c + int_linear +
    AGE + SEX_M + sev5 + sev6,
  data = data_secondary,
  family = binomial(link = "logit")
)

logit_reduced <- glm(
  composite_outcome ~ TRTP_bin + BMI_c +
    AGE + SEX_M + sev5 + sev6,
  data = data_secondary,
  family = binomial(link = "logit")
)

# LRT
lrt_anova <- anova(logit_reduced, logit_full, test = "Chisq")
lrt_chi <- round(lrt_anova$Deviance[2], 3)
lrt_p   <- signif(lrt_anova$`Pr(>Chi)`[2], 3)

cat("BMI x Treatment Interaction (Likelihood Ratio Test):\n")
cat("  chi-square:", lrt_chi, "\n")
cat("  df: 1\n")
cat("  p-value:", lrt_p, "\n\n")

################################################################################
# COMPUTE OR AT EACH BMI USING DELTA METHOD
################################################################################

compute_or_at_bmi <- function(bmi_value, fit, bmi_mean) {
  bmi_c_val <- bmi_value - bmi_mean
  coef_fit  <- as.numeric(coef(fit))
  vcov_fit  <- vcov(fit)
  
  # Model coefficient order:
  # 1=(Intercept), 2=TRTP_bin, 3=BMI_c, 4=int_linear, 5=AGE, 6=SEX_M, 7=sev5, 8=sev6
  trtp_idx <- 2
  int_idx  <- 4
  
  log_or <- coef_fit[trtp_idx] + coef_fit[int_idx] * bmi_c_val
  
  var_log_or <- vcov_fit[trtp_idx, trtp_idx] +
                bmi_c_val^2 * vcov_fit[int_idx, int_idx] +
                2 * bmi_c_val * vcov_fit[trtp_idx, int_idx]
  
  se_log_or <- sqrt(var_log_or)
  
  c(exp(log_or),
    exp(log_or - 1.96 * se_log_or),
    exp(log_or + 1.96 * se_log_or))
}

# Continuous curve
bmi_grid <- seq(min(data_secondary$BMI_kgm2), max(data_secondary$BMI_kgm2), by = 0.5)
or_curve <- tibble(bmi = bmi_grid, or = NA_real_, lower = NA_real_, upper = NA_real_)

for (i in seq_along(bmi_grid)) {
  result <- compute_or_at_bmi(bmi_grid[i], logit_full, bmi_mean_secondary)
  or_curve$or[i]    <- result[1]
  or_curve$lower[i] <- result[2]
  or_curve$upper[i] <- result[3]
}

# Point estimates
bmi_points <- c(20, 25, 30, 35, 40, 45)
or_points <- tibble(bmi = bmi_points, or = NA_real_, lower = NA_real_, upper = NA_real_)

for (i in seq_along(bmi_points)) {
  result <- compute_or_at_bmi(bmi_points[i], logit_full, bmi_mean_secondary)
  or_points$or[i]    <- result[1]
  or_points$lower[i] <- result[2]
  or_points$upper[i] <- result[3]
}

cat("OR at clinically relevant BMI values:\n")
print(or_points %>% mutate(across(c(or, lower, upper), ~ round(.x, 2))))
cat("\n")

################################################################################
# FIGURE 3: OR ACROSS BMI (BLACK AND WHITE)
################################################################################

p_figure3 <- ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "gray40", linewidth = 0.7) +
  geom_ribbon(data = or_curve,
              aes(x = bmi, ymin = lower, ymax = upper),
              fill = "gray80", alpha = 0.5) +
  geom_line(data = or_curve,
            aes(x = bmi, y = or),
            color = "black", linewidth = 1.0) +
  geom_errorbar(data = or_points,
                aes(x = bmi, ymin = lower, ymax = upper),
                width = 0.5, linewidth = 0.8, color = "black") +
  geom_point(data = or_points,
             aes(x = bmi, y = or),
             size = 3, shape = 21, fill = "white",
             color = "black", stroke = 1.2) +
  scale_x_continuous(breaks = c(20, 25, 30, 35, 40, 45), limits = c(18, 47)) +
  scale_y_continuous(breaks = seq(0, 2.0, 0.25), limits = c(0, 2.0)) +
  labs(x = "Body Mass Index (kg/m\u00B2)",
       y = "Odds Ratio") +
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

write_csv(or_points, "04_Results/tables/figure3_or_by_bmi_points.csv")

ggsave("04_Results/figures/figure3_or_by_bmi.tiff",
       p_figure3, width = 7, height = 5, dpi = 500,
       device = "tiff", compression = "lzw")

cat("Figure 3 saved to 04_Results/figures/figure3_or_by_bmi.tiff\n")
cat("OR points saved to 04_Results/tables/figure3_or_by_bmi_points.csv\n")
