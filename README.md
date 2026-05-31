# Body Mass Index and Baricitinib Treatment Effect in Hospitalized Adults With COVID-19

Analysis code for a secondary analysis of the Adaptive COVID-19 Treatment Trial 2 (ACTT-2) examining whether body mass index (BMI) and metabolic phenotypes modify the treatment effect of baricitinib plus remdesivir.

## Repository structure

```
.
├── 01_Analytic_Sample.R          Create analytic cohort from raw data
├── 02_Primary_Analysis.R         Primary HR by BMI (linear, Fine-Gray)
├── 03_Secondary_Outcome.R        OR by BMI for composite secondary outcome
├── 04_Metabolic_Phenotypes.R     HR by metabolic subgroup
├── 05_Sensitivity_Analyses.R     All pre-specified sensitivity analyses
├── 06_Model_Adequacy.R           PH testing, AIC, influence diagnostics
└── README.md
```

## Data access

Individual-level data from ACTT-2 are available from the National Institute of Allergy and Infectious Diseases Clinical Trials Data Repository (https://accessclinicaldata.niaid.nih.gov) upon execution of a data use agreement. Place the raw ACTT2.csv file in `02_Data/raw/` before running the scripts.

## Software requirements

R version 4.5.1 or later. Required packages:

```r
install.packages(c("tidyverse", "cmprsk", "splines", "survival"))
```

## Running the analyses

Scripts are designed to run sequentially. Set the working directory to the project root before running:

```r
setwd("/path/to/project_root")
source("01_Analytic_Sample.R")     # Creates analytic dataset (n=942)
source("02_Primary_Analysis.R")    # Primary HR analysis
source("03_Secondary_Outcome.R")   # Composite outcome analysis
source("04_Metabolic_Phenotypes.R") # Metabolic subgroup analysis
source("05_Sensitivity_Analyses.R") # Sensitivity analyses
source("06_Model_Adequacy.R")      # Model adequacy diagnostics
```

## Primary specification

Per the revised manuscript:

- Fine-Gray competing risk regression for time to recovery, with death as a competing event
- Linear BMI with BMI × treatment interaction
- Adjusted for age (continuous), sex (male vs female), baseline disease severity (four-level factor with score 4 as reference)
- Heterogeneity tested by likelihood ratio test
- HRs and 95% CIs at BMI values of 20, 25, 30, 35, 40, and 45 kg/m² computed using the delta method

## Citation

If you use this code, please cite:

[Citation to be added upon publication]

## License

MIT License. See LICENSE file for details.

## Contact

Tim Gaulton, MD, MSc
Massachusetts General Hospital, Harvard Medical School
