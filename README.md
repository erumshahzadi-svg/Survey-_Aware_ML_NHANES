# Survey-aware learning for a rare outcome: NHANES application

Compares survey-weighted logistic regression (B) with joint survey- and class-weighted logistic
regression (C, class weight lambda = 2, 5, 10, 20) for **undiagnosed diabetes** in US adults
(NHANES 2015-2016 and 2017-2018). Undiagnosed diabetes = HbA1c >= 6.5% with no reported diagnosis.

## Run
1. Install R (base R is enough; the `foreign` package ships with R).
2. In R: `source("nhanes_analysis.R")`. Data download from the CDC site the first time (internet needed).
3. Outputs appear in `results/`: `nhanes_main_table.csv`, `nhanes_paired_difference.csv`, `nhanes_figure.png`.

## What it does
- Keeps adults 20+ with HbA1c and no diagnosis; uses MEC exam weights (halved for two cycles).
- Predictors: age, sex, race/ethnicity, BMI, waist, blood pressure, income ratio (median imputation).
- 10 x 5-fold cross-validation, keeping whole PSUs in one fold.
- Survey-weighted sensitivity, false-positive rate, AUROC and calibration on held-out folds.
- Tuned threshold (FPR about 5%) chosen on training folds only (the "B + tuned threshold" comparator).
- PSU bootstrap interval for the difference between lambda = 20 and B.

## Limits
Single outcome definition, simple imputation, linear logistic model. NHANES has no known truth,
so performance is estimated by design-weighted cross-validation, not computed exactly as in the simulation.
Related simulation code: `pilot.R`, `Pilot_Simulation_Prevalence.Rmd`.
