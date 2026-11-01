Survey-Aware Machine Learning for Rare-Outcome Prediction: NHANES Application

This repository contains the R code and results for the preliminary analysis in my PhD proposal, Survey-Aware Machine Learning for Rare-Outcome Prediction under Complex Sampling Designs.

The question

Health surveys such as NHANES select people with unequal probabilities, in clusters, and many outcomes of interest are rare. Two kinds of weights are used in training a prediction model:

Survey weights say how many people in the population each respondent represents.
Class weights (λ) make the model pay more attention to the rare outcome.

They do different jobs, so it is not obvious how they should be combined. This analysis asks: does adding a class weight λ to a survey-weighted model improve detection of a rare outcome, compared with simply keeping the survey-weighted model and tuning its classification threshold?

Data and outcome
Data: NHANES 2015–2016 and 2017–2018 (US National Center for Health Statistics, CDC). The script downloads the public files itself (demographics, body measures, blood pressure, glycohemoglobin, diabetes questionnaire).
Sample: adults aged 20 and over with an HbA1c measurement and no reported diabetes diagnosis (n = 8,621).
Outcome: undiagnosed diabetes, defined as HbA1c ≥ 6.5% with no reported diagnosis (341 cases; weighted prevalence 2.5%).
Predictors: age, sex, race/ethnicity, BMI, waist circumference, systolic and diastolic blood pressure, income-to-poverty ratio. Missing predictor values are filled with the median.
Weights: the two-year MEC examination weight (WTMEC2YR) divided by 2 to give a combined four-year weight, following NHANES guidance. Strata and PSU variables (SDMVSTRA, SDMVPSU) are kept for PSU-grouped folds and the PSU bootstrap. The weight-based Kish effective sample size is about 4,173; it reflects weight variability only, not the full effect of clustering and stratification.
Methods compared
Label	Model
B	Survey-weighted logistic regression (λ = 1)
C, λ = 2, 5, 10, 20	Survey-weighted logistic regression in which each positive case's weight is multiplied by λ
B + tuned threshold	B with its classification threshold chosen to give a false-positive rate of about 5% (the key comparator; here, "tuned threshold" columns)

Evaluation. Ten repeats of five-fold cross-validation, where all people from the same primary sampling unit (PSU) are kept in the same fold, so the models are always tested on clusters they have not seen. All measures are survey-weighted and computed on held-out data. Tuned thresholds are chosen on the training folds only. A PSU-bootstrap (resampling PSUs within strata, 500 draws) gives the interval for the difference between λ = 20 and B.

Results
Model	Sensitivity at t = 0.5	False-positive rate at t = 0.5	Sensitivity, tuned threshold	Calibration intercept	Calibration slope
B (λ = 1)	0.000	0.000	0.229	−0.26	0.92
C, λ = 2	0.005	0.001	0.227	−0.95	0.90
C, λ = 5	0.085	0.011	0.225	−1.84	0.86
C, λ = 10	0.222	0.051	0.226	−2.49	0.83
C, λ = 20	0.470	0.142	0.231	−3.09	0.78

AUROC was about 0.79 for all models. With the tuned threshold, the false-positive rate was 5.1–5.2% for every model.

Reading the table

At a fixed threshold (0.5), a larger λ raises sensitivity, but only because it also raises the false-positive rate. The model's ranking of people has not improved; its operating point has moved.
Calibration: the calibration intercept (from regressing the observed outcome on the logit of the predicted probability) falls by 2.82 between B and λ = 20, close to −log(20) = −3.00. This matches a simple analytic result: in a correctly specified logistic model, class weighting leaves the slopes unchanged and shifts the model's fitted intercept by log(λ), which inflates the predicted probabilities. The calibration intercept is a different quantity: it measures the agreement between predicted and observed outcomes, and it reflects the fitted-intercept shift indirectly. The calibration slope also worsens, so class-weighted probabilities need recalibration before they are used as risks.
With the threshold tuned to the same false-positive rate, sensitivity is essentially the same for all λ (0.225–0.231). The average difference between λ = 20 and B is +0.002 (95% PSU-bootstrap interval −0.008 to 0.012), which is indistinguishable from zero.

Conclusion. In this setting, class weighting adds nothing beyond a survey-weighted model with a tuned threshold. The research question for the PhD is where this equivalence breaks (model misspecification, penalized or flexible learners, very few positive cases, outcome-related selection).

Files
File	What it is
nhanes_analysis.R	Complete NHANES analysis: data download, cleaning, models, cross-validation, bootstrap, table, figure
results/nhanes_main_table.csv	Table of results above (unrounded)
results/nhanes_paired_difference.csv	Difference in tuned-threshold sensitivity, λ = 20 minus B, with 95% interval
results/nhanes_figure.png	Sensitivity, false-positive rate, and calibration against λ
results/sessionInfo.txt	R version and platform used
pilot.R, Pilot_Simulation_Prevalence.Rmd	Earlier simulation with a known population (prevalence 3%, 4%, 5%) showing the same pattern
How to run
Install R (base R is enough; the foreign package ships with R). An internet connection is needed the first time.
Download this repository and open nhanes_analysis.R in RStudio.
Set the working directory to the folder that holds the script (Session → Set Working Directory → To Source File Location).
Click Source. The analysis takes a few minutes. Downloaded data are saved in a data folder, and outputs are written to results/.

The random seed is fixed, so results should reproduce exactly.

Limitations
One outcome definition (HbA1c only; fasting glucose is not used) and one set of predictors.
A linear logistic model and simple median imputation.
NHANES has no known truth, so performance is estimated by design-weighted cross-validation. The simulation, where the whole population is known, provides the exact comparison.
This is preliminary work to support a research proposal, not a finished study, and it is not clinical guidance.
Data source

National Center for Health Statistics. National Health and Nutrition Examination Survey (NHANES). Centers for Disease Control and Prevention. https://www.cdc.gov/nhanes/

Author

Erum Shahzadi, Department of Statistics, University of Gujrat, Pakistan.

National Center for Health Statistics. National Health and Nutrition Examination Survey (NHANES). Centers for Disease Control and Prevention. https://www.cdc.gov/nhanes/

Author

Erum Shahzadi, Department of Statistics, University of Gujrat, Pakistan.
