# Survey-Aware Machine Learning for Rare-Outcome Prediction: NHANES Application

This repository holds the R code and results for the preliminary analysis in my PhD proposal,
*Survey-Aware Machine Learning for Rare-Outcome Prediction under Complex Sampling Designs: A Multistage Approach*.

## The question

Health surveys such as NHANES select people with unequal probabilities and in clusters. Many outcomes
of interest are also rare. Two kinds of weights can be used when training a prediction model:

- **Survey weights** say how many people in the population each respondent represents.
- **Class weights** (λ) make the model pay more attention to the rare outcome.

They do different jobs, so it is not obvious how they should be combined. This analysis asks one question:
**does adding a class weight λ to a survey-weighted model detect more rare cases than simply keeping the
survey-weighted model and tuning its classification threshold?**

## NHANES data files

The script downloads the public files from the CDC website on first run and saves them in a `data/` folder.
No data files are stored in this repository.

| NHANES component | 2015–2016 file | 2017–2018 file | Variables used |
|---|---|---|---|
| Demographics and design | `DEMO_I.xpt` | `DEMO_J.xpt` | `SEQN`, `RIAGENDR` (sex), `RIDAGEYR` (age), `RIDRETH3` (race/ethnicity), `INDFMPIR` (income-to-poverty ratio), `SDMVSTRA` (strata), `SDMVPSU` (PSU), `WTMEC2YR` (MEC exam weight) |
| Body measures | `BMX_I.xpt` | `BMX_J.xpt` | `BMXBMI` (BMI), `BMXWAIST` (waist circumference) |
| Blood pressure | `BPX_I.xpt` | `BPX_J.xpt` | `BPXSY1` (systolic), `BPXDI1` (diastolic) |
| Glycohemoglobin | `GHB_I.xpt` | `GHB_J.xpt` | `LBXGH` (HbA1c, %) |
| Diabetes questionnaire | `DIQ_I.xpt` | `DIQ_J.xpt` | `DIQ010` (ever told you have diabetes) |

Files are merged within each cycle by `SEQN`, and the two cycles are then stacked.

## Sample, outcome and weights

- **Sample:** adults aged 20 and over with an HbA1c result and a positive MEC weight, who did not report a
  diabetes diagnosis (`DIQ010` = 2 "No" or 3 "Borderline"). n = 8,621.
- **Outcome:** undiagnosed diabetes, defined as HbA1c ≥ 6.5% (341 cases; weighted prevalence 2.5%).
- **Predictors:** age, sex, race/ethnicity, BMI, waist circumference, systolic and diastolic blood pressure,
  and income-to-poverty ratio. Missing predictor values are filled with the median.
- **Weights:** following NHANES guidance, the two-year MEC weight (`WTMEC2YR`) is divided by 2 to give a
  combined four-year weight. Strata and PSU variables are kept for PSU-grouped folds and the PSU bootstrap.
- **Effective sample size:** the weight-based Kish effective sample size is about 4,173. It reflects weight
  variability only, not the full effect of clustering and stratification.

## Methods compared

| Label | Model |
|---|---|
| B | Survey-weighted logistic regression (λ = 1) |
| C, λ = 2, 5, 10, 20 | Survey-weighted logistic regression in which each positive case's weight is multiplied by λ |
| B + tuned threshold | B with its threshold chosen to give a false-positive rate of about 5% (the key comparator; the "tuned threshold" column) |

**Evaluation.** Ten repeats of five-fold cross-validation. All people from the same PSU stay in the same fold,
so models are always tested on clusters they have not seen. All measures are survey-weighted and computed on
held-out data. Tuned thresholds are chosen on the training folds only. A PSU bootstrap (PSUs resampled within
strata, 500 draws) gives the interval for the difference between λ = 20 and B.

## Results

| Model | Sensitivity at t = 0.5 | False-positive rate at t = 0.5 | Sensitivity, tuned threshold | Calibration intercept | Calibration slope | AUROC |
|---|---|---|---|---|---|---|
| B (λ = 1) | 0.000 | 0.000 | 0.229 | −0.26 | 0.92 | 0.786 |
| C, λ = 2 | 0.005 | 0.001 | 0.227 | −0.95 | 0.90 | 0.786 |
| C, λ = 5 | 0.085 | 0.011 | 0.225 | −1.84 | 0.86 | 0.787 |
| C, λ = 10 | 0.222 | 0.051 | 0.226 | −2.49 | 0.83 | 0.788 |
| C, λ = 20 | 0.470 | 0.142 | 0.231 | −3.09 | 0.78 | 0.788 |

With the tuned threshold, the false-positive rate was 5.1–5.2% for every model.

**Reading the table**

1. **At the default threshold (0.5),** a larger λ raises sensitivity. But the false-positive rate rises with it.
   The model's ranking of people has not improved. Only its operating point has moved.
2. **Ranking:** AUROC stays at about 0.79 for every λ. Class weighting does not help the model separate cases
   from non-cases.
3. **Calibration:** the calibration intercept falls by 2.82 between B and λ = 20, close to log(20) = 3.00.
   This is expected. In a correctly specified logistic model, class weighting mainly shifts the fitted intercept
   by about log λ, which inflates the predicted probabilities. The calibration slope also worsens, so
   class-weighted probabilities need recalibration before they are used as risks.
4. **With the threshold tuned to the same false-positive rate,** sensitivity is about the same for every λ
   (0.225–0.231). The difference between λ = 20 and B is **+0.002** (95% PSU-bootstrap interval −0.008 to 0.012).
   This is indistinguishable from zero.

**Conclusion.** In this setting, class weighting adds nothing beyond a survey-weighted model with a tuned
threshold. This is what theory predicts for a near-correctly specified logistic model. The PhD research asks
where this equivalence breaks: under model misspecification, penalized estimation, flexible learners such as
random forests and XGBoost, very few positive cases, or highly variable survey weights.

## Files

| File | What it is |
|---|---|
| `nhanes_analysis.R` | Complete NHANES analysis: data download, cleaning, models, cross-validation, bootstrap, table and figure |
| `nhanes_main_table.csv` | Results table above (unrounded) |
| `nhanes_paired_difference.csv` | Difference in tuned-threshold sensitivity, λ = 20 minus B, with 95% interval |
| `nhanes_figure.png` | Sensitivity, false-positive rate and calibration plotted against λ |
| `sessionInfo.txt` | R version and platform used |
| `NHANES_Analysis.Rmd` | Same NHANES analysis as a step-by-step R Markdown report |
| `make_nhanes_data.R` | Builds the analysis file `nhanes_2015_2018_analysis.csv` from the NHANES public files |
| `pilot.R`, `Pilot_Simulation_Prevalence.Rmd` | Pilot simulation on a known finite population (prevalence 3%, 4% and 5%) showing the same pattern |

The output files are stored at the top level of this repository. Running the script writes them to a `results/` folder.

## How to run

1. Install R. Base R is enough; the `foreign` package comes with R. An internet connection is needed the first time.
2. Download this repository and open `nhanes_analysis.R` in RStudio.
3. Set the working directory to the script's folder (Session → Set Working Directory → To Source File Location).
4. Click **Source**. The analysis takes a few minutes. Data are saved in `data/` and outputs in `results/`.
5. To test the code without downloading anything, set `USE_SYNTHETIC <- TRUE` at the top of the script.

The random seed is fixed, so the results should reproduce exactly.

## Limitations

- One outcome definition (HbA1c only; fasting glucose is not used) and one set of predictors.
- A linear logistic model and simple median imputation.
- NHANES has no known truth, so performance is estimated by design-weighted cross-validation. The pilot
  simulation, where the whole population is known, gives the exact comparison.
- This is preliminary work for a research proposal. It is not a finished study and not clinical guidance.

## Data source

National Center for Health Statistics. National Health and Nutrition Examination Survey (NHANES), 2015–2016
and 2017–2018. Centers for Disease Control and Prevention. https://www.cdc.gov/nhanes/

## Author

Erum Shahzadi, Department of Statistics, University of Gujrat, Pakistan.
