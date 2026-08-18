# Residential Air Pollution and Oculomotor Markers

This repository contains the analysis code accompanying the manuscript:

**Residential Early-Life Air Pollution Exposure and Oculomotor Markers of Cognitive Workload in Children: A Prospective Birth Cohort Study**

The study was conducted within the ENVIRONAGE (ENVIRonmental influence ON early AGEing) birth cohort and investigates associations between residential air pollution exposure and eye-tracking-derived markers of cognitive workload and emotional attention in children.

## Repository contents

The `code/` directory contains the R scripts used for the statistical analyses, sensitivity analyses, model diagnostics, stratified analyses, Bayesian mixture analyses, and figure generation reported in the manuscript and supplementary material.

### Cognitive eye-tracking analyses

- `00_setup_cognitive_analysis.R`  
  Loads the analysis dataset, defines exposure variables, cognitive eye-tracking outcomes and covariates, and contains shared functions used throughout the cognitive analyses.

- `01_primary_cognitive_analysis.R`  
  Performs the primary multiple linear regression analyses for cognitive eye-tracking outcomes. Effect estimates and 95% confidence intervals are expressed per interquartile range (IQR) increase in exposure. Multiple testing is addressed using the effective number of independent tests (Meff).

- `02_cognitive_model_diagnostics.R`  
  Generates the residuals-versus-fitted and normal Q–Q plots reported in the supplementary material. The script also compares the primary linear exposure-response specification with a natural spline model using three degrees of freedom for the selected exposure-outcome associations reported in Supplemental Figure 1 and Supplemental Table 1.

- `03_cognitive_sensitivity_analyses.R`  
  Performs sensitivity analyses for primary cognitive associations with nominal `p < 0.05`. These include additional adjustment for time of examination, reported sleep duration, cumulative lifetime exposure to the corresponding pollutant, and mutual adjustment for exposure during the two remaining trimesters.

- `04_cognitive_sex_stratified_analysis.R`  
  Performs sex-stratified analyses for the cognitive eye-tracking outcomes using the same model specification as the primary analyses, except that child sex is omitted as a covariate within each stratum.

- `05_cognitive_figures.R`  
  Generates the forest plots for the primary cognitive eye-tracking analyses using the output from `01_primary_cognitive_analysis.R`.

### Emotional attention-bias analyses

- `06_setup_emotional_analysis.R`  
  Loads the analysis dataset, defines exposure variables, emotional attention-bias outcomes and covariates, and contains shared functions used throughout the emotional analyses.

- `07_primary_emotional_analysis.R`  
  Performs the primary linear regression analyses for the happy, anger, and fear attention-bias scores. Effect estimates and 95% confidence intervals are expressed per IQR increase in exposure. Multiple testing is addressed using the effective number of independent tests (Meff).

- `08_emotional_sex_stratified_analysis.R`  
  Performs sex-stratified analyses for the emotional attention-bias outcomes using the same model specification as the primary analyses, except that child sex is omitted as a covariate within each stratum.

- `09_emotional_figures.R`  
  Generates the forest plots for the happy, anger, and fear attention-bias outcomes using the output from `07_primary_emotional_analysis.R`.

### BKMR-DLM analyses

- `10_BKMR_DLM_lastYear.R`  
  Performs the Bayesian kernel machine regression distributed lag model (BKMR-DLM) analysis for air pollution exposure during the year preceding the eye-tracking assessment.

- `11_BKMR_DLM_First1000Days.R`  
  Performs the BKMR-DLM analysis for air pollution exposure during the first 1000 days of life.

## Analysis dataset

The R scripts expect the individual-level analysis dataset to be available locally as:

`data/analysis_dataset.xlsx`

The participant-level ENVIRONAGE dataset is **not included in this repository**.

The individual-level data underlying this study cannot currently be made publicly available because they are being used in several ongoing and planned studies that are at different stages of the research and publication process. In addition, access to individual-level cohort data is subject to applicable ethical, privacy, and institutional data-governance requirements.

No individual-level participant data are included in this repository.

## Software

Statistical analyses were performed using:

**R version 4.4.1**

Required R packages are specified within the corresponding analysis scripts.

## Analysis workflow

The cognitive and emotional analyses are organized as separate workflows.

### Cognitive analyses

The intended workflow is:

`00_setup_cognitive_analysis.R`  
→ `01_primary_cognitive_analysis.R`  
→ `02_cognitive_model_diagnostics.R`  
→ `03_cognitive_sensitivity_analyses.R`  
→ `04_cognitive_sex_stratified_analysis.R`  
→ `05_cognitive_figures.R`

The setup script contains the shared variable definitions and analysis functions used by the downstream cognitive scripts.

### Emotional attention-bias analyses

The intended workflow is:

`06_setup_emotional_analysis.R`  
→ `07_primary_emotional_analysis.R`  
→ `08_emotional_sex_stratified_analysis.R`  
→ `09_emotional_figures.R`

The setup script contains the shared variable definitions and analysis functions used by the downstream emotional scripts.

### BKMR-DLM analyses

The BKMR-DLM analyses are provided separately in:

- `10_BKMR_DLM_lastYear.R`
- `11_BKMR_DLM_First1000Days.R`

## Multiple-testing correction

For the primary eye-tracking analyses, multiple testing was addressed using the effective number of independent exposure tests (Meff), estimated using the Li and Ji method based on the correlation structure of the exposure variables.

Early-life and recent exposure windows were treated as separate exposure hypothesis families. The corresponding Meff-adjusted significance threshold is calculated within the analysis scripts.

## Cognitive model diagnostics

For selected primary cognitive exposure-outcome associations, model assumptions were evaluated using residuals-versus-fitted plots and normal Q–Q plots.

Potential non-linearity was evaluated by comparing the primary linear model with a natural spline model using three degrees of freedom for the exposure term.

These analyses correspond to the diagnostic results reported in Supplemental Figure 1 and Supplemental Table 1.

## Sensitivity analyses

For selected cognitive associations with nominal `p < 0.05` in the primary analysis, sensitivity analyses evaluated additional adjustment for:

- time of examination;
- reported sleep duration;
- cumulative lifetime exposure to the corresponding pollutant;
- exposure during the two remaining trimesters for trimester-specific models.

For the trimester-specific sensitivity analysis, the exposure estimates for the two remaining trimesters were entered separately as additional covariates.

## Sex-stratified analyses

Sex-stratified analyses were conducted separately for male and female participants.

The same model specification as in the corresponding primary analysis was used, except that child sex was omitted from the covariate set within each sex stratum.

## Reproducibility

The analysis code is made publicly available to facilitate transparency and reproducibility of the statistical workflow.

Because the underlying participant-level data cannot be publicly distributed, the analyses cannot be reproduced directly from this repository alone. However, the provided scripts document the data-processing steps, variable definitions, statistical models, model diagnostics, sensitivity analyses, multiple-testing procedures, sex-stratified analyses, Bayesian mixture analyses, and figure-generation procedures used for the manuscript.

Generated analysis outputs are written to the `results/` directory and are not included in the repository.

## Data availability

The individual-level ENVIRONAGE cohort data underlying this study are not publicly available.

The data are currently being used in several ongoing and planned studies that are at different stages of the research and publication process. Access to individual-level cohort data is additionally subject to applicable ethical, privacy, and institutional data-governance requirements.

The analysis code used for the present study is publicly available through this repository.

## Citation

If you use or refer to the code provided in this repository, please cite the associated publication once available.

A formal citation for this repository will be added upon publication of the manuscript.

## Contact

For questions regarding the study, analysis code, or data-access procedures, please contact:

**Prof. Michelle Plusquin**  
Corresponding author  
Hasselt University  
Centre for Environmental Sciences  
Email: michelle.plusquin@uhasselt.be
