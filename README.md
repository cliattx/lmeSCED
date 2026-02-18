# lmeSCED: Generalized Least Squares Transformation for Single-Case Experimental Design

[![R-CMD-check](https://github.com/cliattx/lmeSCED/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cliattx/lmeSCED/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
**lmeSCED** is an R package designed for the analysis of Single-Case Experimental Design (SCED) data using Multilevel Modeling (MLM). It addresses two critical methodological challenges inherent in SCED data: **autocorrelation** and **small sample sizes**.

## Key Features

* **GLS Transformation:** Implements a Generalized Least Squares (GLS) transformation to remove first-order autoregressive (AR(1)) autocorrelation from the residuals, satisfying the independence assumption required by standard MLM software.
* **Small Sample Adjustments:** Combines the GLS transformation with Satterthwaite’s degrees of freedom adjustment (via `lmerTest`) to provide accurate fixed effects inference and control Type I error rates.
* **Variance Component Testing:** Provides the Boundary-Corrected Restricted Likelihood Ratio Test (BC-RLRT) for testing random effects variance components and a parametric bootstrap for testing covariance components.

## Installation

You can install the development version of `lmeSCED` directly from GitHub:

```r
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install lmeSCED
devtools::install_github("cliattx/lmeSCED")
```
## Usage Example

This example demonstrates how to use `lmeSCED` to analyze a multiple baseline design. We use the dataset from *English et al. (1997)* (available via the `SingleCaseES` package).

### 1. Load Data and Package

```r
library(lmeSCED)
# Install SingleCaseES for the example data if needed
if (!requireNamespace("SingleCaseES", quietly = TRUE)) install.packages("SingleCaseES")

# Prepare the data
data_analysis <- SingleCaseES::English1997

# Recode phase: Baseline (A) = 0, Intervention (B) = 1
data_analysis$phase <- ifelse(data_analysis$phase == 'A', 0, 1)

head(data_analysis)
```

### 2. Fit the Model (`lme_sced`)

The `lme_sced` function fits the linear mixed-effects model. You can specify the fixed effects, random effects, and the AR(1) structure.

```r
# Fit the model with AR(1) correction
model_english <- lme_sced(
  data = data_analysis,
  fixed = score ~ session + phase,  # Fixed effects: Intercept, Trend, Shift in Level
  random = ~ phase | case,          # Random effects: Intercept and Phase by Case
  ar_1 = ~ session | case           # AR(1) structure: Session within Case
)

# View the AR(1) autocorrelation estimate
print(model_english$ar1)
# Output: 0.109 (approx)

# View the full model summary (Satterthwaite's method applied)
summary(model_english$model)
```

### 3. Test Variance Components (`rlrt`)

The `rlrt` function performs the Boundary-Corrected Restricted Likelihood Ratio Test (BC-RLRT) for variance components and parametric bootstrapping for covariance.

```r
# Run the test (1000 bootstrap replicates by default)
results_rlrt <- rlrt(model_english)

# View P-values for random effects
print(results_rlrt$p_values_random_effects)
# $intercept: p-value for random baseline level
# $phase: p-value for random treatment effect

# View P-value for covariance between intercept and phase
print(results_rlrt$p_value_covariance)
```

## Methodology

This package implements the methods described in Li, Baek, & Luo (2025). The workflow consists of two stages:

1.  **Transformation:** The function estimates the autocorrelation parameter $\rho$ using `nlme`, constructs a block-diagonal transformation matrix, and transforms the data ($Y$, $X$, $Z$) to whiten the residuals.
2.  **Estimation & Inference:** The transformed data is passed to `lmer` (from `lmerTest`), which estimates parameters using REML and computes Satterthwaite-adjusted p-values for fixed effects.

## Citation

If you use this package in your research, please cite the following manuscript:

> **Li, C., Baek, E., & Luo, W. (2026).** Generalized Least Squares Transformation for Single-case Experimental Design: Introducing the R Package lmeSCED. *Behavior Research Methods*. (Manuscript Accepted).

## Authors

* **Chendong Li** (Maintainer) - Texas A&M University
* **Eunkyeng Baek** - Texas A&M University
* **Wen Luo** - Texas A&M University

## License

This project is licensed under the GPL-3 License.

