# lmeSCED

`lmeSCED` is an R package designed to fit mixed-effects models for single-case experimental design (SCED) data.

## Features

- Mixed-effects model with serial correlated residuals.
- Kenward Roger/Satterthwaite's degrees of freedom for fixed inference adjustment.
- Generalized least squares transformation.
- Designed for SCED data. 

## Installation

You can install the package from the source:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install lmeSCED package
devtools::install_github("cliattx/lmeSCED")

```

## Example

```r
lme_sced(data = data_1,
         formula = y ~ time + phase + interaction,
         random_formula = ~ 1 + time + phase + interaction|id,
         ar1_formula = ~ time|id,
         time_col = "time",
         phase_col = "phase",
         interaction_col = "interaction",
         id_col = "id")

```
