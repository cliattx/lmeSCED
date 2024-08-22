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

## Syntax

This package is designed using the same syntax from the `nlme` package.

## Example

```r
lme_sced(data = data_1,
         fixed = y ~ time + phase,
         random = ~ time | id,
         ar_1 = ~ time | id)
```

## Limitations 

The package is limited to the balanced-design right now (i.e., same number of observations across individuals). 
