# lmeSCED

`FitTransformModel` is an R package that fits an initial linear mixed-effects model using the `nlme` package, applies a transformation based on the fitted model, and then fits a transformed model using the `lmerTest` package.

## Installation

You can install the development version of `lmeSCED` from GitHub with:

```R
# Install devtools if you haven't already
install.packages("devtools")

# Install the package from GitHub
devtools::install_github("cliattx/lmeSCED")
