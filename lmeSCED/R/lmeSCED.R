#' Fit and Transform Model
#'
#' This function fits an initial linear mixed-effects model using the `nlme` package, applies a transformation based on the fitted model, and then fits a transformed model using the `lmerTest` package.
#'
#' @param data A data frame containing the data to be modeled.
#' @param formula A formula for the fixed effects in the initial model.
#' @param random_formula A formula for the random effects in the initial model.
#' @param ar1_formula A formula for the AR(1) correlation structure in the initial model.
#' @param time_col The name of the time variable.
#' @param phase_col The name of the phase variable.
#' @param interaction_col The name of the interaction variable.
#' @param id_col The name of the ID variable.
#' @param nmeasurement_list The number of measurements per subject.
#' @param nsubject_list The number of subjects in the dataset.
#'
#' @return A fitted mixed-effects model from the `lmerTest` package.
#' @import nlme
#' @import Matrix
#' @import lmerTest
#' @export
#'
#' @examples
#' \dontrun{
#' fit_and_transform_model(data, y ~ time + phase + interaction, ~ 1 + time + phase + interaction | id, ~ time | id, "time", "phase", "interaction", "id")
#' }
fit_and_transform_model <- function(data, formula, random_formula, ar1_formula,
                                    time_col, phase_col, interaction_col, id_col) {

  # Check and install packages if necessary
  packages <- c("nlme", "Matrix", "lmerTest")
  sapply(packages, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  })

  nmeasurement_list <- table(data[[id_col]])[1]  # assuming equal number of measurements per subject
  nsubject_list <- length(unique(data[[id_col]]))

  results_model <- list()

  # Fit the initial model
  model_fit <- try(
    nlme::lme(
      formula,
      random = random_formula,
      data = data,
      correlation = corAR1(form = ar1_formula),
      control = list(msMaxIter = 10000, opt = "optim")
    ),
    silent = TRUE
  )
  results_model$model_fit <- model_fit

  if (inherits(model_fit, "try-error")) {
    stop("Initial model fitting failed.")
  }

  # Transformation
  a <- unname((1 - coef(model_fit$modelStruct$corStruct, unconstrained = "F")^2)^-1/2)
  b <- a * -coef(model_fit$modelStruct$corStruct, unconstrained = "F")

  matrix_block <- sparseMatrix(
    c(1:nmeasurement_list, 2:nmeasurement_list),
    c(1:nmeasurement_list, 1:(nmeasurement_list - 1)),
    x = c(1, rep(a, times = nmeasurement_list - 1), rep(b, times = nmeasurement_list - 1))
  )

  matrix_R <- Diagonal(x = rep(1, nmeasurement_list * nsubject_list))
  for (k in seq_len(nsubject_list)) {
    start <- (k - 1) * nmeasurement_list + 1
    end <- k * nmeasurement_list
    matrix_R[start:end, start:end] <- matrix_block
  }

  data$trans_y <- as.vector(matrix_R %*% data[[all.vars(formula)[1]]])
  data$trans_intercept <- as.vector(matrix_R %*% rep(1, times = nmeasurement_list * nsubject_list))
  data$trans_phase <- as.vector(matrix_R %*% data[[phase_col]])
  data$trans_time <- as.vector(matrix_R %*% data[[time_col]])
  data$trans_interaction <- as.vector(matrix_R %*% data[[interaction_col]])

  # Fit the second model
  model_fit1 <- try(
    lmerTest::lmer(
      trans_y ~ 0 + trans_intercept + trans_time + trans_phase + trans_interaction +
        (-1 + trans_intercept + trans_time + trans_phase + trans_interaction | data[[id_col]]),
      data = data
    ),
    silent = TRUE
  )
  results_model$model_fit1 <- model_fit1

  if (inherits(model_fit1, "try-error")) {
    stop("Transformed model fitting failed.")
  }

  return(results_model$model_fit1)
}
