#' Fit and Transform Model
#'
#' This function fits an initial model using nlme::lme, performs transformations, and then fits a second model using lmerTest::lmer.
#'
#' @param data Data frame containing the dataset.
#' @param formula A formula for the fixed effects in the initial model.
#' @param random_formula A formula for the random effects in the initial model.
#' @param ar1_formula A formula for the AR1 correlation structure in the initial model.
#' @param time_col The name of the time column in the dataset.
#' @param phase_col The name of the phase column in the dataset.
#' @param interaction_col The name of the interaction column in the dataset.
#' @param id_col The name of the ID column in the dataset.
#' @return The fitted lmer model after transformations.
#' @import nlme Matrix lmerTest
#' @export
#'
#' @examples
#' # Assuming `data` is a data frame with appropriate columns
#' # lme_sced(data, y ~ time + phase + interaction, ~ 1 + time + phase + interaction | id, ~ time | id, "time", "phase", "interaction", "id")
#' 
#' 
lme_sced <- function(data, formula, random_formula, ar1_formula,
                                    time_col, phase_col, interaction_col, id_col) {

  # Check and install packages if necessary
  packages <- c("nlme", "Matrix", "lmerTest")
  sapply(packages, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  })
  
  data[[id_col]] <- as.factor(data[[id_col]])

  nmeasurement_list <- table(data[[id_col]])[1]  # assuming equal number of measurements per subject
  nsubject_list <- length(unique(data[[id_col]]))


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

  # Create the formula for the second model
  second_model_formula <- as.formula(paste0(
    "trans_y ~ 0 + trans_intercept + trans_time + trans_phase + trans_interaction + ",
    "(0 + trans_intercept + trans_time + trans_phase + trans_interaction | ", id_col, ")"
  ))
  
  # Fit the second model
  model_fit1 <- try(
    lmerTest::lmer(
      second_model_formula,
      data = data
    ),
    silent = TRUE
  )
  
  if (inherits(model_fit1, "try-error")) {
    stop("Transformed model fitting failed.")
  }
  
  return(model_fit1)
}