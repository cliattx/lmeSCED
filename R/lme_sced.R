#' Apply GLS Transformation and Mixed Effects Modeling for SCED
#'
#' Fits a linear mixed-effects model to Single-Case Experimental Design (SCED) data.
#' It handles autocorrelation using a Generalized Least Squares (GLS) transformation
#' and provides fixed effects inference using Satterthwaite's method.
#'
#' @param data A data frame containing the SCED data.
#' @param fixed A formula specifying the fixed effects (e.g., \code{score ~ session + phase}).
#' @param random A formula specifying the random effects (e.g., \code{~ phase | case}).
#' @param control A list of control values for the estimation algorithm. Defaults to \code{list(msMaxIter = 10000, opt = "optim")}.
#' @param ar_1 A formula specifying the AR(1) structure (e.g., \code{~ session | case}). If \code{NULL}, a standard LME is fitted.
#'
#' @return A list containing:
#' \item{model}{The fitted \code{lmer} model object (after GLS transformation).}
#' \item{ar1}{The estimated AR(1) autocorrelation coefficient.}
#' \item{data}{The transformed data frame used for the final model.}
#'
#' @importFrom nlme lme corAR1
#' @importFrom Matrix sparseMatrix bdiag
#' @importFrom lmerTest lmer
#' @importFrom stats as.formula coef model.matrix setNames
#' @export
lme_sced <- function (data, fixed, random, control = NULL, ar_1 = NULL) {
  if (!is.null(ar_1)) {
    # Extract variable names
    id_col <- all.vars(random)[length(all.vars(random))]
    y_col <- all.vars(fixed)[1]
    time_col <- all.vars(ar_1)[1]
    fixed_predictors <- all.vars(fixed)[-1]
    fixed_predictors <- fixed_predictors[fixed_predictors != time_col]
    random_predictors <- all.vars(random)[-length(all.vars(random))]

    # Handle Time Column format
    if (!is.numeric(data[[time_col]])) {
      data[[time_col]] <- as.factor(data[[time_col]])
      time_levels <- levels(data[[time_col]])
      time_mapping <- setNames(seq_along(time_levels), time_levels)
      data[[time_col]] <- as.numeric(time_mapping[data[[time_col]]])
    }
    original_data <- data

    # Handle Fixed Factors
    fixed_factor_predictors <- sapply(fixed_predictors, function(var) {
      is.factor(data[[var]]) || is.character(data[[var]])
    })
    fixed_numeric_predictors <- fixed_predictors[!fixed_factor_predictors]
    fixed_factor_predictors <- fixed_predictors[fixed_factor_predictors]

    if (length(fixed_factor_predictors) > 0) {
      for (var in fixed_factor_predictors) {
        data[[var]] <- as.factor(data[[var]])
        dummies <- model.matrix(~data[[var]] - 1)
        colnames(dummies) <- paste0(var, levels(data[[var]]))
        data <- cbind(data, dummies)
      }
      fixed_predictors <- c(fixed_numeric_predictors, colnames(dummies))
    }

    # Handle Random Factors
    random_factor_predictors <- sapply(random_predictors, function(var) {
      is.factor(data[[var]]) || is.character(data[[var]])
    })
    random_numeric_predictors <- random_predictors[!random_factor_predictors]
    random_factor_predictors <- random_predictors[random_factor_predictors]

    if (length(random_factor_predictors) > 0) {
      for (var in random_factor_predictors) {
        data[[var]] <- as.factor(data[[var]])
        dummies <- model.matrix(~data[[var]] - 1)
        colnames(dummies) <- paste0(var, levels(data[[var]]))
        data <- cbind(data, dummies)
      }
      random_predictors <- c(random_numeric_predictors, colnames(dummies))
    }

    data[[id_col]] <- as.factor(data[[id_col]])
    data <- data[order(data[[id_col]], data[[time_col]]), ]
    nmeasurement_list <- as.numeric(table(data[[id_col]]))

    # Initial AR(1) estimation using nlme
    if (is.null(control)) {
      control = list(msMaxIter = 10000, opt = "optim")
    }
    model_fit <- try(nlme::lme(fixed = fixed, random = random,
                               data = original_data, correlation = nlme::corAR1(form = ar_1),
                               control = control), silent = TRUE)
    if (inherits(model_fit, "try-error")) {
      stop("AR(1) estimation error using nlme.")
    }

    # Transformation Matrix Construction
    ar1_value <- coef(model_fit$modelStruct$corStruct, unconstrained = FALSE)
    a <- unname((1 - ar1_value^2)^(-1/2))
    b <- a * -ar1_value

    matrix_list <- lapply(nmeasurement_list, function(n_i) {
      if (n_i == 1) {
        Matrix::sparseMatrix(i = 1, j = 1, x = 1)
      } else {
        row_indices <- c(1:n_i, 2:n_i)
        col_indices <- c(1:n_i, 1:(n_i - 1))
        values <- c(1, rep(a, times = n_i - 1), rep(b, times = n_i - 1))
        Matrix::sparseMatrix(i = row_indices, j = col_indices, x = values, dims = c(n_i, n_i))
      }
    })
    matrix_R <- Matrix::bdiag(matrix_list)

    # Apply Transformation
    data[[y_col]] <- as.vector(matrix_R %*% data[[y_col]])
    data[[time_col]] <- as.vector(matrix_R %*% data[[time_col]])
    data$intercept <- as.vector(matrix_R %*% rep(1, nrow(data)))

    for (i in seq_along(fixed_predictors)) {
      data[[fixed_predictors[i]]] <- as.vector(matrix_R %*% data[[fixed_predictors[i]]])
    }
    for (i in seq_along(random_predictors)) {
      data[[random_predictors[i]]] <- as.vector(matrix_R %*% data[[random_predictors[i]]])
    }

    # Construct Formulas for Transformed Model
    fixed_formula <- paste0(y_col, " ~ 0 + intercept + ", time_col)
    if (length(fixed_predictors) > 0) {
      fixed_formula <- paste0(fixed_formula, " + ", paste(fixed_predictors, collapse = " + "))
    }
    random_formula <- paste0("( ", paste(random_predictors, collapse = " + "), " | ", id_col, ")")
    second_model_formula <- as.formula(paste0(fixed_formula, " + ", random_formula))

    # Fit Final Model using lmerTest
    model <- try(lmerTest::lmer(second_model_formula, data = data), silent = TRUE)

    if (inherits(model, "try-error")) {
      stop("Transformed model fitting failed.")
    }
    return(list(model = model, ar1 = ar1_value, data = data))

  } else {
    # Standard LME if no AR(1) requested
    model <- try(nlme::lme(fixed = fixed, random = random,
                           data = data, control = list(msMaxIter = 10000, opt = "optim")),
                 silent = TRUE)
    if (inherits(model, "try-error")) {
      stop("nlme estimation error")
    }
    return(model)
  }
}
