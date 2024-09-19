#' Fit a Linear Mixed Effects Model for SCED Data
#'
#' This function fits a linear mixed effects model for SCED data, with the option to account for AR(1) correlation.
#'
#' @param data A data frame containing the variables for the model.
#' @param fixed A formula specifying the fixed effects.
#' @param random A formula specifying the random effects.
#' @param control A list of control options for the model (optional).
#' @param ar_1 A formula specifying the AR(1) correlation structure (optional).
#'
#' @return A model object or a list containing the model and AR(1) parameter value.
#' @import nlme Matrix lmerTest pbkrtest
#' @export
#' @examples
#' # Example usage:
#' # model <- lme_sced(data, fixed = y ~ time + phase, random = ~ phase | id, ar_1 = ~ time | id)
#'
#'
lme_sced <- function(data,
                     fixed,
                     random,
                     control = NULL,
                     ar_1 = NULL) {
  if(!is.null(ar_1)){
    # Step 1: Load the necessary packages
    packages <- c("nlme", "Matrix", "lmerTest", "pbkrtest")
    sapply(packages, function(pkg) {
      if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
      }
    })

    # Step 2: From formula to variables
    id_col <- all.vars(random)[length(all.vars(random))]
    y_col <- all.vars(fixed)[1]
    time_col <- all.vars(ar_1)[1]

    # Extract predictors for fixed and random effects separately
    fixed_predictors <- all.vars(fixed)[-1]  # Exclude the response variable
    fixed_predictors <- fixed_predictors[fixed_predictors != time_col]  # Exclude the variable in ar_1

    random_predictors <- all.vars(random)[-length(all.vars(random))]  # Exclude the grouping variable

    # Step 3: Get the number of individuals and the measurements
    data[[id_col]] <- as.factor(data[[id_col]])
    nmeasurement_list <- table(data[[id_col]])[1]  # assuming equal number of measurements per subject
    nsubject_list <- length(unique(data[[id_col]]))

    if(is.null(control)){
      control = list(msMaxIter = 10000, opt = "optim")
    }

    # Step 4: Fit the initial AR(1) model
    model_fit <- try(
      nlme::lme(
        fixed = fixed,
        random = random,
        data = data,
        correlation = corAR1(form = ar_1),
        control = control
      ),
      silent = TRUE
    )

    if (inherits(model_fit, "try-error")) {
      stop("AR_1 estimation error")
    }

    # Step 5: Transformation

    ar1_value <- coef(model_fit$modelStruct$corStruct, unconstrained = "F")

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

    # Step 6: Transforming the outcome and predictors
    data[[y_col]] <- as.vector(matrix_R %*% data[[y_col]])
    data[[time_col]] <- as.vector(matrix_R %*% data[[time_col]])

    # intercept
    data$intercept <- as.vector(matrix_R %*% rep(1, nmeasurement_list * nsubject_list))

    # Transform the fixed predictors but keep their original names
    for (i in seq_along(fixed_predictors)) {
      data[[fixed_predictors[i]]] <- as.vector(matrix_R %*% data[[fixed_predictors[i]]])
    }

    # Transform the random predictors but keep their original names
    for (i in seq_along(random_predictors)) {
      data[[random_predictors[i]]] <- as.vector(matrix_R %*% data[[random_predictors[i]]])
    }

    # Step 7: Create the new formulas dynamically for fixed and random parts
    # Now use the original variable names for the transformed data
    fixed_formula <- paste0(y_col, " ~ 0 + intercept + ", time_col, " + ", paste(fixed_predictors, collapse = " + "))

    random_formula <- paste0("( ", paste(random_predictors, collapse = " + "), " | ", id_col, ")")

    # Combine the fixed and random parts into a full model formula
    second_model_formula <- as.formula(paste0(fixed_formula, " + ", random_formula))

    # Fit the second model
    model <- try(
      lmerTest::lmer(
        second_model_formula,
        data = data
      ),
      silent = TRUE
    )

    if (inherits(model, "try-error")) {
      stop("Transformed model fitting failed.")
    }

    return(list(
      model = model,
      ar1 = ar1_value
    ))

  } else {

    model <- try(
      nlme::lme(
        fixed = fixed,
        random = random,
        data = data,
        control = list(msMaxIter = 10000, opt = "optim")
      ),
      silent = TRUE
    )

    if (inherits(model, "try-error")) {
      stop("nlme estimation error")
    }

    return(model)

  }
}
