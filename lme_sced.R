#' Fit Linear Mixed Effects Models for SCED Data
#'
#' This function fits linear mixed effects models tailored for Single-Case Experimental Designs (SCED) data,
#' with options to account for AR(1) correlation structures.
#'
#' @param data A data frame containing the variables in the model.
#' @param fixed A formula specifying the fixed effects.
#' @param random A formula specifying the random effects.
#' @param control An optional list of control parameters. Defaults to `list(msMaxIter = 10000, opt = "optim")`.
#' @param ar_1 An optional formula specifying the AR(1) correlation structure.
#'
#' @return A list containing the fitted model, AR(1) coefficient, and transformed data (if AR(1) is specified).
#' @export
#'
#' @examples
#' # Example usage of lme_sced()
#' model <- lme_sced(data = my_data, fixed = Y ~ TIME + PHASE, random = ~ TIME | ID, ar_1 = ~ TIME | ID)

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

    # **NEW STEP**: Convert time variable to numeric if it's not
    if (!is.numeric(data[[time_col]])) {
      # Ensure the time variable is a factor to maintain order if it's character
      data[[time_col]] <- as.factor(data[[time_col]])

      # Get the levels of the time variable in the order they appear in the data
      time_levels <- levels(data[[time_col]])

      # Create a mapping from time levels to numeric values
      time_mapping <- setNames(seq_along(time_levels), time_levels)

      # Map the time variable to numeric
      data[[time_col]] <- as.numeric(time_mapping[data[[time_col]]])
    }

    # **NEW STEP**: Handle factor predictors in fixed effects
    # Create a copy of the data to store the original factor variables
    original_data <- data

    # For fixed predictors
    fixed_factor_predictors <- sapply(fixed_predictors, function(var) {
      is.factor(data[[var]]) || is.character(data[[var]])
    })

    fixed_numeric_predictors <- fixed_predictors[!fixed_factor_predictors]
    fixed_factor_predictors <- fixed_predictors[fixed_factor_predictors]

    # Convert factor predictors to dummy variables
    if (length(fixed_factor_predictors) > 0) {
      for (var in fixed_factor_predictors) {
        # Convert to factor if not already
        data[[var]] <- as.factor(data[[var]])
        # Create dummy variables excluding the reference level
        dummies <- model.matrix(~ data[[var]] - 1)
        # Add dummy variables to data with appropriate names
        colnames(dummies) <- paste0(var, levels(data[[var]]))
        data <- cbind(data, dummies)
      }
    }

    # Update fixed predictors list to include dummy variables and remove original factor predictors
    if (length(fixed_factor_predictors) > 0) {
      fixed_predictors <- c(fixed_numeric_predictors, colnames(dummies))
    }

    # For random predictors (handle similarly if needed)
    random_factor_predictors <- sapply(random_predictors, function(var) {
      is.factor(data[[var]]) || is.character(data[[var]])
    })

    random_numeric_predictors <- random_predictors[!random_factor_predictors]
    random_factor_predictors <- random_predictors[random_factor_predictors]

    # Convert factor predictors to dummy variables
    if (length(random_factor_predictors) > 0) {
      for (var in random_factor_predictors) {
        # Convert to factor if not already
        data[[var]] <- as.factor(data[[var]])
        # Create dummy variables excluding the reference level
        dummies <- model.matrix(~ data[[var]] - 1)
        # Add dummy variables to data with appropriate names
        colnames(dummies) <- paste0(var, levels(data[[var]]))
        data <- cbind(data, dummies)
      }
    }

    # Update random predictors list to include dummy variables and remove original factor predictors
    if (length(random_factor_predictors) > 0) {
      random_predictors <- c(random_numeric_predictors, colnames(dummies))
    }

    # Step 3: Prepare data and get number of measurements per subject
    data[[id_col]] <- as.factor(data[[id_col]])

    # Ensure data is ordered by subject and time
    data <- data[order(data[[id_col]], data[[time_col]]), ]

    nmeasurement_list <- as.numeric(table(data[[id_col]]))  # Number of measurements per subject
    nsubject_list <- length(nmeasurement_list)

    if(is.null(control)){
      control = list(msMaxIter = 10000, opt = "optim")
    }

    # Step 4: Fit the initial AR(1) model
    model_fit <- try(
      nlme::lme(
        fixed = fixed,
        random = random,
        data = original_data,  # Use original data with factor variables
        correlation = corAR1(form = ar_1),
        control = control
      ),
      silent = TRUE
    )

    if (inherits(model_fit, "try-error")) {
      stop("AR(1) estimation error")
    }

    # Step 5: Transformation

    ar1_value <- coef(model_fit$modelStruct$corStruct, unconstrained = FALSE)

    a <- unname((1 - ar1_value^2)^(-1/2))
    b <- a * -ar1_value

    # Create transformation matrices for each subject
    matrix_list <- lapply(nmeasurement_list, function(n_i) {
      if (n_i == 1) {
        # For subjects with only one measurement, the transformation matrix is 1x1 identity
        sparseMatrix(i = 1, j = 1, x = 1)
      } else {
        # Create the index vectors
        row_indices <- c(1:n_i, 2:n_i)
        col_indices <- c(1:n_i, 1:(n_i -1))
        values <- c(1, rep(a, times = n_i -1), rep(b, times = n_i -1))
        sparseMatrix(i = row_indices, j = col_indices, x = values, dims = c(n_i, n_i))
      }
    })

    # Create the block diagonal matrix
    matrix_R <- bdiag(matrix_list)

    # Step 6: Transforming the outcome and predictors
    data[[y_col]] <- as.vector(matrix_R %*% data[[y_col]])
    data[[time_col]] <- as.vector(matrix_R %*% data[[time_col]])

    # Intercept
    data$intercept <- as.vector(matrix_R %*% rep(1, nrow(data)))

    # Transform the fixed predictors
    for (i in seq_along(fixed_predictors)) {
      data[[fixed_predictors[i]]] <- as.vector(matrix_R %*% data[[fixed_predictors[i]]])
    }

    # Transform the random predictors
    for (i in seq_along(random_predictors)) {
      data[[random_predictors[i]]] <- as.vector(matrix_R %*% data[[random_predictors[i]]])
    }

    # Step 7: Create the new formulas dynamically for fixed and random parts
    # Now use the transformed variable names for the transformed data
    fixed_formula <- paste0(y_col, " ~ 0 + intercept + ", time_col)
    if (length(fixed_predictors) > 0) {
      fixed_formula <- paste0(fixed_formula, " + ", paste(fixed_predictors, collapse = " + "))
    }

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

    # Return the model, AR(1) coefficient, and transformed data
    return(list(
      model = model,
      ar1 = ar1_value,
      data = data  # Include the transformed dataset in the output
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
