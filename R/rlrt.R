#' Restricted Likelihood Ratio Test for Variance Components
#'
#' Performs a restricted likelihood ratio test (RLRT) to test the significance of variance components
#' in a linear mixed-effects model.
#'
#' @param data A data frame containing the variables in the model.
#' @param id_var A string specifying the name of the grouping variable.
#' @param fixed_vars A character vector of fixed effect variable names.
#' @param random_vars A character vector of random effect variable names.
#' @param control An optional list of control parameters for model fitting.
#'
#' @return A list containing the full model, diagonal covariance model, p-values for random effects, and p-value for covariance between random effects.
#' @export
#'
#' @examples
#' # Example usage of rlrt()
#' results <- rlrt(data = my_data, id_var = "ID", fixed_vars = c("TIME", "PHASE"), random_vars = c("TIME"))


rlrt <- function(data, id_var, fixed_vars, random_vars, control = NULL) {

  # Rename the grouping variable to 'IDLEVEL2' if it's different
  if (id_var != "IDLEVEL2") {
    names(data)[names(data) == id_var] <- "IDLEVEL2"
  }

  if (is.null(control)) {
    control <- list(msMaxIter = 10000, opt = "optim")
  }

  # Prepare formulas for fixed and random effects
  fixed_formula <- as.formula(paste("Y ~ -1 +", paste(fixed_vars, collapse = " + ")))
  random_formula_full <- as.formula(paste("~ 0 +", paste(random_vars, collapse = " + "), "|", "IDLEVEL2"))

  # Full model with covariance between random effects using nlme::lme()
  mA <- try(
    nlme::lme(
      fixed = fixed_formula,
      random = random_formula_full,
      data = data,
      control = control
    ),
    silent = TRUE
  )

  # Model with diagonal covariance matrix using nlme::lme()
  m.vc <- try(
    nlme::lme(
      fixed = fixed_formula,
      random = list(IDLEVEL2 = pdDiag(as.formula(paste("~ 0 +", paste(random_vars, collapse = " + "))))),
      data = data,
      control = control
    ),
    silent = TRUE
  )

  # Check if both models fitted successfully
  if (!inherits(mA, "try-error") && !inherits(m.vc, "try-error")) {

    # Proceed with RLRT as before

    # Models with only one random effect at a time
    p_values <- list()
    for (var in random_vars) {
      random_formula_single <- as.formula(paste("~ 0 +", var, "|", "IDLEVEL2"))
      m_single <- try(
        nlme::lme(
          fixed = fixed_formula,
          random = random_formula_single,
          data = data,
          control = control
        ),
        silent = TRUE
      )

      if (inherits(m_single, "try-error")) {
        p_values[[var]] <- NA
        next
      }

      # Calculate the p-value
      df_diff <- length(random_vars) - 1  # Degrees of freedom difference
      lr_stat <- -2 * (m_single$logLik - m.vc$logLik)
      p_value <- 0.5 * (1 - pchisq(lr_stat, df = df_diff - 1)) +
        0.5 * (1 - pchisq(lr_stat, df = df_diff))
      p_values[[var]] <- p_value
    }

    # Test for covariance between random effects
    Sim <- simulate.lme(m.vc, nsim = 1000, m2 = mA, method = 'REML')
    lTR.null <- -2 * Sim$null$REML[, 2] - (-2 * Sim$alt$REML[, 2])
    lr_stat_covar <- -2 * m.vc$logLik - (-2 * mA$logLik)
    p_covar_PB <- mean(lTR.null > lr_stat_covar, na.rm = TRUE)

    return(list(
      full_model = mA,
      diag_model = m.vc,
      p_values_random_effects = p_values,
      p_value_covariance = p_covar_PB
    ))

  } else {
    # If model fitting failed, use lmerTest::ranova()
    message("Model fitting failed with nlme::lme(). Using lmerTest::lmer() and ranova() instead.")

    # Prepare lmer formulas
    fixed_formula_lmer <- as.formula(paste("Y ~ -1 +", paste(fixed_vars, collapse = " + ")))
    random_formula_full_lmer <- paste("(0 +", paste(random_vars, collapse = " + "), "| IDLEVEL2)")
    random_formula_diag_lmer <- paste("(0 +", paste(random_vars, collapse = " + "), "|| IDLEVEL2)")
    lmer_formula_full <- as.formula(paste(deparse(fixed_formula_lmer), "+", random_formula_full_lmer))
    lmer_formula_diag <- as.formula(paste(deparse(fixed_formula_lmer), "+", random_formula_diag_lmer))

    # Fit the full model using lmerTest::lmer()
    lmer_full_model <- try(
      lmerTest::lmer(
        formula = lmer_formula_full,
        data = data,
        REML = TRUE
      ),
      silent = TRUE
    )

    # Fit the diagonal covariance model using lmerTest::lmer()
    lmer_diag_model <- try(
      lmerTest::lmer(
        formula = lmer_formula_diag,
        data = data,
        REML = TRUE
      ),
      silent = TRUE
    )

    if (inherits(lmer_full_model, "try-error") || inherits(lmer_diag_model, "try-error")) {
      stop("Model fitting failed with both nlme::lme() and lmerTest::lmer().")
    }

    # Perform ranova() on the full model
    ranova_results <- try(
      lmerTest::ranova(lmer_full_model, reduce.terms = TRUE),
      silent = TRUE
    )

    if (inherits(ranova_results, "try-error")) {
      stop("ranova() failed.")
    }

    # Compare full and diagonal models to test covariance between random effects
    anova_covar <- anova(lmer_full_model, lmer_diag_model, refit = FALSE)
    p_value_covar <- anova_covar$`Pr(>Chisq)`[2]

    return(list(
      lmer_full_model = lmer_full_model,
      lmer_diag_model = lmer_diag_model,
      ranova_results = ranova_results,
      p_value_covariance = p_value_covar
    ))
  }
}
