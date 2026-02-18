#' Restricted Likelihood Ratio Test for Variance Components
#'
#' Performs a boundary-corrected restricted likelihood ratio test (BC-RLRT) for variance components
#' and a parametric bootstrap for covariance components in SCED models.
#'
#' @param lme_sced_results A list object returned by the \code{lme_sced} function.
#'
#' @return A list containing:
#' \item{full_model_nlme}{The full model fitted with nlme.}
#' \item{diag_model_nlme}{The reduced model (diagonal covariance) fitted with nlme.}
#' \item{p_values_random_effects}{P-values for the variance components (BC-RLRT).}
#' \item{p_value_covariance}{P-value for the covariance components (Parametric Bootstrap).}
#'
#' @importFrom lme4 findbars ranef nobars
#' @importFrom nlme lme pdDiag lmeControl gls
#' @importFrom stats pchisq simulate formula logLik
#' @export
rlrt <- function (lme_sced_results) {
  if (!is.list(lme_sced_results) || !all(c("model", "data") %in% names(lme_sced_results))) {
    stop("Input must be a list returned by lme_sced(), containing 'model' and 'data'.")
  }

  model <- lme_sced_results$model
  data <- lme_sced_results$data
  full_formula <- formula(model)

  # Extract formula components using lme4 utilities
  random_part <- lme4::findbars(full_formula)[[1]]
  id_var <- names(lme4::ranef(model))
  random_vars <- c("intercept", all.vars(random_part[[2]]))
  fixed_vars <- all.vars(lme4::nobars(full_formula))[-1]
  y_var <- all.vars(full_formula[[2]])

  if (id_var != "IDLEVEL2") {
    names(data)[names(data) == id_var] <- "IDLEVEL2"
  }

  # Construct formulas for nlme
  fixed_formula_nlme <- as.formula(paste(y_var, "~ -1 +", paste(fixed_vars, collapse = " + ")))
  random_formula_nlme <- as.formula(paste("~ 0 +", paste(random_vars, collapse = " + ")))

  # Fit Models for testing
  mA <- try(nlme::lme(fixed = fixed_formula_nlme, random = list(IDLEVEL2 = random_formula_nlme),
                      data = data, control = nlme::lmeControl(msMaxIter = 1000, opt = "optim")), silent = TRUE)

  m.vc <- try(nlme::lme(fixed = fixed_formula_nlme, random = list(IDLEVEL2 = nlme::pdDiag(random_formula_nlme)),
                        data = data, control = nlme::lmeControl(msMaxIter = 1000, opt = "optim")), silent = TRUE)

  if (inherits(mA, "try-error") || inherits(m.vc, "try-error")) {
    stop("Model fitting with nlme::lme failed. This often happens with singular fits.")
  }

  # Variance Components Test
  p_values <- list()
  for (var in random_vars) {
    reduced_vars <- setdiff(random_vars, var)
    if (length(reduced_vars) > 0) {
      random_formula_reduced <- as.formula(paste("~ 0 +", paste(reduced_vars, collapse = " + ")))
      m_reduced <- try(nlme::lme(fixed = fixed_formula_nlme,
                                 random = list(IDLEVEL2 = nlme::pdDiag(random_formula_reduced)),
                                 data = data, control = nlme::lmeControl(msMaxIter = 1000, opt = "optim")), silent = TRUE)
    } else {
      m_reduced <- try(nlme::gls(model = fixed_formula_nlme, data = data), silent = TRUE)
    }

    if (inherits(m_reduced, "try-error")) {
      p_values[[var]] <- NA
      next
    }

    lr_stat <- -2 * (logLik(m_reduced)[1] - logLik(m.vc)[1])
    lr_stat <- max(0, lr_stat)
    p_values[[var]] <- 0.5 * pchisq(lr_stat, df = 1, lower.tail = FALSE)
  }

  # Covariance Test (Parametric Bootstrap)
  sim <- stats::simulate(m.vc, nsim = 1000, m2 = mA, method = "REML")
  obs_lr_stat_covar <- -2 * (logLik(m.vc)[1] - logLik(mA)[1])
  obs_lr_stat_covar <- max(0, obs_lr_stat_covar)

  # Extract null distribution from simulation
  null_distribution <- -2 * sim$null$REML[, 2] - (-2 * sim$alt$REML[, 2])
  p_value_covariance <- mean(null_distribution >= obs_lr_stat_covar, na.rm = TRUE)

  return(list(full_model_nlme = mA, diag_model_nlme = m.vc,
              p_values_random_effects = p_values, p_value_covariance = p_value_covariance))
}
