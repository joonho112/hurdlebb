# ============================================================================
# methods.R --- Standard S3 Methods for hbb_fit Objects
#
# Provides the standard R model-object interface for Hurdle Beta-Binomial
# fits: nobs, coef, vcov, fitted, residuals, summary, and print.summary.
# Enables seamless interoperability with generic model-inspection workflows.
#
# Theory:
#   The Hurdle Beta-Binomial model implies E[Y_i] = n_i * q_i * mu_i,
#   where q_i = logistic(X_i' alpha) is the extensive-margin probability,
#   mu_i = logistic(X_i' beta) is the intensive-margin mean, and kappa
#   governs overdispersion.  The fixed-effect parameter vector is
#   theta = (alpha[1:P], beta[1:P], log_kappa), dimension D = 2P + 1.
#
#   For SVC (state-varying coefficient) models, the linear predictor
#   includes state-level random effects delta[s, 1:2P]:
#     eta_ext[i] = X_i' alpha + X_i' delta_ext[s(i)]
#     eta_int[i] = X_i' beta  + X_i' delta_int[s(i)]
#
#   Pearson residuals use the hurdle variance decomposition:
#     Var[Y_i] = q_i * V_ZT + q_i * (1 - q_i) * E_ZT^2,
#   where E_ZT and V_ZT are the zero-truncated Beta-Binomial mean and
#   variance.  See hurdle_variance() in hurdle.R.
#
#   When a sandwich variance is provided, standard errors and confidence
#   intervals are design-consistent Wald inference:
#     SE_p = sqrt(V_sand[p,p]),
#     CI = theta_hat +/- z_{(1+level)/2} * SE_p.
#   Otherwise, posterior standard deviations and quantile-based intervals
#   are reported.
#
# Existing helpers reused (NOT redefined here):
#   .validate_level(level)     from cholesky-transform.R
#   .make_stan_param_names(P)  from cholesky-transform.R
#   .build_param_labels(hbb_data) from sandwich.R
#   .extract_delta_means(cmdstan_fit, S, K) from sandwich.R
#   hurdle_variance(n, q, mu, kappa) from hurdle.R
#
# Contents:
#   Exported:
#     1. nobs.hbb_fit            --- Number of observations
#     2. coef.hbb_fit            --- Extract fixed-effect coefficients
#     3. vcov.hbb_fit            --- Variance-covariance matrix
#     4. fitted.hbb_fit          --- Fitted (predicted) values
#     5. residuals.hbb_fit       --- Residuals (response or Pearson)
#     6. summary.hbb_fit         --- Tabular summary with CIs
#     7. print.summary.hbb_fit   --- Formatted summary output
#   Internal:
#     8. .validate_hbb_fit_methods --- Comprehensive input validation
#     9. .extract_fixed_effects    --- Fixed-effects data.frame
#    10. .extract_diagnostics      --- MCMC diagnostics list
#    11. .compute_fitted_values    --- Vectorised fitted-value computation
#    12. .print_fe_table           --- Aligned table printing helper
# ============================================================================


# ============================================================================
# 1. nobs.hbb_fit --- Number of observations (exported)
# ============================================================================

#' Number of Observations in an hbb_fit
#'
#' Extracts the number of provider-level observations \eqn{N} used to
#' fit the hurdle Beta-Binomial model.
#'
#' @param object An object of class \code{"hbb_fit"} returned by
#'   \code{\link{hbb}}.
#' @param ... Currently unused; included for S3 method consistency.
#'
#' @return An integer scalar giving the number of observations.
#'
#' @seealso \code{\link{hbb}}, \code{\link{coef.hbb_fit}},
#'   \code{\link{summary.hbb_fit}}
#'
#' @examples
#' \dontrun{
#' fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
#' nobs(fit)
#' }
#'
#' @method nobs hbb_fit
#' @export
nobs.hbb_fit <- function(object, ...) {
    .validate_hbb_fit_methods(object)
    as.integer(object$hbb_data$N)
}


# ============================================================================
# 2. coef.hbb_fit --- Extract fixed-effect coefficients (exported)
# ============================================================================

#' Extract Fixed-Effect Coefficients from an hbb_fit
#'
#' Returns the posterior mean of the fixed-effect parameter vector
#' \eqn{\hat\theta = (\hat\alpha_1, \ldots, \hat\alpha_P,\;
#' \hat\beta_1, \ldots, \hat\beta_P,\; \widehat{\log\kappa})},
#' optionally restricted to a single margin.
#'
#' @param object An object of class \code{"hbb_fit"} returned by
#'   \code{\link{hbb}}.
#' @param margin Character string specifying which parameters to return:
#'   \describe{
#'     \item{\code{"both"}}{(Default.) The full \eqn{D = 2P + 1}
#'       parameter vector.}
#'     \item{\code{"extensive"}}{Only \eqn{\hat\alpha_{1:P}} (P-vector).}
#'     \item{\code{"intensive"}}{Only \eqn{\hat\beta_{1:P}} (P-vector).}
#'   }
#' @param ... Currently unused; included for S3 method consistency.
#'
#' @return A named numeric vector of posterior means.  Length depends on
#'   \code{margin}: \eqn{2P + 1} for \code{"both"}, \eqn{P} for
#'   \code{"extensive"} or \code{"intensive"}.
#'
#' @section Mathematical note:
#' The posterior mean is the Bayes estimator under squared-error loss.
#' For survey-weighted pseudo-posteriors, the posterior mean converges to
#' the pseudo-true parameter under standard regularity conditions
#' (Williams and Savitsky, 2021).
#'
#' @seealso \code{\link{vcov.hbb_fit}}, \code{\link{summary.hbb_fit}}
#'
#' @examples
#' \dontrun{
#' fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
#' coef(fit)                       # full D-vector
#' coef(fit, margin = "extensive") # alpha only
#' coef(fit, margin = "intensive") # beta only
#' }
#'
#' @method coef hbb_fit
#' @export
coef.hbb_fit <- function(object, margin = c("both", "extensive", "intensive"),
                          ...) {

    .validate_hbb_fit_methods(object)
    margin <- match.arg(margin)

    P <- object$hbb_data$P
    D <- 2L * P + 1L
    param_names <- .make_stan_param_names(P)

    # Extract posterior draws with tryCatch
    draws <- tryCatch(
        object$fit$draws(variables = param_names, format = "matrix"),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "Failed to extract posterior draws.",
                "i" = conditionMessage(e)
            ))
        }
    )
    if (!is.matrix(draws)) draws <- as.matrix(draws)

    # Validate draw dimensions
    if (ncol(draws) != D) {
        cli::cli_abort(c(
            "x" = "Draw matrix has {ncol(draws)} columns but expected {D}.",
            "i" = "Parameters: {.val {param_names}}"
        ))
    }

    # Check sufficient MCMC draws
    M <- nrow(draws)
    if (M < 2L) {
        cli::cli_abort(c(
            "x" = "Only {M} MCMC draw{?s} available; need at least 2.",
            "i" = "Check your CmdStan fit for sampling failures."
        ))
    }

    theta_hat <- colMeans(draws)

    # NaN/Inf guard
    if (any(!is.finite(theta_hat))) {
        n_bad <- sum(!is.finite(theta_hat))
        cli::cli_warn(c(
            "!" = "{n_bad} non-finite posterior mean{?s} detected in \\
                   fixed effects.",
            "i" = "Check CmdStan diagnostics for sampling issues."
        ))
    }

    # Apply human-readable labels
    param_labels <- .build_param_labels(object$hbb_data)
    names(theta_hat) <- param_labels

    # Subset by margin
    switch(margin,
        both      = theta_hat,
        extensive = theta_hat[seq_len(P)],
        intensive = theta_hat[(P + 1L):(2L * P)]
    )
}


# ============================================================================
# 3. vcov.hbb_fit --- Variance-covariance matrix (exported)
# ============================================================================

#' Variance-Covariance Matrix of Fixed Effects
#'
#' Returns the \eqn{D \times D} variance-covariance matrix of the
#' fixed-effect parameter vector.
#'
#' @param object An object of class \code{"hbb_fit"} returned by
#'   \code{\link{hbb}}.
#' @param sandwich An optional object of class \code{"hbb_sandwich"}
#'   returned by \code{\link{sandwich_variance}}, or \code{NULL}
#'   (default).
#'   \describe{
#'     \item{If provided:}{The design-consistent sandwich variance
#'       \eqn{V_{\mathrm{sand}} = H_{\mathrm{obs}}^{-1}\,
#'       J_{\mathrm{cluster}}\, H_{\mathrm{obs}}^{-1}} is returned.
#'       This is the recommended choice for survey-weighted models.}
#'     \item{If \code{NULL}:}{The MCMC posterior covariance
#'       \eqn{\Sigma_{\mathrm{MCMC}} = \mathrm{Cov}(\theta^{(m)})}
#'       is returned.  For hierarchical models this is prior-dominated
#'       and should be used only for unweighted base models.}
#'   }
#' @param ... Currently unused; included for S3 method consistency.
#'
#' @return A named numeric matrix of dimension \eqn{D \times D} where
#'   \eqn{D = 2P + 1}.  Row and column names are human-readable
#'   parameter labels (e.g., \code{alpha_intercept}, \code{beta_poverty},
#'   \code{log_kappa}).
#'
#' @section Warning -- prior domination:
#' For hierarchical models with state-varying coefficients,
#' \eqn{\Sigma_{\mathrm{MCMC}}} absorbs prior and random-effect
#' variance, producing Prior Inflation ratios of 600--6500 for
#' fixed effects (see the sandwich variance documentation).  Using
#' \eqn{\Sigma_{\mathrm{MCMC}}} as a variance estimate will
#' substantially over-cover.  The sandwich variance is the
#' appropriate measure of frequentist uncertainty.
#'
#' @seealso \code{\link{coef.hbb_fit}}, \code{\link{sandwich_variance}}
#'
#' @examples
#' \dontrun{
#' fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
#' vcov(fit)  # MCMC posterior covariance
#'
#' # With sandwich variance (recommended for survey data)
#' sand <- sandwich_variance(fit)
#' vcov(fit, sandwich = sand)
#' }
#'
#' @method vcov hbb_fit
#' @export
vcov.hbb_fit <- function(object, sandwich = NULL, ...) {

    .validate_hbb_fit_methods(object)

    P <- object$hbb_data$P
    D <- 2L * P + 1L
    param_labels <- .build_param_labels(object$hbb_data)

    # Validate and dispatch: sandwich V_sand if available
    if (!is.null(sandwich)) {
        if (!inherits(sandwich, "hbb_sandwich")) {
            cli::cli_abort(c(
                "x" = "{.arg sandwich} must be class {.cls hbb_sandwich}.",
                "i" = "Got {.cls {class(sandwich)[1]}}.",
                "i" = "Use {.fn sandwich_variance} to create this object."
            ))
        }
        D_expected <- D
        if (nrow(sandwich$V_sand) != D_expected ||
            ncol(sandwich$V_sand) != D_expected) {
            cli::cli_abort(c(
                "x" = "Dimension mismatch: {.field V_sand} is \\
                       {nrow(sandwich$V_sand)}x{ncol(sandwich$V_sand)}, \\
                       expected {D_expected}x{D_expected}.",
                "i" = "Ensure the sandwich was computed from the same model."
            ))
        }
        V <- sandwich$V_sand
        # Ensure consistent labels
        rownames(V) <- colnames(V) <- param_labels
        return(V)
    }

    # Fallback: MCMC posterior covariance
    param_names <- .make_stan_param_names(P)
    draws <- tryCatch(
        object$fit$draws(variables = param_names, format = "matrix"),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "Failed to extract posterior draws for vcov.",
                "i" = conditionMessage(e)
            ))
        }
    )
    if (!is.matrix(draws)) draws <- as.matrix(draws)

    # Validate dimensions
    if (ncol(draws) != D) {
        cli::cli_abort(c(
            "x" = "Draw matrix has {ncol(draws)} columns but expected {D}.",
            "i" = "Parameters: {.val {param_names}}"
        ))
    }

    # Need M > 1 for cov()
    M <- nrow(draws)
    if (M < 2L) {
        cli::cli_abort(c(
            "x" = "Only {M} MCMC draw{?s} available; need at least 2 \\
                   for covariance estimation."
        ))
    }

    V <- cov(draws)
    rownames(V) <- colnames(V) <- param_labels
    V
}


# ============================================================================
# 4. fitted.hbb_fit --- Fitted (predicted) values (exported)
# ============================================================================

#' Fitted Values from an hbb_fit
#'
#' Computes observation-level fitted values from the posterior means of
#' the fixed effects (and state random effects for SVC models).
#'
#' Three types are available:
#' \describe{
#'   \item{\code{"response"}}{(Default.) The predicted enrollment
#'     proportion \eqn{\hat{q}_i \cdot \hat{\mu}_i}, i.e. the
#'     probability of a positive enrollment times the conditional
#'     enrollment share.  All values lie in \eqn{[0, 1]}.}
#'   \item{\code{"extensive"}}{The predicted participation probability
#'     \eqn{\hat{q}_i = \mathrm{logistic}(X_i' \hat\alpha)}.}
#'   \item{\code{"intensive"}}{The predicted conditional enrollment
#'     share \eqn{\hat\mu_i = \mathrm{logistic}(X_i' \hat\beta)}.}
#' }
#'
#' @param object An object of class \code{"hbb_fit"} returned by
#'   \code{\link{hbb}}.
#' @param type Character string: \code{"response"} (default),
#'   \code{"extensive"}, or \code{"intensive"}.
#' @param ... Currently unused; included for S3 method consistency.
#'
#' @return Numeric vector of length \eqn{N}.
#'
#' @section SVC adjustment:
#' For state-varying coefficient models (\code{model_type \%in\%
#' c("svc", "svc_weighted")}), the linear predictors include
#' state-level random effects:
#' \deqn{
#'   \eta_{\mathrm{ext},i} = X_i' \hat\alpha
#'     + X_i' \hat\delta_{\mathrm{ext}}[s(i)],
#' }
#' where \eqn{\hat\delta} is the posterior mean of the state random
#' effects extracted via \code{.extract_delta_means()}.  If delta
#' extraction fails, a warning is issued and fitted values use fixed
#' effects only.
#'
#' @seealso \code{\link{residuals.hbb_fit}}, \code{\link{coef.hbb_fit}}
#'
#' @examples
#' \dontrun{
#' fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
#' yhat   <- fitted(fit)                      # q * g(mu) (ZT-corrected)
#' q_hat  <- fitted(fit, type = "extensive")  # participation prob
#' mu_hat <- fitted(fit, type = "intensive")  # conditional share
#' }
#'
#' @method fitted hbb_fit
#' @export
fitted.hbb_fit <- function(object, type = c("response", "extensive",
                                              "intensive"), ...) {

    .validate_hbb_fit_methods(object)
    type <- match.arg(type)

    fv <- .compute_fitted_values(object)

    # ZT correction: g(mu) = mu / (1 - p0) for exact unconditional mean
    n_trial  <- as.integer(object$hbb_data$n_trial)
    p0       <- compute_p0(n_trial, fv$mu_hat, fv$kappa_hat)
    mu_trunc <- fv$mu_hat / pmax(1 - p0, .Machine$double.eps)

    switch(type,
        response  = fv$q_hat * mu_trunc,
        extensive = fv$q_hat,
        intensive = mu_trunc
    )
}


# ============================================================================
# 5. residuals.hbb_fit --- Residuals (exported)
# ============================================================================

#' Residuals from an hbb_fit
#'
#' Computes residuals for the hurdle Beta-Binomial model.  Two types
#' are supported:
#'
#' \describe{
#'   \item{\code{"response"}}{Raw residuals on the proportion scale:
#'     \eqn{r_i = y_i / n_i - \hat{q}_i g(\hat{\mu}_i)},
#'     where \eqn{g(\mu) = \mu/(1-p_0)} is the ZT intensity.}
#'   \item{\code{"pearson"}}{Pearson residuals standardised by the
#'     model-implied standard deviation (count scale):
#'     \deqn{
#'       r_i^P = \frac{y_i - n_i \hat{q}_i g(\hat{\mu}_i)}
#'                    {\sqrt{\widehat{\mathrm{Var}}(Y_i)}},
#'     }
#'     where \eqn{\widehat{\mathrm{Var}}(Y_i)} is computed via
#'     \code{\link{hurdle_variance}} using \eqn{(n_i, \hat{q}_i,
#'     \hat\mu_i, \hat\kappa)}.
#'     Under correct model specification, Pearson residuals have
#'     approximately mean zero and variance one.}
#' }
#'
#' @param object An object of class \code{"hbb_fit"} returned by
#'   \code{\link{hbb}}.
#' @param type Character string: \code{"response"} (default) or
#'   \code{"pearson"}.
#' @param ... Currently unused; included for S3 method consistency.
#'
#' @return Numeric vector of length \eqn{N}.
#'
#' @section Numerical stability (Pearson residuals):
#' Several guards are applied to prevent division by zero or NaN:
#' \itemize{
#'   \item \eqn{\hat\mu_i} is clamped to
#'     \eqn{(\varepsilon,\; 1 - \varepsilon)} where
#'     \eqn{\varepsilon = \sqrt{\mathtt{.Machine\$double.eps}}}
#'     (\eqn{\approx 1.49 \times 10^{-8}}).
#'   \item \eqn{\hat{q}_i} is clamped to
#'     \eqn{[\varepsilon,\; 1]}.  (\eqn{q = 0} yields
#'     \eqn{\mathrm{Var}[Y] = 0}, making Pearson undefined.)
#'   \item \eqn{\hat{\kappa}} is clamped to \eqn{\le 10^{15}}.
#'   \item The hurdle variance is floored at
#'     \code{.Machine$double.eps} to prevent division by zero.
#'   \item Non-finite residuals are replaced with \code{NA} and a
#'     warning is issued.
#' }
#'
#' @seealso \code{\link{fitted.hbb_fit}}, \code{\link{hurdle_variance}}
#'
#' @examples
#' \dontrun{
#' fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
#' r_raw     <- residuals(fit)
#' r_pearson <- residuals(fit, type = "pearson")
#' hist(r_pearson, breaks = 50)
#' }
#'
#' @method residuals hbb_fit
#' @export
residuals.hbb_fit <- function(object, type = c("response", "pearson"), ...) {

    .validate_hbb_fit_methods(object)
    type <- match.arg(type)

    hd      <- object$hbb_data
    y       <- hd$y
    n_trial <- as.integer(hd$n_trial)

    fv <- .compute_fitted_values(object)

    # ZT correction: g(mu) = mu / (1 - p0)
    p0       <- compute_p0(n_trial, fv$mu_hat, fv$kappa_hat)
    mu_trunc <- fv$mu_hat / pmax(1 - p0, .Machine$double.eps)

    # Proportion-scale residual: y_i/n_i - q_hat_i * g(mu_hat_i)
    r_response <- y / as.numeric(n_trial) - fv$q_hat * mu_trunc

    if (type == "response") {
        return(r_response)
    }

    # Expected count (ZT-corrected)
    E_y <- as.numeric(n_trial) * fv$q_hat * mu_trunc
    r_count <- y - E_y

    # ---- Pearson residuals with comprehensive guards ---------------

    eps <- .Machine$double.eps^0.5   # ~1.49e-8

    # Clamp mu to (eps, 1-eps) to satisfy hurdle_variance domain
    mu_clamped <- pmin(pmax(fv$mu_hat, eps), 1 - eps)

    # Clamp q to [eps, 1] -- q=0 gives V=0, Pearson undefined
    q_clamped <- pmin(pmax(fv$q_hat, eps), 1)

    # Clamp kappa to prevent overflow in hurdle_variance internals
    kappa_safe <- pmin(fv$kappa_hat, 1e15)

    # Hurdle variance: Var[Y_i] under the two-part model
    var_hat <- hurdle_variance(
        n     = n_trial,
        q     = q_clamped,
        mu    = mu_clamped,
        kappa = kappa_safe
    )

    # Guard against zero/negative variance
    var_hat <- pmax(var_hat, .Machine$double.eps)
    n_zero_var <- sum(var_hat <= .Machine$double.eps)
    if (n_zero_var > 0L) {
        cli::cli_warn(c(
            "!" = "{n_zero_var} observation{?s} with near-zero hurdle \\
                   variance.",
            "i" = "Pearson residuals are clamped for these observations."
        ))
    }

    # Pearson residual on the count scale: (y - n*q*mu) / sqrt(V)
    r_pearson <- r_count / sqrt(var_hat)

    # Final NaN/Inf guard
    n_bad_resid <- sum(!is.finite(r_pearson))
    if (n_bad_resid > 0L) {
        cli::cli_warn(c(
            "!" = "{n_bad_resid} non-finite Pearson residual{?s} detected.",
            "i" = "Setting non-finite residuals to NA."
        ))
        r_pearson[!is.finite(r_pearson)] <- NA_real_
    }

    r_pearson
}


# ============================================================================
# 6. summary.hbb_fit --- Tabular summary with CIs (exported)
# ============================================================================

#' Summary Method for hbb_fit Objects
#'
#' Produces a comprehensive summary of a fitted hurdle Beta-Binomial
#' model, including fixed-effect estimates with standard errors and
#' confidence intervals, dispersion parameter summary, MCMC diagnostics,
#' and optionally design effect ratios from the sandwich variance.
#'
#' @section Inference mode:
#' When \code{sandwich} is supplied, standard errors and confidence
#' intervals are computed from the sandwich variance
#' \eqn{V_{\mathrm{sand}}} (Wald-type intervals):
#' \deqn{
#'   \mathrm{CI}_{1-\alpha}(\theta_p)
#'     = \hat\theta_p \pm z_{(1+\mathrm{level})/2}\,
#'       \sqrt{V_{\mathrm{sand},pp}}.
#' }
#' When \code{sandwich} is \code{NULL}, standard errors are posterior
#' standard deviations and intervals are quantile-based credible
#' intervals from the MCMC draws.
#'
#' @section Dispersion:
#' The dispersion parameter \eqn{\kappa} controls overdispersion
#' relative to the Binomial.  It is estimated on the log scale
#' (\code{log_kappa}) and back-transformed via the delta method:
#' \deqn{\mathrm{SE}(\hat\kappa) = \hat\kappa \cdot
#'   \mathrm{SE}(\widehat{\log\kappa}).}
#' The confidence interval on the \eqn{\kappa} scale is obtained by
#' exponentiating the interval for \eqn{\log\kappa}.
#'
#' @param object An object of class \code{"hbb_fit"} returned by
#'   \code{\link{hbb}}.
#' @param sandwich An optional object of class \code{"hbb_sandwich"}
#'   returned by \code{\link{sandwich_variance}}.  If provided,
#'   produces sandwich-based Wald inference.  If \code{NULL} (default),
#'   produces posterior-based inference.
#' @param level Numeric scalar in \eqn{(0, 1)}.  Confidence/credible
#'   level.  Default is \code{0.95}.
#' @param ... Currently unused; included for S3 method consistency.
#'
#' @return An S3 object of class \code{"summary.hbb_fit"} containing:
#' \describe{
#'   \item{\code{fixed_effects}}{Data frame with columns:
#'     \code{parameter}, \code{estimate}, \code{se}, \code{ci_lower},
#'     \code{ci_upper}, \code{rhat}, \code{ess_bulk}.  One row per
#'     fixed-effect parameter (\eqn{D = 2P + 1}).}
#'   \item{\code{dispersion}}{Named list with \code{kappa_hat}
#'     (point estimate on the natural scale), \code{kappa_se}
#'     (delta-method SE or posterior SD), and \code{kappa_ci}
#'     (2-vector confidence limits on the \eqn{\kappa} scale).}
#'   \item{\code{random_effects}}{\code{NULL} for base/weighted
#'     models.  For SVC models, a list with element \code{tau}
#'     (posterior mean of the random-effect scale parameters).}
#'   \item{\code{diagnostics}}{Named list with \code{n_divergent},
#'     \code{n_max_treedepth}, \code{ebfmi}, \code{max_rhat},
#'     \code{min_ess_bulk}, \code{min_ess_tail}.}
#'   \item{\code{model_info}}{Named list with \code{N}, \code{P},
#'     \code{S} (or \code{NULL}), \code{zero_rate}, \code{model_type},
#'     \code{formula}.}
#'   \item{\code{sandwich_used}}{Logical: whether sandwich SEs were
#'     used.}
#'   \item{\code{DER}}{Named numeric vector of Design Effect Ratios
#'     if \code{sandwich} was provided, \code{NULL} otherwise.}
#'   \item{\code{level}}{The confidence level used.}
#'   \item{\code{call}}{The original model call.}
#' }
#'
#' @seealso
#' \code{\link{print.summary.hbb_fit}} for the print method,
#' \code{\link{coef.hbb_fit}}, \code{\link{vcov.hbb_fit}},
#' \code{\link{sandwich_variance}}
#'
#' @references
#' Williams, M. R. and Savitsky, T. D. (2021).
#' Uncertainty estimation for pseudo-Bayesian inference under complex
#' sampling.
#' \emph{International Statistical Review}, \strong{89}(1), 72--107.
#' \doi{10.1111/insr.12376}
#'
#' @examples
#' \dontrun{
#' fit  <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
#' sand <- sandwich_variance(fit)
#'
#' # Summary with sandwich SEs (recommended for survey data)
#' s <- summary(fit, sandwich = sand, level = 0.95)
#' print(s)
#'
#' # Summary with posterior SDs (for unweighted models)
#' s0 <- summary(fit)
#' print(s0)
#' }
#'
#' @method summary hbb_fit
#' @export
summary.hbb_fit <- function(object, sandwich = NULL, level = 0.95, ...) {

    .validate_hbb_fit_methods(object)
    .validate_level(level)

    hd <- object$hbb_data
    P  <- hd$P
    D  <- 2L * P + 1L

    # -- Sandwich validation -----------------------------------------
    if (!is.null(sandwich)) {
        if (!inherits(sandwich, "hbb_sandwich")) {
            cli::cli_abort(c(
                "x" = "{.arg sandwich} must be class {.cls hbb_sandwich}.",
                "i" = "Got {.cls {class(sandwich)[1]}}.",
                "i" = "Use {.fn sandwich_variance} to create this object."
            ))
        }
        D_expected <- D
        if (nrow(sandwich$V_sand) != D_expected ||
            ncol(sandwich$V_sand) != D_expected) {
            cli::cli_abort(c(
                "x" = "Dimension mismatch: {.field V_sand} is \\
                       {nrow(sandwich$V_sand)}x{ncol(sandwich$V_sand)}, \\
                       expected {D_expected}x{D_expected}.",
                "i" = "Ensure the sandwich was computed from the same model."
            ))
        }
        # Validate DER if present
        if (!is.null(sandwich$DER) && length(sandwich$DER) != D) {
            cli::cli_warn(c(
                "!" = "DER vector length ({length(sandwich$DER)}) does not \\
                       match D ({D}).",
                "i" = "DER will not be reported in the summary."
            ))
        }
    }

    sandwich_used <- !is.null(sandwich)
    param_labels  <- .build_param_labels(hd)

    # ---- Fixed effects table -----------------------------------------------
    fixed_effects <- .extract_fixed_effects(object, sandwich, level,
                                             param_labels)

    # ---- Dispersion (kappa) ------------------------------------------------
    # log_kappa is the last fixed-effect parameter (index D)
    lk_idx <- D
    lk_hat <- fixed_effects$estimate[lk_idx]
    lk_se  <- fixed_effects$se[lk_idx]

    kappa_hat <- exp(lk_hat)
    # Delta method: SE(kappa) = kappa * SE(log_kappa)
    kappa_se  <- kappa_hat * lk_se
    # CI on kappa scale via exp transform of log_kappa CI
    kappa_ci  <- exp(c(fixed_effects$ci_lower[lk_idx],
                       fixed_effects$ci_upper[lk_idx]))

    dispersion <- list(
        kappa_hat = kappa_hat,
        kappa_se  = kappa_se,
        kappa_ci  = kappa_ci
    )

    # ---- Random effects (SVC only) -----------------------------------------
    random_effects <- NULL
    is_svc <- isTRUE(object$model_type %in% c("svc", "svc_weighted"))
    if (is_svc) {
        tau <- tryCatch({
            tau_names <- paste0("tau[", seq_len(2L * P), "]")
            tau_summ <- object$fit$summary(variables = tau_names)
            tau_vec <- tau_summ$mean
            cov_names <- colnames(hd$X)
            if (is.null(cov_names)) cov_names <- paste0("x", seq_len(P))
            names(tau_vec) <- c(
                paste0("tau_alpha_", cov_names),
                paste0("tau_beta_",  cov_names)
            )
            tau_vec
        }, error = function(e) NULL)
        if (!is.null(tau)) {
            random_effects <- list(tau = tau)
        }
    }

    # ---- MCMC diagnostics --------------------------------------------------
    diagnostics <- .extract_diagnostics(object)

    # ---- Model info --------------------------------------------------------
    zero_rate <- 1 - mean(hd$z)

    formula_text <- tryCatch(
        deparse(object$formula$formula, width.cutoff = 80),
        error = function(e) "(formula unavailable)"
    )
    formula_text <- paste(formula_text, collapse = " ")

    model_info <- list(
        N          = hd$N,
        P          = P,
        S          = if (is_svc) hd$S else NULL,
        zero_rate  = zero_rate,
        model_type = object$model_type %||% "(unknown)",
        formula    = formula_text
    )

    # ---- DER ---------------------------------------------------------------
    DER <- NULL
    if (sandwich_used && !is.null(sandwich$DER) &&
        length(sandwich$DER) == D) {
        DER <- sandwich$DER
    }

    # ---- Assemble return object --------------------------------------------
    structure(
        list(
            fixed_effects  = fixed_effects,
            dispersion     = dispersion,
            random_effects = random_effects,
            diagnostics    = diagnostics,
            model_info     = model_info,
            sandwich_used  = sandwich_used,
            DER            = DER,
            level          = level,
            call           = object$call
        ),
        class = "summary.hbb_fit"
    )
}


# ============================================================================
# 7. print.summary.hbb_fit --- Formatted summary output (exported)
# ============================================================================

#' Print Method for summary.hbb_fit Objects
#'
#' Displays a structured summary of the hurdle Beta-Binomial model,
#' including model information, fixed-effect estimates separated by
#' margin with confidence intervals, dispersion parameter, MCMC
#' diagnostics, and design effect ratios (when sandwich SEs are used).
#'
#' The entire method is wrapped in a top-level \code{tryCatch} so that
#' it never crashes, even if the summary object is corrupted or
#' incomplete.
#'
#' @param x An object of class \code{"summary.hbb_fit"} returned by
#'   \code{\link{summary.hbb_fit}}.
#' @param digits Integer; number of significant digits for numeric
#'   output.  Default is 3.
#' @param ... Additional arguments (currently unused).
#'
#' @return The object \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.hbb_fit}}
#'
#' @examples
#' \dontrun{
#' fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
#' s <- summary(fit)
#' print(s, digits = 4)
#' }
#'
#' @method print summary.hbb_fit
#' @export
print.summary.hbb_fit <- function(x, digits = 3, ...) {

    # Top-level tryCatch: NEVER crash from print
    tryCatch({

        mi <- x$model_info
        fe <- x$fixed_effects
        P  <- mi$P

        # ---- Header -------------------------------------------------------
        cat("\n")
        cat(strrep("=", 60), "\n")
        cat("  Hurdle Beta-Binomial Model Summary\n")
        cat(strrep("=", 60), "\n\n")

        # ---- Model info ----------------------------------------------------
        type_labels <- c(
            base         = "Base (unweighted)",
            weighted     = "Weighted (survey-weighted)",
            svc          = "SVC (state-varying coefficients)",
            svc_weighted = "SVC (state-varying, survey-weighted)"
        )
        type_label <- type_labels[mi$model_type]
        if (length(type_label) == 0L || is.na(type_label)) {
            type_label <- mi$model_type %||% "(unknown)"
        }

        cat("  Model type   :", type_label, "\n")
        cat("  Formula      :", mi$formula, "\n")
        cat("  Observations :", mi$N, "\n")
        if (!is.null(mi$S)) {
            cat("  States (S)   :", mi$S, "\n")
        }
        cat("  Zero rate    :", sprintf("%.1f%%", 100 * mi$zero_rate), "\n")

        if (isTRUE(x$sandwich_used)) {
            cat("  Inference    : Wald (sandwich SE)\n")
        } else {
            cat("  Inference    : Posterior (MCMC)\n")
        }
        cat("  Level        :", x$level, "\n")
        cat("\n")

        # ---- Extensive Margin (alpha) --------------------------------------
        cat(strrep("-", 50), "\n")
        cat("  Extensive Margin (alpha)\n")
        cat(strrep("-", 50), "\n")

        alpha_rows <- seq_len(P)
        .print_fe_table(fe[alpha_rows, , drop = FALSE], digits)
        cat("\n")

        # ---- Intensive Margin (beta) ---------------------------------------
        cat(strrep("-", 50), "\n")
        cat("  Intensive Margin (beta)\n")
        cat(strrep("-", 50), "\n")

        beta_rows <- (P + 1L):(2L * P)
        .print_fe_table(fe[beta_rows, , drop = FALSE], digits)
        cat("\n")

        # ---- Dispersion ----------------------------------------------------
        cat(strrep("-", 50), "\n")
        cat("  Dispersion\n")
        cat(strrep("-", 50), "\n")

        disp <- x$dispersion
        if (!is.null(disp) && is.finite(disp$kappa_hat)) {
            lk_se_val <- disp$kappa_se / disp$kappa_hat  # SE of log_kappa
            cat(sprintf("  kappa = %s (log_kappa = %s, SE = %s)\n",
                        format(round(disp$kappa_hat, digits), nsmall = digits),
                        format(round(log(disp$kappa_hat), digits),
                               nsmall = digits),
                        format(round(lk_se_val, digits), nsmall = digits)))
            pct <- round(100 * x$level, 0)
            cat(sprintf("  kappa %d%% CI: [%s, %s]\n",
                        pct,
                        format(round(disp$kappa_ci[1], digits),
                               nsmall = digits),
                        format(round(disp$kappa_ci[2], digits),
                               nsmall = digits)))
        } else {
            cat("  (dispersion unavailable)\n")
        }
        cat("\n")

        # ---- Random effects (SVC) -----------------------------------------
        if (!is.null(x$random_effects)) {
            cat(strrep("-", 50), "\n")
            cat("  Random Effects (tau)\n")
            cat(strrep("-", 50), "\n")
            tau <- x$random_effects$tau
            if (!is.null(tau)) {
                for (i in seq_along(tau)) {
                    cat(sprintf("  %-25s  %s\n",
                                names(tau)[i],
                                format(round(tau[i], digits),
                                       nsmall = digits)))
                }
            }
            cat("\n")
        }

        # ---- MCMC Diagnostics ----------------------------------------------
        cat(strrep("-", 50), "\n")
        cat("  MCMC Diagnostics\n")
        cat(strrep("-", 50), "\n")

        dx <- x$diagnostics
        if (!is.null(dx)) {
            # Helper for safe field printing
            .pf <- function(label, value, warn_fn = NULL) {
                if (is.null(value) ||
                    (length(value) == 1L && is.na(value))) {
                    cat(sprintf("  %-25s: (unavailable)\n", label))
                } else {
                    flag <- ""
                    if (!is.null(warn_fn)) flag <- warn_fn(value)
                    cat(sprintf("  %-25s: %s %s\n", label,
                                format(value, digits = 4), flag))
                }
            }

            .pf("Divergent transitions",
                dx$n_divergent,
                function(v) if (v > 0L) "[WARNING]" else "[OK]")
            .pf("Max treedepth hits",
                dx$n_max_treedepth,
                function(v) if (v > 0L) "[WARNING]" else "[OK]")

            # E-BFMI (may be a vector of per-chain values)
            if (!is.null(dx$ebfmi) && is.numeric(dx$ebfmi) &&
                any(!is.na(dx$ebfmi))) {
                ebfmi_str <- paste(round(dx$ebfmi, 3), collapse = ", ")
                low_flag <- if (any(dx$ebfmi < 0.2, na.rm = TRUE)) {
                    "[WARNING]"
                } else {
                    "[OK]"
                }
                cat(sprintf("  %-25s: %s %s\n", "E-BFMI",
                            ebfmi_str, low_flag))
            }

            .pf("Max Rhat",
                dx$max_rhat,
                function(v) if (v > 1.01) "[WARNING]" else "[OK]")
            .pf("Min ESS (bulk)",
                dx$min_ess_bulk,
                function(v) if (v < 400) "[LOW]" else "[OK]")
            .pf("Min ESS (tail)",
                dx$min_ess_tail,
                function(v) if (v < 400) "[LOW]" else "[OK]")
        } else {
            cat("  (diagnostics unavailable)\n")
        }
        cat("\n")

        # ---- Design Effect Ratios ------------------------------------------
        if (!is.null(x$DER) && length(x$DER) > 0L) {
            cat(strrep("-", 50), "\n")
            cat("  Design Effect Ratios (DER)\n")
            cat(strrep("-", 50), "\n")

            der_labels <- if (!is.null(names(x$DER))) {
                names(x$DER)
            } else if (!is.null(fe)) {
                fe$parameter[seq_along(x$DER)]
            } else {
                paste0("param_", seq_along(x$DER))
            }

            for (i in seq_along(x$DER)) {
                cat(sprintf("  %-25s  %s\n",
                            der_labels[i],
                            format(round(x$DER[i], digits),
                                   nsmall = digits)))
            }
            cat(sprintf("\n  DER range: [%s, %s], mean: %s\n",
                        format(round(min(x$DER), digits), nsmall = digits),
                        format(round(max(x$DER), digits), nsmall = digits),
                        format(round(mean(x$DER), digits),
                               nsmall = digits)))
            cat("\n")
        }

    }, error = function(e) {
        # Defensive fallback: print what we can without crashing
        cat("Error printing summary: ", conditionMessage(e), "\n")
        cat("Raw summary object fields: ",
            paste(names(x), collapse = ", "), "\n")
    })

    invisible(x)
}


# ============================================================================
# Internal helpers
# ============================================================================


# ============================================================================
# 8. .validate_hbb_fit_methods --- Comprehensive input validation (internal)
# ============================================================================

#' Validate that an object is a well-formed hbb_fit for S3 methods
#'
#' Performs comprehensive checks: class verification, presence of
#' \code{hbb_data} with required fields (\code{N}, \code{P}, \code{X}),
#' dimensional consistency of the design matrix, type checks on \code{N}
#' and \code{P}, and presence of the CmdStanMCMC \code{fit} slot.
#' Aborts with structured \code{cli_abort} messages on failure.
#'
#' @param x An object to validate.
#' @return Invisible \code{TRUE}; aborts on failure.
#' @keywords internal
.validate_hbb_fit_methods <- function(x) {

    # 1. Class check
    if (!inherits(x, "hbb_fit")) {
        cli::cli_abort(c(
            "x" = "Expected {.cls hbb_fit} object.",
            "i" = "Got {.cls {class(x)[1]}}."
        ))
    }

    # 2. hbb_data check (needed by almost all methods)
    if (is.null(x$hbb_data)) {
        cli::cli_abort(
            "{.field hbb_data} is NULL in {.cls hbb_fit} object."
        )
    }

    # 3. fit (CmdStanMCMC) check
    if (is.null(x$fit)) {
        cli::cli_abort(c(
            "x" = "{.field fit} is NULL.",
            "i" = "CmdStanMCMC object required for parameter extraction."
        ))
    }

    # 4. Core dimensions with type checks
    N <- x$hbb_data$N
    P <- x$hbb_data$P
    if (is.null(N) || is.null(P)) {
        cli::cli_abort(
            "{.field hbb_data} missing {.field N} or {.field P}."
        )
    }
    if (!is.numeric(N) || N < 1L) {
        cli::cli_abort(
            "{.field N} must be a positive integer; got {.val {N}}."
        )
    }
    if (!is.numeric(P) || P < 1L) {
        cli::cli_abort(
            "{.field P} must be a positive integer; got {.val {P}}."
        )
    }

    # 5. Design matrix X: must exist, be a matrix, match N x P
    if (is.null(x$hbb_data$X) || !is.matrix(x$hbb_data$X)) {
        cli::cli_abort(
            "{.field hbb_data$X} must be a numeric matrix."
        )
    }
    if (nrow(x$hbb_data$X) != N || ncol(x$hbb_data$X) != P) {
        cli::cli_abort(
            "Design matrix X dimensions ({nrow(x$hbb_data$X)} x \\
             {ncol(x$hbb_data$X)}) do not match N={N}, P={P}."
        )
    }

    invisible(TRUE)
}


# ============================================================================
# 9. .extract_fixed_effects --- Fixed-effects data.frame (internal)
# ============================================================================

#' Build the fixed-effects summary table
#'
#' Constructs a data frame with \eqn{D} rows (one per fixed-effect
#' parameter) containing the posterior mean, standard error,
#' confidence/credible interval bounds, Rhat, and bulk ESS.
#'
#' When \code{sandwich} is provided, SEs come from
#' \eqn{\sqrt{\mathrm{diag}(V_{\mathrm{sand}})}} and CIs are Wald
#' intervals: \eqn{\hat\theta_p \pm z_{(1+\mathrm{level})/2} \cdot
#' \mathrm{SE}_p}.  Otherwise, SEs are posterior standard deviations
#' and CIs are quantile-based credible intervals from the MCMC draws.
#'
#' Rhat and ESS are always extracted from the CmdStanR
#' \code{$summary()} method, wrapped in \code{tryCatch}.
#'
#' @param object An \code{hbb_fit} object.
#' @param sandwich An optional \code{hbb_sandwich} object (or
#'   \code{NULL}).
#' @param level Numeric confidence level in (0, 1).
#' @param param_labels Character vector of length D.
#' @return Data frame with columns: parameter, estimate, se, ci_lower,
#'   ci_upper, rhat, ess_bulk.
#' @keywords internal
.extract_fixed_effects <- function(object, sandwich, level, param_labels) {

    P <- object$hbb_data$P
    D <- 2L * P + 1L
    param_names <- .make_stan_param_names(P)

    # -- Extract draws with tryCatch --------------------------------
    draws <- tryCatch(
        object$fit$draws(variables = param_names, format = "matrix"),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "Failed to extract posterior draws for summary.",
                "i" = conditionMessage(e)
            ))
        }
    )
    if (!is.matrix(draws)) draws <- as.matrix(draws)

    # Validate column count
    if (ncol(draws) != D) {
        cli::cli_abort(c(
            "x" = "Draw matrix has {ncol(draws)} columns but expected {D}.",
            "i" = "Parameters: {.val {param_names}}"
        ))
    }

    # Check sufficient MCMC draws
    M <- nrow(draws)
    if (M < 2L) {
        cli::cli_abort(c(
            "x" = "Only {M} MCMC draw{?s} available; need at least 2 \\
                   for SE estimation."
        ))
    }

    theta_hat <- colMeans(draws)

    # NaN/Inf guard on posterior means
    if (any(!is.finite(theta_hat))) {
        n_bad <- sum(!is.finite(theta_hat))
        cli::cli_warn(c(
            "!" = "{n_bad} non-finite posterior mean{?s} detected.",
            "i" = "Summary table may contain NA values."
        ))
    }

    # -- CmdStan summary for Rhat and ESS ------------------------------------
    param_summ <- tryCatch(
        object$fit$summary(variables = param_names),
        error = function(e) NULL
    )

    rhat_vec     <- rep(NA_real_, D)
    ess_bulk_vec <- rep(NA_real_, D)

    if (!is.null(param_summ)) {
        if ("rhat" %in% names(param_summ) &&
            length(param_summ$rhat) >= D) {
            rhat_vec <- param_summ$rhat[seq_len(D)]
        }
        if ("ess_bulk" %in% names(param_summ) &&
            length(param_summ$ess_bulk) >= D) {
            ess_bulk_vec <- param_summ$ess_bulk[seq_len(D)]
        }
    }

    # -- SE and CI: sandwich vs posterior -----------------
    if (!is.null(sandwich) && inherits(sandwich, "hbb_sandwich")) {
        # Sandwich-based Wald inference
        V_sand <- sandwich$V_sand
        if (nrow(V_sand) != D || ncol(V_sand) != D) {
            cli::cli_abort(
                "Sandwich V_sand is {nrow(V_sand)} x {ncol(V_sand)} \\
                 but expected {D} x {D}."
            )
        }

        # Clamp negative diagonals to 0
        v_diag <- diag(V_sand)
        n_neg <- sum(v_diag < 0)
        if (n_neg > 0L) {
            cli::cli_warn(c(
                "!" = "{n_neg} negative diagonal element{?s} in V_sand.",
                "i" = "Clamping to zero before computing SE."
            ))
        }

        se_vec   <- sqrt(pmax(v_diag, 0))
        z_crit   <- qnorm((1 + level) / 2)
        ci_lower <- theta_hat - z_crit * se_vec
        ci_upper <- theta_hat + z_crit * se_vec
    } else {
        # Posterior-based inference
        se_vec <- apply(draws, 2, sd)

        alpha_half <- (1 - level) / 2
        probs <- c(alpha_half, 1 - alpha_half)
        ci_mat   <- apply(draws, 2, quantile, probs = probs)
        ci_lower <- ci_mat[1, ]
        ci_upper <- ci_mat[2, ]
    }

    data.frame(
        parameter = param_labels,
        estimate  = unname(as.numeric(theta_hat)),
        se        = unname(as.numeric(se_vec)),
        ci_lower  = unname(as.numeric(ci_lower)),
        ci_upper  = unname(as.numeric(ci_upper)),
        rhat      = unname(as.numeric(rhat_vec)),
        ess_bulk  = unname(as.numeric(ess_bulk_vec)),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}


# ============================================================================
# 10. .extract_diagnostics --- MCMC diagnostics list (internal)
# ============================================================================

#' Extract MCMC diagnostics from an hbb_fit
#'
#' Safely extracts divergent transitions, max treedepth hits, E-BFMI,
#' max Rhat, and min ESS (bulk and tail) from the CmdStanMCMC object.
#' All extractions are wrapped in \code{tryCatch} so that this function
#' never fails, even if the fit object is corrupted or output files
#' have been moved.
#'
#' @param object An \code{hbb_fit} object.
#' @return Named list with elements: \code{n_divergent},
#'   \code{n_max_treedepth}, \code{ebfmi}, \code{max_rhat},
#'   \code{min_ess_bulk}, \code{min_ess_tail}.  Missing diagnostics
#'   are set to \code{NA}.
#' @keywords internal
.extract_diagnostics <- function(object) {

    result <- list(
        n_divergent     = NA_integer_,
        n_max_treedepth = NA_integer_,
        ebfmi           = NULL,
        max_rhat        = NA_real_,
        min_ess_bulk    = NA_real_,
        min_ess_tail    = NA_real_
    )

    # Diagnostic summary (divergences, treedepth, E-BFMI)
    diag_summ <- tryCatch(
        object$fit$diagnostic_summary(quiet = TRUE),
        error = function(e) NULL
    )

    if (!is.null(diag_summ)) {
        result$n_divergent <- tryCatch(
            sum(diag_summ$num_divergent),
            error = function(e) NA_integer_
        )
        result$n_max_treedepth <- tryCatch(
            sum(diag_summ$num_max_treedepth),
            error = function(e) NA_integer_
        )
        result$ebfmi <- tryCatch(
            diag_summ$ebfmi,
            error = function(e) NULL
        )
    }

    # Rhat and ESS from full summary
    param_summ <- tryCatch(
        object$fit$summary(),
        error = function(e) NULL
    )

    if (!is.null(param_summ)) {
        if ("rhat" %in% names(param_summ)) {
            result$max_rhat <- tryCatch(
                max(param_summ$rhat, na.rm = TRUE),
                error = function(e) NA_real_
            )
        }
        if ("ess_bulk" %in% names(param_summ)) {
            result$min_ess_bulk <- tryCatch(
                min(param_summ$ess_bulk, na.rm = TRUE),
                error = function(e) NA_real_
            )
        }
        if ("ess_tail" %in% names(param_summ)) {
            result$min_ess_tail <- tryCatch(
                min(param_summ$ess_tail, na.rm = TRUE),
                error = function(e) NA_real_
            )
        }
    }

    result
}


# ============================================================================
# 11. .compute_fitted_values --- Vectorised fitted-value computation (internal)
# ============================================================================

#' Compute fitted values from posterior means
#'
#' Uses BLAS-optimised matrix multiplication to compute linear predictors
#' and inverse-logit transformations at the posterior mean of the
#' fixed-effect parameter vector.  For SVC models, adds state-level
#' random effects via \code{.extract_delta_means()}.
#'
#' @section Edge cases handled:
#' \itemize{
#'   \item Extreme \code{log_kappa}: \code{kappa_hat} is clamped via
#'     \code{pmin(exp(log_kappa), 1e15)} to prevent overflow.
#'   \item Non-finite posterior means: a warning is issued but
#'     computation proceeds (\code{plogis} handles +/-Inf gracefully).
#'   \item SVC delta extraction failure: a warning is issued and
#'     fitted values use fixed effects only.
#'   \item All \code{z = 0} (all structural zeros): \code{q_hat}
#'     will be small but not exactly zero (logistic never reaches 0);
#'     \code{mu_hat} is still computed but irrelevant.
#'   \item Draw matrix column mismatch: aborts with informative error.
#'   \item Insufficient MCMC draws (M < 2): aborts with informative error.
#' }
#'
#' @param object An \code{hbb_fit} object (already validated).
#' @return Named list with elements:
#'   \describe{
#'     \item{\code{q_hat}}{N-vector of participation probabilities.}
#'     \item{\code{mu_hat}}{N-vector of conditional intensities.}
#'     \item{\code{kappa_hat}}{Scalar dispersion parameter
#'       (\eqn{\hat\kappa = \exp(\widehat{\log\kappa})}, clamped to
#'       \eqn{[10^{-15}, 10^{15}]}).}
#'     \item{\code{theta_hat}}{D-vector of posterior means.}
#'   }
#' @keywords internal
.compute_fitted_values <- function(object) {

    hd <- object$hbb_data
    X  <- hd$X
    P  <- hd$P
    N  <- hd$N
    D  <- 2L * P + 1L

    param_names <- .make_stan_param_names(P)

    # -- Extract posterior means of fixed effects with tryCatch -----
    draws <- tryCatch(
        object$fit$draws(variables = param_names, format = "matrix"),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "Failed to extract posterior draws for fitted values.",
                "i" = conditionMessage(e)
            ))
        }
    )
    if (!is.matrix(draws)) draws <- as.matrix(draws)

    # Validate column count
    if (ncol(draws) != D) {
        cli::cli_abort(c(
            "x" = "Draw matrix has {ncol(draws)} columns but expected {D}.",
            "i" = "Parameters: {.val {param_names}}"
        ))
    }

    # Check sufficient MCMC draws
    M <- nrow(draws)
    if (M < 2L) {
        cli::cli_abort(c(
            "x" = "Only {M} MCMC draw{?s} available; need at least 2.",
            "i" = "Check your CmdStan fit for sampling failures."
        ))
    }

    theta_hat <- colMeans(draws)

    # NaN/Inf guard on theta_hat
    if (any(!is.finite(theta_hat))) {
        n_bad <- sum(!is.finite(theta_hat))
        cli::cli_warn(c(
            "!" = "{n_bad} non-finite posterior mean{?s} detected.",
            "i" = "Fitted values may be unreliable."
        ))
    }

    alpha_hat     <- theta_hat[seq_len(P)]
    beta_hat      <- theta_hat[(P + 1L):(2L * P)]
    log_kappa_hat <- theta_hat[D]

    # Extreme kappa guard: clamp exp(log_kappa) to [1e-15, 1e15]
    kappa_hat <- pmin(pmax(exp(log_kappa_hat), 1e-15), 1e15)

    # -- Linear predictors via BLAS %*% (no R-level loops) -------------------
    eta_ext <- as.numeric(X %*% alpha_hat)
    eta_int <- as.numeric(X %*% beta_hat)

    # -- SVC adjustment: add delta random effects ----------------------------
    is_svc <- isTRUE(object$model_type %in% c("svc", "svc_weighted"))
    if (is_svc) {
        S <- hd$S
        K <- 2L * P
        delta_hat <- tryCatch(
            .extract_delta_means(object$fit, S, K),
            error = function(e) {
                cli::cli_warn(c(
                    "!" = "Failed to extract delta random effects \\
                           for fitted values.",
                    "i" = conditionMessage(e),
                    "i" = "Returning fitted values without state \\
                           random effects."
                ))
                NULL
            }
        )

        if (!is.null(delta_hat)) {
            state_vec <- hd$state

            # delta_ext: first P columns; delta_int: next P columns
            delta_ext <- delta_hat[state_vec, seq_len(P), drop = FALSE]
            delta_int <- delta_hat[state_vec, (P + 1L):K, drop = FALSE]

            # Vectorised: rowSums(X * delta) adds observation-level RE
            eta_ext <- eta_ext + rowSums(X * delta_ext)
            eta_int <- eta_int + rowSums(X * delta_int)
        }
    }

    # -- Inverse logit (plogis handles overflow/underflow safely) ------------
    q_hat  <- plogis(eta_ext)
    mu_hat <- plogis(eta_int)

    list(
        q_hat     = q_hat,
        mu_hat    = mu_hat,
        kappa_hat = kappa_hat,
        theta_hat = theta_hat
    )
}


# ============================================================================
# 12. .print_fe_table --- Aligned table printing helper (internal)
# ============================================================================

#' Print a subset of the fixed-effects table with aligned columns
#'
#' Formats and prints a fixed-effects data.frame with right-aligned
#' numeric columns using \code{sprintf}.  Used internally by
#' \code{print.summary.hbb_fit} to display the extensive and intensive
#' margin tables separately.
#'
#' @param fe_sub A data.frame (subset of fixed_effects from
#'   \code{summary.hbb_fit}).
#' @param digits Integer; number of decimal places for numeric output.
#' @return Invisible \code{NULL} (side effect: prints to console).
#' @keywords internal
.print_fe_table <- function(fe_sub, digits) {

    # Format numeric columns
    fmt_est  <- format(round(fe_sub$estimate,  digits), nsmall = digits)
    fmt_se   <- format(round(fe_sub$se,        digits), nsmall = digits)
    fmt_lo   <- format(round(fe_sub$ci_lower,  digits), nsmall = digits)
    fmt_hi   <- format(round(fe_sub$ci_upper,  digits), nsmall = digits)
    fmt_rhat <- ifelse(is.na(fe_sub$rhat), "  NA",
                       format(round(fe_sub$rhat, 4), nsmall = 4))
    fmt_ess  <- ifelse(is.na(fe_sub$ess_bulk), "  NA",
                       format(round(fe_sub$ess_bulk, 0)))

    # Compute column widths for alignment
    w_par  <- max(nchar(as.character(fe_sub$parameter)),
                  nchar("Parameter"))
    w_est  <- max(nchar(fmt_est),  nchar("Estimate"))
    w_se   <- max(nchar(fmt_se),   nchar("SE"))
    w_lo   <- max(nchar(fmt_lo),   nchar("CI_lower"))
    w_hi   <- max(nchar(fmt_hi),   nchar("CI_upper"))
    w_rhat <- max(nchar(fmt_rhat), nchar("Rhat"))
    w_ess  <- max(nchar(fmt_ess),  nchar("ESS"))

    # Header row
    cat(sprintf("  %-*s  %*s  %*s  %*s  %*s  %*s  %*s\n",
                w_par, "Parameter",
                w_est, "Estimate",
                w_se,  "SE",
                w_lo,  "CI_lower",
                w_hi,  "CI_upper",
                w_rhat, "Rhat",
                w_ess,  "ESS"))

    # Data rows
    for (i in seq_len(nrow(fe_sub))) {
        cat(sprintf("  %-*s  %*s  %*s  %*s  %*s  %*s  %*s\n",
                    w_par, as.character(fe_sub$parameter[i]),
                    w_est, fmt_est[i],
                    w_se,  fmt_se[i],
                    w_lo,  fmt_lo[i],
                    w_hi,  fmt_hi[i],
                    w_rhat, fmt_rhat[i],
                    w_ess,  fmt_ess[i]))
    }

    invisible(NULL)
}
