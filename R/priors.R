# ============================================================================
# priors.R --- Prior Specification for Hurdle Beta-Binomial Models
#
# User-facing functions for specifying and inspecting priors. The prior
# object is an S3 class ("hbb_prior") that stores distribution families
# and hyperparameters for each model parameter block.
#
# Currently supported distributions:
#   - Normal (the only option for now; `dist` field allows future extension)
#
# Stan parameter mapping:
#   alpha     ~ Normal(mean, sd)           [extensive-margin fixed effects]
#   beta      ~ Normal(mean, sd)           [intensive-margin fixed effects]
#   log_kappa ~ Normal(mean, sd)           [log concentration]
#   gamma     ~ Normal(mean, sd)           [policy moderator coefficients]
#   tau       ~ Normal(mean, sd) [lower=0] [random effect scales; half-normal]
#   L_Omega   ~ LKJ(eta)                  [Cholesky correlation factor]
#
# Contents:
#   1. hbb_prior      --- User-facing prior specification (exported)
#   2. default_prior   --- Default prior specification (exported)
#   3. print.hbb_prior --- Print method
#   4. .validate_prior_component --- Internal: validate a single prior list
# ============================================================================


# ============================================================================
# 1. hbb_prior --- User-facing prior specification (exported)
# ============================================================================

#' Specify Priors for a Hurdle Beta-Binomial Model
#'
#' Creates a prior specification object that controls the prior
#' distributions used in the Stan model. All arguments are optional;
#' unspecified components receive the defaults from [default_prior()].
#'
#' @param alpha Prior on the extensive-margin fixed effects
#'   \eqn{\alpha} (logit scale). A named list with elements `dist`
#'   (character), `mean` (numeric scalar), and `sd` (positive numeric
#'   scalar). Currently only `dist = "normal"` is supported.
#' @param beta Prior on the intensive-margin fixed effects \eqn{\beta}
#'   (logit scale). Same structure as `alpha`.
#' @param log_kappa Prior on the log-concentration parameter
#'   \eqn{\log \kappa}. Same structure as `alpha`.
#' @param gamma Prior on the policy moderator coefficients
#'   \eqn{\Gamma}. Same structure as `alpha`.
#' @param tau Prior on the random effect scale parameters \eqn{\tau}.
#'   Same structure as `alpha`. Note: in Stan, \eqn{\tau} is declared
#'   with `lower = 0`, so a `Normal(mean, sd)` prior becomes a
#'   **half-normal** prior in practice. The `mean` component is
#'   typically set to 0 for half-normal priors.
#' @param lkj_eta Concentration parameter \eqn{\eta > 0} for the
#'   LKJ prior on the Cholesky factor of the correlation matrix
#'   \eqn{L_\Omega}. Higher values shrink the correlation toward the
#'   identity matrix. A scalar \eqn{\eta = 1} gives a uniform prior
#'   over correlation matrices; \eqn{\eta = 2} (default) mildly favours
#'   the identity.
#'
#' @return An S3 object of class `"hbb_prior"` with components:
#' \describe{
#'   \item{`alpha`}{Named list with `dist`, `mean`, `sd`.}
#'   \item{`beta`}{Named list with `dist`, `mean`, `sd`.}
#'   \item{`log_kappa`}{Named list with `dist`, `mean`, `sd`.}
#'   \item{`gamma`}{Named list with `dist`, `mean`, `sd`.}
#'   \item{`tau`}{Named list with `dist`, `mean`, `sd`.}
#'   \item{`lkj_eta`}{Positive numeric scalar.}
#' }
#'
#' @details
#' ## Extensibility
#'
#' The `dist` field in each component currently accepts only `"normal"`.
#' The structure is designed to accommodate future distributions (e.g.,
#' Student-t) without breaking the API.
#'
#' ## Model-specific usage
#'
#' Not all parameters are present in every model variant:
#' \itemize{
#'   \item **hbb_base** and **hbb_weighted**: use `alpha`, `beta`,
#'     `log_kappa` only.
#'   \item **hbb_svc** and **hbb_svc_weighted**: additionally use
#'     `gamma`, `tau`, and `lkj_eta`.
#' }
#' Unused prior components are silently ignored by the fitting function.
#'
#' @examples
#' # Default priors
#' default_prior()
#'
#' # Tighter prior on fixed effects
#' hbb_prior(alpha = list(dist = "normal", mean = 0, sd = 1))
#'
#' # Custom prior on dispersion and LKJ
#' hbb_prior(
#'   log_kappa = list(dist = "normal", mean = 3, sd = 1),
#'   lkj_eta   = 4
#' )
#'
#' # Override everything
#' p <- hbb_prior(
#'   alpha     = list(dist = "normal", mean = 0, sd = 1),
#'   beta      = list(dist = "normal", mean = 0, sd = 1),
#'   log_kappa = list(dist = "normal", mean = 2, sd = 1),
#'   gamma     = list(dist = "normal", mean = 0, sd = 0.5),
#'   tau       = list(dist = "normal", mean = 0, sd = 0.5),
#'   lkj_eta   = 5
#' )
#' print(p)
#'
#' @seealso [default_prior()] for the package defaults.
#' @family priors
#' @export
hbb_prior <- function(alpha     = NULL,
                       beta      = NULL,
                       log_kappa = NULL,
                       gamma     = NULL,
                       tau       = NULL,
                       lkj_eta   = NULL) {

    # -- Retrieve defaults for any unspecified component -----------------------
    defaults <- .default_prior_raw()

    if (is.null(alpha))     alpha     <- defaults$alpha
    if (is.null(beta))      beta      <- defaults$beta
    if (is.null(log_kappa)) log_kappa <- defaults$log_kappa
    if (is.null(gamma))     gamma     <- defaults$gamma
    if (is.null(tau))       tau       <- defaults$tau
    if (is.null(lkj_eta))   lkj_eta   <- defaults$lkj_eta

    # -- Validate each component ----------------------------------------------
    .validate_prior_component(alpha,     "alpha")
    .validate_prior_component(beta,      "beta")
    .validate_prior_component(log_kappa, "log_kappa")
    .validate_prior_component(gamma,     "gamma")
    .validate_prior_component(tau,       "tau")

    # Validate lkj_eta
    assert_number(lkj_eta, lower = .Machine$double.eps,
                  .var.name = "lkj_eta")

    # -- Construct S3 object --------------------------------------------------
    structure(
        list(
            alpha     = alpha,
            beta      = beta,
            log_kappa = log_kappa,
            gamma     = gamma,
            tau       = tau,
            lkj_eta   = lkj_eta
        ),
        class = "hbb_prior"
    )
}


# ============================================================================
# 2. default_prior --- Default prior specification (exported)
# ============================================================================

#' Default Prior Specification for Hurdle Beta-Binomial Models
#'
#' Returns the package default priors, which are weakly informative and
#' match the Stan model code:
#'
#' \tabular{ll}{
#'   **Parameter** \tab **Prior** \cr
#'   `alpha`     \tab Normal(0, 2) \cr
#'   `beta`      \tab Normal(0, 2) \cr
#'   `log_kappa` \tab Normal(2, 1.5) \cr
#'   `gamma`     \tab Normal(0, 1) \cr
#'   `tau`       \tab HalfNormal(0, 1) \cr
#'   `L_Omega`   \tab LKJ(2) \cr
#' }
#'
#' @return An S3 object of class `"hbb_prior"`. See [hbb_prior()] for
#'   the full structure.
#'
#' @examples
#' dp <- default_prior()
#' dp
#'
#' # Inspect a specific component
#' dp$alpha
#' dp$lkj_eta
#'
#' @seealso [hbb_prior()] for custom prior specification.
#' @family priors
#' @export
default_prior <- function() {
    hbb_prior(
        alpha     = list(dist = "normal", mean = 0,   sd = 2),
        beta      = list(dist = "normal", mean = 0,   sd = 2),
        log_kappa = list(dist = "normal", mean = 2,   sd = 1.5),
        gamma     = list(dist = "normal", mean = 0,   sd = 1),
        tau       = list(dist = "normal", mean = 0,   sd = 1),
        lkj_eta   = 2
    )
}


# ============================================================================
# 3. print.hbb_prior --- Print method
# ============================================================================

#' @export
print.hbb_prior <- function(x, ...) {

    cat("Hurdle Beta-Binomial Prior Specification\n")
    cat("----------------------------------------\n")

    # Helper: format a Normal prior component
    fmt <- function(comp, label, half = FALSE) {
        prefix <- if (half) "HalfNormal" else "Normal"
        sprintf("  %-9s ~ %s(%s, %s)",
                label, prefix,
                format(comp$mean, digits = 4, drop0trailing = TRUE),
                format(comp$sd, digits = 4, drop0trailing = TRUE))
    }

    cat(fmt(x$alpha,     "alpha"),     "\n")
    cat(fmt(x$beta,      "beta"),      "\n")
    cat(fmt(x$log_kappa, "log_kappa"), "\n")
    cat(fmt(x$gamma,     "gamma"),     "\n")
    cat(fmt(x$tau,       "tau", half = TRUE), "\n")
    cat(sprintf("  %-9s ~ LKJ(%s)",
                "L_Omega",
                format(x$lkj_eta, digits = 4, drop0trailing = TRUE)),
        "\n")

    invisible(x)
}


# ============================================================================
# 4. .validate_prior_component --- Internal validation helper
# ============================================================================

#' Validate a single prior component
#'
#' Checks that a prior component is a named list with `dist`, `mean`,
#' and `sd` entries, and that each has the correct type and range.
#'
#' @param comp A named list with elements `dist`, `mean`, `sd`.
#' @param name Character. The parameter name (for error messages).
#' @return Invisibly returns `TRUE` if all checks pass. Throws an
#'   error otherwise.
#' @noRd
.validate_prior_component <- function(comp, name) {

    # Must be a list
    if (!is.list(comp)) {
        cli_abort(c(
            "Prior for {.val {name}} must be a named list.",
            "x" = "Got {.cls {class(comp)}}.",
            "i" = 'Expected: {.code list(dist = "normal", mean = 0, sd = 2)}.'
        ))
    }

    # Must have the required names
    required <- c("dist", "mean", "sd")
    missing_names <- setdiff(required, names(comp))
    if (length(missing_names) > 0L) {
        cli_abort(c(
            "Prior for {.val {name}} is missing required element{?s}.",
            "x" = "Missing: {.val {missing_names}}.",
            "i" = 'Required elements: {.val {required}}.'
        ))
    }

    # dist: character, currently only "normal"
    assert_string(comp$dist, .var.name = paste0(name, "$dist"))
    supported_dists <- c("normal")
    if (!comp$dist %in% supported_dists) {
        cli_abort(c(
            'Unsupported distribution for {.val {name}}.',
            "x" = "Got {.val {comp$dist}}.",
            "i" = "Currently supported: {.val {supported_dists}}."
        ))
    }

    # mean: finite numeric scalar (no NA, NaN, Inf)
    assert_number(comp$mean, finite = TRUE,
                  .var.name = paste0(name, "$mean"))

    # sd: positive finite numeric scalar
    assert_number(comp$sd, lower = .Machine$double.eps, finite = TRUE,
                  .var.name = paste0(name, "$sd"))

    invisible(TRUE)
}


# ============================================================================
# 5. .default_prior_raw --- Internal: raw defaults (no validation)
# ============================================================================

#' Raw default prior values (avoids infinite recursion)
#'
#' Returns the default prior hyperparameters as a plain list.
#' This is used internally by [hbb_prior()] to fill in unspecified
#' components before validation. It avoids the circular call
#' `hbb_prior() -> default_prior() -> hbb_prior()`.
#'
#' @return A plain list (not an S3 object) with default prior values.
#' @noRd
.default_prior_raw <- function() {
    list(
        alpha     = list(dist = "normal", mean = 0,   sd = 2),
        beta      = list(dist = "normal", mean = 0,   sd = 2),
        log_kappa = list(dist = "normal", mean = 2,   sd = 1.5),
        gamma     = list(dist = "normal", mean = 0,   sd = 1),
        tau       = list(dist = "normal", mean = 0,   sd = 1),
        lkj_eta   = 2
    )
}
