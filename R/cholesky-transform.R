# ============================================================================
# cholesky-transform.R --- Cholesky Posterior Recalibration for HBB Models
#
# Implements the Williams-Savitsky (2021, Theorem 4.1) affine Cholesky
# correction that maps prior-dominated MCMC draws to draws whose
# covariance equals the design-consistent sandwich variance V_sand.
#
# Part of the hurdlebb package.
#
# Contents:
#   1. cholesky_correct        --- Top-level function (exported)
#   2. compute_wald_ci         --- Standalone Wald CIs (exported)
#   3. print.hbb_cholesky      --- S3 print method (exported)
#   4. .validate_cholesky_inputs --- Input validation (internal)
#   5. .validate_level          --- Level validation (internal)
#   6. .make_stan_param_names   --- CmdStanR param names (internal)
#   7. .ensure_cholesky_pd      --- PD enforcement (internal)
#   8. .safe_chol_lower         --- Defensive Cholesky (internal)
#   9. .extract_cholesky_draws  --- Draw extraction (internal)
#  10. .verify_cholesky_transform --- Algebraic verification (internal)
#  11. .build_comparison_table  --- CI comparison table (internal)
# ============================================================================


# ============================================================================
# 1. cholesky_correct --- Top-level function (exported)
# ============================================================================

#' Cholesky Posterior Recalibration for Survey-Weighted Inference
#'
#' Applies the affine Cholesky correction of Williams and Savitsky (2021,
#' Theorem 4.1) to MCMC posterior draws, replacing the prior-dominated
#' posterior covariance with the design-consistent sandwich variance while
#' preserving the posterior mean.
#'
#' @description
#' In survey-weighted Bayesian models the pseudo-posterior covariance
#' \eqn{\Sigma_{\mathrm{MCMC}}} typically over-covers because the prior
#' contributes substantially to the curvature (Prior Inflation ratios of
#' 600--6500 are common for fixed effects).
#' The Cholesky correction constructs an affine map
#' \deqn{
#'   \theta^{*(m)} = \hat\theta + A\bigl(\theta^{(m)} - \hat\theta\bigr),
#'   \quad A = L_{\mathrm{sand}}\, L_{\mathrm{MCMC}}^{-1},
#' }
#' where \eqn{L_{\mathrm{sand}}} and \eqn{L_{\mathrm{MCMC}}} are the
#' lower-triangular Cholesky factors of the sandwich variance
#' \eqn{V_{\mathrm{sand}}} and the MCMC covariance
#' \eqn{\Sigma_{\mathrm{MCMC}}}, respectively.
#'
#' @section Mathematical properties:
#' The transformation satisfies three exact algebraic identities:
#' \enumerate{
#'   \item \strong{Mean preservation.}
#'     \eqn{E[\theta^*] = \hat\theta} because the affine map is centred
#'     at the posterior mean.
#'   \item \strong{Covariance recovery.}
#'     \deqn{
#'       \mathrm{Cov}(\theta^*) = A\,\Sigma_{\mathrm{MCMC}}\,A^\top
#'         = L_{\mathrm{sand}}\,L_{\mathrm{MCMC}}^{-1}\,
#'           \Sigma_{\mathrm{MCMC}}\,
#'           L_{\mathrm{MCMC}}^{-\top}\,L_{\mathrm{sand}}^\top
#'         = V_{\mathrm{sand}}.
#'     }
#'     The middle cancellation uses
#'     \eqn{\Sigma_{\mathrm{MCMC}} = L_{\mathrm{MCMC}}\,L_{\mathrm{MCMC}}^\top}.
#'   \item \strong{Affine equivariance.}
#'     Quantile ordering is preserved under monotone marginal
#'     transformations, so credible intervals from the corrected draws
#'     are valid frequentist confidence intervals with correct nominal
#'     coverage.
#' }
#'
#' @section Prior domination and shrinkage:
#' For parameters where the prior dominates the likelihood (Prior Inflation
#' \eqn{\mathrm{PI}_p = \Sigma_{\mathrm{MCMC},pp} / H_{\mathrm{obs},pp}^{-1}
#' \gg 1}), the diagonal of \eqn{A} is much less than unity
#' (typically 0.01--0.06).
#' This \emph{shrinks} the over-dispersed MCMC draws towards
#' \eqn{\hat\theta}, replacing prior-inflated credible intervals with
#' data-driven confidence intervals.
#' This is the expected and correct behaviour: the sandwich variance
#' \eqn{V_{\mathrm{sand}}} is smaller than both \eqn{\Sigma_{\mathrm{MCMC}}}
#' and \eqn{H_{\mathrm{obs}}^{-1}} when design effects are moderate.
#'
#' For parameters where the prior contributes little (e.g., \code{log_kappa}
#' with \eqn{\mathrm{PI} \approx 2.3}), \eqn{A_{pp} \approx 1.06},
#' producing slight inflation consistent with design-effect adjustment.
#'
#' @section Design Effect Ratio:
#' The DER is the key diagnostic for the correction magnitude:
#' \deqn{
#'   \mathrm{DER}_p = \frac{V_{\mathrm{sand},pp}}{H_{\mathrm{obs},pp}^{-1}}.
#' }
#' Values in the range 1--5 are typical for complex survey designs.
#'
#' @section Relationship to Wald inference:
#' Wald confidence intervals
#' \eqn{\hat\theta_p \pm z_{1-\alpha/2}\,\sqrt{V_{\mathrm{sand},pp}}}
#' are algebraically equivalent to the marginal quantiles of the corrected
#' draws in large samples.  They are recommended as the primary reporting
#' device because they do not depend on MCMC sampling variability.
#'
#' @param fit An object of class \code{"hbb_fit"} returned by
#'   \code{\link{hbb}} or similar model-fitting function.  Must contain
#'   a CmdStanR fit accessible via \code{fit$fit$draws()}.
#' @param sandwich An object of class \code{"hbb_sandwich"} returned by
#'   \code{\link{sandwich_variance}}.  Must contain \code{V_sand},
#'   \code{Sigma_MCMC}, \code{H_obs_inv}, \code{DER}, and
#'   \code{param_labels}.
#' @param level Numeric scalar in \eqn{(0,1)}.  Confidence level for
#'   interval construction.  Default is \code{0.95}.
#'
#' @return An S3 object of class \code{"hbb_cholesky"} containing:
#' \describe{
#'   \item{\code{theta_corrected}}{Numeric matrix of dimension M by D
#'     containing the Cholesky-corrected posterior draws.}
#'   \item{\code{theta_hat}}{Numeric vector of length D: posterior means
#'     (column means of the original draws).}
#'   \item{\code{A}}{The D by D affine transformation matrix
#'     \eqn{A = L_{\mathrm{sand}}\,L_{\mathrm{MCMC}}^{-1}}.}
#'   \item{\code{L_MCMC}}{Lower-triangular Cholesky factor of the MCMC
#'     posterior covariance.}
#'   \item{\code{L_sand}}{Lower-triangular Cholesky factor of the sandwich
#'     variance.}
#'   \item{\code{comparison_table}}{Data frame with one row per parameter
#'     containing naive, corrected, and Wald confidence intervals, width
#'     ratios, DER, DER-vs-MCMC, prior inflation, and the square root
#'     of DER.}
#'   \item{\code{verification}}{Named list of numerical verification
#'     checks: mean preservation error, relative variance error,
#'     A algebraic identity error, and logical pass/fail flags.}
#'   \item{\code{level}}{The confidence level used.}
#'   \item{\code{D}}{Integer: number of parameters (2P + 1).}
#'   \item{\code{M}}{Integer: number of MCMC draws.}
#'   \item{\code{P}}{Integer: number of covariates per margin.}
#'   \item{\code{pd_corrections}}{List of PD correction details for each
#'     matrix (whether ridge or nearPD was applied, condition numbers).}
#' }
#'
#' @seealso
#' \code{\link{sandwich_variance}} for the sandwich variance computation,
#' \code{\link{compute_wald_ci}} for standalone Wald confidence intervals,
#' \code{\link{print.hbb_cholesky}} for the print method.
#'
#' @references
#' Williams, M. R. and Savitsky, T. D. (2021).
#' Uncertainty estimation for pseudo-Bayesian inference under complex
#' sampling.
#' \emph{International Statistical Review}, \strong{89}(1), 72--107.
#' \doi{10.1111/insr.12376}
#'
#' Ghosal, R., Ghosh, S. K., and Maiti, T. (2020).
#' Two-part regression models for longitudinal zero-inflated count data.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' \strong{183}(4), 1603--1626.
#'
#' @examples
#' \dontrun{
#' # After fitting and computing sandwich variance:
#' sand <- sandwich_variance(fit)
#' chol_obj <- cholesky_correct(fit, sand, level = 0.95)
#' print(chol_obj)
#'
#' # Compare interval widths
#' chol_obj$comparison_table[, c("parameter", "naive_width",
#'                               "corrected_width", "wald_width")]
#' }
#'
#' @export
cholesky_correct <- function(fit, sandwich, level = 0.95) {

    # ---- 0. Input validation -----------------------------------
    .validate_cholesky_inputs(fit, sandwich)
    .validate_level(level)

    # ---- 1. Extract dimensions and parameter names -----------------------
    P <- sandwich$P
    D <- sandwich$D
    param_names  <- .make_stan_param_names(P)
    param_labels <- sandwich$param_labels

    # Fallback labels if missing or wrong length
    if (is.null(param_labels) || length(param_labels) != D) {
        cli::cli_alert_warning(
            "param_labels missing or wrong length; using CmdStanR names."
        )
        param_labels <- param_names
    }

    cli::cli_alert_info(
        "Cholesky correction: D = {D}, P = {P}, level = {level}"
    )

    # ---- 2. Extract MCMC draws (M x D) ------------------------
    theta_draws <- .extract_cholesky_draws(fit, param_names)
    M <- nrow(theta_draws)
    cli::cli_alert_info("Extracted {M} MCMC draws")

    # ---- 3. Posterior mean -----------------------------------------------
    theta_hat <- colMeans(theta_draws)
    names(theta_hat) <- param_labels

    # ---- 4. Retrieve and validate covariance matrices --------------------
    Sigma_MCMC <- sandwich$Sigma_MCMC
    V_sand     <- sandwich$V_sand
    H_obs_inv  <- sandwich$H_obs_inv

    # ---- 5. Ensure positive definiteness -------------
    pd_corrections <- list()

    # Sigma_MCMC: ridge method (preserves matrix structure)
    sigma_pd <- .ensure_cholesky_pd(Sigma_MCMC, "Sigma_MCMC", method = "ridge")
    Sigma_MCMC <- sigma_pd$mat
    pd_corrections$Sigma_MCMC <- sigma_pd[c("corrected", "details",
                                              "min_eig", "cond_number")]

    # V_sand: nearPD method (handles complex non-PD from clustering)
    vsand_pd <- .ensure_cholesky_pd(V_sand, "V_sand", method = "nearpd")
    V_sand <- vsand_pd$mat
    pd_corrections$V_sand <- vsand_pd[c("corrected", "details",
                                          "min_eig", "cond_number")]

    # ---- 6. Cholesky decomposition (safe wrapper) ----------------
    L_MCMC <- .safe_chol_lower(Sigma_MCMC, "Sigma_MCMC")
    L_sand <- .safe_chol_lower(V_sand, "V_sand")

    # ---- 7. A = L_sand * L_MCMC^{-1} via forwardsolve ---------
    # forwardsolve exploits known lower-triangular structure: O(D^2)
    L_MCMC_inv <- forwardsolve(L_MCMC, diag(D))
    A <- L_sand %*% L_MCMC_inv

    # ---- 8. Extreme shrinkage detection ------------------------
    A_diag <- diag(A)
    names(A_diag) <- param_labels
    extreme_idx <- which(abs(A_diag) < 0.001)
    if (length(extreme_idx) > 0L) {
        extreme_labs <- param_labels[extreme_idx]
        extreme_vals <- formatC(A_diag[extreme_idx], format = "e", digits = 3)
        cli::cli_alert_warning(
            "Extreme shrinkage (|A_dd| < 0.001) for {length(extreme_idx)} parameter(s)."
        )
        for (k in seq_along(extreme_idx)) {
            cli::cli_alert_warning(
                "  {extreme_labs[k]}: A[{extreme_idx[k]},{extreme_idx[k]}] = {extreme_vals[k]}"
            )
        }
    }

    n_shrink <- sum(A_diag < 1)
    n_expand <- sum(A_diag > 1)
    cli::cli_alert_info(
        "A diagonal: {n_shrink} shrink (< 1), {n_expand} expand (> 1). Range [{signif(min(A_diag), 3)}, {signif(max(A_diag), 3)}]."
    )

    # ---- 9. Vectorised affine transformation -------------------
    # theta_corrected[m,] = theta_hat + A %*% (theta[m,] - theta_hat)
    # sweep + %*% is O(M*D^2) with no R-level loops over M
    theta_centered  <- sweep(theta_draws, 2, theta_hat, "-")
    theta_corrected <- sweep(theta_centered %*% t(A), 2, theta_hat, "+")
    colnames(theta_corrected) <- param_labels

    # ---- 10. Verification checks -------------------------------
    verification <- .verify_cholesky_transform(
        theta_corrected, theta_hat, A, Sigma_MCMC, V_sand
    )

    if (verification$mean_pass) {
        cli::cli_alert_success(
            "Mean preservation: PASS (error = {format(verification$mean_preservation_error, digits = 3)})"
        )
    } else {
        cli::cli_alert_warning(
            "Mean preservation: FAIL (error = {format(verification$mean_preservation_error, digits = 3)})"
        )
    }
    if (verification$A_pass) {
        cli::cli_alert_success(
            "Algebraic identity A*Sigma*A' = V_sand: PASS (error = {format(verification$A_verification_error, digits = 3)})"
        )
    } else {
        cli::cli_alert_warning(
            "Algebraic identity A*Sigma*A': FAIL (error = {format(verification$A_verification_error, digits = 3)})"
        )
    }

    # ---- 11. Comparison table ----------------------------------
    comparison_table <- .build_comparison_table(
        theta_hat       = theta_hat,
        theta_draws     = theta_draws,
        theta_corrected = theta_corrected,
        V_sand          = V_sand,
        Sigma_MCMC      = Sigma_MCMC,
        H_obs_inv       = H_obs_inv,
        param_labels    = param_labels,
        level           = level
    )

    # ---- 12. Construct return object -------------------------------------
    cli::cli_alert_success("Cholesky correction complete")

    structure(
        list(
            theta_corrected  = theta_corrected,
            theta_hat        = theta_hat,
            A                = A,
            L_MCMC           = L_MCMC,
            L_sand           = L_sand,
            comparison_table = comparison_table,
            verification     = verification,
            level            = level,
            D                = D,
            M                = M,
            P                = P,
            pd_corrections   = pd_corrections
        ),
        class = "hbb_cholesky"
    )
}


# ============================================================================
# 2. compute_wald_ci --- Standalone Wald CIs (exported)
# ============================================================================

#' Wald Confidence Intervals from Sandwich Variance
#'
#' Computes Wald confidence intervals using the sandwich variance
#' \eqn{V_{\mathrm{sand}}}.  This is a standalone function that does not
#' require a fitted model object or MCMC draws.
#'
#' @description
#' The Wald confidence interval for parameter \eqn{\theta_p} is
#' \deqn{
#'   \mathrm{CI}_{1-\alpha}(\theta_p)
#'     = \hat\theta_p \pm z_{1-\alpha/2}\,
#'       \sqrt{V_{\mathrm{sand},pp}},
#' }
#' where \eqn{z_{1-\alpha/2}} is the standard normal quantile.
#' Because the sandwich variance is design-consistent under mild regularity
#' conditions (Williams and Savitsky, 2021, Theorem 3.2), Wald intervals
#' achieve correct frequentist coverage for the pseudo-true parameter.
#'
#' @section Recommendation:
#' Wald intervals are recommended as the primary inference device because:
#' \enumerate{
#'   \item They depend only on the point estimate and
#'     \eqn{V_{\mathrm{sand}}}, not on MCMC sampling variability.
#'   \item They are algebraically equivalent to the marginal quantiles of
#'     the Cholesky-corrected draws in large samples.
#'   \item They avoid the non-trivial Monte Carlo error in tail quantile
#'     estimation from finite MCMC samples.
#' }
#'
#' @section Significance testing:
#' The Wald z-statistic is
#' \eqn{z_p = \hat\theta_p / \sqrt{V_{\mathrm{sand},pp}}},
#' with two-sided p-value \eqn{2\,\Phi(-|z_p|)}.
#' A parameter is flagged significant when the p-value falls below
#' \eqn{1 - \mathrm{level}}.
#'
#' @param theta_hat Numeric vector of length D containing point estimates
#'   (typically posterior means).
#' @param V_sand Numeric symmetric positive-(semi)definite matrix of
#'   dimension D by D, the sandwich variance.  Only the diagonal entries
#'   are used for marginal Wald intervals.
#' @param level Numeric scalar in \eqn{(0,1)}.
#'   Confidence level.  Default is \code{0.95}.
#' @param param_labels Optional character vector of length D giving
#'   parameter names.  If \code{NULL} (default), names are taken from
#'   \code{names(theta_hat)} or generated as \code{param_1}, \code{param_2},
#'   etc.
#'
#' @return A data frame with one row per parameter and columns:
#' \describe{
#'   \item{\code{parameter}}{Character: parameter label.}
#'   \item{\code{post_mean}}{Numeric: point estimate.}
#'   \item{\code{se}}{Numeric: sandwich standard error.}
#'   \item{\code{z_stat}}{Numeric: Wald z-statistic.}
#'   \item{\code{p_value}}{Numeric: two-sided p-value.}
#'   \item{\code{ci_lo}}{Numeric: lower confidence limit.}
#'   \item{\code{ci_hi}}{Numeric: upper confidence limit.}
#'   \item{\code{ci_width}}{Numeric: interval width.}
#'   \item{\code{significant}}{Logical: TRUE if p-value is below
#'     the nominal significance level.}
#' }
#'
#' @seealso
#' \code{\link{cholesky_correct}} for the full Cholesky recalibration,
#' \code{\link{sandwich_variance}} for obtaining \code{V_sand}.
#'
#' @references
#' Williams, M. R. and Savitsky, T. D. (2021).
#' Uncertainty estimation for pseudo-Bayesian inference under complex
#' sampling.
#' \emph{International Statistical Review}, \strong{89}(1), 72--107.
#' \doi{10.1111/insr.12376}
#'
#' @examples
#' # Minimal example with known values
#' theta_hat <- c(alpha_1 = -0.324, beta_1 = 0.090, log_kappa = 1.92)
#' V_sand <- diag(c(0.0027, 0.00035, 0.0041))
#' compute_wald_ci(theta_hat, V_sand, level = 0.95)
#'
#' @export
compute_wald_ci <- function(theta_hat, V_sand, level = 0.95,
                            param_labels = NULL) {

    # ---- Validation --------------------------------------------
    if (!is.numeric(theta_hat) || !is.vector(theta_hat)) {
        cli::cli_abort("{.arg theta_hat} must be a numeric vector.")
    }
    D <- length(theta_hat)

    if (!is.matrix(V_sand) || !is.numeric(V_sand)) {
        cli::cli_abort("{.arg V_sand} must be a numeric matrix.")
    }
    if (nrow(V_sand) != D || ncol(V_sand) != D) {
        cli::cli_abort(
            "{.arg V_sand} must be {D} x {D} to match {.arg theta_hat}."
        )
    }
    .validate_level(level)

    # Parameter labels
    if (is.null(param_labels)) {
        param_labels <- if (!is.null(names(theta_hat))) {
            names(theta_hat)
        } else {
            paste0("param_", seq_len(D))
        }
    }
    if (length(param_labels) != D) {
        cli::cli_abort(
            "{.arg param_labels} must have length {D} to match {.arg theta_hat}."
        )
    }

    # ---- Diagonal safety ---------------------------------------
    v_diag <- diag(V_sand)
    bad_idx <- which(!is.na(v_diag) & v_diag <= 0)
    if (length(bad_idx) > 0L) {
        bad_labs <- param_labels[bad_idx]
        cli::cli_alert_warning(
            "V_sand diagonal has {length(bad_idx)} non-positive value(s): {.val {bad_labs}}. SEs clamped to 0."
        )
    }

    # ---- Computation -------------------------------------------
    z_crit  <- qnorm((1 + level) / 2)
    se      <- sqrt(pmax(v_diag, 0))
    z_stat  <- theta_hat / se
    p_value <- 2 * pnorm(-abs(z_stat))

    ci_lo <- theta_hat - z_crit * se
    ci_hi <- theta_hat + z_crit * se

    data.frame(
        parameter   = param_labels,
        post_mean   = as.numeric(theta_hat),
        se          = se,
        z_stat      = z_stat,
        p_value     = p_value,
        ci_lo       = ci_lo,
        ci_hi       = ci_hi,
        ci_width    = ci_hi - ci_lo,
        significant = p_value < (1 - level),
        stringsAsFactors = FALSE,
        row.names   = NULL
    )
}


# ============================================================================
# 3. print.hbb_cholesky --- S3 print method (exported)
# ============================================================================

#' Print Method for Cholesky-Corrected Posterior Objects
#'
#' Displays a structured summary of the Cholesky posterior recalibration,
#' including the transformation matrix diagonal (with shrinkage/inflation
#' labels), verification checks, DER summary, and an abbreviated comparison
#' table.
#'
#' @param x An object of class \code{"hbb_cholesky"} returned by
#'   \code{\link{cholesky_correct}}.
#' @param digits Integer; number of significant digits for numeric output.
#'   Default is 4.
#' @param ... Additional arguments passed to \code{\link[base]{format}}.
#'
#' @return The object \code{x}, invisibly.
#'
#' @seealso \code{\link{cholesky_correct}}
#'
#' @export
print.hbb_cholesky <- function(x, digits = 4, ...) {

    # Full tryCatch for defensive operation on corrupted objects
    tryCatch({

        cat("\n")
        cat(strrep("=", 60), "\n")
        cat("  Cholesky Posterior Recalibration (Williams-Savitsky 2021)\n")
        cat(strrep("=", 60), "\n\n")

        # Dimensions
        cat("  Parameters (D):", x$D, "\n")
        cat("  MCMC draws (M):", x$M, "\n")
        cat("  Covariates (P):", x$P, "\n")
        cat("  Confidence level:", x$level, "\n")

        # PD corrections
        if (!is.null(x$pd_corrections)) {
            sigma_corr <- isTRUE(x$pd_corrections$Sigma_MCMC$corrected)
            vsand_corr <- isTRUE(x$pd_corrections$V_sand$corrected)
            if (sigma_corr || vsand_corr) {
                cat("\n  PD Corrections Applied:\n")
                if (sigma_corr) {
                    cat("    Sigma_MCMC:", x$pd_corrections$Sigma_MCMC$details, "\n")
                }
                if (vsand_corr) {
                    cat("    V_sand:", x$pd_corrections$V_sand$details, "\n")
                }
            }
        }
        cat("\n")

        # A diagonal with shrink/inflate labels
        cat(strrep("-", 50), "\n")
        cat("  Transformation Matrix A (diagonal)\n")
        cat(strrep("-", 50), "\n")
        a_diag <- diag(x$A)
        labels <- if (!is.null(x$comparison_table)) {
            x$comparison_table$parameter
        } else {
            paste0("param_", seq_along(a_diag))
        }
        for (i in seq_along(a_diag)) {
            direction <- if (a_diag[i] < 0.99) {
                "SHRINK"
            } else if (a_diag[i] > 1.01) {
                "inflate"
            } else {
                "~unity"
            }
            cat(sprintf("    %-20s  %s  [%s]\n",
                        labels[i],
                        format(a_diag[i], digits = digits),
                        direction))
        }
        cat("\n")

        # Verification checks
        cat(strrep("-", 50), "\n")
        cat("  Verification Checks\n")
        cat(strrep("-", 50), "\n")
        v <- x$verification
        .print_cholesky_check("Mean preservation",
                               v$mean_pass,
                               format(v$mean_preservation_error, digits = 3))
        .print_cholesky_check("Variance recovery",
                               v$variance_pass,
                               format(v$variance_relative_error, digits = 3))
        .print_cholesky_check("A algebraic identity",
                               v$A_pass,
                               format(v$A_verification_error, digits = 3))
        cat("\n")

        # Abbreviated comparison table
        cat(strrep("-", 50), "\n")
        cat("  CI Comparison (width_ratio = corrected / naive)\n")
        cat(strrep("-", 50), "\n")
        if (!is.null(x$comparison_table)) {
            tbl <- x$comparison_table[, c("parameter", "post_mean",
                                           "naive_width", "corrected_width",
                                           "wald_width", "width_ratio", "DER")]
            tbl[, -1] <- lapply(tbl[, -1], function(col) {
                format(round(col, digits), nsmall = min(digits, 4))
            })
            print(tbl, row.names = FALSE, right = TRUE)
        }
        cat("\n")

        # DER summary
        if (!is.null(x$comparison_table)) {
            DER <- x$comparison_table$DER
            cat(sprintf("  DER range: [%.3f, %.3f], mean: %.3f\n",
                        min(DER), max(DER), mean(DER)))
        }
        cat("\n")

    }, error = function(e) {
        # Defensive fallback
        cat("Error printing hbb_cholesky object:", conditionMessage(e), "\n")
        cat("Object class:", paste(class(x), collapse = ", "), "\n")
        cat("Object names:", paste(names(x), collapse = ", "), "\n")
    })

    invisible(x)
}


# ============================================================================
# Internal helpers
# ============================================================================


# ---- 4. .validate_cholesky_inputs ----------------------------------------

#' Validate inputs for cholesky_correct
#'
#' @param fit An hbb_fit object.
#' @param sandwich An hbb_sandwich object.
#' @return Invisible NULL; aborts on validation failure.
#' @keywords internal
.validate_cholesky_inputs <- function(fit, sandwich) {

    # Class checks
    if (!inherits(fit, "hbb_fit")) {
        cli::cli_abort(c(
            "x" = "{.arg fit} must be an object of class {.cls hbb_fit}.",
            "i" = "Received class: {.cls {class(fit)}}.",
            "i" = "Use {.fn hbb} to create a fitted model object."
        ))
    }
    if (!inherits(sandwich, "hbb_sandwich")) {
        cli::cli_abort(c(
            "x" = "{.arg sandwich} must be an object of class {.cls hbb_sandwich}.",
            "i" = "Received class: {.cls {class(sandwich)}}.",
            "i" = "Use {.fn sandwich_variance} to create this object."
        ))
    }

    # Required fields in sandwich
    required <- c("V_sand", "Sigma_MCMC", "H_obs_inv", "DER", "D", "P")
    missing_fields <- setdiff(required, names(sandwich))
    if (length(missing_fields) > 0L) {
        cli::cli_abort(
            "{.arg sandwich} is missing required field{?s}: {.field {missing_fields}}."
        )
    }

    # D = 2P + 1 consistency
    expected_D <- 2L * sandwich$P + 1L
    if (sandwich$D != expected_D) {
        cli::cli_abort(c(
            "x" = "Sandwich dimension mismatch: D = {sandwich$D} but 2*P+1 = {expected_D}.",
            "i" = "P = {sandwich$P}."
        ))
    }

    # P consistency between fit and sandwich
    P_fit <- fit$hbb_data$P
    P_sand <- sandwich$P
    if (!is.null(P_fit) && P_fit != P_sand) {
        cli::cli_abort(
            "Dimension mismatch: fit has P = {P_fit} but sandwich has P = {P_sand}."
        )
    }

    invisible(NULL)
}


# ---- 5. .validate_level -------------------------------------------------

#' Validate confidence level parameter
#'
#' @param level Confidence level to validate.
#' @return Invisible TRUE; aborts on failure.
#' @keywords internal
.validate_level <- function(level) {
    if (!is.numeric(level) || length(level) != 1L || is.na(level)) {
        cli::cli_abort(c(
            "x" = "{.arg level} must be a single numeric value in (0, 1).",
            "i" = "Received: {.val {level}} of type {.cls {class(level)}}."
        ))
    }
    if (level <= 0 || level >= 1) {
        cli::cli_abort(c(
            "x" = "{.arg level} must be strictly in (0, 1).",
            "i" = "Received: {.val {level}}.",
            "i" = "Common choices: 0.90, 0.95, 0.99."
        ))
    }
    invisible(TRUE)
}


# ---- 6. .make_stan_param_names ------------------------------------------

#' Generate CmdStanR parameter names
#'
#' Produces \code{c("alpha[1]", ..., "alpha[P]", "beta[1]", ...,
#' "beta[P]", "log_kappa")} matching CmdStanR naming conventions.
#'
#' @param P Integer; number of covariates per margin.
#' @return Character vector of length 2P + 1.
#' @keywords internal
.make_stan_param_names <- function(P) {
    c(
        paste0("alpha[", seq_len(P), "]"),
        paste0("beta[", seq_len(P), "]"),
        "log_kappa"
    )
}


# ---- 7. .ensure_cholesky_pd ---------------------------------------------

#' Ensure a matrix is positive definite for Cholesky decomposition
#'
#' Checks positive definiteness and applies corrections if needed.
#' For Sigma_MCMC, uses ridge regularisation (additive diagonal perturbation).
#' For V_sand, uses nearPD projection (Higham 2002).
#'
#' @param mat Symmetric matrix.
#' @param mat_name Character label for diagnostic messages.
#' @param method Character: \code{"ridge"} adds a diagonal ridge,
#'   \code{"nearpd"} projects via \code{Matrix::nearPD}.
#' @return A list with elements: \code{mat} (corrected matrix),
#'   \code{corrected} (logical), \code{details} (character or NULL),
#'   \code{min_eig}, \code{max_eig}, \code{cond_number}.
#' @keywords internal
.ensure_cholesky_pd <- function(mat, mat_name = "matrix",
                                 method = c("ridge", "nearpd")) {
    method <- match.arg(method)

    # Force exact symmetry
    mat <- (mat + t(mat)) / 2

    eig <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
    min_eig <- min(eig)
    max_eig <- max(eig)
    cond <- if (min_eig > 0) max_eig / min_eig else Inf

    result <- list(mat = mat, corrected = FALSE, details = NULL,
                   min_eig = min_eig, max_eig = max_eig, cond_number = cond)

    # PD check
    if (min_eig <= 0) {
        if (method == "ridge") {
            ridge <- max(2 * abs(min_eig),
                         .Machine$double.eps * norm(mat, "F"))
            mat <- mat + ridge * diag(nrow(mat))
            detail_msg <- sprintf(
                "%s not PD (min eig = %.3e). Ridge = %.3e applied.",
                mat_name, min_eig, ridge
            )
        } else {
            mat <- as.matrix(Matrix::nearPD(mat, corr = FALSE)$mat)
            detail_msg <- sprintf(
                "%s not PD (min eig = %.3e). nearPD projection applied.",
                mat_name, min_eig
            )
        }
        cli::cli_alert_warning(detail_msg)

        # Recompute eigenvalues
        eig_new <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
        result$min_eig <- min(eig_new)
        result$max_eig <- max(eig_new)
        result$cond_number <- max(eig_new) / min(eig_new)
        result$mat <- mat
        result$corrected <- TRUE
        result$details <- detail_msg
    }

    # Condition number warning
    if (is.finite(result$cond_number) && result$cond_number > 1e10) {
        cli::cli_alert_warning(
            "{mat_name} condition number = {formatC(result$cond_number, format = 'e', digits = 2)}. Numerical instability possible."
        )
    }

    result
}


# ---- 8. .safe_chol_lower ------------------------------------------------

#' Safely compute lower-triangular Cholesky factor
#'
#' Returns \code{t(chol(mat))}, i.e., the lower-triangular factor L
#' such that \code{mat = L \%*\% t(L)}.  Wraps in tryCatch for
#' informative error messages.
#'
#' @param mat A positive-definite matrix.
#' @param mat_name Name for error messages.
#' @return Lower-triangular Cholesky factor.
#' @keywords internal
.safe_chol_lower <- function(mat, mat_name = "matrix") {
    tryCatch(
        t(chol(mat)),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "Cholesky decomposition of {.arg {mat_name}} failed.",
                "i" = "The matrix may not be numerically positive definite even after PD safeguards.",
                "i" = "Original error: {conditionMessage(e)}"
            ))
        }
    )
}


# ---- 9. .extract_cholesky_draws -----------------------------------------

#' Extract posterior draws from an hbb_fit for Cholesky correction
#'
#' Validates that \code{fit$fit} exists and draws match expected dimensions.
#'
#' @param fit An hbb_fit object.
#' @param param_names Character vector of CmdStanR parameter names.
#' @return Numeric matrix of draws (M x D).
#' @keywords internal
.extract_cholesky_draws <- function(fit, param_names) {
    if (is.null(fit$fit)) {
        cli::cli_abort(c(
            "x" = "{.code fit$fit} is {.val NULL}.",
            "i" = "The {.cls hbb_fit} object must contain a valid CmdStanR fit",
            "i" = "in its {.field $fit} slot."
        ))
    }

    draws <- tryCatch(
        fit$fit$draws(variables = param_names, format = "matrix"),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "Failed to extract draws from {.cls hbb_fit} object.",
                "i" = "Requested parameters: {.val {param_names}}",
                "i" = "Error: {conditionMessage(e)}"
            ))
        }
    )

    if (!is.matrix(draws)) draws <- as.matrix(draws)

    if (ncol(draws) != length(param_names)) {
        cli::cli_abort(c(
            "x" = "Draw matrix has {ncol(draws)} columns but expected {length(param_names)}.",
            "i" = "Parameters: {.val {param_names}}"
        ))
    }

    draws
}


# ---- 10. .verify_cholesky_transform -------------------------------------

#' Verify the Cholesky transform satisfies algebraic identities
#'
#' Checks three properties of the affine correction:
#' mean preservation, finite-sample variance recovery, and the
#' algebraic identity A Sigma A' = V_sand.
#'
#' @param theta_corrected Corrected M by D draw matrix.
#' @param theta_hat D-vector of posterior means.
#' @param A D by D transformation matrix.
#' @param Sigma_MCMC D by D MCMC covariance.
#' @param V_sand D by D sandwich variance.
#' @return Named list with numeric diagnostics and logical pass flags.
#' @keywords internal
.verify_cholesky_transform <- function(theta_corrected, theta_hat,
                                        A, Sigma_MCMC, V_sand) {

    # 1. Mean preservation (exact up to floating point)
    mean_corrected <- colMeans(theta_corrected)
    mean_error <- max(abs(mean_corrected - theta_hat))

    # 2. Finite-sample variance recovery
    Cov_corrected <- cov(theta_corrected)
    safe_denom <- pmax(abs(diag(V_sand)), .Machine$double.eps)
    var_rel_error <- max(abs(diag(Cov_corrected) - diag(V_sand)) / safe_denom)

    # 3. Algebraic identity: A Sigma_MCMC A' = V_sand
    A_Sigma_At <- A %*% Sigma_MCMC %*% t(A)
    A_error <- norm(A_Sigma_At - V_sand, "F") /
        max(norm(V_sand, "F"), .Machine$double.eps)

    list(
        mean_preservation_error = mean_error,
        variance_relative_error = var_rel_error,
        A_verification_error    = A_error,
        mean_pass     = mean_error < 1e-10,
        variance_pass = var_rel_error < 0.10,
        A_pass        = A_error < 1e-10
    )
}


# ---- 11. .build_comparison_table ----------------------------------------

#' Build a comparison table of naive, corrected, and Wald CIs
#'
#' @param theta_hat Numeric vector length D: posterior means.
#' @param theta_draws M x D matrix: original MCMC draws.
#' @param theta_corrected M x D matrix: corrected draws.
#' @param V_sand D x D sandwich variance.
#' @param Sigma_MCMC D x D MCMC posterior covariance.
#' @param H_obs_inv D x D inverse observed Fisher information.
#' @param param_labels Character vector length D.
#' @param level Confidence level.
#' @return A data.frame with one row per parameter and 16 columns.
#' @keywords internal
.build_comparison_table <- function(theta_hat, theta_draws, theta_corrected,
                                     V_sand, Sigma_MCMC, H_obs_inv,
                                     param_labels, level) {
    D <- length(theta_hat)
    alpha <- 1 - level
    z_crit <- qnorm(1 - alpha / 2)

    # Quantile probabilities
    probs <- c(alpha / 2, 1 - alpha / 2)

    # Naive CIs from raw MCMC quantiles (vectorized via apply)
    naive_q <- apply(theta_draws, 2, quantile, probs = probs)
    naive_lo    <- naive_q[1, ]
    naive_hi    <- naive_q[2, ]
    naive_width <- naive_hi - naive_lo

    # Corrected CIs
    corr_q <- apply(theta_corrected, 2, quantile, probs = probs)
    corrected_lo    <- corr_q[1, ]
    corrected_hi    <- corr_q[2, ]
    corrected_width <- corrected_hi - corrected_lo

    # Wald CIs (analytic)
    se_sand    <- sqrt(diag(V_sand))
    wald_lo    <- theta_hat - z_crit * se_sand
    wald_hi    <- theta_hat + z_crit * se_sand
    wald_width <- 2 * z_crit * se_sand

    # Diagnostic ratios
    diag_V    <- diag(V_sand)
    diag_S    <- diag(Sigma_MCMC)
    diag_H    <- diag(H_obs_inv)

    width_ratio     <- corrected_width / naive_width
    DER             <- diag_V / diag_H
    DER_vs_MCMC     <- diag_V / diag_S
    prior_inflation <- diag_S / diag_H
    sqrt_DER        <- sqrt(DER)

    data.frame(
        parameter       = param_labels,
        post_mean       = unname(as.numeric(theta_hat)),
        naive_lo        = unname(naive_lo),
        naive_hi        = unname(naive_hi),
        naive_width     = unname(naive_width),
        corrected_lo    = unname(corrected_lo),
        corrected_hi    = unname(corrected_hi),
        corrected_width = unname(corrected_width),
        wald_lo         = unname(wald_lo),
        wald_hi         = unname(wald_hi),
        wald_width      = unname(wald_width),
        width_ratio     = unname(width_ratio),
        DER             = unname(DER),
        DER_vs_MCMC     = unname(DER_vs_MCMC),
        prior_inflation = unname(prior_inflation),
        sqrt_DER        = unname(sqrt_DER),
        stringsAsFactors = FALSE,
        row.names        = NULL
    )
}


# ---- 12. .print_cholesky_check ------------------------------------------

#' Print a verification check line (internal helper for print method)
#'
#' @param label Check label.
#' @param pass Logical pass/fail.
#' @param value Formatted numeric value.
#' @keywords internal
.print_cholesky_check <- function(label, pass, value) {
    status <- if (isTRUE(pass)) "PASS" else "FAIL"
    cat(sprintf("    [%s] %s (error = %s)\n", status, label, value))
}
