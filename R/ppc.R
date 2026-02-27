# ============================================================================
# ppc.R --- Posterior Predictive Checks for Hurdle Beta-Binomial Models
#
# Implements graphical and numerical posterior predictive checks (PPC) for
# the HBB model, following the framework of Gabry et al. (2019). Two
# complementary test statistics are computed:
#
#   1. Zero rate: T_0(y) = mean(y == 0)
#      The proportion of structural zeros in the observed or replicated
#      data.  Under the model, E[T_0(y_rep)] = 1 - mean(q_i).
#
#   2. IT share: T_s(y) = mean(y[y > 0] / n[y > 0])
#      The mean enrollment share conditional on participation. Under
#      the model, E[T_s(y_rep) | y > 0] approx mean(mu_i).
#
# These two statistics are motivated by the decomposition of the HBB
# likelihood into an extensive margin and an intensive margin: T_0 targets
# the extensive margin and T_s targets the intensive margin.  Both must
# fall within the posterior predictive distribution for the model to be
# considered adequately calibrated.
#
# Two extraction methods are supported:
#   "stan"     : Extracts y_rep[N] draws from the Stan GQ block. This is
#                the primary method and uses the actual joint posterior.
#                If y_rep draws are unavailable, emits a cli_warn and
#                automatically falls back to "simulate".
#   "simulate" : Regenerates y_rep in R using the hurdle Beta-Binomial
#                hierarchy, composing MCMC draws of (alpha, beta, log_kappa)
#                with rbetabinom()/rztbetabinom(). Guards against NaN from
#                plogis, Inf from exp(log_kappa), and non-finite parameters.
#
# References:
#   Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A.
#   (2019). Visualisation in Bayesian workflow. J. R. Stat. Soc. A,
#   182(2), 389-402.
#
#   Ghosal, R., Ghosh, S. K., and Maiti, T. (2020). Two-part regression
#   models for longitudinal zero-inflated count data. Journal of the Royal
#   Statistical Society: Series A, 183(4), 1603-1626.
#
# Contents:
#   1. ppc                    --- Main PPC function (exported)
#   2. print.hbb_ppc          --- Defensive S3 print method (exported)
#   3. plot.hbb_ppc           --- S3 plot method (exported)
#   4. .ppc_validate_inputs   --- Comprehensive validation + fallback (internal)
#   5. .ppc_extract_yrep      --- Extract y_rep from Stan GQ (internal)
#   6. .ppc_simulate_yrep     --- Simulate y_rep in R (internal)
#   7. .ppc_compute_stats     --- Compute T_0 and T_s per draw (internal)
#   8. .ppc_summarize_draws   --- Summary statistics for draws (internal)
# ============================================================================


# ============================================================================
# 1. ppc --- Main PPC function (exported)
# ============================================================================

#' Posterior Predictive Checks for Hurdle Beta-Binomial Models
#'
#' Computes numerical posterior predictive checks (PPC) for a fitted
#' \code{hbb_fit} model by comparing observed test statistics to their
#' posterior predictive distributions.
#'
#' @description
#' The function evaluates two test statistics designed to probe the two
#' structural components of the hurdle Beta-Binomial model:
#'
#' \describe{
#'   \item{Zero rate}{
#'     \deqn{T_0(\mathbf{y}) = \frac{1}{N} \sum_{i=1}^{N}
#'       \mathbf{1}(y_i = 0).}
#'     Targets the extensive margin: whether a center serves
#'     infant/toddler children at all.  Under the fitted model,
#'     \eqn{E[T_0(\mathbf{y}^{\mathrm{rep}})] = 1 - N^{-1}\sum_i q_i}.
#'   }
#'   \item{IT share}{
#'     \deqn{T_s(\mathbf{y}) = \frac{1}{N^+}
#'       \sum_{i:\, y_i > 0} \frac{y_i}{n_i},}
#'     where \eqn{N^+ = \#\{i : y_i > 0\}} is the count of participating
#'     providers.  Targets the intensive margin: the enrollment share
#'     conditional on participation.  Under the fitted model,
#'     \eqn{E[T_s(\mathbf{y}^{\mathrm{rep}}) \mid y^{\mathrm{rep}} > 0]
#'     \approx N^{-1}\sum_i \mu_i}.
#'   }
#' }
#'
#' @section Theory and Motivation:
#' A well-calibrated Bayesian model should produce posterior predictive
#' distributions \eqn{p(\mathbf{y}^{\mathrm{rep}} \mid \mathbf{y})} that
#' are consistent with the observed data.  Formally, the Bayesian p-value
#' for a test statistic \eqn{T} is:
#' \deqn{
#'   p_B = \Pr\bigl(T(\mathbf{y}^{\mathrm{rep}}) \ge T(\mathbf{y})
#'         \mid \mathbf{y}\bigr),
#' }
#' which should be near \eqn{0.5} for a well-specified model and near
#' \eqn{0} or \eqn{1} for systematic misspecification (Gelman et al.,
#' 1996).
#'
#' The zero rate \eqn{T_0} and the IT share \eqn{T_s} are chosen because
#' they directly correspond to the two structural parts of the hurdle
#' model.  A defect in \eqn{T_0} indicates misspecification of the
#' Bernoulli participation equation; a defect in \eqn{T_s} indicates
#' misspecification of the zero-truncated Beta-Binomial intensity
#' equation.  Together they provide a minimal but targeted diagnostic
#' toolkit, following the visualisation philosophy of Gabry et al. (2019).
#'
#' @section Coverage Criterion:
#' A test statistic "passes" the PPC at level \eqn{\alpha} if the
#' observed value falls within the \eqn{(1-\alpha)} central posterior
#' predictive interval:
#' \deqn{
#'   T(\mathbf{y}) \;\in\;
#'   \bigl[Q_{\alpha/2}\bigl(T(\mathbf{y}^{\mathrm{rep}})\bigr),\;
#'          Q_{1-\alpha/2}\bigl(T(\mathbf{y}^{\mathrm{rep}})\bigr)
#'   \bigr].
#' }
#' Coverage is reported in the \code{coverage} element of the returned
#' object.
#'
#' @section Extraction Methods:
#' \describe{
#'   \item{\code{"stan"}}{Extracts the posterior predictive replications
#'     \eqn{\mathbf{y}^{\mathrm{rep}}} directly from the \code{y_rep[N]}
#'     array in the Stan generated quantities block.  This is the
#'     preferred method: it uses the exact joint posterior and incurs no
#'     additional R-side simulation cost.}
#'   \item{\code{"simulate"}}{Reconstructs
#'     \eqn{\mathbf{y}^{\mathrm{rep}}} in R by composing posterior draws
#'     of \eqn{(\alpha, \beta, \log\kappa)} with the hurdle
#'     Beta-Binomial PMF via \code{\link{rhurdle_betabinom}}.  Useful
#'     for cross-validating the Stan GQ block and for models where GQ
#'     draws are unavailable.}
#' }
#'
#' @param fit An object of class \code{"hbb_fit"} returned by
#'   \code{\link{hbb}}.  Must contain a CmdStanMCMC fit with generated
#'   quantities \code{log_lik[N]} and \code{y_rep[N]}.
#' @param type Character string; which test statistics to compute.
#'   One of \code{"both"} (default), \code{"zero_rate"}, or
#'   \code{"it_share"}.
#' @param method Character string; how to obtain posterior predictive
#'   draws.  One of \code{"stan"} (default, extracts \code{y_rep} from
#'   the Stan GQ block) or \code{"simulate"} (generates draws in R via
#'   \code{\link{rhurdle_betabinom}}).
#' @param n_draws Integer or \code{NULL}.  Number of posterior draws to
#'   use.  If \code{NULL} (default), all available draws are used.  If
#'   specified, draws are subsampled by systematic thinning.
#' @param level Numeric scalar in \eqn{(0, 1)}.  Coverage level for the
#'   posterior predictive interval.  Default \code{0.95}.
#'
#' @return An S3 object of class \code{"hbb_ppc"} containing:
#' \describe{
#'   \item{\code{observed}}{Named list: \code{zero_rate} (scalar),
#'     \code{it_share} (scalar or \code{NA} if no positive observations),
#'     \code{N} (total observations), \code{N_pos} (count of positives).}
#'   \item{\code{predicted}}{Named list with elements \code{zero_rate}
#'     and \code{it_share} (each \code{NULL} if not requested), where
#'     each non-null element is a list containing:
#'     \code{draws} (numeric vector of length \eqn{M}),
#'     \code{mean}, \code{median}, \code{ci} (named length-2 vector with
#'     entries \code{"lower"} and \code{"upper"}), \code{sd}.}
#'   \item{\code{coverage}}{Named list with elements \code{zero_rate}
#'     and \code{it_share} (each \code{NULL} if not requested), where
#'     each non-null element contains \code{in_ci} (logical) and
#'     \code{ci} (named length-2 vector).}
#'   \item{\code{n_draws}}{Integer; number of draws used.}
#'   \item{\code{n_total_draws}}{Integer; total available draws before
#'     thinning.}
#'   \item{\code{type}}{Character; the \code{type} argument used.}
#'   \item{\code{method}}{Character; the \code{method} argument used.}
#'   \item{\code{level}}{Numeric; the \code{level} argument used.}
#'   \item{\code{model_type}}{Character; from \code{fit$model_type}.}
#' }
#'
#' @examples
#' \dontrun{
#' # Fit a model
#' fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data,
#'            weights = "weight")
#'
#' # Full PPC (both statistics, Stan draws)
#' ppc_result <- ppc(fit)
#' print(ppc_result)
#' plot(ppc_result)
#'
#' # Zero-rate only, using 500 draws
#' ppc_zr <- ppc(fit, type = "zero_rate", n_draws = 500)
#'
#' # Cross-check via R simulation
#' ppc_sim <- ppc(fit, method = "simulate", n_draws = 200)
#' }
#'
#' @seealso
#' \code{\link{loo.hbb_fit}} for LOO-CV model comparison,
#' \code{\link{hbb}} for model fitting,
#' \code{\link{rhurdle_betabinom}} for the hurdle BetaBinomial sampler.
#'
#' @references
#' Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A.
#' (2019). Visualisation in Bayesian workflow. \emph{Journal of the Royal
#' Statistical Society: Series A}, \strong{182}(2), 389--402.
#'
#' Gelman, A., Meng, X.-L., and Stern, H. (1996). Posterior predictive
#' assessment of model fitness via realized discrepancies.
#' \emph{Statistica Sinica}, \strong{6}(4), 733--760.
#'
#' Ghosal, R., Ghosh, S. K., and Maiti, T. (2020). Two-part regression
#' models for longitudinal zero-inflated count data. \emph{Journal of the
#' Royal Statistical Society: Series A}, \strong{183}(4), 1603--1626.
#'
#' @family model-checking
#' @export
ppc <- function(fit,
                type    = c("both", "zero_rate", "it_share"),
                method  = c("stan", "simulate"),
                n_draws = NULL,
                level   = 0.95) {

    # =========================================================================
    # 0. INPUT VALIDATION
    # =========================================================================

    type   <- match.arg(type)
    method <- match.arg(method)

    # .ppc_validate_inputs may return list(method = "simulate") if y_rep is
    # unavailable and automatic fallback is triggered.
    val_result <- .ppc_validate_inputs(fit, type, method, n_draws, level)
    method     <- val_result$method   # may be updated to "simulate"

    N       <- fit$hbb_data$N
    y       <- fit$hbb_data$y
    n_trial <- fit$hbb_data$n_trial

    # =========================================================================
    # 1. OBSERVED TEST STATISTICS
    # =========================================================================

    obs_zero_rate <- mean(y == 0L)
    pos_idx_obs   <- which(y > 0L)
    N_pos         <- length(pos_idx_obs)

    obs_it_share <- if (N_pos > 0L) {
        mean(y[pos_idx_obs] / n_trial[pos_idx_obs])
    } else {
        cli::cli_warn(c(
            "!" = "No positive observations ({.code y > 0}) in the data.",
            "i" = "Observed IT share is {.val NA}."
        ))
        NA_real_
    }

    # =========================================================================
    # 2. EXTRACT / SIMULATE POSTERIOR PREDICTIVE DRAWS
    # =========================================================================

    if (method == "stan") {
        y_rep_mat <- .ppc_extract_yrep(fit, n_draws)
    } else {
        y_rep_mat <- .ppc_simulate_yrep(fit, n_draws)
    }

    M_used        <- nrow(y_rep_mat)
    M_total_draws <- attr(y_rep_mat, "M_total") %||% M_used

    cli::cli_alert_info(
        "PPC: N = {N}, M draws = {M_used}, type = '{type}', method = '{method}'"
    )

    # =========================================================================
    # 3. COMPUTE PPC STATISTICS PER DRAW
    # =========================================================================

    stats_list <- .ppc_compute_stats(y_rep_mat, n_trial, type)

    # =========================================================================
    # 4. SUMMARIZE DRAWS
    # =========================================================================

    pred_zero_rate <- NULL
    pred_it_share  <- NULL

    if (!is.null(stats_list$zero_rate)) {
        pred_zero_rate <- .ppc_summarize_draws(stats_list$zero_rate, level)
    }
    if (!is.null(stats_list$it_share)) {
        pred_it_share <- .ppc_summarize_draws(stats_list$it_share, level)
    }

    # =========================================================================
    # 5. COVERAGE CHECK
    # =========================================================================

    cov_zero_rate <- NULL
    cov_it_share  <- NULL

    if (!is.null(pred_zero_rate)) {
        ci <- pred_zero_rate$ci
        cov_zero_rate <- list(
            in_ci = (obs_zero_rate >= ci[[1L]] && obs_zero_rate <= ci[[2L]]),
            ci    = ci
        )
    }

    if (!is.null(pred_it_share) && !is.na(obs_it_share)) {
        ci <- pred_it_share$ci
        cov_it_share <- list(
            in_ci = (obs_it_share >= ci[[1L]] && obs_it_share <= ci[[2L]]),
            ci    = ci
        )
    }

    # =========================================================================
    # 6. CLI SUMMARY
    # =========================================================================

    if (!is.null(cov_zero_rate)) {
        status <- if (cov_zero_rate$in_ci) "PASS" else "FAIL"
        cli::cli_alert_info(
            "Zero-rate PPC: obs = {round(obs_zero_rate, 4)}, pred mean = {round(pred_zero_rate$mean, 4)} [{status}]"
        )
    }
    if (!is.null(pred_it_share) && !is.na(obs_it_share)) {
        status <- if (!is.null(cov_it_share) && cov_it_share$in_ci) "PASS" else "FAIL"
        cli::cli_alert_info(
            "IT-share  PPC: obs = {round(obs_it_share, 4)}, pred mean = {round(pred_it_share$mean, 4)} [{status}]"
        )
    }

    cli::cli_alert_success("PPC complete")

    # =========================================================================
    # 7. CONSTRUCT RETURN OBJECT
    # =========================================================================

    structure(
        list(
            observed = list(
                zero_rate = obs_zero_rate,
                it_share  = obs_it_share,
                N         = N,
                N_pos     = N_pos
            ),
            predicted = list(
                zero_rate = pred_zero_rate,
                it_share  = pred_it_share
            ),
            coverage = list(
                zero_rate = cov_zero_rate,
                it_share  = cov_it_share
            ),
            n_draws       = M_used,
            n_total_draws = M_total_draws,
            type          = type,
            method        = method,
            level         = level,
            model_type    = fit$model_type %||% "(unknown)"
        ),
        class = "hbb_ppc"
    )
}


# ============================================================================
# 2. print.hbb_ppc --- S3 print method (exported)
# ============================================================================

#' Print Method for hbb_ppc Objects
#'
#' Displays a compact numerical summary of the posterior predictive check,
#' showing the observed test statistic, the posterior predictive mean and
#' median, the credible interval, and whether the observed value falls
#' within the interval.
#'
#' @param x An object of class \code{"hbb_ppc"} returned by
#'   \code{\link{ppc}}.
#' @param digits Integer; number of decimal places for numeric output.
#'   Default \code{4}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The object \code{x}, invisibly.
#'
#' @seealso \code{\link{ppc}}, \code{\link{plot.hbb_ppc}}
#'
#' @family model-checking
#' @method print hbb_ppc
#' @export
print.hbb_ppc <- function(x, digits = 4L, ...) {

    tryCatch({

        cat("\n")
        cat(strrep("=", 65), "\n")
        cat("  Posterior Predictive Check (HBB Model)\n")
        cat(strrep("=", 65), "\n\n")

        # -- Safe field helper (avoids crashes on corrupted objects) -----------
        .sf <- function(expr_fn, fallback) {
            tryCatch(expr_fn(), error = function(e) fallback)
        }

        # -- Metadata (all field accesses via .sf) -----------------------------
        model_type <- .sf(function() x$model_type %||% "(unknown)", "(unknown)")
        x_method   <- .sf(function() x$method     %||% "(unknown)", "(unknown)")
        x_level    <- .sf(function() {
            l <- x$level %||% 0.95
            if (!is.numeric(l) || !is.finite(l)) 0.95 else l
        }, 0.95)
        n_drw  <- .sf(function() x$n_draws        %||% NA_integer_, NA_integer_)
        n_tot  <- .sf(function() x$n_total_draws  %||% NA_integer_, NA_integer_)
        obs_N  <- .sf(function() x$observed$N     %||% NA_integer_, NA_integer_)
        obs_Np <- .sf(function() x$observed$N_pos %||% NA_integer_, NA_integer_)

        cat(sprintf("  Model type    : %s\n", model_type))
        cat(sprintf("  Observations  : %s\n",
                    if (is.na(obs_N)) "?" else as.character(obs_N)))
        if (!is.na(obs_N) && !is.na(obs_Np)) {
            cat(sprintf("  Positive obs  : %d (%.1f%% participation)\n",
                        obs_Np, 100 * obs_Np / max(obs_N, 1L)))
        }
        cat(sprintf("  Draws used    : %s (of %s total)\n",
                    if (is.na(n_drw)) "?" else as.character(n_drw),
                    if (is.na(n_tot)) "?" else as.character(n_tot)))
        cat(sprintf("  Method        : %s\n", x_method))
        cat(sprintf("  CI level      : %.2f\n", x_level))
        cat("\n")

        # -- Column header -----------------------------------------------------
        ci_pct <- sprintf("%.0f%% PI", 100 * x_level)
        fmt <- "  %-12s  %8s  %8s  %8s  %-22s  %s\n"
        cat(sprintf(fmt,
                    "Statistic", "Observed", "Pred.Mean", "Pred.Med",
                    ci_pct, "Coverage"))
        cat(strrep("-", 72), "\n")

        # -- Helper for one row ------------------------------------------------
        .print_row <- function(stat_name, obs_val, pred) {
            tryCatch({
                if (is.null(pred)) return(invisible(NULL))

                cv     <- .sf(function() x$coverage[[stat_name]], NULL)
                in_ci  <- .sf(function() cv$in_ci %||% NA, NA)

                status <- if (!is.na(in_ci) && !is.null(obs_val) &&
                               !is.na(obs_val)) {
                    if (isTRUE(in_ci)) "[PASS]" else "[FAIL]"
                } else if (is.null(obs_val) || is.na(obs_val)) {
                    "[NA: no positives]"
                } else {
                    ""
                }

                p_mean   <- .sf(function() pred$mean   %||% NA_real_, NA_real_)
                p_median <- .sf(function() pred$median %||% NA_real_, NA_real_)
                ci_lo    <- .sf(function() pred$ci[[1L]] %||% NA_real_, NA_real_)
                ci_hi    <- .sf(function() pred$ci[[2L]] %||% NA_real_, NA_real_)
                n_finite <- .sf(function() pred$n_finite %||% NA_integer_,
                                NA_integer_)
                n_na_d   <- .sf(function() pred$n_na %||% 0L, 0L)

                ci_str <- if (!is.na(ci_lo) && !is.na(ci_hi)) {
                    sprintf("[%.*f, %.*f]", digits, ci_lo, digits, ci_hi)
                } else {
                    "(unavailable)"
                }

                obs_str <- if (!is.null(obs_val) && !is.na(obs_val)) {
                    sprintf("%8.*f", digits, obs_val)
                } else {
                    sprintf("%8s", "NA")
                }
                pm_str  <- if (!is.na(p_mean))   sprintf("%8.*f", digits, p_mean)
                           else sprintf("%8s", "NA")
                pmd_str <- if (!is.na(p_median)) sprintf("%8.*f", digits, p_median)
                           else sprintf("%8s", "NA")

                cat(sprintf("  %-12s  %s  %s  %s  %-22s  %s\n",
                            stat_name, obs_str, pm_str, pmd_str,
                            ci_str, status))

                # Note: all-NA draws
                if (!is.na(n_na_d) && n_na_d > 0L) {
                    cat(sprintf(
                        "    (Note: %d draw(s) had all-zero replications -- "
                        , n_na_d))
                    cat("IT share was NA for those draws.)\n")
                }
                if (!is.na(n_finite) && !is.na(n_drw) && n_finite < n_drw) {
                    cat(sprintf(
                        "    (Note: %d/%d draws had finite %s values.)\n",
                        n_finite, n_drw, stat_name))
                }

            }, error = function(e) {
                cat(sprintf("  %-12s  (print error: %s)\n",
                            stat_name, conditionMessage(e)))
            })
        }

        pred_lst <- .sf(function() x$predicted, list())
        .print_row("zero_rate",
                   .sf(function() x$observed$zero_rate, NULL),
                   .sf(function() pred_lst$zero_rate, NULL))
        .print_row("it_share",
                   .sf(function() x$observed$it_share, NULL),
                   .sf(function() pred_lst$it_share, NULL))

        cat(strrep("-", 72), "\n")
        cat("\n  Note: 'Coverage' indicates whether observed is inside the\n")
        cat("        posterior predictive interval (PASS/FAIL).\n\n")

    }, error = function(e) {
        # Last-resort: bare-minimum output so print() never crashes
        cat("Error printing hbb_ppc object:", conditionMessage(e), "\n")
        cat("Object class:",
            paste(tryCatch(class(x), error = function(e2) "?"),
                  collapse = ", "), "\n")
        cat("Object names:",
            paste(tryCatch(names(x), error = function(e2) "?"),
                  collapse = ", "), "\n")
    })

    invisible(x)
}


# ============================================================================
# 3. plot.hbb_ppc --- S3 plot method (exported)
# ============================================================================

#' Plot Method for hbb_ppc Objects
#'
#' Produces a histogram of the posterior predictive draws for each
#' requested test statistic, overlaid with a vertical line at the
#' observed value.  When \code{type = "both"}, the two panels are
#' combined using \pkg{patchwork}.
#'
#' @description
#' Each panel shows:
#' \itemize{
#'   \item A histogram of the \eqn{M} posterior predictive draws of
#'     the test statistic \eqn{T(\mathbf{y}^{\mathrm{rep}})}.
#'   \item A solid red vertical line at the observed value
#'     \eqn{T(\mathbf{y})}, annotated with its numeric value.
#'   \item Dashed green vertical lines at the lower and upper bounds of
#'     the \eqn{(1-\alpha)} posterior predictive interval.
#' }
#' A pass/fail indicator appears in the panel subtitle.
#'
#' @section ggplot2 Requirement:
#' This function requires \pkg{ggplot2} (and \pkg{patchwork} for
#' \code{type = "both"}).  Both are in \code{Suggests}; they will be
#' checked via \code{rlang::check_installed()} and the user will be
#' prompted to install them if absent.
#'
#' @param x An object of class \code{"hbb_ppc"} returned by
#'   \code{\link{ppc}}.
#' @param type Character string or \code{NULL}; which statistic to plot.
#'   If \code{NULL} (default), uses \code{x$type}.  Otherwise one of
#'   \code{"both"}, \code{"zero_rate"}, \code{"it_share"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{ggplot} object (single-statistic) or a
#'   \code{patchwork} object (\code{type = "both"}), returned invisibly.
#'
#' @seealso \code{\link{ppc}}, \code{\link{print.hbb_ppc}}
#'
#' @references
#' Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A.
#' (2019). Visualisation in Bayesian workflow. \emph{Journal of the Royal
#' Statistical Society: Series A}, \strong{182}(2), 389--402.
#'
#' @family model-checking
#' @method plot hbb_ppc
#' @export
plot.hbb_ppc <- function(x, type = NULL, ...) {

    rlang::check_installed(
        "ggplot2",
        reason = "to plot hbb_ppc objects"
    )

    # -- Resolve type ----------------------------------------------------------
    if (is.null(type)) {
        type <- x$type
    } else {
        type <- match.arg(type, c("both", "zero_rate", "it_share"))
    }

    # -- Internal panel builder ------------------------------------------------
    .make_ppc_panel <- function(draws, obs_val, ci, level, stat_label) {

        ci_pct <- sprintf("%.0f%%", 100 * level)

        # Subtitle: pass/fail indicator
        if (!is.null(ci) && !is.na(obs_val)) {
            in_ci    <- (obs_val >= ci[[1L]] && obs_val <= ci[[2L]])
            subtitle <- if (in_ci) {
                sprintf("PASS: observed within %s predictive interval", ci_pct)
            } else {
                sprintf("FAIL: observed OUTSIDE %s predictive interval", ci_pct)
            }
        } else if (is.na(obs_val)) {
            subtitle <- "Observed: NA (no positive observations)"
        } else {
            subtitle <- ""
        }

        df_draws <- data.frame(stat = draws)

        p <- ggplot2::ggplot(
            df_draws,
            ggplot2::aes(x = .data[["stat"]])
        ) +
            ggplot2::geom_histogram(
                bins   = 40L,
                fill   = "#4393C3",
                colour = "white",
                alpha  = 0.75
            ) +
            ggplot2::labs(
                title    = sprintf("PPC: %s", stat_label),
                subtitle = subtitle,
                x        = sprintf("Posterior predictive %s", stat_label),
                y        = "Count"
            ) +
            ggplot2::theme_bw(base_size = 11L)

        # Dashed green lines for CI bounds
        if (!is.null(ci)) {
            p <- p +
                ggplot2::geom_vline(
                    xintercept = ci[[1L]],
                    linetype   = "dashed",
                    colour     = "#1B7837",
                    linewidth  = 0.8
                ) +
                ggplot2::geom_vline(
                    xintercept = ci[[2L]],
                    linetype   = "dashed",
                    colour     = "#1B7837",
                    linewidth  = 0.8
                )
        }

        # Solid red line + annotation for observed value
        if (!is.na(obs_val)) {
            p <- p +
                ggplot2::geom_vline(
                    xintercept = obs_val,
                    colour     = "#D6604D",
                    linewidth  = 1.1
                ) +
                ggplot2::annotate(
                    "text",
                    x      = obs_val,
                    y      = Inf,
                    label  = sprintf("obs = %.4f", obs_val),
                    colour = "#D6604D",
                    hjust  = -0.1,
                    vjust  = 1.5,
                    size   = 3.5
                )
        }

        p
    }

    # -- Build individual panels -----------------------------------------------

    p_zr <- NULL
    p_it <- NULL

    if (type %in% c("both", "zero_rate")) {
        pz <- x$predicted$zero_rate
        if (!is.null(pz)) {
            p_zr <- .make_ppc_panel(
                draws      = pz$draws,
                obs_val    = x$observed$zero_rate,
                ci         = pz$ci,
                level      = x$level,
                stat_label = "Zero Rate"
            )
        } else {
            cli::cli_alert_warning(
                "Zero-rate draws not available; cannot plot zero_rate panel."
            )
        }
    }

    if (type %in% c("both", "it_share")) {
        pi_ <- x$predicted$it_share
        if (!is.null(pi_)) {
            p_it <- .make_ppc_panel(
                draws      = pi_$draws,
                obs_val    = x$observed$it_share,
                ci         = pi_$ci,
                level      = x$level,
                stat_label = "IT Share (conditional on positive)"
            )
        } else {
            cli::cli_alert_warning(
                "IT-share draws not available; cannot plot it_share panel."
            )
        }
    }

    # -- Combine panels --------------------------------------------------------

    if (type == "both" && !is.null(p_zr) && !is.null(p_it)) {
        rlang::check_installed(
            "patchwork",
            reason = "to combine both PPC panels side-by-side"
        )
        out <- p_zr + p_it +
            patchwork::plot_annotation(
                title   = "Posterior Predictive Checks (HBB Model)",
                caption = sprintf(
                    "Dashed green = %d%% PI  |  Red = observed value  |  M = %d draws",
                    as.integer(100 * x$level),
                    x$n_draws
                )
            )
    } else if (!is.null(p_zr)) {
        out <- p_zr
    } else if (!is.null(p_it)) {
        out <- p_it
    } else {
        cli::cli_abort("No PPC panels available to plot.")
    }

    print(out)
    invisible(out)
}


# ============================================================================
# Internal helpers
# ============================================================================


# ---- 4. .ppc_validate_inputs -----------------------------------------------

#' Comprehensive input validation for ppc()
#'
#' Validates all arguments.  For \code{method = "stan"}, attempts a trial
#' extraction of \code{y_rep} draws from the Stan fit.  If extraction
#' fails (e.g., model compiled without generated quantities), emits a
#' \code{cli_warn} and returns \code{list(method = "simulate")} to trigger
#' automatic fallback in \code{ppc()}.
#'
#' @param fit An object to validate as \code{hbb_fit}.
#' @param type Character: already matched by \code{match.arg}.
#' @param method Character: already matched by \code{match.arg}.
#' @param n_draws \code{NULL} or positive integer.
#' @param level Numeric in \eqn{(0, 1)}.
#' @return Named list with element \code{method} (possibly updated to
#'   \code{"simulate"} after fallback).  Aborts on unrecoverable errors.
#' @noRd
.ppc_validate_inputs <- function(fit, type, method, n_draws, level) {

    # 1. fit class
    if (!inherits(fit, "hbb_fit")) {
        cli::cli_abort(c(
            "x" = "{.arg fit} must be an {.cls hbb_fit} object, not {.cls {class(fit)}}.",
            "i" = "Use {.fn hbb} to fit a hurdle Beta-Binomial model."
        ))
    }

    # 2. fit$fit not NULL
    if (is.null(fit$fit)) {
        cli::cli_abort(c(
            "x" = "{.arg fit} does not contain a CmdStan fit object.",
            "i" = "{.code fit$fit} is {.val NULL}.",
            "i" = "Was the model successfully sampled with {.fn hbb}?"
        ))
    }

    # 3. hbb_data required fields
    hd <- fit$hbb_data
    if (is.null(hd)) {
        cli::cli_abort(c(
            "x" = "{.code fit$hbb_data} is {.val NULL}.",
            "i" = "The {.cls hbb_fit} object appears to be incomplete."
        ))
    }
    for (fld in c("y", "n_trial", "N")) {
        if (is.null(hd[[fld]])) {
            cli::cli_abort(c(
                "x" = "{.code fit$hbb_data${fld}} is {.val NULL}; required for PPC.",
                "i" = "Refit the model to regenerate {.code hbb_data}."
            ))
        }
    }
    if (length(hd$y) != hd$N || length(hd$n_trial) != hd$N) {
        cli::cli_abort(c(
            "x" = "Length mismatch in {.code fit$hbb_data}.",
            "i" = "{.code y} has length {length(hd$y)}, {.code n_trial} has \\
                   length {length(hd$n_trial)}, but N = {hd$N}."
        ))
    }

    # 4. level in (0, 1) exclusive
    .validate_level(level)

    # 5. n_draws: NULL or positive integer (checkmate for clear messages)
    if (!is.null(n_draws)) {
        checkmate::assert_count(n_draws,
                                positive  = TRUE,
                                null.ok   = FALSE,
                                .var.name = "n_draws")
    }

    # 6. For method = "stan": trial extraction of y_rep to verify availability.
    #    On failure: emit cli_warn and return list(method = "simulate") to
    #    signal automatic fallback to the caller.
    if (method == "stan") {
        test_ok <- tryCatch({
            test_mat <- fit$fit$draws(variables = "y_rep", format = "matrix")
            !is.null(test_mat) && length(test_mat) > 0L
        }, error = function(e) FALSE)

        if (!test_ok) {
            cli::cli_warn(c(
                "!" = "{.code y_rep} draws are not available in the Stan fit.",
                "i" = "This occurs when the model was compiled without a \\
                       generated quantities block producing {.code y_rep[N]}.",
                "i" = "Automatically falling back to {.arg method} = {.val simulate}.",
                "i" = "To suppress this warning, call \\
                       {.code ppc(fit, method = \"simulate\")} directly."
            ))
            return(list(method = "simulate"))
        }
    }

    list(method = method)
}


# ---- 5. .ppc_extract_yrep --------------------------------------------------

#' Extract y_rep from the Stan GQ block
#'
#' Extracts the \code{y_rep[N]} array produced by the Stan generated
#' quantities block and returns it as an \eqn{M \times N} integer matrix.
#' Systematic thinning is applied if \code{n_draws < M_total}.
#'
#' @details
#' The Stan GQ block is assumed to produce the variable \code{y_rep}
#' (an array of length \code{N}), which CmdStanR exposes as
#' \code{y_rep[1]}, \ldots, \code{y_rep[N]}.  The extraction uses
#' \code{fit$fit$draws(variables = "y_rep", format = "matrix")} which
#' returns an \eqn{M \times N} matrix directly.
#'
#' @param fit An \code{hbb_fit} object.
#' @param n_draws \code{NULL} or positive integer; target draw count
#'   after thinning.
#' @return An \eqn{M \times N} integer matrix.  The attribute
#'   \code{"M_total"} records the pre-thinning draw count.
#' @noRd
.ppc_extract_yrep <- function(fit, n_draws) {

    N <- fit$hbb_data$N

    # -- Extract draws matrix (M_total x N) -----------------------------------
    y_rep_mat <- tryCatch(
        fit$fit$draws(variables = "y_rep", format = "matrix"),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "Failed to extract {.code y_rep} from the CmdStanMCMC fit.",
                "i" = "Ensure the Stan model was compiled and run with generated quantities.",
                "i" = "Error: {conditionMessage(e)}"
            ))
        }
    )

    if (!is.matrix(y_rep_mat)) y_rep_mat <- as.matrix(y_rep_mat)

    # -- Empty matrix check (edge case: zero rows) ----------------------------
    M_total <- nrow(y_rep_mat)
    if (M_total < 1L) {
        cli::cli_abort(c(
            "x" = "Extracted {.code y_rep} matrix has zero rows.",
            "i" = "The Stan fit may be empty or corrupted."
        ))
    }

    # -- Dimension check -------------------------------------------------------
    if (ncol(y_rep_mat) != N) {
        cli::cli_abort(c(
            "x" = "Extracted {.code y_rep} has {ncol(y_rep_mat)} columns but N = {N}.",
            "i" = "Ensure the Stan model was fit with the same data (same N)."
        ))
    }

    # -- NA/NaN detection (warn; do not abort --downstream handles NA) --------
    n_na <- sum(is.na(y_rep_mat))
    if (n_na > 0L) {
        pct_na <- round(100 * n_na / length(y_rep_mat), 2)
        cli::cli_warn(c(
            "!" = "{n_na} NA/NaN value{?s} detected in {.code y_rep} draws ({pct_na}%).",
            "i" = "These are treated as {.val NA} in PPC statistics."
        ))
    }

    # -- Infinite value detection and replacement -----------------------------
    n_inf <- sum(is.infinite(y_rep_mat))
    if (n_inf > 0L) {
        cli::cli_warn(c(
            "!" = "{n_inf} infinite value{?s} in {.code y_rep} draws.",
            "i" = "These will be replaced with {.val NA}."
        ))
        y_rep_mat[is.infinite(y_rep_mat)] <- NA_integer_
    }

    # -- Systematic thinning --------------------------------------------------
    if (!is.null(n_draws)) {
        if (n_draws >= M_total) {
            cli::cli_alert_info(
                "{.arg n_draws} = {n_draws} >= available draws ({M_total}). Using all draws."
            )
        } else {
            thin_idx  <- round(seq(1L, M_total, length.out = as.integer(n_draws)))
            thin_idx  <- unique(thin_idx)
            y_rep_mat <- y_rep_mat[thin_idx, , drop = FALSE]
            cli::cli_alert_info(
                "Thinned y_rep from {M_total} to {nrow(y_rep_mat)} draws"
            )
        }
    }

    attr(y_rep_mat, "M_total") <- M_total
    storage.mode(y_rep_mat)    <- "integer"
    y_rep_mat
}


# ---- 6. .ppc_simulate_yrep -------------------------------------------------

#' Simulate y_rep in R via rhurdle_betabinom
#'
#' For each posterior draw of \eqn{(\alpha, \beta, \log\kappa)},
#' computes the linear predictors, applies the inverse logit transform
#' to obtain \eqn{(q_i, \mu_i)}, exponentiates \eqn{\log\kappa} to
#' obtain \eqn{\kappa}, then calls \code{\link{rhurdle_betabinom}} to
#' generate a replicated data vector
#' \eqn{\mathbf{y}^{\mathrm{rep},(m)} \in \{0,\ldots,n_i\}^N}.
#'
#' @details
#' The composing steps for draw \eqn{m} are:
#' \deqn{
#'   q_i^{(m)}  = \mathrm{logistic}(X_i' \alpha^{(m)}), \quad
#'   \mu_i^{(m)} = \mathrm{logistic}(X_i' \beta^{(m)}), \quad
#'   \kappa^{(m)} = \exp(\log\kappa^{(m)}).
#' }
#' All values are clamped to valid ranges before calling
#' \code{rhurdle_betabinom} to prevent edge-case failures in the
#' Beta-Binomial sampler.
#'
#' @param fit An \code{hbb_fit} object.
#' @param n_draws \code{NULL} or positive integer; target draw count
#'   after thinning.
#' @return An \eqn{M \times N} integer matrix with attribute
#'   \code{"M_total"}.
#' @noRd
.ppc_simulate_yrep <- function(fit, n_draws) {

    X       <- fit$hbb_data$X
    P       <- fit$hbb_data$P
    N       <- fit$hbb_data$N
    n_trial <- fit$hbb_data$n_trial

    # -- Validate X availability -----------------------------------------------
    if (is.null(X) || !is.matrix(X)) {
        cli::cli_abort(c(
            "x" = "Design matrix {.code fit$hbb_data$X} is missing or not a matrix.",
            "i" = "Cannot simulate {.code y_rep} without the covariate matrix."
        ))
    }
    if (nrow(X) != N) {
        cli::cli_abort(c(
            "x" = "Dimension mismatch: {.code nrow(X)} = {nrow(X)} but N = {N}.",
            "i" = "Check that {.code hbb_data} is internally consistent."
        ))
    }

    # -- Extract parameter draws with name validation -------------------------
    alpha_names <- paste0("alpha[", seq_len(P), "]")
    beta_names  <- paste0("beta[",  seq_len(P), "]")
    param_names <- c(alpha_names, beta_names, "log_kappa")

    draws_mat <- tryCatch(
        fit$fit$draws(variables = param_names, format = "matrix"),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "Failed to extract parameter draws for R-side simulation.",
                "i" = "Error: {conditionMessage(e)}"
            ))
        }
    )
    if (!is.matrix(draws_mat)) draws_mat <- as.matrix(draws_mat)

    # Validate parameter names exist in the extracted matrix
    missing_pars <- setdiff(param_names, colnames(draws_mat))
    if (length(missing_pars) > 0L) {
        cli::cli_abort(c(
            "x" = "Expected parameter{?s} not found in draws: {.val {missing_pars}}.",
            "i" = "Model may have a different parameterisation than expected."
        ))
    }

    M_total <- nrow(draws_mat)

    # -- Systematic thinning --------------------------------------------------
    if (!is.null(n_draws)) {
        if (n_draws >= M_total) {
            cli::cli_alert_info(
                "{.arg n_draws} = {n_draws} >= available draws ({M_total}). Using all draws."
            )
        } else {
            thin_idx  <- round(seq(1L, M_total, length.out = as.integer(n_draws)))
            thin_idx  <- unique(thin_idx)
            draws_mat <- draws_mat[thin_idx, , drop = FALSE]
            cli::cli_alert_info(
                "Thinned parameter draws from {M_total} to {nrow(draws_mat)} draws"
            )
        }
    }

    M_use <- nrow(draws_mat)
    eps   <- .Machine$double.eps

    # -- Bulk check for non-finite log_kappa -----------------------------------
    lk_all   <- as.numeric(draws_mat[, "log_kappa"])
    n_bad_lk <- sum(!is.finite(lk_all))
    if (n_bad_lk > 0L) {
        cli::cli_warn(c(
            "!" = "{n_bad_lk} non-finite {.code log_kappa} draw{?s} \\
                   ({round(100*n_bad_lk/M_use, 1)}% of draws).",
            "i" = "Affected draws will yield all-NA rows in {.code y_rep}."
        ))
    }

    # -- Allocate output -------------------------------------------------------
    y_rep_mat <- matrix(NA_integer_, nrow = M_use, ncol = N)

    n_skip      <- 0L
    n_nan_total <- 0L

    # -- Simulate per draw -----------------------------------------------------
    for (m in seq_len(M_use)) {

        alpha_m     <- as.numeric(draws_mat[m, alpha_names])
        beta_m      <- as.numeric(draws_mat[m, beta_names])
        log_kappa_m <- as.numeric(draws_mat[m, "log_kappa"])

        # Guard: non-finite parameters ->skip this draw (leave row as NA)
        if (!is.finite(log_kappa_m) || any(!is.finite(alpha_m)) ||
            any(!is.finite(beta_m))) {
            n_skip <- n_skip + 1L
            next
        }

        kappa_m <- exp(log_kappa_m)

        # Guard: Inf kappa (extreme log_kappa)
        if (!is.finite(kappa_m) || kappa_m <= 0) {
            n_skip <- n_skip + 1L
            next
        }

        q_m  <- plogis(as.numeric(X %*% alpha_m))
        mu_m <- plogis(as.numeric(X %*% beta_m))

        # Guard against NaN from plogis (extreme linear predictors)
        n_nan_q  <- sum(is.nan(q_m))
        n_nan_mu <- sum(is.nan(mu_m))
        if (n_nan_q + n_nan_mu > 0L) {
            n_nan_total <- n_nan_total + n_nan_q + n_nan_mu
            # Replace NaN with boundary-neutral value (0.5)
            q_m[is.nan(q_m)]   <- 0.5
            mu_m[is.nan(mu_m)] <- 0.5
        }

        # Clamp to valid ranges for sampler stability
        q_m     <- pmin(pmax(q_m,  eps), 1 - eps)
        mu_m    <- pmin(pmax(mu_m, eps), 1 - eps)
        kappa_m <- max(kappa_m, eps)

        y_rep_mat[m, ] <- rhurdle_betabinom(
            nn    = N,
            n     = n_trial,
            q     = q_m,
            mu    = mu_m,
            kappa = kappa_m
        )
    }

    # -- Report skipped draws once (avoid per-draw spam) ----------------------
    if (n_skip > 0L) {
        cli::cli_warn(c(
            "!" = "{n_skip} draw{?s} skipped due to non-finite parameters.",
            "i" = "Those rows in {.code y_rep} remain all {.val NA}."
        ))
    }
    if (n_nan_total > 0L) {
        cli::cli_warn(c(
            "!" = "{n_nan_total} NaN value{?s} in linear predictors (replaced with 0.5).",
            "i" = "Consider checking the posterior for extreme parameter draws."
        ))
    }

    attr(y_rep_mat, "M_total") <- M_total
    storage.mode(y_rep_mat)    <- "integer"
    y_rep_mat
}


# ---- 7. .ppc_compute_stats -------------------------------------------------

#' Compute PPC test statistics for each posterior predictive draw
#'
#' For each row of the \eqn{M \times N} replicated data matrix, computes
#' the zero rate \eqn{T_0} and/or the conditional IT share \eqn{T_s}.
#'
#' @details
#' \describe{
#'   \item{Zero rate (vectorised)}{
#'     \deqn{T_0(\mathbf{y}^{\mathrm{rep},m})
#'       = \frac{1}{N} \sum_{i=1}^{N}
#'         \mathbf{1}(y_i^{\mathrm{rep},m} = 0).}
#'     Computed as \code{rowMeans(y_rep_mat == 0L)}, which avoids an
#'     explicit loop over draws.
#'   }
#'   \item{IT share (per-draw loop)}{
#'     \deqn{T_s(\mathbf{y}^{\mathrm{rep},m})
#'       = \frac{1}{N_m^+}
#'         \sum_{i:\, y_i^{\mathrm{rep},m} > 0}
#'         \frac{y_i^{\mathrm{rep},m}}{n_i},}
#'     where \eqn{N_m^+} varies by draw and cannot be vectorised
#'     across draws without memory overhead.  Draws where
#'     \eqn{N_m^+ = 0} contribute \code{NA}.
#'   }
#' }
#'
#' @param y_rep_mat \eqn{M \times N} integer matrix of replicated data.
#' @param n_trial Integer vector of length \eqn{N}: trial sizes.
#' @param type Character: \code{"both"}, \code{"zero_rate"}, or
#'   \code{"it_share"}.
#' @return Named list with elements \code{zero_rate} and
#'   \code{it_share}; each is a numeric vector of length \eqn{M}, or
#'   \code{NULL} if not requested.
#' @noRd
.ppc_compute_stats <- function(y_rep_mat, n_trial, type) {

    M <- nrow(y_rep_mat)
    N <- ncol(y_rep_mat)

    zr_vec <- NULL
    it_vec <- NULL

    # -- Zero rate: vectorised row means of binary indicator ------------------
    if (type %in% c("both", "zero_rate")) {
        zr_vec <- rowMeans(y_rep_mat == 0L, na.rm = TRUE)
    }

    # -- IT share: per-draw mean over positive subset -------------------------
    if (type %in% c("both", "it_share")) {

        # Detect and clamp negative values (simulation artifact)
        neg_flag <- !is.na(y_rep_mat) & (y_rep_mat < 0L)
        n_neg    <- sum(neg_flag, na.rm = TRUE)
        if (n_neg > 0L) {
            cli::cli_warn(c(
                "!" = "{n_neg} negative {.code y_rep} value{?s} detected (clamped to 0).",
                "i" = "This may indicate a simulation implementation issue."
            ))
            y_rep_mat[neg_flag] <- 0L
        }

        # Detect values exceeding n_trial (unexpected from ZT-BB)
        n_trial_mat <- matrix(rep(as.numeric(n_trial), each = M),
                              nrow = M, ncol = N)
        exceed_flag <- !is.na(y_rep_mat) & (y_rep_mat > n_trial_mat)
        n_exceed    <- sum(exceed_flag, na.rm = TRUE)
        if (n_exceed > 0L) {
            cli::cli_warn(c(
                "!" = "{n_exceed} {.code y_rep} value{?s} exceed {.code n_trial}.",
                "i" = "This may indicate a simulation inconsistency."
            ))
        }

        # Pre-compute share matrix (M x N); integer division ->numeric
        share_mat <- sweep(
            matrix(as.numeric(y_rep_mat), nrow = M, ncol = N),
            MARGIN = 2L,
            STATS  = as.numeric(n_trial),
            FUN    = "/"
        )

        n_all_zero <- 0L
        it_vec     <- numeric(M)

        for (m in seq_len(M)) {
            pos_mask_m <- !is.na(y_rep_mat[m, ]) &
                          (y_rep_mat[m, ] > 0L) &
                          (n_trial > 0L)
            if (any(pos_mask_m, na.rm = TRUE)) {
                it_vec[m] <- mean(share_mat[m, pos_mask_m], na.rm = TRUE)
            } else {
                it_vec[m] <- NA_real_
                n_all_zero <- n_all_zero + 1L
            }
        }

        if (n_all_zero > 0L) {
            cli::cli_warn(c(
                "!" = "{n_all_zero} draw{?s} had all-zero replications (IT share = NA).",
                "i" = "This is normal if the zero rate is high.",
                "i" = "NA draws are excluded from predictive summaries."
            ))
        }
        if (all(is.na(it_vec))) {
            cli::cli_warn(c(
                "!" = "All {M} draws resulted in NA IT share.",
                "i" = "The IT share PPC cannot be computed.",
                "i" = "Check that {.code n_trial > 0} for at least some observations."
            ))
        }

        attr(it_vec, "n_na") <- n_all_zero
    }

    list(zero_rate = zr_vec, it_share = it_vec)
}


# ---- 8. .ppc_summarize_draws -----------------------------------------------

#' Summary statistics for a vector of posterior predictive draws
#'
#' Returns the posterior predictive mean, median, credible interval,
#' and standard deviation for a single test statistic.  \code{NA} values
#' (possible for \code{it_share} when all replicates are zero) are silently
#' excluded.  Non-finite values are also excluded, and their count is
#' reported in \code{n_finite}.  If all draws are non-finite, returns a
#' list with all-\code{NA} summary fields (no error or abort).  Notes
#' zero-variance distributions via \code{cli_alert_info} without aborting.
#'
#' @param draws Numeric vector of test-statistic draws.
#' @param level Numeric in \eqn{(0, 1)}: CI level.
#' @return Named list: \code{draws}, \code{mean}, \code{median},
#'   \code{ci} (named numeric(2) with entries \code{"lower"} and
#'   \code{"upper"}), \code{sd}, \code{n_finite} (integer count of finite
#'   draws used), \code{n_na} (integer from \code{attr(draws, "n_na")}).
#' @noRd
.ppc_summarize_draws <- function(draws, level) {

    n_na_attr    <- attr(draws, "n_na") %||% 0L
    finite_mask  <- is.finite(draws)
    n_finite     <- sum(finite_mask)
    finite_draws <- draws[finite_mask]

    # All non-finite: return NA summary without error
    if (n_finite == 0L) {
        ci_empty <- c(lower = NA_real_, upper = NA_real_)
        return(list(
            draws    = draws,
            mean     = NA_real_,
            median   = NA_real_,
            ci       = ci_empty,
            sd       = NA_real_,
            n_finite = 0L,
            n_na     = n_na_attr
        ))
    }

    alpha_half <- (1 - level) / 2
    probs      <- c(alpha_half, 1 - alpha_half)

    ci         <- quantile(finite_draws, probs = probs, na.rm = TRUE)
    names(ci)  <- c("lower", "upper")

    sd_val <- if (n_finite >= 2L) sd(finite_draws, na.rm = TRUE) else NA_real_

    # Zero-variance note (degenerate predictive --informative, not fatal)
    if (!is.na(sd_val) && sd_val == 0) {
        cli::cli_alert_info(
            "Posterior predictive draws are constant (all values = {finite_draws[1L]})."
        )
    }

    list(
        draws    = draws,
        mean     = mean(finite_draws,   na.rm = TRUE),
        median   = median(finite_draws, na.rm = TRUE),
        ci       = ci,
        sd       = sd_val,
        n_finite = n_finite,
        n_na     = n_na_attr
    )
}
