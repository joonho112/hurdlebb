# ============================================================================
# loo.R --- Leave-One-Out Cross-Validation for Hurdle Beta-Binomial Models
#
# Implements LOO-CV via Pareto-smoothed importance sampling (PSIS-LOO)
# following Vehtari, Gelman, and Gabry (2017), and model comparison
# utilities for hbb_fit objects.
#
# Theory:
#   PSIS-LOO approximates the leave-one-out predictive density
#
#     p(y_i | y_{-i}) = integral p(y_i | theta) p(theta | y_{-i}) dtheta
#
#   using importance sampling with weights proportional to
#   1 / p(y_i | theta).  Raw IS weights have infinite variance for
#   heavy-tailed posteriors; PSIS stabilises them by fitting a
#   generalised Pareto distribution to the largest weights and
#   replacing the upper tail.  The Pareto shape parameter k provides
#   a diagnostic: k < 0.5 indicates reliable estimates, k in [0.5, 0.7)
#   indicates moderate issues, and k >= 0.7 indicates that LOO estimates
#   may be unreliable.
#
#   The expected log predictive density (ELPD) is:
#
#     elpd_LOO = sum_{i=1}^N log p(y_i | y_{-i})
#
#   Models with higher elpd_LOO predict held-out data better.
#
# Log-likelihood note:
#   The Stan GQ block produces UNWEIGHTED pointwise log-likelihoods
#   log p(y_i | theta^(m)).  For survey-weighted models, these are the
#   nominal (design-ignorant) likelihoods.  LOO is used for model
#   selection among HBB specifications, not for design-based inference;
#   for the latter, use sandwich-corrected Wald statistics.
#
# Relative effective sample size:
#   The accuracy of PSIS is improved by accounting for within-chain
#   autocorrelation via the relative effective sample size (r_eff):
#
#     r_eff_i = ESS_i / M
#
#   where ESS_i is the effective sample size for observation i.
#   Passing r_eff to loo::loo() corrects the Pareto k computation.
#
# References:
#   Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian
#   model evaluation using leave-one-out cross-validation and WAIC.
#   Statistics and Computing, 27(5), 1413-1432.
#
#   Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A.
#   (2019). Visualisation in Bayesian workflow. J. R. Stat. Soc. A,
#   182(2), 389-402.
#
#   Ghosal, R., Ghosh, S. K., and Maiti, T. (2020). Two-part regression
#   models for longitudinal zero-inflated count data. Journal of the
#   Royal Statistical Society: Series A, 183(4), 1603-1626.
#
# Contents:
#   1. loo.hbb_fit           --- S3 loo method (registered via @rawNamespace)
#   2. hbb_loo_compare       --- Multi-model LOO comparison (exported)
#   3. print.hbb_loo_compare --- S3 print method (exported)
#   4. .loo_extract_loglik   --- Extract M x N log_lik matrix (internal)
#   5. .loo_validate_inputs  --- Input validation (internal)
#   6. .loo_chain_id         --- Build chain_id vector from metadata (internal)
#   7. .loo_pareto_diagnostics --- Emit Pareto-k diagnostics via cli (internal)
#   8. .loo_pairwise_z       --- Pairwise z-ratios for consecutive models (internal)
# ============================================================================


# ============================================================================
# 1. loo.hbb_fit --- S3 loo method
# ============================================================================

#' Leave-One-Out Cross-Validation for Hurdle Beta-Binomial Models
#'
#' Computes Pareto-smoothed importance sampling LOO-CV (PSIS-LOO) for a
#' fitted \code{hbb_fit} object using the pointwise log-likelihood
#' \code{log_lik[N]} computed in the Stan generated quantities block.
#'
#' @description
#' PSIS-LOO approximates the leave-one-out predictive density:
#' \deqn{
#'   p(y_i \mid \mathbf{y}_{-i})
#'   = \int p(y_i \mid \theta) \, p(\theta \mid \mathbf{y}_{-i}) \, d\theta,
#' }
#' using importance sampling with weights proportional to
#' \eqn{1/p(y_i \mid \theta)}.  Raw importance weights have heavy
#' tails; PSIS stabilises them by fitting a generalised Pareto
#' distribution to the largest weights and replacing the upper tail.
#'
#' The expected log predictive density (ELPD) is:
#' \deqn{
#'   \widehat{\mathrm{elpd}}_{\mathrm{LOO}}
#'   = \sum_{i=1}^{N} \log p(y_i \mid \mathbf{y}_{-i}),
#' }
#' and the LOO information criterion is
#' \eqn{\mathrm{LOOIC} = -2 \cdot \widehat{\mathrm{elpd}}_{\mathrm{LOO}}}.
#' Models with higher (less negative) ELPD predict held-out data better.
#'
#' @section Pareto-k Diagnostic:
#' The Pareto shape parameter \eqn{k} for each observation provides a
#' diagnostic of LOO reliability:
#' \itemize{
#'   \item \eqn{k < 0.5}: Good --LOO estimate is reliable.
#'   \item \eqn{0.5 \le k < 0.7}: OK --LOO estimate is moderately
#'     reliable; the finite-moment condition is met.
#'   \item \eqn{k \ge 0.7}: Bad --LOO estimate may be unreliable; the
#'     importance weight distribution has infinite variance. Consider
#'     moment matching (\code{loo::loo_moment_match}) or \eqn{K}-fold CV.
#' }
#' Pareto-k values are automatically reported via \code{cli} upon
#' completion.
#'
#' @section Relative Effective Sample Size:
#' When \code{r_eff = TRUE} (default), the relative effective sample
#' size \eqn{r_{\mathrm{eff},i} = \mathrm{ESS}_i / M} is computed via
#' \code{\link[loo:relative_eff]{loo::relative_eff()}} and passed to
#' \code{\link[loo:loo]{loo::loo()}} to correct for within-chain
#' autocorrelation.  The correction uses the chain membership vector
#' constructed from CmdStanR metadata (draws are in row-major chain
#' order: chain 1 first, then chain 2, etc.).
#'
#' @section Survey-Weighted Models:
#' For survey-weighted models the \code{log_lik[N]} array in the Stan
#' GQ block contains the \emph{unweighted} nominal log-likelihoods
#' \eqn{\log p(y_i \mid \theta)}, not the pseudo-log-likelihoods.
#' LOO therefore provides model selection among alternative HBB
#' specifications on a common likelihood scale, but does not
#' account for survey design weighting.  For design-based inference
#' use the sandwich-corrected Wald confidence intervals from
#' \code{\link{sandwich_variance}}.
#'
#' @param x An object of class \code{"hbb_fit"} returned by
#'   \code{\link{hbb}}.  Must contain a CmdStanMCMC fit with the
#'   generated quantity \code{log_lik[N]}.
#' @param ... Additional arguments passed to
#'   \code{\link[loo:loo]{loo::loo()}} (e.g., \code{save_psis = TRUE}).
#' @param r_eff Logical; if \code{TRUE} (default), compute relative
#'   effective sample sizes to improve PSIS accuracy.  If
#'   \code{FALSE}, sets all \code{r_eff = 1}.
#' @param cores Integer; number of cores for parallel PSIS computation.
#'   Default is \code{1}.
#'
#' @return An object of class \code{"psis_loo"} as returned by
#'   \code{\link[loo:loo]{loo::loo()}}.  Key elements:
#' \describe{
#'   \item{\code{estimates}}{3 x 2 matrix with rows \code{elpd_loo},
#'     \code{p_loo}, \code{looic} and columns \code{Estimate},
#'     \code{SE}.}
#'   \item{\code{diagnostics}}{List with \code{pareto_k} (N-vector)
#'     and \code{n_eff} (N-vector).}
#'   \item{\code{pointwise}}{N x 5 matrix of per-observation LOO
#'     summaries.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
#' loo_result <- loo(fit)
#' print(loo_result)
#' }
#'
#' @seealso
#' \code{\link{hbb_loo_compare}} for multi-model LOO comparison,
#' \code{\link{ppc}} for posterior predictive checks,
#' \code{\link[loo:loo]{loo::loo}} for the underlying PSIS-LOO engine.
#'
#' @references
#' Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian
#' model evaluation using leave-one-out cross-validation and WAIC.
#' \emph{Statistics and Computing}, \strong{27}(5), 1413--1432.
#'
#' Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A.
#' (2019). Visualisation in Bayesian workflow. \emph{Journal of the
#' Royal Statistical Society: Series A}, \strong{182}(2), 389--402.
#'
#' @family model-checking
#' @rawNamespace S3method(loo::loo, hbb_fit)
#' @importFrom loo loo relative_eff loo_compare
loo.hbb_fit <- function(x, ..., r_eff = TRUE, cores = 1L) {

    # =========================================================================
    # 0. VALIDATE
    # =========================================================================

    .loo_validate_inputs(x)

    # -- Validate cores ----------------------------------------------------------
    checkmate::assert_count(cores, positive = TRUE)

    # -- Validate r_eff ----------------------------------------------------------
    if (!is.logical(r_eff) || length(r_eff) != 1L || is.na(r_eff)) {
        cli::cli_abort(c(
            "x" = "{.arg r_eff} must be a single non-NA logical value.",
            "i" = "Received: {.val {r_eff}}. Use {.code TRUE} or {.code FALSE}."
        ))
    }

    N <- x$hbb_data$N

    cli::cli_alert_info("Computing PSIS-LOO: N = {N}")

    # =========================================================================
    # 1. EXTRACT M x N LOG-LIKELIHOOD MATRIX
    # =========================================================================

    log_lik <- .loo_extract_loglik(x)
    M       <- nrow(log_lik)

    cli::cli_alert_info(
        "Extracted log_lik: {M} draws x {N} observations"
    )

    # =========================================================================
    # 2. RELATIVE EFFECTIVE SAMPLE SIZES
    # =========================================================================

    r_eff_val <- NULL

    if (isTRUE(r_eff)) {

        # Pass M so .loo_chain_id can validate length
        chain_id <- tryCatch(
            .loo_chain_id(x, M_expected = M),
            error = function(e) {
                cli::cli_warn(c(
                    "!" = "Could not construct chain_id: {conditionMessage(e)}",
                    "i" = "Proceeding without r_eff correction (r_eff = 1)."
                ))
                NULL
            }
        )

        if (!is.null(chain_id)) {

            if (length(chain_id) != M) {
                cli::cli_warn(c(
                    "!" = "chain_id length ({length(chain_id)}) != M ({M}); skipping r_eff."
                ))
                chain_id <- NULL
            }
        }

        if (!is.null(chain_id)) {
            r_eff_val <- tryCatch(
                loo::relative_eff(
                    exp(log_lik),
                    chain_id = chain_id,
                    cores    = cores
                ),
                error = function(e) {
                    cli::cli_warn(c(
                        "!" = "relative_eff() failed: {conditionMessage(e)}",
                        "i" = "Proceeding with r_eff = 1."
                    ))
                    NULL
                }
            )
        }
    }

    # =========================================================================
    # 3. COMPUTE PSIS-LOO
    # =========================================================================

    loo_obj <- loo::loo(
        log_lik,
        r_eff = r_eff_val,
        cores = cores,
        ...
    )

    # =========================================================================
    # 4. PARETO-K DIAGNOSTICS
    # =========================================================================

    model_name <- x$model_name %||% x$model_type %||% "hbb_fit"
    .loo_pareto_diagnostics(loo_obj, model_name)

    elpd_est <- tryCatch(
        loo_obj$estimates["elpd_loo", "Estimate"],
        error = function(e) NA_real_
    )
    cli::cli_alert_success(
        "LOO-CV complete (ELPD = {round(elpd_est, 1)})"
    )

    loo_obj
}


# ============================================================================
# 2. hbb_loo_compare --- Multi-model LOO comparison (exported)
# ============================================================================

#' Compare Hurdle Beta-Binomial Models via LOO-CV
#'
#' Computes LOO-CV for each supplied \code{hbb_fit} object and calls
#' \code{\link[loo:loo_compare]{loo::loo_compare()}} to rank models by
#' expected log predictive density (ELPD).  Also reports pairwise
#' z-statistics for consecutive model pairs.
#'
#' @description
#' Models are ranked from highest (best) to lowest (worst) ELPD.
#' The ELPD difference \eqn{\Delta \mathrm{elpd}_{ij} =
#' \mathrm{elpd}_i - \mathrm{elpd}_j} for the best model against each
#' other model is reported together with its standard error.
#'
#' For consecutive pairs in the ranked list, the pairwise z-ratio is:
#' \deqn{
#'   z_{AB} = \frac{\Delta\,\mathrm{elpd}_{AB}}
#'                 {\mathrm{SE}(\Delta\,\mathrm{elpd}_{AB})},
#' }
#' where \eqn{\mathrm{SE}} is estimated from the pointwise ELPD
#' differences via:
#' \deqn{
#'   \mathrm{SE}(\Delta\,\mathrm{elpd}_{AB})
#'   = \sqrt{N} \cdot \mathrm{SD}\bigl(\hat{\ell}_i^{(A)}
#'       - \hat{\ell}_i^{(B)}\bigr),
#' }
#' where \eqn{\hat{\ell}_i^{(A)}} is the per-observation LOO log-density
#' for model \eqn{A}.  A magnitude \eqn{|z| > 2} is taken as
#' substantial evidence for the higher-ranked model (Vehtari et al., 2017).
#'
#' @section Input Format:
#' Pass models as named arguments (e.g., \code{m0 = fit0, m1 = fit1, ...})
#' or as a single named list.  All models must share the same number of
#' observations \eqn{N}.
#'
#' @param ... Two or more named \code{hbb_fit} objects, or a single named
#'   list of \code{hbb_fit} objects.
#' @param cores Integer; number of cores for parallel PSIS computation.
#'   Default \code{1}.
#'
#' @return An S3 object of class \code{"hbb_loo_compare"} containing:
#' \describe{
#'   \item{\code{comparison}}{Data frame from
#'     \code{\link[loo:loo_compare]{loo::loo_compare()}}: models ranked
#'     best to worst by ELPD, with columns \code{model},
#'     \code{elpd_loo}, \code{elpd_diff}, \code{se_diff}, \code{looic}.}
#'   \item{\code{loo_list}}{Named list of \code{"psis_loo"} objects,
#'     one per model.}
#'   \item{\code{pairwise}}{Data frame of pairwise z-ratios for
#'     consecutive model pairs (in ELPD-ranked order), with columns
#'     \code{comparison}, \code{delta_elpd}, \code{se_diff},
#'     \code{z_ratio}.}
#'   \item{\code{pareto_k}}{Named list of Pareto-k diagnostic vectors
#'     (one per model, length \eqn{N}).}
#'   \item{\code{model_names}}{Character vector of model names (as
#'     supplied by the user).}
#'   \item{\code{n_models}}{Integer: number of models compared.}
#' }
#'
#' @examples
#' \dontrun{
#' fit0 <- hbb(y | trials(n) ~ 1,              data = my_data)
#' fit1 <- hbb(y | trials(n) ~ poverty,        data = my_data)
#' fit2 <- hbb(y | trials(n) ~ poverty + urban, data = my_data)
#' comp <- hbb_loo_compare(m0 = fit0, m1 = fit1, m2 = fit2)
#' print(comp)
#' }
#'
#' @seealso
#' \code{\link{loo.hbb_fit}} for single-model LOO-CV,
#' \code{\link{ppc}} for posterior predictive checks.
#'
#' @references
#' Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian
#' model evaluation using leave-one-out cross-validation and WAIC.
#' \emph{Statistics and Computing}, \strong{27}(5), 1413--1432.
#'
#' @family model-checking
#' @export
hbb_loo_compare <- function(..., cores = 1L) {

    # =========================================================================
    # 0. PARSE ... INTO A NAMED LIST OF hbb_fit OBJECTS
    # =========================================================================

    # -- Validate cores -----------------------------------------------------------
    checkmate::assert_count(cores, positive = TRUE)

    dots <- list(...)

    # Handle single-list call: hbb_loo_compare(list(m0 = fit0, m1 = fit1))
    if (length(dots) == 1L &&
        is.list(dots[[1L]]) &&
        !inherits(dots[[1L]], "hbb_fit")) {
        fits <- dots[[1L]]
    } else {
        fits <- dots
    }

    # -- Names required ----------------------------------------------------------
    nms <- names(fits)
    if (is.null(nms) || any(nchar(nms) == 0L)) {
        blank_pos <- if (!is.null(nms)) which(nchar(nms) == 0L) else seq_along(fits)
        cli::cli_abort(c(
            "x" = "All models passed to {.fn hbb_loo_compare} must be named.",
            "i" = "Unnamed positions: {.val {blank_pos}}.",
            "i" = "Use: {.code hbb_loo_compare(m0 = fit0, m1 = fit1, ...)}."
        ))
    }

    model_names <- nms
    n_models    <- length(fits)

    # -- At least two models --------------------------------------------------
    if (n_models < 2L) {
        cli::cli_abort(c(
            "x" = "{.fn hbb_loo_compare} requires at least 2 models.",
            "i" = "Received {n_models} model{?s}."
        ))
    }

    # -- All must be hbb_fit --------------------------------------------------
    not_hbb <- !vapply(fits, inherits, logical(1L), what = "hbb_fit")
    if (any(not_hbb)) {
        bad_names <- model_names[not_hbb]
        bad_pos   <- which(not_hbb)
        cli::cli_abort(c(
            "x" = "All arguments must be {.cls hbb_fit} objects.",
            "i" = "Non-{.cls hbb_fit} arguments: {.val {bad_names}} at positions {.val {bad_pos}}.",
            "i" = "Use {.fn hbb} to create {.cls hbb_fit} objects."
        ))
    }

    # -- All must share the same N ------------------------------------------------
    N_vals <- vapply(fits, function(f) {
        n <- tryCatch(f$hbb_data$N, error = function(e) NA_integer_)
        if (is.null(n)) NA_integer_ else as.integer(n)
    }, integer(1L))

    if (anyNA(N_vals)) {
        na_names <- model_names[is.na(N_vals)]
        cli::cli_abort(c(
            "x" = "Could not retrieve {.code hbb_data$N} for {length(na_names)} model{?s}.",
            "i" = "Affected model{?s}: {.val {na_names}}.",
            "i" = "The {.cls hbb_fit} object{?s} may be incomplete."
        ))
    }

    if (length(unique(N_vals)) > 1L) {
        N_report <- paste(model_names, N_vals, sep = " (N=", collapse = "), ")
        N_report <- paste0(N_report, ")")
        cli::cli_abort(c(
            "x" = "All models must have the same number of observations N.",
            "i" = "Observed: {N_report}.",
            "i" = "Ensure all models were fit on the same dataset."
        ))
    }

    cli::cli_alert_info(
        "Comparing {n_models} models: {paste(model_names, collapse = ', ')}"
    )

    # =========================================================================
    # 1. COMPUTE LOO FOR EACH MODEL
    # =========================================================================

    loo_list <- vector("list", length = n_models)
    names(loo_list) <- model_names

    for (nm in model_names) {
        cli::cli_alert_info("Computing LOO for model {.val {nm}}...")
        loo_list[[nm]] <- loo.hbb_fit(fits[[nm]], cores = cores)
    }

    # =========================================================================
    # 2. LOO COMPARISON
    # =========================================================================

    comparison_raw <- loo::loo_compare(loo_list)
    # loo_compare() returns a matrix; rows are ranked model names
    comparison_df  <- as.data.frame(comparison_raw)
    comparison_df  <- cbind(
        model            = rownames(comparison_df),
        comparison_df,
        stringsAsFactors = FALSE,
        row.names        = NULL
    )
    ranked_names <- rownames(comparison_raw)   # best ->worst

    # =========================================================================
    # 3. PAIRWISE Z-RATIOS
    # =========================================================================

    pairwise <- .loo_pairwise_z(loo_list, ranked_names)

    # =========================================================================
    # 4. PARETO-K SUMMARY
    # =========================================================================

    pareto_k <- lapply(loo_list, function(lo) {
        lo$diagnostics$pareto_k
    })
    names(pareto_k) <- model_names

    cli::cli_alert_success("LOO comparison complete")

    # =========================================================================
    # 5. RETURN OBJECT
    # =========================================================================

    structure(
        list(
            comparison  = comparison_df,
            loo_list    = loo_list,
            pairwise    = pairwise,
            pareto_k    = pareto_k,
            model_names = model_names,
            n_models    = n_models
        ),
        class = "hbb_loo_compare"
    )
}


# ============================================================================
# 3. print.hbb_loo_compare --- S3 print method (exported)
# ============================================================================

#' Print Method for hbb_loo_compare Objects
#'
#' Displays a formatted LOO comparison table including ELPD differences,
#' standard errors, and pairwise z-statistics.
#'
#' @description
#' The display is structured into three sections:
#' \enumerate{
#'   \item A model ranking table with \eqn{\mathrm{ELPD}},
#'     \eqn{\Delta\mathrm{ELPD}}, \eqn{\mathrm{SE}(\Delta)}, and
#'     \eqn{\mathrm{LOOIC}}.
#'   \item A pairwise z-ratio table for consecutive model pairs.
#'   \item A Pareto-k diagnostic summary (counts by category).
#' }
#'
#' @param x An object of class \code{"hbb_loo_compare"} returned by
#'   \code{\link{hbb_loo_compare}}.
#' @param digits Integer; decimal places for numeric values.
#'   Default is \code{1}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The object \code{x}, invisibly.
#'
#' @seealso \code{\link{hbb_loo_compare}}, \code{\link{loo.hbb_fit}}
#'
#' @family model-checking
#' @method print hbb_loo_compare
#' @export
print.hbb_loo_compare <- function(x, digits = 1L, ...) {

    # -- Safe-field helper: returns fallback if any expression errors -----------
    .sf <- function(expr_fn, fallback) {
        tryCatch(expr_fn(), error = function(e) fallback)
    }

    # -- Outer guard: last-resort fallback if everything else fails ------------
    tryCatch({

        cat("\n")
        cat(strrep("=", 70), "\n")
        cat("  LOO-CV Model Comparison (HBB Models)\n")
        cat(strrep("=", 70), "\n\n")

        n_models    <- .sf(function() as.integer(x$n_models),     NA_integer_)
        model_names <- .sf(function() as.character(x$model_names), character(0L))

        if (!is.na(n_models)) {
            cat(sprintf("  Models compared: %d\n\n", n_models))
        }

        # -- Main comparison table (ELPD ranking) ------------------------------
        df <- .sf(function() x$comparison, NULL)

        if (!is.null(df) && is.data.frame(df) && nrow(df) > 0L) {

            has_model     <- "model"     %in% names(df)
            has_elpd_loo  <- "elpd_loo"  %in% names(df)
            has_elpd_diff <- "elpd_diff" %in% names(df)
            has_se_diff   <- "se_diff"   %in% names(df)
            has_looic     <- "looic"     %in% names(df)

            cat(strrep("-", 70), "\n")
            cat(sprintf("  %-14s  %10s  %10s  %10s  %10s\n",
                        "Model", "ELPD_LOO", "ELPD_diff", "SE_diff", "LOOIC"))
            cat(strrep("-", 70), "\n")

            for (i in seq_len(nrow(df))) {
                tryCatch({
                    r <- df[i, ]
                    model_nm  <- if (has_model)     .sf(function() as.character(r[["model"]]),    "(?)") else "(?)"
                    elpd_loo  <- if (has_elpd_loo)  .sf(function() as.numeric(r[["elpd_loo"]]),  NA_real_) else NA_real_
                    elpd_diff <- if (has_elpd_diff) .sf(function() as.numeric(r[["elpd_diff"]]), NA_real_) else NA_real_
                    se_diff   <- if (has_se_diff)   .sf(function() as.numeric(r[["se_diff"]]),   NA_real_) else NA_real_
                    looic_val <- if (has_looic)     .sf(function() as.numeric(r[["looic"]]),     NA_real_) else NA_real_

                    # Format NA as "NA" (right-aligned within the field width)
                    fmt_num <- function(v, d) {
                        if (is.na(v)) formatC("NA", width = 10, flag = " ") else
                            sprintf("%10.*f", d, v)
                    }
                    cat(sprintf("  %-14s  %s  %s  %s  %s\n",
                                model_nm,
                                fmt_num(elpd_loo,  digits),
                                fmt_num(elpd_diff, digits),
                                fmt_num(se_diff,   digits),
                                fmt_num(looic_val, digits)))
                }, error = function(e) {
                    cat(sprintf("  (row %d: error formatting --%s)\n",
                                i, conditionMessage(e)))
                })
            }
            cat(strrep("-", 70), "\n\n")

        } else {
            cat("  (comparison table unavailable)\n\n")
        }

        # -- Pairwise z-ratios ------------------------------------------------
        pw <- .sf(function() x$pairwise, NULL)
        if (!is.null(pw) && is.data.frame(pw) && nrow(pw) > 0L) {
            cat("  Pairwise z-ratios (consecutive pairs, ranked by ELPD):\n")
            cat(strrep("-", 70), "\n")
            cat(sprintf("  %-26s  %10s  %10s  %10s\n",
                        "Comparison", "Delta_ELPD", "SE", "z-ratio"))
            cat(strrep("-", 70), "\n")

            fmt_num <- function(v, d) {
                if (is.na(v)) formatC("NA", width = 10, flag = " ") else
                    sprintf("%10.*f", d, v)
            }

            for (i in seq_len(nrow(pw))) {
                tryCatch({
                    r <- pw[i, ]
                    cmp_lbl  <- .sf(function() as.character(r[["comparison"]]), "(?)")
                    d_elpd   <- .sf(function() as.numeric(r[["delta_elpd"]]),  NA_real_)
                    se_d     <- .sf(function() as.numeric(r[["se_diff"]]),     NA_real_)
                    zr       <- .sf(function() as.numeric(r[["z_ratio"]]),     NA_real_)
                    sig_flag <- if (!is.na(zr) && abs(zr) > 2) " *" else "  "
                    cat(sprintf("  %-26s  %s  %s  %s%s\n",
                                cmp_lbl,
                                fmt_num(d_elpd, digits),
                                fmt_num(se_d,   digits),
                                fmt_num(zr,     digits),
                                sig_flag))
                }, error = function(e) {
                    cat(sprintf("  (row %d: error formatting --%s)\n",
                                i, conditionMessage(e)))
                })
            }
            cat(strrep("-", 70), "\n")
            cat("  * |z| > 2 indicates substantial preference for the higher-ranked model.\n\n")
        }

        # -- Pareto-k summary -------------------------------------------------
        cat("  Pareto-k diagnostics (k < 0.5 good; 0.5-0.7 ok; >= 0.7 bad):\n")
        cat(strrep("-", 70), "\n")

        pk_list <- .sf(function() x$pareto_k, list())
        if (length(model_names) == 0L) model_names <- names(pk_list)

        for (nm in model_names) {
            tryCatch({
                pk <- .sf(function() pk_list[[nm]], NULL)
                if (!is.null(pk) && length(pk) > 0L) {
                    n_tot  <- length(pk)
                    n_good <- sum(pk <  0.5,              na.rm = TRUE)
                    n_ok   <- sum(pk >= 0.5 & pk < 0.7,  na.rm = TRUE)
                    n_bad  <- sum(pk >= 0.7,              na.rm = TRUE)
                    pct_bad <- round(100 * n_bad / n_tot, 1)
                    cat(sprintf(
                        "  %-14s  good: %d/%d  ok: %d/%d  bad: %d/%d (%.1f%%)\n",
                        nm, n_good, n_tot, n_ok, n_tot, n_bad, n_tot, pct_bad
                    ))
                } else {
                    cat(sprintf("  %-14s  (no Pareto-k data)\n", nm))
                }
            }, error = function(e) {
                cat(sprintf("  %-14s  (error: %s)\n", nm, conditionMessage(e)))
            })
        }
        cat(strrep("-", 70), "\n\n")

    }, error = function(e) {
        # Last-resort fallback: never crash, always show something
        cat("Error printing hbb_loo_compare:", conditionMessage(e), "\n")
        cat("Object class:", paste(
            tryCatch(class(x),  error = function(e2) "?"), collapse = ", "), "\n")
        cat("Object names:", paste(
            tryCatch(names(x),  error = function(e2) "?"), collapse = ", "), "\n")
    })

    invisible(x)
}


# ============================================================================
# Internal helpers
# ============================================================================


# ---- 4. .loo_extract_loglik ------------------------------------------------

#' Extract M x N log_lik matrix from a CmdStanMCMC fit
#'
#' Reads the \code{log_lik[1:N]} generated quantity from the CmdStanR
#' object and returns it as an \eqn{M \times N} numeric matrix, where
#' \eqn{M} is the total number of post-warmup draws across all chains.
#'
#' @details
#' The log-likelihood extracted here is the \emph{unweighted} nominal
#' log-likelihood \eqn{\log p(y_i \mid \theta^{(m)})}, not the
#' survey-weighted pseudo-likelihood.  This is appropriate for LOO as
#' a model-selection tool among HBB specifications.
#'
#' Pathological values (NA, NaN, Inf) are detected and reported via
#' \code{cli_warn()}; they are not replaced and will propagate to
#' \code{loo::loo()}, where they manifest as very high Pareto-k values.
#'
#' @param fit An \code{hbb_fit} object with a valid CmdStanMCMC fit.
#' @return A numeric matrix of dimension \eqn{M \times N}.
#' @noRd
.loo_extract_loglik <- function(fit) {

    N <- fit$hbb_data$N

    log_lik <- tryCatch(
        fit$fit$draws(variables = "log_lik", format = "matrix"),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "Failed to extract {.code log_lik} from CmdStanMCMC fit.",
                "i" = "Error: {conditionMessage(e)}",
                "i" = "Ensure the Stan model includes {.code log_lik[N]} in generated quantities."
            ))
        }
    )

    if (!is.matrix(log_lik)) log_lik <- as.matrix(log_lik)

    # -- Dimension check ---------------------------------------------------------
    if (ncol(log_lik) != N) {
        cli::cli_abort(c(
            "x" = "Extracted {.code log_lik} has {ncol(log_lik)} column{?s} but {.code hbb_data$N} = {N}.",
            "i" = "Ensure the Stan model was fit with the same dataset (N observations).",
            "i" = "Actual dimensions: {nrow(log_lik)} draws x {ncol(log_lik)} columns."
        ))
    }

    # -- Pathological value checks ------------------------------------------------
    # Note: is.nan() must precede is.na() since NaN is also NA in R.
    n_total    <- length(log_lik)
    n_nan      <- sum(is.nan(log_lik))
    n_na_only  <- sum(is.na(log_lik) & !is.nan(log_lik))  # NA but not NaN
    n_pos_inf  <- sum(is.infinite(log_lik) & log_lik > 0)
    n_neg_inf  <- sum(is.infinite(log_lik) & log_lik < 0)

    if (n_nan > 0L) {
        pct <- round(100 * n_nan / n_total, 2)
        cli::cli_warn(c(
            "!" = "{n_nan} NaN value{?s} ({pct}% of all {.code log_lik} entries).",
            "i" = "NaN log-likelihoods produce undefined PSIS importance weights.",
            "i" = "These observations will receive very high Pareto-k values.",
            "i" = "Check your model specification and data for numerical issues."
        ))
    }

    if (n_na_only > 0L) {
        pct <- round(100 * n_na_only / n_total, 2)
        cli::cli_warn(c(
            "!" = "{n_na_only} NA value{?s} ({pct}% of all {.code log_lik} entries).",
            "i" = "NA log-likelihoods will be treated as {.code -Inf} by {.pkg loo}."
        ))
    }

    if (n_pos_inf > 0L) {
        pct <- round(100 * n_pos_inf / n_total, 2)
        cli::cli_warn(c(
            "!" = "{n_pos_inf} +Inf value{?s} ({pct}% of all {.code log_lik} entries).",
            "i" = "+Inf log-likelihoods are physically impossible and indicate a model error.",
            "i" = "These entries have infinite importance weights and will cause LOO to fail."
        ))
    }

    if (n_neg_inf > 0L) {
        pct <- round(100 * n_neg_inf / n_total, 2)
        # -Inf is less alarming (zero-probability events) but still warrants notice
        if (pct > 1) {
            cli::cli_warn(c(
                "!" = "{n_neg_inf} -Inf value{?s} ({pct}% of all {.code log_lik} entries).",
                "i" = "-Inf log-likelihoods correspond to zero-probability observations.",
                "i" = "A large proportion may indicate a model-data mismatch."
            ))
        }
    }

    log_lik
}


# ---- 5. .loo_validate_inputs -----------------------------------------------

#' Validate inputs for loo.hbb_fit()
#'
#' Checks class, CmdStanMCMC availability, \code{hbb_data$N},
#' and (via metadata inspection) whether \code{log_lik} exists in the
#' Stan output.
#'
#' @param fit Object to validate as \code{hbb_fit}.
#' @return Invisible \code{NULL}; aborts on failure.
#' @noRd
.loo_validate_inputs <- function(fit) {

    # -- Class check -----------------------------------------------------------
    if (!inherits(fit, "hbb_fit")) {
        cli::cli_abort(c(
            "x" = "{.arg fit} must be of class {.cls hbb_fit}.",
            "i" = "Received class: {.cls {class(fit)}}.",
            "i" = "Use {.fn hbb} to fit a hurdle Beta-Binomial model."
        ))
    }

    # -- fit$fit presence ------------------------------------------------------
    if (is.null(fit$fit)) {
        cli::cli_abort(c(
            "x" = "{.code fit$fit} is {.val NULL}.",
            "i" = "The {.cls hbb_fit} must contain a valid CmdStanMCMC object."
        ))
    }

    # -- hbb_data$N presence ---------------------------------------------------
    if (is.null(fit$hbb_data) || is.null(fit$hbb_data$N)) {
        cli::cli_abort(c(
            "x" = "{.code fit$hbb_data$N} is missing.",
            "i" = "The {.cls hbb_fit} object appears incomplete."
        ))
    }

    # -- log_lik presence: first try metadata, then trial draw ----------------
    #    Strategy: metadata() is cheap but may be unavailable (e.g. loaded from
    #    RDS); trial draw is definitive but slow.  Try metadata first; only fall
    #    back to trial draw if metadata is unavailable.

    metadata_checked <- FALSE

    tryCatch({
        vars <- fit$fit$metadata()$stan_variables
        if (!is.null(vars)) {
            metadata_checked <- TRUE
            if (!any(grepl("^log_lik", vars))) {
                cli::cli_abort(c(
                    "x" = "No {.code log_lik} variable found in Stan output.",
                    "i" = "Ensure the Stan model includes {.code log_lik[N]} in generated quantities.",
                    "i" = "The UNWEIGHTED pointwise log-likelihood is required for PSIS-LOO."
                ))
            }
        }
    }, error = function(e) {
        # Re-throw rlang errors from our own cli_abort
        if (inherits(e, "rlang_error")) stop(e)
        # Metadata unavailable --warn and fall through to trial draw below
        cli::cli_warn(c(
            "!" = "Could not verify {.code log_lik} via Stan metadata.",
            "i" = "Metadata error: {conditionMessage(e)}",
            "i" = "Attempting trial extraction to verify {.code log_lik} availability."
        ))
    })

    # If metadata check was not conclusive, do a trial draw
    if (!metadata_checked) {
        trial_ok <- tryCatch({
            test_mat <- fit$fit$draws(variables = "log_lik", format = "matrix")
            !is.null(test_mat) && length(test_mat) > 0L
        }, error = function(e) FALSE)

        if (!trial_ok) {
            cli::cli_abort(c(
                "x" = "{.code log_lik} draws are not available in the Stan fit.",
                "i" = "Ensure the Stan model includes {.code log_lik[N]} in the generated quantities block.",
                "i" = "The UNWEIGHTED pointwise log-likelihood is required for PSIS-LOO."
            ))
        }
    }

    invisible(NULL)
}


# ---- 6. .loo_chain_id ------------------------------------------------------

#' Build the chain_id vector for loo::relative_eff()
#'
#' Constructs an integer vector of length
#' \eqn{M = n\_\mathrm{chains} \times \mathrm{iter\_sampling}}
#' where entry \eqn{m} identifies the chain that produced draw \eqn{m}.
#'
#' @details
#' CmdStanR stores draws in row-major order: all \code{iter_sampling}
#' draws from chain 1 first, then chain 2, etc.  The chain_id vector
#' is therefore \code{rep(1:n_chains, each = iter_sampling)}, which
#' matches the row ordering of the matrix returned by
#' \code{fit$fit$draws(format = "matrix")}.
#'
#' If CmdStanR metadata is unavailable (e.g. the fit was loaded from an
#' RDS file that did not preserve the CmdStan output directory), the
#' function falls back to assigning all \eqn{M} draws to chain 1 and
#' emits a \code{cli_warn}.  This is equivalent to \code{r_eff = 1} for
#' the ESS correction and is conservative but always valid.
#'
#' @param fit An \code{hbb_fit} object.
#' @param M_expected Integer; the number of rows in the log_lik matrix.
#'   If provided, the constructed chain_id is validated against this
#'   length and a \code{cli_warn} is emitted on mismatch.
#' @return Integer vector of length \eqn{M = n\_\mathrm{chains} \times
#'   \mathrm{iter\_sampling}} (or \eqn{M\_\mathrm{expected}} in the
#'   single-chain fallback).
#' @noRd
.loo_chain_id <- function(fit, M_expected = NULL) {

    # -- Try to extract n_chains and iter_sampling from CmdStanR metadata -----
    n_chains      <- tryCatch(fit$fit$num_chains(), error = function(e) NULL)
    iter_sampling <- tryCatch({
        meta <- fit$fit$metadata()
        meta$iter_sampling
    }, error = function(e) NULL)

    # -- Fallback: single chain of length M_expected --------------------------
    if (is.null(n_chains) || is.null(iter_sampling) ||
        is.na(iter_sampling) || is.na(n_chains)) {

        if (!is.null(M_expected) && M_expected > 0L) {
            cli::cli_warn(c(
                "!" = "CmdStanR metadata unavailable; falling back to single-chain assignment.",
                "i" = "All {M_expected} draw{?s} are assigned to chain 1.",
                "i" = "This is equivalent to {.code r_eff = 1} (conservative but valid).",
                "i" = "Metadata details: n_chains={.val {n_chains}}, ",
                      "iter_sampling={.val {iter_sampling}}."
            ))
            return(rep(1L, M_expected))
        } else {
            cli::cli_abort(c(
                "x" = "Could not determine chain structure from CmdStanR metadata.",
                "i" = "n_chains: {.val {n_chains}}, iter_sampling: {.val {iter_sampling}}.",
                "i" = "Provide {.arg M_expected} to enable the single-chain fallback."
            ))
        }
    }

    chain_id <- rep(seq_len(n_chains), each = iter_sampling)

    # -- Length validation against M_expected ---------------------------------
    if (!is.null(M_expected) && length(chain_id) != M_expected) {
        cli::cli_warn(c(
            "!" = "chain_id length ({length(chain_id)}) != M_expected ({M_expected}).",
            "i" = "Derived from n_chains={n_chains}, iter_sampling={iter_sampling}.",
            "i" = "chain_id will be returned as-is; caller should validate length."
        ))
    }

    chain_id
}


# ---- 7. .loo_pareto_diagnostics --------------------------------------------

#' Emit Pareto-k diagnostic messages via cli
#'
#' Reports counts of Pareto-k values in three quality tiers (good,
#' ok, bad) and issues a \code{cli_warn} for models with any bad
#' observations.
#'
#' @section Pareto-k Thresholds:
#' The thresholds follow Vehtari et al. (2017):
#' \describe{
#'   \item{\eqn{k < 0.5}}{Good: LOO estimate is reliable.}
#'   \item{\eqn{0.5 \le k < 0.7}}{OK: moderate reliability;
#'     finite-variance condition is met.}
#'   \item{\eqn{k \ge 0.7}}{Bad: LOO may be unreliable; the
#'     importance weight distribution has infinite variance.
#'     Consider \code{loo::loo_moment_match} or \eqn{K}-fold CV.}
#' }
#'
#' @param loo_obj A \code{"psis_loo"} object.
#' @param model_name Character string: name for display purposes.
#' @return Invisible \code{NULL}.
#' @noRd
.loo_pareto_diagnostics <- function(loo_obj, model_name) {

    # -- NULL / empty guard ------------------------------------------------------
    if (is.null(loo_obj)) {
        cli::cli_warn(c(
            "!" = "Cannot report Pareto-k diagnostics: {.arg loo_obj} is {.val NULL}."
        ))
        return(invisible(NULL))
    }

    pk <- tryCatch(
        loo_obj$diagnostics$pareto_k,
        error = function(e) NULL
    )
    if (is.null(pk) || length(pk) == 0L) {
        cli::cli_warn(c(
            "!" = "{model_name}: Pareto-k diagnostics are unavailable.",
            "i" = "{.code loo_obj$diagnostics$pareto_k} is NULL or empty."
        ))
        return(invisible(NULL))
    }

    n_total  <- length(pk)
    n_good   <- sum(pk <  0.5,              na.rm = TRUE)
    n_ok     <- sum(pk >= 0.5 & pk < 0.7,  na.rm = TRUE)
    n_bad    <- sum(pk >= 0.7,              na.rm = TRUE)
    pct_good <- round(100 * n_good  / n_total, 1)
    pct_ok   <- round(100 * n_ok   / n_total, 1)
    pct_bad  <- round(100 * n_bad  / n_total, 1)

    # Summary line with both counts and percentages
    cli::cli_inform(c(
        "i" = paste0("{model_name} Pareto-k: ",
                     "{n_good}/{n_total} ({pct_good}%) good (<0.5), ",
                     "{n_ok}/{n_total} ({pct_ok}%) ok (0.5-0.7), ",
                     "{n_bad}/{n_total} ({pct_bad}%) bad (>=0.7)")
    ))

    if (n_bad > 0L) {
        bad_idx  <- which(pk >= 0.7)
        n_show   <- min(n_bad, 10L)
        idx_str  <- paste(utils::head(bad_idx, n_show), collapse = ", ")
        tail_str <- if (n_bad > 10L) ", ..." else ""
        cli::cli_warn(c(
            "!" = "{model_name}: {n_bad} ({pct_bad}%) observation{?s} with Pareto-k >= 0.7.",
            "i" = "First {n_show} bad index{?es}: {idx_str}{tail_str}.",
            "i" = "LOO estimates may be unreliable for these observations.",
            "i" = "Consider {.fn loo::loo_moment_match} or K-fold CV."
        ))
    } else if (n_ok > 0L) {
        cli::cli_inform(c(
            "v" = "{model_name}: LOO is reliable (no bad Pareto-k values).",
            "i" = "{n_ok} ({pct_ok}%) observation{?s} in the moderate range (0.5-0.7)."
        ))
    } else {
        cli::cli_inform(c(
            "v" = "{model_name}: All Pareto-k < 0.5 ({pct_good}%). LOO is reliable."
        ))
    }

    invisible(NULL)
}


# ---- 8. .loo_pairwise_z ----------------------------------------------------

#' Compute pairwise z-ratios for consecutive LOO model pairs
#'
#' For each consecutive pair in the ELPD-ranked model list, computes
#' the z-ratio of the ELPD difference divided by its standard error.
#' The SE is estimated from the pointwise ELPD differences using the
#' standard formula:
#' \deqn{
#'   \mathrm{SE}(\Delta\,\mathrm{elpd})
#'   = \sqrt{N} \cdot \mathrm{SD}\bigl(\hat{\ell}_i^{(A)}
#'       - \hat{\ell}_i^{(B)}\bigr).
#' }
#'
#' @param loo_list Named list of \code{"psis_loo"} objects.
#' @param ranked_names Character vector of model names ordered from
#'   best to worst ELPD (as returned by
#'   \code{rownames(loo::loo_compare(...))}).
#' @return Data frame with columns: \code{comparison} (character),
#'   \code{delta_elpd} (numeric), \code{se_diff} (numeric),
#'   \code{z_ratio} (numeric).
#' @noRd
.loo_pairwise_z <- function(loo_list, ranked_names) {

    n <- length(ranked_names)

    if (n < 2L) {
        return(data.frame(
            comparison  = character(0L),
            delta_elpd  = numeric(0L),
            se_diff     = numeric(0L),
            z_ratio     = numeric(0L),
            stringsAsFactors = FALSE
        ))
    }

    n_pairs     <- n - 1L
    comparison  <- character(n_pairs)
    delta_elpd  <- numeric(n_pairs)
    se_diff_vec <- numeric(n_pairs)
    z_ratio_vec <- numeric(n_pairs)

    for (i in seq_len(n_pairs)) {

        nm_a <- ranked_names[i]          # better model
        nm_b <- ranked_names[i + 1L]     # worse model

        comparison[i] <- paste0(nm_a, " vs ", nm_b)

        elpd_a <- tryCatch(
            loo_list[[nm_a]]$estimates["elpd_loo", "Estimate"],
            error = function(e) NA_real_
        )
        elpd_b <- tryCatch(
            loo_list[[nm_b]]$estimates["elpd_loo", "Estimate"],
            error = function(e) NA_real_
        )
        delta_elpd[i] <- elpd_a - elpd_b

        # Per-observation ELPD difference for SE estimation
        ll_a <- tryCatch(
            loo_list[[nm_a]]$pointwise[, "elpd_loo"],
            error = function(e) NULL
        )
        ll_b <- tryCatch(
            loo_list[[nm_b]]$pointwise[, "elpd_loo"],
            error = function(e) NULL
        )

        if (!is.null(ll_a) && !is.null(ll_b) &&
            length(ll_a) == length(ll_b) && length(ll_a) > 1L) {
            N_obs           <- length(ll_a)
            diff_pw         <- ll_a - ll_b
            se_diff_vec[i]  <- sqrt(N_obs) * sd(diff_pw)
        } else {
            se_diff_vec[i] <- NA_real_
        }

        z_ratio_vec[i] <- if (!is.na(se_diff_vec[i]) &&
                               se_diff_vec[i] > 0) {
            delta_elpd[i] / se_diff_vec[i]
        } else {
            NA_real_
        }
    }

    data.frame(
        comparison  = comparison,
        delta_elpd  = delta_elpd,
        se_diff     = se_diff_vec,
        z_ratio     = z_ratio_vec,
        stringsAsFactors = FALSE,
        row.names   = NULL
    )
}
