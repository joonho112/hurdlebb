# ============================================================================
# hbb-fit-class.R --- S3 Class Definition for hbb_fit Objects
#
# Defines the hbb_fit class returned by hbb(). Provides type-checking
# and a compact print method that summarises the fitted model without
# printing full parameter estimates. All diagnostic extractions are
# wrapped in tryCatch so that print() never fails, even if the
# CmdStanMCMC object is corrupted or output files have been moved.
#
# Contents:
#   1. is.hbb_fit    --- Check if an object is an hbb_fit (exported)
#   2. print.hbb_fit --- Compact summary print method (S3 method)
# ============================================================================


# ============================================================================
# 1. is.hbb_fit --- Check if an object is an hbb_fit (exported)
# ============================================================================

#' Test if an Object is an hbb_fit
#'
#' Returns `TRUE` if `x` inherits from class `"hbb_fit"`, `FALSE`
#' otherwise. This is the recommended way to check whether an object
#' was produced by [hbb()].
#'
#' @param x Any R object.
#' @return Logical scalar.
#'
#' @examples
#' is.hbb_fit(1)          # FALSE
#' is.hbb_fit(list())     # FALSE
#'
#' \dontrun{
#' fit <- hbb(y | trials(n_trial) ~ poverty, data = my_data)
#' is.hbb_fit(fit)        # TRUE
#' }
#'
#' @seealso [hbb()] for fitting models.
#' @family fitting
#' @export
is.hbb_fit <- function(x) {
    inherits(x, "hbb_fit")
}


# ============================================================================
# 2. print.hbb_fit --- Compact summary print method (S3 method)
# ============================================================================

#' Print Method for hbb_fit Objects
#'
#' Displays a compact summary of a fitted hurdle Beta-Binomial model,
#' including the model type, formula, data dimensions, MCMC
#' configuration, elapsed time, and basic diagnostic indicators.
#' Does **not** print parameter estimates; use `summary()` for that.
#'
#' All diagnostic extractions are defensive (wrapped in `tryCatch`)
#' so printing never fails, even if the underlying CmdStan output
#' files have been moved or the fit object is incomplete.
#'
#' @param x An object of class `"hbb_fit"`, as returned by [hbb()].
#' @param ... Currently unused; included for S3 method consistency.
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' \dontrun{
#' fit <- hbb(y | trials(n_trial) ~ poverty, data = my_data)
#' print(fit)
#' }
#'
#' @seealso [hbb()], [is.hbb_fit()]
#' @family fitting
#' @method print hbb_fit
#' @export
print.hbb_fit <- function(x, ...) {

    # -- Header ---------------------------------------------------------------
    cat("Hurdle Beta-Binomial Model Fit\n")
    cat("==============================\n\n")

    # -- Model type -----------------------------------------------------------
    type_labels <- c(
        base         = "Base (no random slopes, unweighted)",
        weighted     = "Weighted (no random slopes, survey-weighted)",
        svc          = "SVC (state-varying coefficients, unweighted)",
        svc_weighted = "SVC (state-varying coefficients, survey-weighted)"
    )
    type_label <- if (!is.null(x$model_type)) type_labels[x$model_type] else NA_character_
    if (length(type_label) == 0L || is.na(type_label)) type_label <- x$model_type %||% "(unknown)"
    cat("  Model type   :", type_label, "\n")
    cat("  Stan model   :", x$model_name %||% "(unknown)", "\n")

    # -- Formula --------------------------------------------------------------
    formula_text <- tryCatch(
        deparse(x$formula$formula, width.cutoff = 80),
        error = function(e) "(formula unavailable)"
    )
    formula_text <- paste(formula_text, collapse = "\n                  ")
    cat("  Formula      :", formula_text, "\n")

    # -- Dimensions -----------------------------------------------------------
    hd <- x$hbb_data
    if (!is.null(hd)) {
        cat("\n  Observations (N) :", hd$N, "\n")
        cat("  Covariates   (P) :", hd$P,
            paste0("(intercept + ", hd$P - 1L, " predictors)"), "\n")

        if (x$model_type %in% c("svc", "svc_weighted")) {
            cat("  Groups       (S) :", hd$S, "\n")
            if (!is.null(hd$Q)) cat("  Policy vars  (Q) :", hd$Q, "\n")
            if (!is.null(hd$K)) cat("  RE dimension (K) :", hd$K, "\n")
        }

        zero_rate <- round(1 - mean(hd$z), 3)
        cat("  Zero rate        :", zero_rate, "\n")

        has_weights <- !is.null(hd$w_tilde)
        cat("  Survey weights   :", if (has_weights) "yes" else "no", "\n")
    } else {
        cat("\n  (data summary unavailable)\n")
    }

    # -- MCMC configuration ---------------------------------------------------
    cat("\nMCMC:\n")

    n_chains <- tryCatch(
        x$fit$num_chains(),
        error = function(e) NA_integer_
    )
    metadata <- tryCatch(
        x$fit$metadata(),
        error = function(e) NULL
    )
    n_warmup   <- if (!is.null(metadata)) metadata$iter_warmup   else NA_integer_
    n_sampling <- if (!is.null(metadata)) metadata$iter_sampling else NA_integer_

    cat("  Chains           :",
        if (is.na(n_chains)) "(unknown)" else n_chains, "\n")
    cat("  Warmup           :",
        if (is.na(n_warmup)) "(unknown)" else n_warmup, "\n")
    cat("  Sampling         :",
        if (is.na(n_sampling)) "(unknown)" else n_sampling, "\n")

    total_draws <- if (!is.na(n_chains) && !is.na(n_sampling)) {
        n_chains * n_sampling
    } else {
        NA_integer_
    }
    cat("  Total draws      :",
        if (is.na(total_draws)) "(unknown)" else total_draws, "\n")

    # -- Elapsed time (formatted: seconds / minutes / hours) ------------------
    elapsed <- x$elapsed
    if (!is.null(elapsed) && is.numeric(elapsed) && !is.na(elapsed)) {
        if (elapsed < 60) {
            time_str <- sprintf("%.1f seconds", elapsed)
        } else if (elapsed < 3600) {
            time_str <- sprintf("%.1f minutes", elapsed / 60)
        } else {
            time_str <- sprintf("%.1f hours", elapsed / 3600)
        }
        cat("  Elapsed          :", time_str, "\n")
    }

    # -- MCMC diagnostics -----------------------------------------------------
    diag <- tryCatch(
        x$fit$diagnostic_summary(quiet = TRUE),
        error = function(e) NULL
    )

    if (!is.null(diag)) {
        cat("\nDiagnostics:\n")

        # Divergences
        n_div <- tryCatch(sum(diag$num_divergent), error = function(e) NA)
        if (!is.na(n_div)) {
            cat("  Divergent transitions :", n_div,
                if (n_div > 0L) " [WARNING]" else " [OK]", "\n")
        }

        # Max treedepth
        n_max_td <- tryCatch(sum(diag$num_max_treedepth), error = function(e) NA)
        if (!is.na(n_max_td)) {
            cat("  Max treedepth hits    :", n_max_td,
                if (n_max_td > 0L) " [WARNING]" else " [OK]", "\n")
        }

        # E-BFMI
        ebfmi <- tryCatch(diag$ebfmi, error = function(e) NULL)
        if (!is.null(ebfmi) && is.numeric(ebfmi)) {
            cat("  E-BFMI                :",
                paste(round(ebfmi, 3), collapse = ", "), "\n")
            if (any(ebfmi < 0.2)) {
                cat("    ** WARNING: Low E-BFMI detected (< 0.2) **\n")
            }
        }
    }

    # -- Rhat and ESS from summary -------------------------------------------
    summ <- tryCatch(x$fit$summary(), error = function(e) NULL)

    if (!is.null(summ) &&
        "rhat" %in% names(summ) && "ess_bulk" %in% names(summ)) {

        max_rhat     <- max(summ$rhat, na.rm = TRUE)
        min_ess_bulk <- min(summ$ess_bulk, na.rm = TRUE)
        min_ess_tail <- if ("ess_tail" %in% names(summ)) {
            min(summ$ess_tail, na.rm = TRUE)
        } else {
            NA_real_
        }

        cat("  Max Rhat              :", round(max_rhat, 4),
            if (!is.na(max_rhat) && max_rhat > 1.01) " [WARNING]" else " [OK]",
            "\n")
        cat("  Min bulk ESS          :", round(min_ess_bulk, 0),
            if (!is.na(min_ess_bulk) && min_ess_bulk < 400) " [LOW]" else " [OK]",
            "\n")
        if (!is.na(min_ess_tail)) {
            cat("  Min tail ESS          :", round(min_ess_tail, 0),
                if (min_ess_tail < 400) " [LOW]" else " [OK]",
                "\n")
        }
    }

    cat("\nUse summary() for parameter estimates.\n")

    invisible(x)
}
