# ============================================================================
# plot.R --- Plot Methods for hbb_fit Objects
#
# Provides diagnostic and inferential plots for Hurdle Beta-Binomial
# model fits.  Four plot types are supported:
#
#   1. "coefficients"    --- Forest plot of fixed-effect estimates with CIs
#   2. "trace"           --- MCMC trace plots for selected parameters
#   3. "random_effects"  --- Caterpillar plots of state random effects (SVC)
#   4. "residuals"       --- Residuals-vs-fitted and QQ diagnostic panels
#
#
# Existing helpers reused (NOT redefined here):
#   .validate_hbb_fit_methods(x)      from methods.R
#   .validate_level(level)            from cholesky-transform.R
#   .make_stan_param_names(P)         from cholesky-transform.R
#   .build_param_labels(hbb_data)     from sandwich.R
#   .extract_delta_means(fit, S, K)   from sandwich.R
#   summary.hbb_fit(object, ...)      from methods.R
#   fitted.hbb_fit(object, ...)       from methods.R
#   residuals.hbb_fit(object, ...)    from methods.R
#
# Project colour palette:
#   #4393C3 (blue)   --- extensive margin
#   #D6604D (red)    --- intensive margin
#   #1B7837 (green)  --- reference lines / intervals
#   #762A83 (purple) --- dispersion
#
# Contents:
#   Exported:
#     1. plot.hbb_fit          --- S3 plot dispatcher
#   Internal:
#     2. .plot_coefficients    --- Forest plot of fixed effects
#     3. .plot_trace           --- MCMC trace plots
#     4. .plot_random_effects  --- Caterpillar plot of state REs
#     5. .plot_residuals       --- Two-panel residual diagnostics
# ============================================================================


# ============================================================================
# 1. plot.hbb_fit --- S3 plot dispatcher (exported)
# ============================================================================

#' Plot Diagnostics and Inferential Summaries for hbb_fit Objects
#'
#' Produces diagnostic and summary plots for fitted Hurdle Beta-Binomial
#' models.  Four plot types are available, selected via the \code{type}
#' argument.
#'
#' @description
#' This method dispatches to four specialised internal helpers:
#' \describe{
#'   \item{\code{"coefficients"}}{A forest plot of fixed-effect point
#'     estimates with confidence/credible intervals, faceted by model
#'     margin (extensive, intensive, dispersion).  Uses
#'     \code{\link{summary.hbb_fit}} internally, so the same inference
#'     mode (Wald or posterior) applies.}
#'   \item{\code{"trace"}}{MCMC trace plots showing the evolution of
#'     sampled parameter values across iterations, coloured by chain.
#'     Useful for diagnosing convergence and mixing.}
#'   \item{\code{"random_effects"}}{Caterpillar plots of the posterior
#'     mean and credible interval for each state-level random effect
#'     \eqn{\delta_{s,k}}, faceted by covariate.  Only available for
#'     SVC models (\code{model_type \%in\% c("svc", "svc_weighted")}).}
#'   \item{\code{"residuals"}}{A two-panel diagnostic combining
#'     residuals-versus-fitted and a normal QQ plot of Pearson
#'     residuals.  Panels are arranged via \pkg{patchwork}.}
#' }
#'
#' @section ggplot2 requirement:
#' This function requires \pkg{ggplot2} (and \pkg{patchwork} for
#' \code{type = "residuals"}).  Both are in \code{Suggests}; they are
#' checked via \code{rlang::check_installed()} and the user will be
#' prompted to install them if absent.
#'
#' @section Inference mode (coefficients plot):
#' When \code{sandwich} is supplied, the forest plot shows Wald
#' confidence intervals derived from the sandwich variance:
#' \deqn{
#'   \mathrm{CI}_{1-\alpha}(\theta_p)
#'     = \hat\theta_p \pm z_{(1+\mathrm{level})/2}\,
#'       \sqrt{V_{\mathrm{sand},pp}}.
#' }
#' When \code{sandwich = NULL}, quantile-based credible intervals from
#' the MCMC posterior are shown instead.
#'
#' @section Colour palette:
#' The project palette is used throughout:
#' \itemize{
#'   \item \code{#4393C3} (blue) --- extensive margin
#'   \item \code{#D6604D} (red) --- intensive margin
#'   \item \code{#1B7837} (green) --- reference lines, intervals
#'   \item \code{#762A83} (purple) --- dispersion parameter
#' }
#'
#' @param x An object of class \code{"hbb_fit"} returned by
#'   \code{\link{hbb}}.
#' @param type Character string: one of \code{"coefficients"} (default),
#'   \code{"trace"}, \code{"random_effects"}, or \code{"residuals"}.
#' @param sandwich An optional object of class \code{"hbb_sandwich"}
#'   returned by \code{\link{sandwich_variance}}.  Used only when
#'   \code{type = "coefficients"} to produce Wald-based intervals.
#'   If \code{NULL} (default), posterior-based intervals are shown.
#' @param level Numeric scalar in \eqn{(0, 1)}.  Confidence/credible
#'   level for interval construction.  Default is \code{0.95}.
#' @param margin Character string controlling which coefficients appear
#'   in the forest plot.  One of \code{"both"} (default),
#'   \code{"extensive"}, or \code{"intensive"}.  Ignored for other
#'   plot types.
#' @param pars Optional character vector of parameter names to include
#'   in the trace plot.  If \code{NULL} (default), all \eqn{D = 2P+1}
#'   fixed-effect parameters are plotted.  Ignored for other plot types.
#' @param ... Additional arguments passed to internal helpers (currently
#'   unused).
#'
#' @return A \code{ggplot} or \code{patchwork} object, returned
#'   invisibly.  The plot is printed as a side effect.
#'
#' @seealso
#' \code{\link{summary.hbb_fit}} for the tabular summary reused by the
#' coefficients plot,
#' \code{\link{residuals.hbb_fit}} for residual computation,
#' \code{\link{fitted.hbb_fit}} for fitted values,
#' \code{\link{ppc}} and \code{\link{plot.hbb_ppc}} for posterior
#' predictive checks.
#'
#' @references
#' Ghosal, R., Ghosh, S. K., and Maiti, T. (2020).
#' Two-part regression models for longitudinal zero-inflated count data.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' \strong{183}(4), 1603--1626.
#'
#' Williams, M. R. and Savitsky, T. D. (2021).
#' Uncertainty estimation for pseudo-Bayesian inference under complex
#' sampling.
#' \emph{International Statistical Review}, \strong{89}(1), 72--107.
#' \doi{10.1111/insr.12376}
#'
#' @examples
#' \dontrun{
#' fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
#'
#' # Forest plot with posterior intervals
#' plot(fit)
#'
#' # Forest plot with sandwich (Wald) intervals
#' sand <- sandwich_variance(fit)
#' plot(fit, type = "coefficients", sandwich = sand)
#'
#' # Trace plots for a subset of parameters
#' plot(fit, type = "trace", pars = c("alpha[1]", "beta[1]"))
#'
#' # Residual diagnostics
#' plot(fit, type = "residuals")
#'
#' # Random effects caterpillar (SVC models only)
#' plot(fit, type = "random_effects")
#' }
#'
#' @method plot hbb_fit
#' @export
plot.hbb_fit <- function(x,
                          type = c("coefficients", "trace",
                                   "random_effects", "residuals"),
                          sandwich = NULL,
                          level = 0.95,
                          margin = c("both", "extensive", "intensive"),
                          pars = NULL,
                          ...) {

    # -- Top-level tryCatch: NEVER crash from plot ------------------------------
    tryCatch({

        # -- Dependency guard --------------------------------------------------
        rlang::check_installed(
            "ggplot2",
            reason = "to plot hbb_fit objects"
        )

        # -- Input validation --------------------------------------------------
        .validate_hbb_fit_methods(x)
        type   <- match.arg(type)
        margin <- match.arg(margin)

        # Validate level (reuse existing helper from cholesky-transform.R)
        .validate_level(level)

        # -- Dispatch ----------------------------------------------------------
        out <- switch(type,
            coefficients   = .plot_coefficients(x, sandwich = sandwich,
                                                 level = level,
                                                 margin = margin),
            trace          = .plot_trace(x, pars = pars),
            random_effects = .plot_random_effects(x, level = level),
            residuals      = .plot_residuals(x)
        )

        print(out)
        invisible(out)

    }, error = function(e) {
        # Defensive fallback: NEVER crash from plot
        cli::cli_warn(c(
            "!" = "plot.hbb_fit encountered an error.",
            "x" = conditionMessage(e),
            "i" = "Returning NULL invisibly."
        ))
        invisible(NULL)
    })
}


# ============================================================================
# 2. .plot_coefficients --- Forest plot of fixed effects (internal)
# ============================================================================

#' Forest plot of fixed-effect estimates
#'
#' Constructs a ggplot forest plot by calling \code{summary.hbb_fit()}
#' and extracting the \code{fixed_effects} data frame.  Parameters are
#' grouped into three facets: Extensive (\eqn{\alpha}), Intensive
#' (\eqn{\beta}), and Dispersion (\eqn{\log\kappa}).  A vertical
#' reference line at zero facilitates significance assessment.
#'
#' @param object An \code{hbb_fit} object.
#' @param sandwich An \code{hbb_sandwich} or NULL.
#' @param level Confidence level in (0,1).
#' @param margin One of "both", "extensive", "intensive".
#' @return A ggplot object.
#' @keywords internal
.plot_coefficients <- function(object, sandwich, level, margin) {

    P <- object$hbb_data$P
    D <- 2L * P + 1L

    # Use summary() for consistent inference logic
    summ <- tryCatch(
        summary(object, sandwich = sandwich, level = level),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "Failed to compute summary for coefficient plot.",
                "i" = conditionMessage(e)
            ))
        }
    )

    fe <- summ$fixed_effects

    # Build margin group labels based on row position
    margin_group <- character(D)
    margin_group[seq_len(P)]        <- "Extensive"
    margin_group[(P + 1L):(2L * P)] <- "Intensive"
    margin_group[D]                 <- "Dispersion"

    fe$margin_group <- margin_group

    # Filter by requested margin
    if (margin == "extensive") {
        fe <- fe[fe$margin_group == "Extensive", , drop = FALSE]
    } else if (margin == "intensive") {
        fe <- fe[fe$margin_group == "Intensive", , drop = FALSE]
    }

    if (nrow(fe) == 0L) {
        cli::cli_abort("No parameters to plot for margin = {.val {margin}}.")
    }

    # Factor parameter with reversed levels for top-to-bottom reading
    fe$parameter <- factor(fe$parameter, levels = rev(fe$parameter))

    # Factor margin_group for facet ordering
    fe$margin_group <- factor(
        fe$margin_group,
        levels = c("Extensive", "Intensive", "Dispersion")
    )

    # Colour palette --- project standard
    margin_colors <- c(
        Extensive  = "#4393C3",
        Intensive  = "#D6604D",
        Dispersion = "#762A83"
    )

    # Inference label for subtitle
    pct <- round(100 * level, 0)
    ci_type <- if (isTRUE(summ$sandwich_used)) {
        paste0(pct, "% Wald CI (sandwich SE)")
    } else {
        paste0(pct, "% Credible Interval (posterior)")
    }

    out <- ggplot2::ggplot(
        fe,
        ggplot2::aes(
            x      = .data[["estimate"]],
            y      = .data[["parameter"]],
            xmin   = .data[["ci_lower"]],
            xmax   = .data[["ci_upper"]],
            colour = .data[["margin_group"]]
        )
    ) +
        ggplot2::geom_vline(
            xintercept = 0,
            linetype   = "dashed",
            colour     = "#1B7837",
            linewidth  = 0.5
        ) +
        ggplot2::geom_pointrange(
            size      = 0.5,
            linewidth = 0.7,
            fatten    = 2.5
        ) +
        ggplot2::scale_colour_manual(
            values = margin_colors,
            guide  = "none"
        ) +
        ggplot2::facet_wrap(
            ~ .data[["margin_group"]],
            scales = "free_y",
            ncol   = 1L
        ) +
        ggplot2::labs(
            title    = "Fixed-Effect Coefficient Estimates",
            subtitle = ci_type,
            x        = "Estimate",
            y        = NULL
        ) +
        ggplot2::theme_bw(base_size = 11L) +
        ggplot2::theme(
            strip.background = ggplot2::element_rect(
                fill = "grey92", colour = NA
            ),
            strip.text = ggplot2::element_text(
                face = "bold", size = 10L
            ),
            panel.grid.minor = ggplot2::element_blank()
        )

    out
}


# ============================================================================
# 3. .plot_trace --- MCMC trace plots (internal)
# ============================================================================

#' MCMC trace plots for fixed-effect parameters
#'
#' Extracts posterior draws from the CmdStanMCMC object and plots the
#' iteration-by-iteration trace for each selected parameter, coloured
#' by chain.  Useful for visual assessment of convergence and mixing.
#'
#' Uses \code{format = "array"} to extract draws with chain structure
#' preserved, avoiding the need for chain-ID reconstruction from flat
#' matrices.
#'
#' @section Chain-ID reconstruction:
#' CmdStanR \code{$draws(format = "array")} returns an
#' \eqn{\mathrm{iter} \times \mathrm{chains} \times \mathrm{params}}
#' array, so chain identity is directly available.
#'
#' @param object An \code{hbb_fit} object.
#' @param pars Character vector of CmdStanR parameter names
#'   (e.g., \code{"alpha[1]"}, \code{"beta[2]"}) or \code{NULL}
#'   for all fixed effects.
#' @return A ggplot object.
#' @keywords internal
.plot_trace <- function(object, pars) {

    P <- object$hbb_data$P
    D <- 2L * P + 1L

    param_names  <- .make_stan_param_names(P)
    param_labels <- .build_param_labels(object$hbb_data)

    # Build a name-to-label lookup
    label_map <- setNames(param_labels, param_names)

    # Subset if pars specified
    if (!is.null(pars)) {
        if (!is.character(pars) || length(pars) == 0L) {
            cli::cli_abort(
                "{.arg pars} must be a non-empty character vector of \\
                 CmdStanR parameter names."
            )
        }
        bad_pars <- setdiff(pars, param_names)
        if (length(bad_pars) > 0L) {
            cli::cli_abort(c(
                "x" = "Unknown parameter name{?s}: {.val {bad_pars}}.",
                "i" = "Available parameters: {.val {param_names}}."
            ))
        }
        sel_stan   <- pars
        sel_labels <- label_map[pars]
    } else {
        sel_stan   <- param_names
        sel_labels <- param_labels
    }

    # Extract draws as 3-D array [iter, chain, param]
    draws_array <- tryCatch(
        object$fit$draws(variables = sel_stan, format = "array"),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "Failed to extract posterior draws for trace plot.",
                "i" = conditionMessage(e)
            ))
        }
    )

    n_iter   <- dim(draws_array)[1L]
    n_chains <- dim(draws_array)[2L]
    n_pars   <- dim(draws_array)[3L]

    # Edge case: fewer than 2 iterations
    if (n_iter < 2L) {
        cli::cli_abort(c(
            "x" = "Only {n_iter} iteration{?s} available; need at least 2 for trace plot.",
            "i" = "Check CmdStan fit for sampling failures."
        ))
    }

    # Build long-format data frame (vectorised approach, adapted for array)
    df_list <- vector("list", n_chains * n_pars)
    idx <- 1L
    for (p in seq_len(n_pars)) {
        for (ch in seq_len(n_chains)) {
            df_list[[idx]] <- data.frame(
                iteration = seq_len(n_iter),
                chain     = factor(ch),
                parameter = sel_labels[p],
                value     = as.numeric(draws_array[, ch, p]),
                stringsAsFactors = FALSE
            )
            idx <- idx + 1L
        }
    }
    df_long <- do.call(rbind, df_list)

    # Preserve parameter ordering
    df_long$parameter <- factor(df_long$parameter,
                                 levels = unique(sel_labels))

    # Chain colours --- project palette extended to 8 chains
    chain_colours <- c("#4393C3", "#D6604D", "#1B7837", "#762A83",
                        "#E08214", "#542788", "#B2182B", "#2166AC")

    out <- ggplot2::ggplot(
        df_long,
        ggplot2::aes(
            x      = .data[["iteration"]],
            y      = .data[["value"]],
            colour = .data[["chain"]]
        )
    ) +
        ggplot2::geom_line(
            linewidth = 0.25,
            alpha     = 0.7
        ) +
        ggplot2::facet_wrap(
            ~ .data[["parameter"]],
            scales = "free_y"
        ) +
        ggplot2::scale_colour_manual(
            values = chain_colours[seq_len(n_chains)]
        ) +
        ggplot2::labs(
            title  = "MCMC Trace Plots",
            x      = "Iteration",
            y      = "Parameter value",
            colour = "Chain"
        ) +
        ggplot2::theme_bw(base_size = 11L) +
        ggplot2::theme(
            strip.background = ggplot2::element_rect(
                fill = "grey92", colour = NA
            ),
            strip.text       = ggplot2::element_text(face = "bold"),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position  = "bottom"
        )

    out
}


# ============================================================================
# 4. .plot_random_effects --- Caterpillar plot of state REs (internal)
# ============================================================================

#' Caterpillar plot of state-varying random effects
#'
#' For SVC models (\code{model_type \%in\% c("svc", "svc_weighted")}),
#' extracts the posterior summaries (mean and quantile-based intervals)
#' of the state random effects \eqn{\delta_{s,k}} and displays them
#' as a caterpillar plot faceted by covariate.
#'
#' Each panel shows all \eqn{S} states ordered by their posterior mean,
#' with horizontal bars spanning the \code{level} credible interval.
#' The vertical reference line at zero aids interpretation: states
#' with intervals that exclude zero exhibit statistically meaningful
#' departures from the fixed-effect mean.
#'
#' @param object An \code{hbb_fit} object (must be SVC model).
#' @param level Numeric in (0,1) for credible interval width.
#' @return A ggplot object.
#' @keywords internal
.plot_random_effects <- function(object, level) {

    # Validate SVC model type
    is_svc <- isTRUE(object$model_type %in% c("svc", "svc_weighted"))
    if (!is_svc) {
        cli::cli_abort(c(
            "x" = "Random effects plot requires an SVC model.",
            "i" = "Current model type: {.val {object$model_type}}.",
            "i" = "Use {.code plot(fit, type = 'coefficients')} instead."
        ))
    }

    hd <- object$hbb_data
    S  <- hd$S
    P  <- hd$P
    K  <- 2L * P   # total RE dimension

    # Build delta parameter names for CmdStanR summary
    delta_names <- character(S * K)
    idx <- 1L
    for (s in seq_len(S)) {
        for (k in seq_len(K)) {
            delta_names[idx] <- sprintf("delta[%d,%d]", s, k)
            idx <- idx + 1L
        }
    }

    # Extract posterior summaries
    delta_summ <- tryCatch(
        object$fit$summary(variables = delta_names),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "Failed to extract delta summaries for random effects plot.",
                "i" = conditionMessage(e)
            ))
        }
    )

    # Build readable covariate labels
    cov_names <- colnames(hd$X)
    if (is.null(cov_names)) cov_names <- paste0("x", seq_len(P))

    re_labels <- c(
        paste0("alpha_", cov_names),
        paste0("beta_",  cov_names)
    )

    # State labels
    state_labels <- if (!is.null(hd$group_levels)) {
        hd$group_levels
    } else {
        as.character(seq_len(S))
    }

    # Compute quantile column names based on level
    alpha_half <- (1 - level) / 2
    lo_pct <- round(alpha_half * 100)
    hi_pct <- round((1 - alpha_half) * 100)
    lo_col <- paste0("q", lo_pct)
    hi_col <- paste0("q", hi_pct)

    # Check if the needed quantile columns exist; fall back to q5/q95
    if (!lo_col %in% names(delta_summ)) {
        lo_col <- "q5"
        hi_col <- "q95"
        cli::cli_warn(c(
            "!" = "Quantile columns for level {level} not found in CmdStan summary.",
            "i" = "Falling back to q5/q95 (90% interval)."
        ))
    }

    # Build data frame: S*K rows
    # Row order from CmdStanR summary: (s=1,k=1), (s=1,k=2), ..., (s=1,k=K),
    # (s=2,k=1), ... => state varies slowly, covariate varies fast
    df_list <- vector("list", K)
    for (k in seq_len(K)) {
        row_indices <- seq(k, S * K, by = K)
        df_list[[k]] <- data.frame(
            state     = state_labels,
            covariate = re_labels[k],
            mean      = delta_summ$mean[row_indices],
            ci_lower  = delta_summ[[lo_col]][row_indices],
            ci_upper  = delta_summ[[hi_col]][row_indices],
            stringsAsFactors = FALSE
        )
    }
    df_re <- do.call(rbind, df_list)

    # Factor covariate to preserve order
    df_re$covariate <- factor(df_re$covariate, levels = re_labels)

    # Assign margin colouring for point-range
    df_re$margin <- ifelse(
        grepl("^alpha_", df_re$covariate),
        "Extensive", "Intensive"
    )

    margin_colors <- c(
        Extensive = "#4393C3",
        Intensive = "#D6604D"
    )

    pct <- round(level * 100)

    out <- ggplot2::ggplot(
        df_re,
        ggplot2::aes(
            x      = .data[["mean"]],
            y      = stats::reorder(.data[["state"]], .data[["mean"]]),
            xmin   = .data[["ci_lower"]],
            xmax   = .data[["ci_upper"]],
            colour = .data[["margin"]]
        )
    ) +
        ggplot2::geom_vline(
            xintercept = 0,
            linetype   = "dashed",
            colour     = "#1B7837",
            linewidth  = 0.5
        ) +
        ggplot2::geom_pointrange(
            size      = 0.25,
            linewidth = 0.4,
            fatten    = 1.5
        ) +
        ggplot2::scale_colour_manual(
            values = margin_colors,
            name   = "Margin"
        ) +
        ggplot2::facet_wrap(
            ~ .data[["covariate"]],
            scales = "free",
            ncol   = 2L
        ) +
        ggplot2::labs(
            title    = sprintf("State Random Effects (%d%% Credible Intervals)",
                               pct),
            subtitle = "Posterior mean, sorted by mean within each facet",
            x        = "Random effect",
            y        = "State"
        ) +
        ggplot2::theme_bw(base_size = 10L) +
        ggplot2::theme(
            strip.background = ggplot2::element_rect(
                fill = "grey92", colour = NA
            ),
            strip.text       = ggplot2::element_text(face = "bold", size = 8L),
            axis.text.y      = ggplot2::element_text(size = 6L),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position  = "bottom"
        )

    out
}


# ============================================================================
# 5. .plot_residuals --- Residual diagnostic panels (internal)
# ============================================================================

#' Two-panel residual diagnostic plot
#'
#' Constructs two diagnostic panels using existing S3 methods:
#' \enumerate{
#'   \item \strong{Residuals vs fitted}: Pearson residuals (vertical)
#'     against fitted values \eqn{\hat{q}_i \hat\mu_i} (horizontal).
#'     Under correct model specification, the residuals should scatter
#'     randomly around zero with no systematic pattern.  A LOESS
#'     smoother (red) is overlaid to highlight trends.
#'   \item \strong{Normal QQ plot}: Sorted Pearson residuals against
#'     theoretical standard normal quantiles
#'     \eqn{\Phi^{-1}((i - 0.375)/(N + 0.25))}.  Departures from the
#'     diagonal indicate non-normality of the residual distribution.
#' }
#'
#' The panels are combined via \pkg{patchwork}.
#'
#' @section Pearson residuals:
#' The Pearson residual for observation \eqn{i} is
#' \deqn{
#'   r_i^P = \frac{y_i - n_i \hat{q}_i \hat\mu_i}
#'                {\sqrt{\widehat{\mathrm{Var}}(Y_i)}},
#' }
#' where \eqn{\widehat{\mathrm{Var}}(Y_i)} is the hurdle variance
#' computed from \eqn{(n_i, \hat{q}_i, \hat\mu_i, \hat\kappa)}.
#' See \code{\link{hurdle_variance}} for the variance decomposition.
#'
#' @param object An \code{hbb_fit} object.
#' @return A patchwork object combining two ggplot panels.
#' @keywords internal
.plot_residuals <- function(object) {

    rlang::check_installed(
        c("ggplot2", "patchwork"),
        reason = "to create residual diagnostic plots"
    )

    # Reuse existing S3 methods for consistency
    fv <- fitted(object, type = "response")

    # Pearson residuals with graceful fallback
    resid_pearson <- tryCatch(
        residuals(object, type = "pearson"),
        error = function(e) {
            cli::cli_warn(c(
                "!" = "Failed to compute Pearson residuals; QQ panel will be empty.",
                "i" = conditionMessage(e)
            ))
            NULL
        }
    )

    # ---- Panel 1: Residuals vs Fitted ----------------------------------------
    df_resid <- data.frame(
        fitted   = fv,
        residual = if (!is.null(resid_pearson)) resid_pearson else rep(NA_real_, length(fv)),
        stringsAsFactors = FALSE
    )

    # Guard NaN/Inf
    df_resid$residual[!is.finite(df_resid$residual)] <- NA_real_

    n_valid <- sum(is.finite(df_resid$residual) & is.finite(df_resid$fitted))

    if (n_valid < 2L) {
        cli::cli_abort(c(
            "x" = "Fewer than 2 valid residuals available for plotting.",
            "i" = "Check model fit and diagnostics."
        ))
    }

    p1 <- ggplot2::ggplot(
        df_resid,
        ggplot2::aes(
            x = .data[["fitted"]],
            y = .data[["residual"]]
        )
    ) +
        ggplot2::geom_hline(
            yintercept = 0,
            linetype   = "dashed",
            colour     = "#1B7837",
            linewidth  = 0.5
        ) +
        ggplot2::geom_point(
            colour = "#4393C3",
            alpha  = 0.3,
            size   = 0.8
        )

    # Only add smoother if enough non-NA points
    if (n_valid >= 5L) {
        p1 <- p1 + ggplot2::geom_smooth(
            method    = "loess",
            formula   = y ~ x,
            se        = FALSE,
            colour    = "#D6604D",
            linewidth = 0.8
        )
    }

    p1 <- p1 +
        ggplot2::labs(
            title = "Residuals vs Fitted",
            x     = expression(hat(q)[i] %.% hat(mu)[i]),
            y     = "Pearson residual"
        ) +
        ggplot2::theme_bw(base_size = 10L) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank()
        )

    # ---- Panel 2: Normal QQ Plot ---------------------------------------------
    if (!is.null(resid_pearson)) {
        rp_clean <- resid_pearson[is.finite(resid_pearson)]

        if (length(rp_clean) >= 2L) {
            n_clean <- length(rp_clean)
            theoretical <- stats::qnorm(stats::ppoints(n_clean))
            sample_q    <- sort(rp_clean)

            df_qq <- data.frame(
                theoretical = theoretical,
                sample      = sample_q,
                stringsAsFactors = FALSE
            )

            p2 <- ggplot2::ggplot(
                df_qq,
                ggplot2::aes(
                    x = .data[["theoretical"]],
                    y = .data[["sample"]]
                )
            ) +
                ggplot2::geom_abline(
                    slope     = 1,
                    intercept = 0,
                    linetype  = "dashed",
                    colour    = "#1B7837",
                    linewidth = 0.5
                ) +
                ggplot2::geom_point(
                    colour = "#762A83",
                    alpha  = 0.3,
                    size   = 0.8
                ) +
                ggplot2::labs(
                    title = "Normal QQ Plot",
                    x     = "Theoretical quantiles",
                    y     = "Sample quantiles (Pearson)"
                ) +
                ggplot2::theme_bw(base_size = 10L) +
                ggplot2::theme(
                    panel.grid.minor = ggplot2::element_blank()
                )
        } else {
            # Fewer than 2 valid Pearson residuals (fallback)
            p2 <- ggplot2::ggplot() +
                ggplot2::labs(
                    title    = "Normal QQ Plot",
                    subtitle = "(insufficient valid residuals)"
                ) +
                ggplot2::theme_void(base_size = 10L)
        }
    } else {
        # Pearson residuals failed entirely (fallback)
        p2 <- ggplot2::ggplot() +
            ggplot2::labs(
                title    = "Normal QQ Plot",
                subtitle = "(Pearson residual computation failed)"
            ) +
            ggplot2::theme_void(base_size = 10L)
    }

    # ---- Combine via patchwork -----------------------------------------------
    out <- p1 + p2 +
        patchwork::plot_annotation(
            title   = "Hurdle Beta-Binomial Residual Diagnostics",
            caption = sprintf("N = %d valid observations", n_valid),
            theme   = ggplot2::theme(
                plot.title = ggplot2::element_text(
                    face = "bold", size = 13L
                )
            )
        )

    out
}
