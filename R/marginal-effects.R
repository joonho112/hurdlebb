# ============================================================================
# marginal-effects.R --- Average Marginal Effects Decomposition for HBB Models
#
# Implements the AME decomposition of the Hurdle Beta-Binomial model into
# extensive-margin (participation) and intensive-margin (enrollment share)
# contributions, with posterior uncertainty quantification and optional
# delta-method Wald inference via the sandwich variance.
#
# Theory:
#   The exact unconditional mean of the HBB model is:
#     E[y_i / n_i] = q_i * g(mu_i, n_i, kappa)
#   where g(mu) = mu / (1 - p_0) >= mu is the intensity function (Eq. 2.5).
#   The AME computation uses this exact formula.
#   For the NSECE data (typical n ~ 50, kappa ~ 7), p_0 < 0.01, so
#   g(mu) ~ mu and the correction is typically < 1%.
#
#   By the product rule on q*mu, the marginal effect of covariate k decomposes:
#     d(q*mu) / dx_k = alpha_k * q_i(1-q_i) * mu_i
#                     + beta_k  * mu_i(1-mu_i) * q_i
#                     = [extensive]  +  [intensive]
#
#   Note: The manuscript's formal decomposition (Proposition 2.5, LAE/LIE)
#   operates on the log scale and exactly incorporates h(mu) via the
#   elasticity factor epsilon_h ~ 0.93.  This level-scale product-rule
#   decomposition is a complementary summary that is standard in the
#   hurdle/two-part model literature.
#
#   The Average Marginal Effect (AME) averages over all N observations:
#     AME_k = (1/N) sum_i [ext_ik + int_ik]
#
#   Key vectorisation insight: since the N-vectors q_deriv*mu and
#   mu_deriv*q are independent of covariate index k, we factor out the scalar
#   means:
#     AME_k^ext = alpha_k * mean(q_deriv * mu)
#     AME_k^int = beta_k  * mean(mu_deriv * q)
#   This reduces the per-draw cost from O(N*P) to O(N).
#
#
# Contents:
#   1. ame                   --- Main AME computation (exported)
#   2. print.hbb_ame         --- S3 print method (exported)
#   3. ame_decomposition     --- Extract decomposition table (exported)
#   4. .ame_at_theta         --- AME at a single theta vector (internal)
#   5. .ame_summarize        --- Posterior summary for M x P matrix (internal)
#   6. .ame_build_decomp     --- Build decomposition data.frame (internal)
#   7. .ame_compute_wald     --- Delta-method Wald CIs (internal)
#   8. .ame_reversal_probs   --- Reversal probability computation (internal)
#   9. .ame_validate_inputs  --- Comprehensive input validation (internal)
#  10. .ame_extract_draws    --- Draw extraction logic (internal)
# ============================================================================


# ============================================================================
# 1. ame --- Main AME computation (exported)
# ============================================================================

#' Average Marginal Effects Decomposition for Hurdle Beta-Binomial Models
#'
#' Computes the Average Marginal Effect (AME) of each covariate on the
#' expected IT enrollment share \eqn{E[y_i / n_i] = q_i \cdot g(\mu_i)}
#' where \eqn{g(\mu) = \mu / (1 - p_0)} is the zero-truncation intensity
#' function, decomposed into extensive-margin
#' and intensive-margin contributions.
#'
#' @description
#' The hurdle Beta-Binomial model implies two channels through which
#' covariates affect the expected outcome.  The extensive margin governs
#' whether a center serves infant/toddler children at all
#' (\eqn{q_i = \mathrm{logistic}(X_i' \alpha)}), while the intensive
#' margin governs the enrollment share conditional on participation
#' (\eqn{\mu_i = \mathrm{logistic}(X_i' \beta)}).
#'
#' By the product rule, the marginal effect of covariate \eqn{k} on
#' \eqn{E[y_i / n_i]} decomposes as:
#' \deqn{
#'   \frac{\partial E[y_i/n_i]}{\partial x_{ik}}
#'     = \underbrace{\alpha_k \, q_i(1-q_i) \, \mu_i}_{\text{extensive}}
#'     + \underbrace{\beta_k \, \mu_i(1-\mu_i) \, q_i}_{\text{intensive}}.
#' }
#'
#' The Average Marginal Effect averages this over all \eqn{N} observations:
#' \deqn{
#'   \mathrm{AME}_k = \frac{1}{N} \sum_{i=1}^{N}
#'     \bigl[\alpha_k \, q_i(1-q_i) \, \mu_i
#'           + \beta_k \, \mu_i(1-\mu_i) \, q_i\bigr]
#'     = \mathrm{AME}_k^{\mathrm{ext}} + \mathrm{AME}_k^{\mathrm{int}}.
#' }
#'
#' @section Poverty Reversal:
#' A covariate exhibits a "reversal" when its extensive and intensive
#' AME components have opposing signs.
#' For the poverty variable in the NSECE application:
#' \eqn{\alpha_{\mathrm{poverty}} < 0} (higher poverty reduces participation)
#' but \eqn{\beta_{\mathrm{poverty}} > 0} (higher poverty increases IT share
#' among servers).  The posterior probability of reversal is
#' \eqn{\Pr(\mathrm{AME}_k^{\mathrm{ext}} < 0 \;\text{AND}\;
#' \mathrm{AME}_k^{\mathrm{int}} > 0)}.
#'
#' @section Population-Average AME:
#' This function computes the population-average AME (PA-AME) using
#' global (population-level) coefficients \eqn{\alpha} and \eqn{\beta}
#' only, without state random effects.  This is consistent with the
#' sandwich-corrected inference framework and answers the question:
#' "what is the average marginal effect for a new (arbitrary) state?"
#'
#' @section Wald Inference:
#' When \code{sandwich} is provided, the function also computes
#' delta-method Wald confidence intervals for each AME component.
#' The numerical gradient of \eqn{\mathrm{AME}_k(\theta)} with respect
#' to \eqn{\theta} is computed via central differences (step size
#' \eqn{\epsilon = 10^{-5}}), and the Wald variance is:
#' \deqn{
#'   \mathrm{Var}(\mathrm{AME}_k)
#'     \approx \nabla_\theta \mathrm{AME}_k(\hat\theta)^\top
#'             V_{\mathrm{sand}}\,
#'             \nabla_\theta \mathrm{AME}_k(\hat\theta).
#' }
#'
#' @section Performance:
#' The computation is vectorised over observations within each MCMC draw.
#' The key insight is that the N-vectors \eqn{q'}, \eqn{\mu}, \eqn{\mu'},
#' \eqn{q} are independent of the covariate index \eqn{k}, so:
#' \deqn{
#'   \mathrm{AME}_k^{\mathrm{ext}} = \alpha_k \cdot \overline{q'(1-q) \mu},
#'   \quad
#'   \mathrm{AME}_k^{\mathrm{int}} = \beta_k \cdot \overline{\mu'(1-\mu) q},
#' }
#' reducing the per-draw cost from \eqn{O(NP)} to \eqn{O(N)}.
#'
#' @param fit An object of class \code{"hbb_fit"} returned by
#'   \code{\link{hbb}}.  Must contain a CmdStanR fit with parameters
#'   \code{alpha[1:P]}, \code{beta[1:P]}, \code{log_kappa}, and an
#'   \code{hbb_data} element with design matrix \code{X}.
#' @param cholesky An optional object of class \code{"hbb_cholesky"}
#'   returned by \code{\link{cholesky_correct}}.  If provided, uses
#'   the Cholesky-corrected draws (\code{cholesky$theta_corrected}) and
#'   posterior mean (\code{cholesky$theta_hat}) for AME computation.
#'   If \code{NULL} (default), uses raw MCMC draws from the fit object.
#' @param sandwich An optional object of class \code{"hbb_sandwich"}
#'   returned by \code{\link{sandwich_variance}}.  If provided,
#'   computes delta-method Wald confidence intervals for each AME
#'   component.  If \code{NULL} (default), Wald results are omitted.
#' @param level Numeric scalar in \eqn{(0, 1)}.  Confidence level for
#'   posterior credible intervals and Wald confidence intervals.
#'   Default is \code{0.95}.
#' @param n_draws Integer or \code{NULL}.  Number of MCMC draws to use.
#'   If \code{NULL} (default), uses all available draws.  If an integer,
#'   subsamples by systematic thinning to approximately \code{n_draws}
#'   draws.
#'
#' @return An S3 object of class \code{"hbb_ame"} containing:
#' \describe{
#'   \item{\code{theta_draws_used}}{Numeric matrix of dimension
#'     \code{M_use x D} containing the (possibly subsampled) draws.}
#'   \item{\code{theta_hat}}{Numeric vector of length D: point estimate
#'     (posterior mean or Cholesky theta_hat).}
#'   \item{\code{ext_ame_draws}}{Numeric matrix \code{M_use x P}:
#'     extensive AME draws for each covariate.}
#'   \item{\code{int_ame_draws}}{Numeric matrix \code{M_use x P}:
#'     intensive AME draws for each covariate.}
#'   \item{\code{total_ame_draws}}{Numeric matrix \code{M_use x P}:
#'     total AME draws for each covariate.}
#'   \item{\code{ext_summary}}{Data frame with posterior summaries of
#'     extensive AME (7 columns: covariate, post_mean, post_median,
#'     ci_lo, ci_hi, post_sd, pr_positive).}
#'   \item{\code{int_summary}}{Data frame with posterior summaries of
#'     intensive AME.}
#'   \item{\code{total_summary}}{Data frame with posterior summaries of
#'     total AME.}
#'   \item{\code{decomp_table}}{Data frame for non-intercept covariates
#'     with columns: covariate, ext_ame, ext_ci_lo, ext_ci_hi,
#'     int_ame, int_ci_lo, int_ci_hi, total_ame, total_ci_lo,
#'     total_ci_hi, ext_share, int_share, sign_pattern.}
#'   \item{\code{reversal_probs}}{Named numeric vector: for each
#'     non-intercept covariate, the posterior probability of opposing
#'     signs between extensive and intensive components.}
#'   \item{\code{wald_summary}}{Data frame of Wald-based AME inference,
#'     or \code{NULL} if \code{sandwich} was not provided.}
#'   \item{\code{ame_ext_hat}}{Numeric vector of length P: extensive
#'     AME point estimates at \code{theta_hat}.}
#'   \item{\code{ame_int_hat}}{Numeric vector of length P: intensive
#'     AME point estimates at \code{theta_hat}.}
#'   \item{\code{ame_total_hat}}{Numeric vector of length P: total
#'     AME point estimates at \code{theta_hat}.}
#'   \item{\code{mean_q}}{Numeric vector of length \code{M_use}:
#'     mean participation probability across observations per draw.}
#'   \item{\code{mean_mu}}{Numeric vector of length \code{M_use}:
#'     mean conditional intensity across observations per draw.}
#'   \item{\code{N, P, D, M_total, M_use, level, cov_labels}}{
#'     Dimensional and metadata scalars.}
#' }
#'
#' @examples
#' \dontrun{
#' # After fitting:
#' fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data,
#'            weights = "weight")
#' sand <- sandwich_variance(fit)
#' chol <- cholesky_correct(fit, sand)
#'
#' # Full AME decomposition with Wald comparison
#' ame_result <- ame(fit, cholesky = chol, sandwich = sand)
#' print(ame_result)
#'
#' # Extract decomposition table
#' ame_decomposition(ame_result)
#' }
#'
#' @seealso
#' \code{\link{ame_decomposition}} for extracting the decomposition table,
#' \code{\link{cholesky_correct}} for Cholesky-corrected draws,
#' \code{\link{sandwich_variance}} for the sandwich variance.
#'
#' @references
#' Ghosal, R., Ghosh, S. K., and Maiti, T. (2020).
#' Two-part regression models for longitudinal zero-inflated count data.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' \strong{183}(4), 1603--1626.
#'
#' @export
ame <- function(fit,
                cholesky = NULL,
                sandwich = NULL,
                level    = 0.95,
                n_draws  = NULL) {

    # =========================================================================
    # 0. INPUT VALIDATION
    # =========================================================================

    .ame_validate_inputs(fit, cholesky, sandwich, level, n_draws)

    # =========================================================================
    # 1. EXTRACT DATA AND DRAWS
    # =========================================================================

    X <- fit$hbb_data$X
    P <- fit$hbb_data$P
    N <- fit$hbb_data$N
    D <- 2L * P + 1L
    n_trial <- as.integer(fit$hbb_data$n_trial)

    # -- Covariate labels ------------------------------------------------------
    if (!is.null(colnames(X))) {
        cov_labels <- colnames(X)
    } else {
        cov_labels <- c("intercept", paste0("x", seq_len(P - 1L)))
    }

    cli::cli_alert_info(
        "AME decomposition: N = {N}, P = {P}, D = {D}"
    )

    # -- Extract draws ---------------------------------------------------------
    draws_info <- .ame_extract_draws(fit, cholesky, P)
    theta_all <- draws_info$draws
    theta_hat <- draws_info$theta_hat
    M_total   <- draws_info$M_total

    # -- Subsample by thinning if requested ------------------------------------
    if (!is.null(n_draws) && n_draws < M_total) {
        thin_step <- max(1L, floor(M_total / n_draws))
        thin_idx <- seq(1L, M_total, by = thin_step)
        if (length(thin_idx) > n_draws) thin_idx <- thin_idx[seq_len(n_draws)]
        theta_sub <- theta_all[thin_idx, , drop = FALSE]
        cli::cli_alert_info(
            "Subsampled to {nrow(theta_sub)} draws (from {M_total}, every {thin_step}-th)"
        )
    } else {
        theta_sub <- theta_all
    }

    M_use <- nrow(theta_sub)

    cli::cli_alert_info("Using {M_use} draws (from {M_total} total)")

    # =========================================================================
    # 2. COMPUTE AME FOR EACH DRAW
    # =========================================================================

    cli::cli_alert_info("Computing AME for {M_use} draws x {N} observations...")

    # Pre-allocate storage
    ext_ame   <- matrix(NA_real_, nrow = M_use, ncol = P)
    int_ame   <- matrix(NA_real_, nrow = M_use, ncol = P)
    total_ame <- matrix(NA_real_, nrow = M_use, ncol = P)
    mean_q    <- numeric(M_use)
    mean_mu   <- numeric(M_use)

    colnames(ext_ame)   <- cov_labels
    colnames(int_ame)   <- cov_labels
    colnames(total_ame) <- cov_labels

    # NaN/Inf tracking
    n_nan_draws <- 0L

    for (m in seq_len(M_use)) {
        alpha_m <- theta_sub[m, seq_len(P)]
        beta_m  <- theta_sub[m, (P + 1L):(2L * P)]
        log_kappa_m <- theta_sub[m, D]
        kappa_m <- pmin(exp(log_kappa_m), 1e15)

        # Linear predictors (N-vectors)
        eta_ext <- as.numeric(X %*% alpha_m)
        eta_int <- as.numeric(X %*% beta_m)

        # Inverse logit (plogis handles overflow/underflow safely)
        q_m  <- plogis(eta_ext)
        mu_m <- plogis(eta_int)

        # Logistic derivatives
        q_deriv  <- q_m * (1 - q_m)
        mu_deriv <- mu_m * (1 - mu_m)

        # --- ZT correction: g(mu) = mu / (1 - p0) ---
        p0_m    <- compute_p0(n_trial, mu_m, kappa_m)
        one_mp0 <- pmax(1 - p0_m, .Machine$double.eps)
        g_m     <- mu_m / one_mp0

        # g'(mu) for intensive kernel
        b_m       <- (1 - mu_m) * kappa_m
        Lambda_m  <- kappa_m * (digamma(b_m + n_trial) - digamma(b_m))
        g_prime_m <- 1 / one_mp0 + mu_m * p0_m * Lambda_m / one_mp0^2

        # Scalar kernels (ZT-corrected)
        mean_ext_kernel <- mean(q_deriv * g_m)
        mean_int_kernel <- mean(g_prime_m * mu_deriv * q_m)

        # NaN/Inf check
        if (!is.finite(mean_ext_kernel) || !is.finite(mean_int_kernel)) {
            n_nan_draws <- n_nan_draws + 1L
            next
        }

        # Vectorised across covariates: AME_k = coeff_k * kernel
        ext_ame[m, ]   <- alpha_m * mean_ext_kernel
        int_ame[m, ]   <- beta_m  * mean_int_kernel
        total_ame[m, ] <- ext_ame[m, ] + int_ame[m, ]

        # Diagnostics
        mean_q[m]  <- mean(q_m)
        mean_mu[m] <- mean(g_m)
    }

    if (n_nan_draws > 0L) {
        cli::cli_alert_warning(
            "{n_nan_draws} draw(s) produced NaN/Inf kernels and were skipped."
        )
    }

    cli::cli_alert_info(
        "Mean P(serve IT) = {round(mean(mean_q, na.rm = TRUE), 4)}, Mean E[IT share | serve] = {round(mean(mean_mu, na.rm = TRUE), 4)}"
    )

    # =========================================================================
    # 3. POSTERIOR SUMMARIES
    # =========================================================================

    ext_summary   <- .ame_summarize(ext_ame,   cov_labels, level)
    int_summary   <- .ame_summarize(int_ame,   cov_labels, level)
    total_summary <- .ame_summarize(total_ame, cov_labels, level)

    # =========================================================================
    # 4. DECOMPOSITION TABLE
    # =========================================================================

    decomp_idx <- if (P >= 2L) seq(2L, P) else integer(0)

    decomp_table <- .ame_build_decomp(
        ext_summary, int_summary, total_summary, cov_labels, decomp_idx
    )

    # =========================================================================
    # 5. REVERSAL PROBABILITIES
    # =========================================================================

    reversal_probs <- .ame_reversal_probs(
        ext_ame, int_ame, cov_labels, decomp_idx
    )

    # =========================================================================
    # 6. POINT ESTIMATES AT THETA_HAT
    # =========================================================================

    point_est <- .ame_at_theta(theta_hat, X, n_trial, P)
    ame_ext_hat   <- point_est$ext
    ame_int_hat   <- point_est$int
    ame_total_hat <- point_est$total
    names(ame_ext_hat)   <- cov_labels
    names(ame_int_hat)   <- cov_labels
    names(ame_total_hat) <- cov_labels

    # =========================================================================
    # 7. WALD INFERENCE
    # =========================================================================

    wald_summary <- NULL
    if (!is.null(sandwich)) {
        # Validate sandwich dimensions
        if (is.null(sandwich$V_sand) ||
            nrow(sandwich$V_sand) != D || ncol(sandwich$V_sand) != D) {
            cli::cli_alert_warning(
                "Sandwich V_sand dimensions do not match D = {D}. Skipping Wald AME."
            )
        } else {
            cli::cli_alert_info("Computing delta-method Wald AME CIs...")
            wald_summary <- .ame_compute_wald(
                theta_hat = theta_hat,
                X         = X,
                n_trial   = n_trial,
                V_sand    = sandwich$V_sand,
                P         = P,
                D         = D,
                level     = level
            )
        }
    }

    # =========================================================================
    # 8. CONSTRUCT RETURN OBJECT
    # =========================================================================

    cli::cli_alert_success("AME decomposition complete")

    structure(
        list(
            theta_draws_used = theta_sub,
            theta_hat        = theta_hat,
            ext_ame_draws    = ext_ame,
            int_ame_draws    = int_ame,
            total_ame_draws  = total_ame,
            ext_summary      = ext_summary,
            int_summary      = int_summary,
            total_summary    = total_summary,
            decomp_table     = decomp_table,
            reversal_probs   = reversal_probs,
            wald_summary     = wald_summary,
            ame_ext_hat      = ame_ext_hat,
            ame_int_hat      = ame_int_hat,
            ame_total_hat    = ame_total_hat,
            mean_q           = mean_q,
            mean_mu          = mean_mu,
            N                = N,
            P                = P,
            D                = D,
            M_total          = M_total,
            M_use            = M_use,
            level            = level,
            cov_labels       = cov_labels
        ),
        class = "hbb_ame"
    )
}


# ============================================================================
# 2. print.hbb_ame --- S3 print method (exported)
# ============================================================================

#' Print Method for hbb_ame Objects
#'
#' Displays a structured summary of the AME decomposition, including
#' the decomposition table with sign patterns, reversal probabilities,
#' and optional Wald comparison.
#'
#' @param x An object of class \code{"hbb_ame"} returned by
#'   \code{\link{ame}}.
#' @param digits Integer; number of significant digits for numeric
#'   output.  Default is 4.
#' @param ... Additional arguments (currently unused).
#'
#' @return The object \code{x}, invisibly.
#'
#' @seealso \code{\link{ame}}, \code{\link{ame_decomposition}}
#'
#' @method print hbb_ame
#' @export
print.hbb_ame <- function(x, digits = 4, ...) {

    # Full tryCatch for defensive operation on corrupted objects
    tryCatch({

        cat("\n")
        cat(strrep("=", 65), "\n")
        cat("  Average Marginal Effects Decomposition (HBB Model)\n")
        cat(strrep("=", 65), "\n\n")

        # -- Dimensions --------------------------------------------------------
        cat(sprintf("  Observations (N):   %d\n", x$N))
        cat(sprintf("  Covariates   (P):   %d (incl. intercept)\n", x$P))
        cat(sprintf("  Parameters   (D):   %d\n", x$D))
        cat(sprintf("  Draws used (M_use): %d (of %d total)\n",
                     x$M_use, x$M_total))
        cat(sprintf("  Confidence level:   %.2f\n", x$level))
        cat("\n")

        # -- Population-level predictions --------------------------------------
        cat(sprintf("  Mean P(serve IT):          %.4f\n",
                    mean(x$mean_q, na.rm = TRUE)))
        cat(sprintf("  Mean E[IT share | serve]:  %.4f\n",
                    mean(x$mean_mu, na.rm = TRUE)))
        cat(sprintf("  Mean E[IT share]:          %.4f\n\n",
                     mean(x$mean_q, na.rm = TRUE) *
                     mean(x$mean_mu, na.rm = TRUE)))

        # -- Decomposition table -----------------------------------------------
        dt <- x$decomp_table
        if (!is.null(dt) && nrow(dt) > 0L) {
            cat(strrep("-", 65), "\n")
            cat("  AME Decomposition (non-intercept covariates)\n")
            cat(strrep("-", 65), "\n")
            cat(sprintf("  %-12s  %10s  %10s  %10s  %6s  %6s  %s\n",
                         "Covariate", "Ext_AME", "Int_AME", "Total_AME",
                         "Ext%", "Int%", "Pattern"))
            cat(strrep("-", 65), "\n")
            for (i in seq_len(nrow(dt))) {
                r <- dt[i, ]
                cat(sprintf("  %-12s  %+10.*f  %+10.*f  %+10.*f  %5.1f%%  %5.1f%%  %s\n",
                             r$covariate,
                             digits, r$ext_ame,
                             digits, r$int_ame,
                             digits, r$total_ame,
                             r$ext_share, r$int_share,
                             r$sign_pattern))
            }
            cat(strrep("-", 65), "\n\n")
        }

        # -- Reversal probabilities --------------------------------------------
        rp <- x$reversal_probs
        if (!is.null(rp) && length(rp) > 0L) {
            cat("  Reversal probabilities (Pr of opposing signs):\n")
            for (nm in names(rp)) {
                cat(sprintf("    %-12s  %.4f\n", nm, rp[[nm]]))
            }
            cat("\n")
        }

        # -- Wald comparison ---------------------------------------------------
        if (!is.null(x$wald_summary)) {
            ws <- x$wald_summary
            cat(strrep("-", 65), "\n")
            cat("  Wald AME (delta method, sandwich SE)\n")
            cat(strrep("-", 65), "\n")
            cat(sprintf("  %-12s  %10s  %10s  %10s  %10s\n",
                         "Covariate", "Total_AME", "SE", "CI_lo", "CI_hi"))
            cat(strrep("-", 65), "\n")
            for (i in seq_len(nrow(ws))) {
                r <- ws[i, ]
                cat(sprintf("  %-12s  %+10.*f  %10.*f  %+10.*f  %+10.*f\n",
                             r$covariate,
                             digits, r$total_point,
                             digits, r$total_se,
                             digits, r$total_lo,
                             digits, r$total_hi))
            }
            cat(strrep("-", 65), "\n\n")

            # Sign agreement check
            decomp_idx <- seq(2L, x$P)
            n_agree <- 0L
            n_check <- 0L
            for (i in decomp_idx) {
                total_post <- x$total_summary$post_mean[i]
                total_wald <- ws$total_point[i]
                if (sign(total_post) == sign(total_wald)) {
                    n_agree <- n_agree + 1L
                }
                n_check <- n_check + 1L
            }
            cat(sprintf("  Posterior-Wald sign agreement: %d/%d covariates\n\n",
                         n_agree, n_check))
        }

    }, error = function(e) {
        cat("Error printing hbb_ame object:", conditionMessage(e), "\n")
        cat("Object class:", paste(class(x), collapse = ", "), "\n")
        cat("Object names:", paste(names(x), collapse = ", "), "\n")
    })

    invisible(x)
}


# ============================================================================
# 3. ame_decomposition --- Extract decomposition table (exported)
# ============================================================================

#' Extract the AME Decomposition Table
#'
#' Extracts the decomposition table from an \code{hbb_ame} object,
#' showing the extensive and intensive components of the Average
#' Marginal Effect for each non-intercept covariate.
#'
#' @param ame_result An object of class \code{"hbb_ame"} returned by
#'   \code{\link{ame}}.
#'
#' @return A data frame with one row per non-intercept covariate and
#'   columns: \code{covariate}, \code{ext_ame}, \code{ext_ci_lo},
#'   \code{ext_ci_hi}, \code{int_ame}, \code{int_ci_lo},
#'   \code{int_ci_hi}, \code{total_ame}, \code{total_ci_lo},
#'   \code{total_ci_hi}, \code{ext_share}, \code{int_share},
#'   \code{sign_pattern}.
#'
#' @seealso \code{\link{ame}}
#'
#' @examples
#' \dontrun{
#' ame_result <- ame(fit, cholesky = chol)
#' decomp <- ame_decomposition(ame_result)
#' print(decomp)
#' }
#'
#' @export
ame_decomposition <- function(ame_result) {
    if (!inherits(ame_result, "hbb_ame")) {
        cli::cli_abort(c(
            "x" = "{.arg ame_result} must be an object of class {.cls hbb_ame}.",
            "i" = "Received class: {.cls {class(ame_result)}}.",
            "i" = "Use {.fn ame} to compute marginal effects first."
        ))
    }
    ame_result$decomp_table
}


# ============================================================================
# Internal helpers
# ============================================================================


# ---- 4. .ame_at_theta --- AME at a single theta vector ---------------------

#' Compute AME at a given theta vector
#'
#' @param theta Numeric vector of length D = 2P + 1.
#' @param X Numeric matrix N x P.
#' @param n_trial Integer vector of length N: trial sizes for ZT correction.
#' @param P Integer: number of covariates per margin.
#' @return Named list with elements: ext (P-vector), int (P-vector),
#'   total (P-vector), mean_q (scalar), mean_mu (scalar).
#' @keywords internal
.ame_at_theta <- function(theta, X, n_trial, P) {

    D <- 2L * P + 1L

    alpha <- theta[seq_len(P)]
    beta  <- theta[(P + 1L):(2L * P)]
    log_kappa <- theta[D]
    kappa <- pmin(exp(log_kappa), 1e15)   # numeric guard

    eta_ext <- as.numeric(X %*% alpha)
    eta_int <- as.numeric(X %*% beta)

    q_v  <- plogis(eta_ext)
    mu_v <- plogis(eta_int)

    q_deriv  <- q_v * (1 - q_v)
    mu_deriv <- mu_v * (1 - mu_v)

    # --- ZT correction: g(mu) = mu / (1 - p0) --------------------------------
    p0_v     <- compute_p0(n_trial, mu_v, kappa)
    one_mp0  <- pmax(1 - p0_v, .Machine$double.eps)
    g_v      <- mu_v / one_mp0             # intensity function h(mu)

    # g'(mu) = 1/(1-p0) + mu * p0 * Lambda / (1-p0)^2
    # where Lambda = kappa * [psi(b+n) - psi(b)], b = (1-mu)*kappa
    b_v       <- (1 - mu_v) * kappa
    Lambda_v  <- kappa * (digamma(b_v + n_trial) - digamma(b_v))
    g_prime_v <- 1 / one_mp0 + mu_v * p0_v * Lambda_v / one_mp0^2

    # Scalar kernels (ZT-corrected)
    ext_kernel <- mean(q_deriv * g_v)
    int_kernel <- mean(g_prime_v * mu_deriv * q_v)

    ext_v   <- alpha * ext_kernel
    int_v   <- beta  * int_kernel
    total_v <- ext_v + int_v

    list(
        ext    = ext_v,
        int    = int_v,
        total  = total_v,
        mean_q  = mean(q_v),
        mean_mu = mean(g_v)        # report ZT-corrected mean intensity
    )
}


# ---- 5. .ame_summarize --- Posterior summary for M x P matrix ---------------

#' Compute posterior summary for an AME draw matrix
#'
#' @param ame_mat M x P numeric matrix of AME draws.
#' @param labels Character vector of length P: covariate labels.
#' @param level Numeric scalar in (0,1): credible interval level.
#' @return Data frame with columns: covariate, post_mean, post_median,
#'   ci_lo, ci_hi, post_sd, pr_positive.
#' @keywords internal
.ame_summarize <- function(ame_mat, labels, level) {

    alpha_half <- (1 - level) / 2
    probs <- c(alpha_half, 1 - alpha_half)

    quantiles <- apply(ame_mat, 2, quantile, probs = probs, na.rm = TRUE)

    data.frame(
        covariate   = labels,
        post_mean   = colMeans(ame_mat, na.rm = TRUE),
        post_median = apply(ame_mat, 2, median, na.rm = TRUE),
        ci_lo       = quantiles[1, ],
        ci_hi       = quantiles[2, ],
        post_sd     = apply(ame_mat, 2, sd, na.rm = TRUE),
        pr_positive = colMeans(ame_mat > 0, na.rm = TRUE),
        stringsAsFactors = FALSE,
        row.names   = NULL
    )
}


# ---- 6. .ame_build_decomp --- Build decomposition data.frame ---------------

#' Build the AME decomposition table
#'
#' @param ext_summ Data frame: extensive summary from .ame_summarize.
#' @param int_summ Data frame: intensive summary from .ame_summarize.
#' @param total_summ Data frame: total summary from .ame_summarize.
#' @param labels Character vector of length P.
#' @param idx Integer vector: indices of non-intercept covariates.
#' @return Data frame with 13 columns.
#' @keywords internal
.ame_build_decomp <- function(ext_summ, int_summ, total_summ, labels, idx) {

    if (length(idx) == 0L) {
        return(data.frame(
            covariate    = character(0),
            ext_ame      = numeric(0),
            ext_ci_lo    = numeric(0),
            ext_ci_hi    = numeric(0),
            int_ame      = numeric(0),
            int_ci_lo    = numeric(0),
            int_ci_hi    = numeric(0),
            total_ame    = numeric(0),
            total_ci_lo  = numeric(0),
            total_ci_hi  = numeric(0),
            ext_share    = numeric(0),
            int_share    = numeric(0),
            sign_pattern = character(0),
            stringsAsFactors = FALSE
        ))
    }

    decomp <- data.frame(
        covariate   = labels[idx],
        ext_ame     = ext_summ$post_mean[idx],
        ext_ci_lo   = ext_summ$ci_lo[idx],
        ext_ci_hi   = ext_summ$ci_hi[idx],
        int_ame     = int_summ$post_mean[idx],
        int_ci_lo   = int_summ$ci_lo[idx],
        int_ci_hi   = int_summ$ci_hi[idx],
        total_ame   = total_summ$post_mean[idx],
        total_ci_lo = total_summ$ci_lo[idx],
        total_ci_hi = total_summ$ci_hi[idx],
        stringsAsFactors = FALSE,
        row.names   = NULL
    )

    # Shares based on absolute magnitudes (division-by-zero guard)
    abs_ext <- abs(decomp$ext_ame)
    abs_int <- abs(decomp$int_ame)
    denom   <- abs_ext + abs_int
    safe_denom <- pmax(denom, .Machine$double.eps)
    decomp$ext_share <- abs_ext / safe_denom * 100
    decomp$int_share <- abs_int / safe_denom * 100

    # Sign pattern
    decomp$sign_pattern <- ifelse(
        sign(decomp$ext_ame) == sign(decomp$int_ame),
        "reinforcing", "opposing"
    )

    decomp
}


# ---- 7. .ame_compute_wald --- Delta-method Wald CIs -----------------------

#' Compute Wald delta-method CIs for AME
#'
#' Uses central-difference numerical gradient of the AME with respect to
#' the D-vector theta, combined with the sandwich variance V_sand, to
#' obtain standard errors via the delta method:
#' SE(AME_k) = sqrt(grad_k' V_sand grad_k).
#'
#' @param theta_hat Numeric D-vector: posterior mean.
#' @param X Numeric N x P matrix.
#' @param n_trial Integer vector of length N: trial sizes for ZT correction.
#' @param V_sand D x D sandwich variance matrix.
#' @param P Integer: number of covariates.
#' @param D Integer: total parameters (2P + 1).
#' @param level Confidence level.
#' @return Data frame with P rows and 13 columns: covariate,
#'   ext_point, ext_se, ext_lo, ext_hi, int_point, int_se,
#'   int_lo, int_hi, total_point, total_se, total_lo, total_hi.
#' @keywords internal
.ame_compute_wald <- function(theta_hat, X, n_trial, V_sand, P, D, level) {

    eps <- 1e-5
    z_crit <- qnorm((1 + level) / 2)

    # Covariate labels
    if (!is.null(colnames(X))) {
        cov_labels <- colnames(X)
    } else {
        cov_labels <- c("intercept", paste0("x", seq_len(P - 1L)))
    }

    # Compute numerical gradient: P x D matrices for ext, int, total
    grad_ext   <- matrix(0, nrow = P, ncol = D)
    grad_int   <- matrix(0, nrow = P, ncol = D)
    grad_total <- matrix(0, nrow = P, ncol = D)

    for (d in seq_len(D)) {
        theta_plus  <- theta_hat
        theta_minus <- theta_hat
        theta_plus[d]  <- theta_plus[d]  + eps
        theta_minus[d] <- theta_minus[d] - eps

        ame_plus  <- .ame_at_theta(theta_plus,  X, n_trial, P)
        ame_minus <- .ame_at_theta(theta_minus, X, n_trial, P)

        grad_ext[, d]   <- (ame_plus$ext   - ame_minus$ext)   / (2 * eps)
        grad_int[, d]   <- (ame_plus$int   - ame_minus$int)   / (2 * eps)
        grad_total[, d] <- (ame_plus$total - ame_minus$total) / (2 * eps)
    }

    # Delta-method variance: Var(AME_k) = grad_k' V_sand grad_k
    var_ext   <- numeric(P)
    var_int   <- numeric(P)
    var_total <- numeric(P)

    for (k in seq_len(P)) {
        var_ext[k]   <- as.numeric(
            grad_ext[k, ] %*% V_sand %*% grad_ext[k, ]
        )
        var_int[k]   <- as.numeric(
            grad_int[k, ] %*% V_sand %*% grad_int[k, ]
        )
        var_total[k] <- as.numeric(
            grad_total[k, ] %*% V_sand %*% grad_total[k, ]
        )
    }

    # Clamp negative variance to 0
    se_ext   <- sqrt(pmax(var_ext, 0))
    se_int   <- sqrt(pmax(var_int, 0))
    se_total <- sqrt(pmax(var_total, 0))

    # Point estimates at theta_hat
    ame_hat <- .ame_at_theta(theta_hat, X, n_trial, P)

    data.frame(
        covariate   = cov_labels,
        ext_point   = ame_hat$ext,
        ext_se      = se_ext,
        ext_lo      = ame_hat$ext - z_crit * se_ext,
        ext_hi      = ame_hat$ext + z_crit * se_ext,
        int_point   = ame_hat$int,
        int_se      = se_int,
        int_lo      = ame_hat$int - z_crit * se_int,
        int_hi      = ame_hat$int + z_crit * se_int,
        total_point = ame_hat$total,
        total_se    = se_total,
        total_lo    = ame_hat$total - z_crit * se_total,
        total_hi    = ame_hat$total + z_crit * se_total,
        stringsAsFactors = FALSE,
        row.names   = NULL
    )
}


# ---- 8. .ame_reversal_probs --- Reversal probability computation ----------

#' Compute reversal probabilities for each non-intercept covariate
#'
#' For each non-intercept covariate, computes the posterior probability
#' that the extensive and intensive AME components have opposing signs.
#'
#' @param ext_draws M by P matrix of extensive AME draws.
#' @param int_draws M by P matrix of intensive AME draws.
#' @param labels Character P-vector of covariate labels.
#' @param idx Integer vector of non-intercept indices.
#' @return Named numeric vector of reversal probabilities.
#' @keywords internal
.ame_reversal_probs <- function(ext_draws, int_draws, labels, idx) {

    if (length(idx) == 0L) {
        return(structure(numeric(0), names = character(0)))
    }

    probs <- numeric(length(idx))
    names(probs) <- labels[idx]

    for (i in seq_along(idx)) {
        k <- idx[i]
        ext_k <- ext_draws[, k]
        int_k <- int_draws[, k]

        # Opposing signs: (ext < 0 AND int > 0) OR (ext > 0 AND int < 0)
        probs[i] <- mean(
            (ext_k < 0 & int_k > 0) | (ext_k > 0 & int_k < 0),
            na.rm = TRUE
        )
    }

    probs
}


# ---- 9. .ame_validate_inputs --- Comprehensive input validation ------------

#' Validate inputs for the ame() function
#'
#' @param fit An object to validate as hbb_fit.
#' @param cholesky An optional hbb_cholesky object.
#' @param sandwich An optional hbb_sandwich object.
#' @param level Confidence level.
#' @param n_draws Optional integer.
#' @return Invisible NULL; aborts on failure.
#' @keywords internal
.ame_validate_inputs <- function(fit, cholesky, sandwich, level, n_draws) {

    # -- fit class
    if (!inherits(fit, "hbb_fit")) {
        cli::cli_abort(c(
            "x" = "{.arg fit} must be an object of class {.cls hbb_fit}.",
            "i" = "Received class: {.cls {class(fit)}}.",
            "i" = "Use {.fn hbb} to create a fitted model object."
        ))
    }

    # -- hbb_data fields
    hd <- fit$hbb_data
    if (is.null(hd)) {
        cli::cli_abort("{.code fit$hbb_data} is {.val NULL}.")
    }
    if (is.null(hd$X) || !is.matrix(hd$X)) {
        cli::cli_abort("{.code fit$hbb_data$X} must be a numeric matrix.")
    }
    if (is.null(hd$N) || is.null(hd$P)) {
        cli::cli_abort(
            "{.code fit$hbb_data} must contain {.field N} and {.field P}."
        )
    }
    if (nrow(hd$X) != hd$N || ncol(hd$X) != hd$P) {
        cli::cli_abort(
            "Design matrix X dimensions ({nrow(hd$X)} x {ncol(hd$X)}) do not match N={hd$N}, P={hd$P}."
        )
    }

    # -- fit$fit check (need draws unless cholesky provided)
    if (is.null(cholesky) && is.null(fit$fit)) {
        cli::cli_abort(c(
            "x" = "{.code fit$fit} is {.val NULL} and no {.arg cholesky} object provided.",
            "i" = "Either provide a {.cls hbb_cholesky} object or ensure the fit contains MCMC draws."
        ))
    }

    # -- cholesky (optional)
    P <- hd$P
    D <- 2L * P + 1L
    if (!is.null(cholesky)) {
        if (!inherits(cholesky, "hbb_cholesky")) {
            cli::cli_abort(c(
                "x" = "{.arg cholesky} must be an {.cls hbb_cholesky} object or {.val NULL}.",
                "i" = "Received class: {.cls {class(cholesky)}}."
            ))
        }
        if (!is.null(cholesky$theta_corrected) &&
            ncol(cholesky$theta_corrected) != D) {
            cli::cli_abort(
                "cholesky$theta_corrected has {ncol(cholesky$theta_corrected)} columns but expected D = {D}."
            )
        }
    }

    # -- sandwich (optional)
    if (!is.null(sandwich)) {
        if (!inherits(sandwich, "hbb_sandwich")) {
            cli::cli_abort(c(
                "x" = "{.arg sandwich} must be an {.cls hbb_sandwich} object or {.val NULL}.",
                "i" = "Received class: {.cls {class(sandwich)}}."
            ))
        }
    }

    # -- level in (0, 1)
    .validate_level(level)

    # -- n_draws must be NULL or positive integer
    if (!is.null(n_draws)) {
        if (!is.numeric(n_draws) || length(n_draws) != 1L ||
            is.na(n_draws) || n_draws < 1 || n_draws != round(n_draws)) {
            cli::cli_abort(c(
                "x" = "{.arg n_draws} must be a positive integer or {.val NULL}.",
                "i" = "Received: {.val {n_draws}}."
            ))
        }
    }

    invisible(NULL)
}


# ---- 10. .ame_extract_draws --- Draw extraction logic ----------------------

#' Extract parameter draws for AME computation
#'
#' If a cholesky object is provided, uses theta_corrected and theta_hat.
#' Otherwise, extracts draws from the CmdStanR fit.
#'
#' @param fit An hbb_fit object.
#' @param cholesky Optional hbb_cholesky object.
#' @param P Integer: number of covariates per margin.
#' @return Named list: draws (M x D matrix), theta_hat (D-vector),
#'   M_total (integer).
#' @keywords internal
.ame_extract_draws <- function(fit, cholesky, P) {

    D <- 2L * P + 1L

    if (!is.null(cholesky)) {
        draws     <- cholesky$theta_corrected
        theta_hat <- cholesky$theta_hat
        M_total   <- nrow(draws)

        cli::cli_alert_info("Using Cholesky-corrected draws")
    } else {
        # Extract from CmdStanR fit
        if (is.null(fit$fit)) {
            cli::cli_abort(c(
                "x" = "{.code fit$fit} is {.val NULL}.",
                "i" = "Either provide a {.cls hbb_cholesky} object via {.arg cholesky}",
                "i" = "or ensure the hbb_fit contains a valid CmdStanR fit."
            ))
        }

        param_names <- c(
            paste0("alpha[", seq_len(P), "]"),
            paste0("beta[", seq_len(P), "]"),
            "log_kappa"
        )

        draws <- tryCatch(
            fit$fit$draws(variables = param_names, format = "matrix"),
            error = function(e) {
                cli::cli_abort(c(
                    "x" = "Failed to extract draws from {.cls hbb_fit} object.",
                    "i" = "Error: {conditionMessage(e)}"
                ))
            }
        )

        # Ensure plain numeric matrix (posterior::draws_matrix subsetting
        # preserves class, which breaks %*% in the AME loop)
        draws <- as.matrix(unclass(draws))
        theta_hat <- colMeans(draws)
        M_total   <- nrow(draws)

        cli::cli_alert_info("Using raw MCMC draws (no Cholesky correction)")
    }

    list(draws = draws, theta_hat = theta_hat, M_total = M_total)
}
