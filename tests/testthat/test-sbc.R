# ============================================================================
# test-sbc.R --- Simulation-Based Calibration (SBC) for Hurdle BB Models
#
# Validates posterior calibration by the Talts et al. (2018) protocol:
#   1. Draw theta* from the prior
#   2. Simulate y ~ p(y | theta*)
#   3. Fit the model to get posterior draws theta^(1), ..., theta^(M)
#   4. Compute rank statistic R_d = #{theta^(m)_d < theta*_d}
#   5. Repeat B times and test H0: R_d ~ Uniform(0, M)
#
# Tests require CmdStan and are skipped on CRAN.
# Expected runtime: ~1-2 minutes (base model, 100 replications, N=100).
#
# Sections:
#   1. Configuration
#   2. SBC Replication Runner
#   3. Uniformity Test Helper
#   4. Run SBC (stores rank matrix)
#   5. Per-Parameter Uniformity Checks
#   6. Aggregate Calibration Checks
# ============================================================================


# ============================================================================
# Section 1: Configuration
# ============================================================================

.sbc_B         <- 100L       # number of replications
.sbc_N         <- 100L       # observations per replication
.sbc_P         <- 2L         # intercept + 1 covariate
.sbc_D         <- 2L * .sbc_P + 1L  # total parameters: alpha[1:P], beta[1:P], log_kappa
.sbc_chains    <- 2L
.sbc_warmup    <- 300L
.sbc_sampling  <- 300L
.sbc_M         <- .sbc_chains * .sbc_sampling  # 600 total posterior draws
.sbc_seed      <- 2026L      # master seed for reproducibility

# Tighter priors (same used in prior draws AND in Stan model)
# This avoids pathological data realizations from extreme prior tails.
.sbc_prior <- hbb_prior(
    alpha     = list(dist = "normal", mean = 0, sd = 1),
    beta      = list(dist = "normal", mean = 0, sd = 1),
    log_kappa = list(dist = "normal", mean = 2, sd = 0.5)
)

# Parameter names (must match Stan model)
.sbc_param_names <- c(
    paste0("alpha[", seq_len(.sbc_P), "]"),
    paste0("beta[", seq_len(.sbc_P), "]"),
    "log_kappa"
)

# Human-readable labels for test output
.sbc_param_labels <- c(
    paste0("alpha_", seq_len(.sbc_P)),
    paste0("beta_", seq_len(.sbc_P)),
    "log_kappa"
)

# Storage for results (shared across test_that blocks)
.sbc_ranks  <- NULL   # B x D matrix of rank statistics
.sbc_n_ok   <- 0L     # number of successful replications


# ============================================================================
# Section 2: SBC Replication Runner
# ============================================================================

#' Run a single SBC replication
#'
#' @param alpha_star Numeric vector of length P (true extensive-margin coefs).
#' @param beta_star  Numeric vector of length P (true intensive-margin coefs).
#' @param lk_star    Numeric scalar (true log_kappa).
#' @param rep_id     Integer. Replication index (used as MCMC seed).
#' @return Named numeric vector of length D with rank statistics,
#'   or NULL if the replication failed.
.sbc_one_rep <- function(alpha_star, beta_star, lk_star, rep_id) {

    kappa_star <- exp(lk_star)
    theta_star <- c(alpha_star, beta_star, lk_star)

    # -- Generate covariates and trial sizes --
    X_cov   <- rnorm(.sbc_N)
    X       <- cbind(1, X_cov)
    n_trial <- sample(10L:30L, .sbc_N, replace = TRUE)

    # -- Linear predictors and link --
    eta_ext <- as.numeric(X %*% alpha_star)
    eta_int <- as.numeric(X %*% beta_star)
    q       <- plogis(eta_ext)
    mu      <- plogis(eta_int)

    # -- Hurdle response --
    z <- rbinom(.sbc_N, 1L, q)
    y <- integer(.sbc_N)

    active <- which(z == 1L)
    if (length(active) > 0L) {
        y[active] <- hurdlebb::rztbetabinom(
            nn    = length(active),
            n     = n_trial[active],
            mu    = mu[active],
            kappa = kappa_star
        )
    }

    # -- Safety: skip if too few or too many zeros --
    n_nonzero <- sum(y > 0L)
    if (n_nonzero < 5L || n_nonzero > .sbc_N - 5L) return(NULL)

    # -- Fit model --
    sim_data <- data.frame(y = y, n_trial = n_trial, x1 = X_cov)
    fit <- tryCatch(
        hbb(
            formula       = y | trials(n_trial) ~ x1,
            data          = sim_data,
            prior         = .sbc_prior,
            chains        = .sbc_chains,
            iter_warmup   = .sbc_warmup,
            iter_sampling = .sbc_sampling,
            seed          = rep_id,
            refresh       = 0L
        ),
        error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)

    # -- Extract posterior draws (M x D matrix) --
    draws <- tryCatch(
        fit$fit$draws(variables = .sbc_param_names, format = "matrix"),
        error = function(e) NULL
    )
    if (is.null(draws) || nrow(draws) != .sbc_M) return(NULL)

    # -- Compute rank statistic --
    ranks <- vapply(seq_len(.sbc_D), function(d) {
        sum(draws[, d] < theta_star[d])
    }, numeric(1L))

    names(ranks) <- .sbc_param_labels
    ranks
}


# ============================================================================
# Section 3: Uniformity Test Helper
# ============================================================================

#' Chi-squared goodness-of-fit test for uniformity
#'
#' Bins the rank statistics into J equal-width bins and tests
#' H0: ranks ~ Uniform(0, M) via chi-squared.
#'
#' @param ranks  Numeric vector of rank statistics (subset of 0:M).
#' @param M      Integer. Number of posterior draws.
#' @param J      Integer. Number of bins (default: 10).
#' @return A list with `stat` (chi-squared statistic), `df` (degrees of
#'   freedom), and `p_value`.
.sbc_chisq_test <- function(ranks, M, J = 10L) {

    # Remove NAs (failed replications)
    ranks <- ranks[!is.na(ranks)]
    B <- length(ranks)
    if (B < 20L) return(list(stat = NA, df = NA, p_value = NA))

    # Bin edges: 0, M/J, 2M/J, ..., M
    breaks <- seq(0, M + 1, length.out = J + 1L)
    obs    <- tabulate(findInterval(ranks, breaks, rightmost.closed = TRUE), J)
    expected <- B / J

    stat <- sum((obs - expected)^2 / expected)
    df   <- J - 1L
    p    <- pchisq(stat, df = df, lower.tail = FALSE)

    list(stat = stat, df = df, p_value = p)
}


# ============================================================================
# Section 4: Run SBC (100 replications)
# ============================================================================

test_that("SBC: all replications complete successfully", {
    skip_if_no_cmdstan()
    skip_on_cran()

    set.seed(.sbc_seed)

    # -- Pre-draw all true parameters from the prior --
    prior_alpha     <- matrix(rnorm(.sbc_B * .sbc_P, 0, 1),
                              nrow = .sbc_B, ncol = .sbc_P)
    prior_beta      <- matrix(rnorm(.sbc_B * .sbc_P, 0, 1),
                              nrow = .sbc_B, ncol = .sbc_P)
    prior_log_kappa <- rnorm(.sbc_B, mean = 2, sd = 0.5)

    # -- Set new RNG state for data generation (independent of prior draws) --
    # This ensures prior draws and data generation use separate streams.
    set.seed(.sbc_seed + 1L)

    rank_mat <- matrix(NA_real_, nrow = .sbc_B, ncol = .sbc_D,
                       dimnames = list(NULL, .sbc_param_labels))
    n_ok <- 0L

    for (b in seq_len(.sbc_B)) {

        result <- .sbc_one_rep(
            alpha_star = prior_alpha[b, ],
            beta_star  = prior_beta[b, ],
            lk_star    = prior_log_kappa[b],
            rep_id     = b
        )

        if (!is.null(result)) {
            rank_mat[b, ] <- result
            n_ok <- n_ok + 1L
        }
    }

    # Store for downstream tests
    .sbc_ranks <<- rank_mat
    .sbc_n_ok  <<- n_ok

    # At least 80 of 100 replications should succeed
    expect_true(
        n_ok >= 80L,
        info = sprintf("Only %d/%d SBC reps succeeded (need >=80)", n_ok, .sbc_B)
    )
})


# ============================================================================
# Section 5: Per-Parameter Uniformity Checks
# ============================================================================

# One test per parameter: chi-squared p-value > 0.001 (very conservative)

test_that("SBC: alpha_1 ranks are uniform", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.sbc_ranks), "SBC ranks not available")

    result <- .sbc_chisq_test(.sbc_ranks[, "alpha_1"], .sbc_M)
    skip_if(is.na(result$p_value), "Too few successful reps for alpha_1")

    expect_true(
        result$p_value > 0.001,
        info = sprintf(
            "alpha_1: chi2=%.2f, df=%d, p=%.4f (need p>0.001)",
            result$stat, result$df, result$p_value
        )
    )
})


test_that("SBC: alpha_2 ranks are uniform", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.sbc_ranks), "SBC ranks not available")

    result <- .sbc_chisq_test(.sbc_ranks[, "alpha_2"], .sbc_M)
    skip_if(is.na(result$p_value), "Too few successful reps for alpha_2")

    expect_true(
        result$p_value > 0.001,
        info = sprintf(
            "alpha_2: chi2=%.2f, df=%d, p=%.4f (need p>0.001)",
            result$stat, result$df, result$p_value
        )
    )
})


test_that("SBC: beta_1 ranks are uniform", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.sbc_ranks), "SBC ranks not available")

    result <- .sbc_chisq_test(.sbc_ranks[, "beta_1"], .sbc_M)
    skip_if(is.na(result$p_value), "Too few successful reps for beta_1")

    expect_true(
        result$p_value > 0.001,
        info = sprintf(
            "beta_1: chi2=%.2f, df=%d, p=%.4f (need p>0.001)",
            result$stat, result$df, result$p_value
        )
    )
})


test_that("SBC: beta_2 ranks are uniform", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.sbc_ranks), "SBC ranks not available")

    result <- .sbc_chisq_test(.sbc_ranks[, "beta_2"], .sbc_M)
    skip_if(is.na(result$p_value), "Too few successful reps for beta_2")

    expect_true(
        result$p_value > 0.001,
        info = sprintf(
            "beta_2: chi2=%.2f, df=%d, p=%.4f (need p>0.001)",
            result$stat, result$df, result$p_value
        )
    )
})


test_that("SBC: log_kappa ranks are uniform", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.sbc_ranks), "SBC ranks not available")

    result <- .sbc_chisq_test(.sbc_ranks[, "log_kappa"], .sbc_M)
    skip_if(is.na(result$p_value), "Too few successful reps for log_kappa")

    expect_true(
        result$p_value > 0.001,
        info = sprintf(
            "log_kappa: chi2=%.2f, df=%d, p=%.4f (need p>0.001)",
            result$stat, result$df, result$p_value
        )
    )
})


# ============================================================================
# Section 6: Aggregate Calibration Checks
# ============================================================================

test_that("SBC: at most 1 parameter shows marginal miscalibration (p < 0.01)", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.sbc_ranks), "SBC ranks not available")

    p_values <- vapply(seq_len(.sbc_D), function(d) {
        result <- .sbc_chisq_test(.sbc_ranks[, d], .sbc_M)
        result$p_value
    }, numeric(1L))

    names(p_values) <- .sbc_param_labels
    n_reject <- sum(!is.na(p_values) & p_values < 0.01)

    expect_true(
        n_reject <= 1L,
        info = sprintf(
            "%d/%d parameters rejected at alpha=0.01: %s",
            n_reject, .sbc_D,
            paste(names(p_values)[!is.na(p_values) & p_values < 0.01],
                  collapse = ", ")
        )
    )
})


test_that("SBC: mean ranks are within 20% of expected midpoint", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.sbc_ranks), "SBC ranks not available")

    # Under uniformity, E[R] = M/2
    expected_mean <- .sbc_M / 2
    tolerance     <- 0.20 * .sbc_M  # 20% of M

    for (d in seq_len(.sbc_D)) {
        ranks_d  <- .sbc_ranks[, d]
        ranks_d  <- ranks_d[!is.na(ranks_d)]
        if (length(ranks_d) < 20L) next

        obs_mean <- mean(ranks_d)
        expect_true(
            abs(obs_mean - expected_mean) < tolerance,
            info = sprintf(
                "%s: mean rank = %.1f, expected = %.1f (tol = %.1f)",
                .sbc_param_labels[d], obs_mean, expected_mean, tolerance
            )
        )
    }
})


test_that("SBC: rank standard deviations are consistent with uniformity", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.sbc_ranks), "SBC ranks not available")

    # Under Uniform(0, M), theoretical SD = M / sqrt(12) ≈ 173.2 for M=600
    expected_sd <- .sbc_M / sqrt(12)

    for (d in seq_len(.sbc_D)) {
        ranks_d <- .sbc_ranks[, d]
        ranks_d <- ranks_d[!is.na(ranks_d)]
        if (length(ranks_d) < 20L) next

        obs_sd <- sd(ranks_d)

        # Allow 50% deviation from theoretical SD
        # (wide tolerance for B=100)
        expect_true(
            obs_sd > 0.5 * expected_sd && obs_sd < 1.5 * expected_sd,
            info = sprintf(
                "%s: rank SD = %.1f, expected = %.1f (50%% tolerance)",
                .sbc_param_labels[d], obs_sd, expected_sd
            )
        )
    }
})
