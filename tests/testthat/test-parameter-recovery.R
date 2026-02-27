# ============================================================================
# test-parameter-recovery.R --- Parameter Recovery Tests for Hurdle BB Models
#
# Tests that the hbb() MCMC sampler recovers known true parameter values
# from simulated data. These tests verify the statistical correctness of
# the full estimation pipeline: data generation -> Stan compilation ->
# MCMC sampling -> posterior summarisation.
#
# Sections:
#   1. Data Generation Helper
#   2. Base Model Parameter Recovery (6 tests)
#   3. Weighted Model Parameter Recovery (3 tests)
#   4. SVC Model with Random Slopes (3 tests)
#   5. Cross-Model Consistency (2 tests)
#
# All tests require CmdStan and are skipped on CRAN.
# ============================================================================


# ============================================================================
# Section 1: Data Generation Helper
# ============================================================================

#' Generate Hurdle Beta-Binomial data with known parameters
#'
#' Creates a data frame from a known DGP for parameter recovery testing.
#' The design matrix X includes an intercept (column 1) plus P-1 standard
#' normal covariates. The response is generated via the two-part hurdle:
#'   z_i ~ Bernoulli(q_i), y_i = ZT-BB(n_i, mu_i, kappa) if z_i=1, else 0.
#'
#' @param N     Integer. Number of observations.
#' @param alpha Numeric vector of length P. Extensive-margin coefficients
#'              (including intercept).
#' @param beta  Numeric vector of length P. Intensive-margin coefficients
#'              (including intercept).
#' @param log_kappa Numeric scalar. Log-dispersion parameter.
#' @param n_trial_range Integer vector of length 2. Range for trial sizes.
#' @param seed  Integer. Random seed.
#' @return A data.frame with columns: y, n_trial, x1, x2, ..., x_{P-1}.
generate_hbb_data <- function(N,
                              alpha,
                              beta,
                              log_kappa,
                              n_trial_range = c(10L, 50L),
                              seed = 42L) {

    stopifnot(length(alpha) == length(beta))
    P <- length(alpha)
    kappa <- exp(log_kappa)

    set.seed(seed)

    # Design matrix: intercept + (P-1) standard normal covariates
    if (P > 1L) {
        X_cov <- matrix(rnorm(N * (P - 1L)), nrow = N, ncol = P - 1L)
        X <- cbind(1, X_cov)
    } else {
        X <- matrix(1, nrow = N, ncol = 1L)
    }

    # Linear predictors
    eta_ext <- as.numeric(X %*% alpha)
    eta_int <- as.numeric(X %*% beta)

    q  <- plogis(eta_ext)
    mu <- plogis(eta_int)

    # Trial sizes
    n_trial <- sample(
        seq(n_trial_range[1], n_trial_range[2]),
        N, replace = TRUE
    )

    # Generate hurdle response
    z <- rbinom(N, size = 1L, prob = q)
    y <- integer(N)

    active <- which(z == 1L)
    if (length(active) > 0L) {
        y[active] <- hurdlebb::rztbetabinom(
            nn    = length(active),
            n     = n_trial[active],
            mu    = mu[active],
            kappa = kappa
        )
    }

    # Assemble data frame
    df <- data.frame(y = y, n_trial = n_trial)
    if (P > 1L) {
        for (j in seq_len(P - 1L)) {
            df[[paste0("x", j)]] <- X[, j + 1L]
        }
    }

    df
}


# ============================================================================
# Section 2: Base Model Parameter Recovery
# ============================================================================

# -- Shared true values and data for Section 2 --------------------------------
.pr_alpha     <- c(0.5, -0.20, 0.25)
.pr_beta      <- c(0.0, 0.12, -0.10)
.pr_log_kappa <- 1.5
.pr_P         <- length(.pr_alpha)    # P = 3 (intercept + 2 covariates)
.pr_N         <- 500L
.pr_seed      <- 13L

# Generate data once for the entire section
.pr_data <- generate_hbb_data(
    N         = .pr_N,
    alpha     = .pr_alpha,
    beta      = .pr_beta,
    log_kappa = .pr_log_kappa,
    seed      = .pr_seed
)

# Fit the base model once (shared across Section 2 tests)
# This is wrapped in a local environment so fit_base is available to all tests
.pr_base_fit <- NULL


test_that("Base model: fits without error", {
    skip_if_no_cmdstan()
    skip_on_cran()

    .pr_base_fit <<- hbb(
        formula       = y | trials(n_trial) ~ x1 + x2,
        data          = .pr_data,
        chains        = 2L,
        iter_warmup   = 500L,
        iter_sampling = 500L,
        seed          = .pr_seed,
        refresh       = 0L
    )

    expect_s3_class(.pr_base_fit, "hbb_fit")
    expect_equal(.pr_base_fit$model_type, "base")
})


test_that("Base model: aggregate parameter coverage >= 85% at 90% CI", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.pr_base_fit), "Base model fit not available")

    s <- summary(.pr_base_fit, level = 0.90)
    fe <- s$fixed_effects

    # Assemble true values: alpha[1..P], beta[1..P], log_kappa
    truth <- c(.pr_alpha, .pr_beta, .pr_log_kappa)
    D <- nrow(fe)
    expect_equal(D, length(truth))

    covered <- sum(truth >= fe$ci_lower & truth <= fe$ci_upper)
    coverage_rate <- covered / D

    # Allow 1 miss out of 7 (85.7%)
    expect_true(
        coverage_rate >= 0.85,
        info = sprintf(
            "Coverage = %.1f%% (%d/%d parameters covered at 90%% level)",
            100 * coverage_rate, covered, D
        )
    )
})


test_that("Base model: log_kappa within 90% CI", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.pr_base_fit), "Base model fit not available")

    s <- summary(.pr_base_fit, level = 0.90)
    fe <- s$fixed_effects

    # log_kappa is the last row (index 2P + 1)
    lk_idx <- 2L * .pr_P + 1L
    expect_true(
        .pr_log_kappa >= fe$ci_lower[lk_idx] &&
            .pr_log_kappa <= fe$ci_upper[lk_idx],
        info = sprintf(
            "log_kappa = %.3f not in 90%% CI [%.3f, %.3f]",
            .pr_log_kappa, fe$ci_lower[lk_idx], fe$ci_upper[lk_idx]
        )
    )
})


test_that("Base model: Rhat < 1.05 for all parameters", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.pr_base_fit), "Base model fit not available")

    s <- summary(.pr_base_fit, level = 0.90)
    fe <- s$fixed_effects

    for (i in seq_len(nrow(fe))) {
        expect_true(
            !is.na(fe$rhat[i]) && fe$rhat[i] < 1.05,
            info = sprintf(
                "%s: Rhat = %.4f (should be < 1.05)",
                fe$parameter[i], fe$rhat[i]
            )
        )
    }
})


test_that("Base model: ESS > 200 for all parameters", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.pr_base_fit), "Base model fit not available")

    s <- summary(.pr_base_fit, level = 0.90)
    fe <- s$fixed_effects

    for (i in seq_len(nrow(fe))) {
        expect_true(
            !is.na(fe$ess_bulk[i]) && fe$ess_bulk[i] > 200,
            info = sprintf(
                "%s: ESS_bulk = %.0f (should be > 200)",
                fe$parameter[i], fe$ess_bulk[i]
            )
        )
    }
})


# ============================================================================
# Section 3: Weighted Model Parameter Recovery
# ============================================================================

# Add weights and survey design columns to the same base data
.pr_wt_data <- .pr_data
set.seed(99L)
.pr_wt_data$weight  <- runif(.pr_N, 0.5, 2.0)
# Rescale weights to sum to N
.pr_wt_data$weight  <- .pr_wt_data$weight / sum(.pr_wt_data$weight) * .pr_N
# Add stratum and PSU columns needed for sandwich variance
.pr_wt_data$stratum <- rep(seq_len(5L), length.out = .pr_N)
.pr_wt_data$psu     <- seq_len(.pr_N)

.pr_wt_fit <- NULL


test_that("Weighted model: fits without error", {
    skip_if_no_cmdstan()
    skip_on_cran()

    .pr_wt_fit <<- hbb(
        formula       = y | trials(n_trial) ~ x1 + x2,
        data          = .pr_wt_data,
        weights       = "weight",
        stratum       = "stratum",
        psu           = "psu",
        chains        = 2L,
        iter_warmup   = 500L,
        iter_sampling = 500L,
        seed          = .pr_seed,
        refresh       = 0L
    )

    expect_s3_class(.pr_wt_fit, "hbb_fit")
    expect_equal(.pr_wt_fit$model_type, "weighted")
})


test_that("Weighted model: fixed effects coverage >= 85%", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.pr_wt_fit), "Weighted model fit not available")

    s <- summary(.pr_wt_fit, level = 0.90)
    fe <- s$fixed_effects

    # Collect all true parameter values
    true_vals <- c(.pr_alpha, .pr_beta, .pr_log_kappa)
    D <- length(true_vals)

    # Count how many are covered
    covered <- vapply(seq_len(D), function(i) {
        true_vals[i] >= fe$ci_lower[i] && true_vals[i] <= fe$ci_upper[i]
    }, logical(1L))

    coverage_rate <- mean(covered)

    # With 7 parameters and 90% CI, we expect ~6.3 to be covered.
    # Allow some slack: require at least 85% (6/7 ~ 85.7%).
    expect_true(
        coverage_rate >= 0.85,
        info = sprintf(
            "Coverage = %.1f%% (%d/%d parameters covered at 90%% level)",
            100 * coverage_rate, sum(covered), D
        )
    )
})


test_that("Weighted model: diagnostics pass (Rhat, ESS)", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.pr_wt_fit), "Weighted model fit not available")

    s <- summary(.pr_wt_fit, level = 0.90)
    fe <- s$fixed_effects

    # All Rhat < 1.05
    expect_true(
        all(!is.na(fe$rhat) & fe$rhat < 1.05),
        info = sprintf("Max Rhat = %.4f", max(fe$rhat, na.rm = TRUE))
    )

    # All ESS > 200
    expect_true(
        all(!is.na(fe$ess_bulk) & fe$ess_bulk > 200),
        info = sprintf("Min ESS = %.0f", min(fe$ess_bulk, na.rm = TRUE))
    )
})


# ============================================================================
# Section 4: SVC Model with Random Slopes
# ============================================================================

# Note: Random-intercept-only models (1 | state_id) map to the "base"
# Stan variant. To test the SVC Stan model, we use (x1 | state_id) which
# includes a random slope and triggers the "svc" model type.

.pr_svc_N <- 300L
.pr_svc_S <- 5L
.pr_svc_alpha <- c(0.5, -0.20)   # intercept + x1
.pr_svc_beta  <- c(0.0, 0.12)    # intercept + x1
.pr_svc_log_kappa <- 1.5

.pr_svc_fit <- NULL


test_that("SVC model: fits without error", {
    skip_if_no_cmdstan()
    skip_on_cran()

    set.seed(77L)

    P <- length(.pr_svc_alpha)

    # Generate covariates
    X_cov <- rnorm(.pr_svc_N)
    X <- cbind(1, X_cov)

    # State membership: uniform random assignment
    state_id <- sample(seq_len(.pr_svc_S), .pr_svc_N, replace = TRUE)

    # State random effects
    delta_ext <- rnorm(.pr_svc_S, mean = 0, sd = 0.3)
    delta_int <- rnorm(.pr_svc_S, mean = 0, sd = 0.2)

    # Linear predictors with state-specific intercept shifts
    eta_ext <- as.numeric(X %*% .pr_svc_alpha) + delta_ext[state_id]
    eta_int <- as.numeric(X %*% .pr_svc_beta) + delta_int[state_id]

    q  <- plogis(eta_ext)
    mu <- plogis(eta_int)
    kappa <- exp(.pr_svc_log_kappa)

    n_trial <- sample(10L:50L, .pr_svc_N, replace = TRUE)
    z <- rbinom(.pr_svc_N, size = 1L, prob = q)
    y <- integer(.pr_svc_N)
    active <- which(z == 1L)
    if (length(active) > 0L) {
        y[active] <- hurdlebb::rztbetabinom(
            nn    = length(active),
            n     = n_trial[active],
            mu    = mu[active],
            kappa = kappa
        )
    }

    svc_data <- data.frame(
        y        = y,
        n_trial  = n_trial,
        x1       = X_cov,
        state_id = state_id
    )

    .pr_svc_fit <<- hbb(
        formula       = y | trials(n_trial) ~ x1 + (x1 | state_id),
        data          = svc_data,
        chains        = 2L,
        iter_warmup   = 500L,
        iter_sampling = 500L,
        seed          = .pr_seed,
        refresh       = 0L
    )

    expect_s3_class(.pr_svc_fit, "hbb_fit")
    expect_equal(.pr_svc_fit$model_type, "svc")
})


test_that("SVC model: global fixed effects (alpha, beta) covered", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.pr_svc_fit), "SVC model fit not available")

    s <- summary(.pr_svc_fit, level = 0.90)
    fe <- s$fixed_effects

    P <- length(.pr_svc_alpha)
    true_vals <- c(.pr_svc_alpha, .pr_svc_beta, .pr_svc_log_kappa)
    D <- length(true_vals)

    covered <- vapply(seq_len(D), function(i) {
        true_vals[i] >= fe$ci_lower[i] && true_vals[i] <= fe$ci_upper[i]
    }, logical(1L))

    coverage_rate <- mean(covered)

    # With random effects absorbing variation, global fixed effects may
    # shift slightly. Require >= 60% coverage (3/5 params).
    expect_true(
        coverage_rate >= 0.60,
        info = sprintf(
            "SVC coverage = %.1f%% (%d/%d parameters covered at 90%% level)",
            100 * coverage_rate, sum(covered), D
        )
    )
})


test_that("SVC model: diagnostics pass", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.pr_svc_fit), "SVC model fit not available")

    s <- summary(.pr_svc_fit, level = 0.90)
    fe <- s$fixed_effects

    # Rhat < 1.10 (more relaxed for SVC models with partial pooling)
    expect_true(
        all(!is.na(fe$rhat) & fe$rhat < 1.10),
        info = sprintf("Max Rhat = %.4f", max(fe$rhat, na.rm = TRUE))
    )

    # ESS > 100 (lower threshold for SVC due to hierarchical structure)
    expect_true(
        all(!is.na(fe$ess_bulk) & fe$ess_bulk > 100),
        info = sprintf("Min ESS = %.0f", min(fe$ess_bulk, na.rm = TRUE))
    )
})


# ============================================================================
# Section 5: Cross-Model Consistency
# ============================================================================

test_that("Cross-model: LOO runs without error on base model", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.pr_base_fit), "Base model fit not available")

    loo_result <- expect_no_error(loo(.pr_base_fit))

    # Basic structure checks
    expect_true(inherits(loo_result, "loo") || inherits(loo_result, "psis_loo"))
    expect_true("estimates" %in% names(loo_result))
})


test_that("Cross-model: sandwich variance computes without error on weighted model", {
    skip_if_no_cmdstan()
    skip_on_cran()
    skip_if(is.null(.pr_wt_fit), "Weighted model fit not available")

    sand <- expect_no_error(sandwich_variance(.pr_wt_fit))

    # Basic structure checks
    expect_true(inherits(sand, "hbb_sandwich"))
    expect_true("V_sand" %in% names(sand))

    D <- 2L * .pr_P + 1L
    expect_equal(nrow(sand$V_sand), D)
    expect_equal(ncol(sand$V_sand), D)

    # V_sand diagonal should be positive (variances)
    expect_true(all(diag(sand$V_sand) > 0))

    # DER should be present and reasonable (> 0)
    if (!is.null(sand$DER)) {
        expect_true(all(sand$DER > 0))
    }
})
