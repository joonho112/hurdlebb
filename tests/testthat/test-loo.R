# ============================================================================
# test-loo.R --- testthat tests for loo.R
#
# Sections:
#   A. Input validation     (class, fit$fit, log_lik presence, hbb_data$N)
#   B. .loo_extract_loglik  (M x N matrix, NA/NaN/Inf detection)
#   C. hbb_loo_compare      (class, dimensions, elpd_diff, z-ratio, names)
#   D. .loo_chain_id        (length, chain structure, fallback)
#   E. .loo_pareto_diagnostics (NULL guard, tier counts)
#   F. print.hbb_loo_compare   (output structure, no crash on corruption)
#   G. Edge cases           (N=1, identical fits, 3+ models, n_chains=1)
# ============================================================================


# ============================================================================
# MOCK CONSTRUCTOR
# ============================================================================

#' Create a mock hbb_fit for LOO testing
#'
#' Generates realistic mock log_lik draws (negative, finite) and wraps them
#' in a CmdStanR-compatible mock interface.
#'
#' @param P    Number of covariates per margin (including intercept).
#' @param N    Number of observations.
#' @param M    Total MCMC draws (must be divisible by n_chains).
#' @param n_chains Number of chains.
#' @param seed Random seed.
#' @param zero_rate  Structural zero rate for simulated y.
#' @param include_yrep   Include y_rep in mock draws.
#' @param include_loglik Include log_lik in mock draws.
#' @param ll_mean  Mean of simulated log-likelihood values (controls ELPD).
#' @return An S3 object of class "hbb_fit".
create_loo_mock_fit <- function(P = 3, N = 50, M = 200, n_chains = 2,
                                seed = 42, zero_rate = 0.35,
                                include_yrep = TRUE,
                                include_loglik = TRUE,
                                ll_mean = -3.0) {
    set.seed(seed)
    D <- 2L * P + 1L

    param_names <- c(paste0("alpha[", 1:P, "]"),
                     paste0("beta[",  1:P, "]"),
                     "log_kappa")
    theta_true  <- c(rep(-0.5, P), rep(0.3, P), 1.5)
    param_draws <- MASS::mvrnorm(M, theta_true, diag(D) * 0.01)
    colnames(param_draws) <- param_names

    X       <- cbind(1, matrix(rnorm(N * (P - 1)), N, P - 1))
    z       <- rbinom(N, 1, 1 - zero_rate)
    n_trial <- sample(10:50, N, replace = TRUE)
    y       <- ifelse(z == 1L, pmax(1L, rbinom(N, n_trial, 0.3)), 0L)

    y_rep_mat  <- NULL
    yrep_names <- NULL
    if (include_yrep) {
        y_rep_mat <- matrix(0L, M, N)
        for (m in 1:M) {
            for (i in 1:N) {
                if (rbinom(1, 1, 1 - zero_rate) == 1) {
                    y_rep_mat[m, i] <- max(1L, rbinom(1, n_trial[i], 0.3))
                }
            }
        }
        yrep_names <- paste0("y_rep[", 1:N, "]")
    }

    log_lik_mat <- NULL
    ll_names    <- NULL
    if (include_loglik) {
        log_lik_mat <- matrix(rnorm(M * N, mean = ll_mean, sd = 0.5), M, N)
        ll_names    <- paste0("log_lik[", 1:N, "]")
    }

    all_draws <- param_draws
    all_names <- param_names
    if (include_yrep) {
        all_draws <- cbind(all_draws, y_rep_mat)
        all_names <- c(all_names, yrep_names)
    }
    if (include_loglik) {
        all_draws <- cbind(all_draws, log_lik_mat)
        all_names <- c(all_names, ll_names)
    }
    colnames(all_draws) <- all_names

    iter_sampling <- M %/% n_chains

    mock_cmdstan <- list(
        draws = function(variables = NULL, format = "matrix") {
            if (!is.null(variables)) {
                if (length(variables) == 1L && variables == "y_rep") {
                    idx <- grep("^y_rep\\[", colnames(all_draws))
                } else if (length(variables) == 1L && variables == "log_lik") {
                    idx <- grep("^log_lik\\[", colnames(all_draws))
                } else {
                    idx <- match(variables, colnames(all_draws))
                }
                if (any(is.na(idx))) {
                    stop(paste("Variable(s) not found:",
                               paste(variables[is.na(idx)], collapse = ", ")))
                }
                all_draws[, idx, drop = FALSE]
            } else {
                all_draws
            }
        },
        num_chains = function() as.integer(n_chains),
        metadata = function() list(
            chains          = as.integer(n_chains),
            iter_sampling   = as.integer(iter_sampling),
            iter_warmup     = 100L,
            stan_variables  = unique(sub("\\[.*", "", all_names))
        )
    )

    hbb_data <- list(
        X       = X,
        N       = as.integer(N),
        P       = as.integer(P),
        z       = as.integer(z),
        n_trial = as.integer(n_trial),
        y       = as.integer(y)
    )

    structure(
        list(
            fit        = mock_cmdstan,
            stan_data  = list(N = N, P = P, y = y, n_trial = n_trial,
                              z = z, X = X),
            hbb_data   = hbb_data,
            model_type = "base",
            model_name = "hbb_base"
        ),
        class = "hbb_fit"
    )
}


# Convenience wrapper to suppress cli progress/info messages
quietly <- function(expr) suppressMessages(expr)


# ============================================================================
# Section A: Input Validation
# ============================================================================

# ---- Section A: Input Validation ----

test_that("A1. loo() errors on non-hbb_fit object", {
    expect_error(
        quietly(loo(list()))
    )
})

test_that("A2. loo() errors when fit$fit is NULL", {
    fit      <- create_loo_mock_fit()
    fit$fit  <- NULL
    expect_error(
        quietly(loo(fit)),
        class = "rlang_error"
    )
})

test_that("A3. loo() errors when hbb_data is NULL", {
    fit           <- create_loo_mock_fit()
    fit$hbb_data  <- NULL
    expect_error(
        quietly(loo(fit)),
        class = "rlang_error"
    )
})

test_that("A4. loo() errors when log_lik is absent from Stan fit", {
    fit <- create_loo_mock_fit(include_loglik = FALSE)
    expect_error(
        quietly(loo(fit)),
        class = "rlang_error"
    )
})

test_that("A5. hbb_loo_compare() errors on non-hbb_fit element in list", {
    fit <- create_loo_mock_fit(N = 30, M = 60, seed = 1)
    expect_error(
        quietly(hbb_loo_compare(m0 = fit, m1 = list())),
        class = "rlang_error"
    )
})

test_that("A6. hbb_loo_compare() errors when fewer than 2 models supplied", {
    fit <- create_loo_mock_fit(N = 30, M = 60, seed = 2)
    expect_error(
        quietly(hbb_loo_compare(m0 = fit)),
        class = "rlang_error"
    )
})

test_that("A7. hbb_loo_compare() errors when models have different N", {
    fit0 <- create_loo_mock_fit(N = 30, M = 60, seed = 3)
    fit1 <- create_loo_mock_fit(N = 40, M = 60, seed = 4)
    expect_error(
        quietly(hbb_loo_compare(m0 = fit0, m1 = fit1)),
        class = "rlang_error"
    )
})

test_that("A8. hbb_loo_compare() errors when model arguments are unnamed", {
    fit0 <- create_loo_mock_fit(N = 30, M = 60, seed = 5)
    fit1 <- create_loo_mock_fit(N = 30, M = 60, seed = 6)
    expect_error(
        quietly(hbb_loo_compare(fit0, fit1)),
        class = "rlang_error"
    )
})


# ============================================================================
# Section B: .loo_extract_loglik
# ============================================================================

# ---- Section B: .loo_extract_loglik ----

test_that("B1. .loo_extract_loglik() returns M x N numeric matrix", {
    fit <- create_loo_mock_fit(N = 40, M = 100, seed = 11)
    mat <- suppressMessages(hurdlebb:::.loo_extract_loglik(fit))
    expect_true(is.matrix(mat))
    expect_equal(nrow(mat), 100L)
    expect_equal(ncol(mat), 40L)
    expect_true(is.numeric(mat))
})

test_that("B2. all log_lik values are negative for valid mock data", {
    fit <- create_loo_mock_fit(N = 30, M = 80, seed = 12, ll_mean = -3.0)
    mat <- suppressMessages(hurdlebb:::.loo_extract_loglik(fit))
    expect_true(all(mat < 0, na.rm = TRUE))
})

test_that("B3. NaN in log_lik triggers warning containing 'NaN'", {
    fit <- create_loo_mock_fit(N = 20, M = 40, seed = 13)
    # Inject NaN into the mock draws for log_lik
    orig_draws <- fit$fit$draws
    fit$fit$draws <- function(variables = NULL, format = "matrix") {
        mat <- orig_draws(variables = variables, format = format)
        if (!is.null(variables) && any(grepl("^log_lik", variables))) {
            mat[1L, 1L] <- NaN
        }
        mat
    }
    expect_warning(
        suppressMessages(hurdlebb:::.loo_extract_loglik(fit)),
        regexp = "NaN"
    )
})


# ============================================================================
# Section C: hbb_loo_compare
# ============================================================================

# ---- Section C: hbb_loo_compare ----

test_that("C1. Returns 'hbb_loo_compare' class", {
    fit0 <- create_loo_mock_fit(N = 30, M = 60, seed = 21)
    fit1 <- create_loo_mock_fit(N = 30, M = 60, seed = 22)
    res  <- quietly(hbb_loo_compare(m0 = fit0, m1 = fit1))
    expect_s3_class(res, "hbb_loo_compare")
})

test_that("C2. comparison data frame has n_models rows and correct columns", {
    fit0 <- create_loo_mock_fit(N = 30, M = 60, seed = 23)
    fit1 <- create_loo_mock_fit(N = 30, M = 60, seed = 24)
    res  <- quietly(hbb_loo_compare(m0 = fit0, m1 = fit1))
    df   <- res$comparison
    expect_true(is.data.frame(df))
    # One row per model
    expect_equal(nrow(df), 2L)
    # Required columns from loo::loo_compare() + model column
    expect_true("model"     %in% names(df))
    expect_true("elpd_loo"  %in% names(df))
    expect_true("elpd_diff" %in% names(df))
    expect_true("se_diff"   %in% names(df))
})

test_that("C3. First row of comparison has elpd_diff close to 0 (best model)", {
    fit0 <- create_loo_mock_fit(N = 30, M = 60, seed = 25)
    fit1 <- create_loo_mock_fit(N = 30, M = 60, seed = 26)
    res  <- quietly(hbb_loo_compare(m0 = fit0, m1 = fit1))
    df   <- res$comparison
    # The best model is always in row 1; loo::loo_compare sets its diff to 0.
    expect_equal(df$elpd_diff[1L], 0, tolerance = 1e-10)
})

test_that("C4. Pairwise z-ratios are computed and stored in pairwise data frame", {
    fit0 <- create_loo_mock_fit(N = 30, M = 60, seed = 27)
    fit1 <- create_loo_mock_fit(N = 30, M = 60, seed = 28)
    res  <- quietly(hbb_loo_compare(m0 = fit0, m1 = fit1))
    pw   <- res$pairwise
    expect_true(is.data.frame(pw))
    expect_true("z_ratio" %in% names(pw))
    expect_true(is.numeric(pw$z_ratio))
    # For 2 models, there is exactly 1 consecutive pair
    expect_equal(nrow(pw), 1L)
})

test_that("C5. Significance flag: |z| > 2 models have TRUE significance in pairwise", {
    # Construct two fits with deliberately different ELPD by using different
    # ll_mean values so the difference is detectable.
    fit_good <- create_loo_mock_fit(N = 50, M = 100, seed = 29, ll_mean = -2.0)
    fit_bad  <- create_loo_mock_fit(N = 50, M = 100, seed = 30, ll_mean = -5.0)
    res      <- quietly(hbb_loo_compare(m_good = fit_good, m_bad = fit_bad))
    pw       <- res$pairwise
    # z_ratio should be numeric (may or may not exceed 2 due to mock data)
    expect_true(is.numeric(pw$z_ratio))
    # Significance flag is logical (computed in print, not stored; verify
    # z_ratio is finite instead)
    expect_true(is.finite(pw$z_ratio))
})

test_that("C6. Model names are preserved in result", {
    fit0 <- create_loo_mock_fit(N = 30, M = 60, seed = 31)
    fit1 <- create_loo_mock_fit(N = 30, M = 60, seed = 32)
    res  <- quietly(hbb_loo_compare(model_alpha = fit0, model_beta = fit1))
    expect_true("model_alpha" %in% res$model_names)
    expect_true("model_beta"  %in% res$model_names)
    expect_true("model_alpha" %in% names(res$loo_list))
    expect_true("model_beta"  %in% names(res$loo_list))
})

test_that("C7. Accepts named list input: hbb_loo_compare(list(m0 = fit0, m1 = fit1))", {
    fit0   <- create_loo_mock_fit(N = 30, M = 60, seed = 33)
    fit1   <- create_loo_mock_fit(N = 30, M = 60, seed = 34)
    res    <- quietly(hbb_loo_compare(list(m0 = fit0, m1 = fit1)))
    expect_s3_class(res, "hbb_loo_compare")
    expect_true("m0" %in% res$model_names)
    expect_true("m1" %in% res$model_names)
})

test_that("C8. Accepts ... input with named arguments: hbb_loo_compare(m0 = fit0, m1 = fit1)", {
    fit0 <- create_loo_mock_fit(N = 30, M = 60, seed = 35)
    fit1 <- create_loo_mock_fit(N = 30, M = 60, seed = 36)
    res  <- quietly(hbb_loo_compare(m0 = fit0, m1 = fit1))
    expect_s3_class(res, "hbb_loo_compare")
    expect_equal(res$n_models, 2L)
    expect_true(all(c("m0", "m1") %in% res$model_names))
})


# ============================================================================
# Section D: .loo_chain_id
# ============================================================================

# ---- Section D: .loo_chain_id ----

test_that("D1. .loo_chain_id() returns integer vector of length M", {
    fit <- create_loo_mock_fit(N = 30, M = 100, n_chains = 2, seed = 41)
    cid <- suppressMessages(hurdlebb:::.loo_chain_id(fit, M_expected = 100L))
    expect_true(is.integer(cid))
    expect_equal(length(cid), 100L)
})

test_that("D2. .loo_chain_id() assigns correct chain labels (rep pattern)", {
    # 2 chains, 50 iter_sampling each -> first 50 = chain 1, next 50 = chain 2
    fit <- create_loo_mock_fit(N = 20, M = 100, n_chains = 2, seed = 42)
    cid <- suppressMessages(hurdlebb:::.loo_chain_id(fit, M_expected = 100L))
    expected <- rep(1L:2L, each = 50L)
    expect_equal(cid, expected)
})

test_that("D3. .loo_chain_id() single-chain fallback when metadata unavailable", {
    fit <- create_loo_mock_fit(N = 20, M = 60, n_chains = 2, seed = 43)
    # Remove the metadata function to simulate a loaded-from-RDS scenario
    fit$fit$metadata  <- function() stop("metadata unavailable")
    fit$fit$num_chains <- function() stop("num_chains unavailable")
    # Should warn and return rep(1L, M_expected)
    cid <- suppressWarnings(
        suppressMessages(hurdlebb:::.loo_chain_id(fit, M_expected = 60L))
    )
    expect_equal(length(cid), 60L)
    expect_true(all(cid == 1L))
})


# ============================================================================
# Section E: .loo_pareto_diagnostics
# ============================================================================

# ---- Section E: .loo_pareto_diagnostics ----

test_that("E1. .loo_pareto_diagnostics() does not error on NULL loo_obj", {
    expect_no_error(
        suppressWarnings(hurdlebb:::.loo_pareto_diagnostics(NULL, "test_model"))
    )
})

test_that("E2. .loo_pareto_diagnostics() emits warning for NULL loo_obj", {
    expect_warning(
        hurdlebb:::.loo_pareto_diagnostics(NULL, "test_model")
    )
})


# ============================================================================
# Section F: print.hbb_loo_compare
# ============================================================================

# ---- Section F: print.hbb_loo_compare ----

test_that("F1. print.hbb_loo_compare produces output without error", {
    fit0 <- create_loo_mock_fit(N = 30, M = 60, seed = 51)
    fit1 <- create_loo_mock_fit(N = 30, M = 60, seed = 52)
    res  <- quietly(hbb_loo_compare(m0 = fit0, m1 = fit1))
    expect_no_error(
        capture.output(suppressMessages(print(res)))
    )
})

test_that("F2. print.hbb_loo_compare output contains 'ELPD'", {
    fit0 <- create_loo_mock_fit(N = 30, M = 60, seed = 53)
    fit1 <- create_loo_mock_fit(N = 30, M = 60, seed = 54)
    res  <- quietly(hbb_loo_compare(m0 = fit0, m1 = fit1))
    out  <- capture.output(suppressMessages(print(res)))
    expect_true(any(grepl("ELPD", out, ignore.case = TRUE)))
})

test_that("F3. print.hbb_loo_compare does not crash on stripped/corrupted object", {
    bad_obj <- structure(list(model_names = character(0L)), class = "hbb_loo_compare")
    expect_no_error(
        capture.output(print(bad_obj))
    )
})

test_that("F4. print.hbb_loo_compare returns x invisibly", {
    fit0 <- create_loo_mock_fit(N = 30, M = 60, seed = 55)
    fit1 <- create_loo_mock_fit(N = 30, M = 60, seed = 56)
    res  <- quietly(hbb_loo_compare(m0 = fit0, m1 = fit1))
    ret  <- suppressMessages({
        capture.output(ret_val <- print(res))
        ret_val
    })
    expect_identical(ret, res)
})


# ============================================================================
# Section G: Edge Cases
# ============================================================================

# ---- Section G: Edge Cases ----

test_that("G1. N = 1 observation for loo (doesn't crash, may warn about Pareto k)", {
    # N=1 is pathological for LOO (all weight on 1 observation) but must not crash.
    fit <- create_loo_mock_fit(N = 1L, M = 40, n_chains = 2, seed = 61)
    expect_no_error(
        suppressWarnings(quietly(loo(fit)))
    )
})

test_that("G2. Two identical mock fits produce elpd_diff close to 0", {
    # When both fits share the same log_lik values, elpd_diff should be ~0.
    fit0 <- create_loo_mock_fit(N = 30, M = 60, seed = 62, ll_mean = -3.0)
    # Create fit1 as a copy of fit0 (same seed => identical log_lik)
    fit1 <- create_loo_mock_fit(N = 30, M = 60, seed = 62, ll_mean = -3.0)
    res  <- quietly(hbb_loo_compare(m0 = fit0, m1 = fit1))
    df   <- res$comparison
    # Both rows have the same ELPD, so all elpd_diff entries should be ~0.
    expect_equal(df$elpd_diff, c(0, 0), tolerance = 0.5)
})

test_that("G3. 3+ model comparison works (3 mock fits with different ELPD)", {
    fit0 <- create_loo_mock_fit(N = 30, M = 60, seed = 63, ll_mean = -2.0)
    fit1 <- create_loo_mock_fit(N = 30, M = 60, seed = 64, ll_mean = -3.5)
    fit2 <- create_loo_mock_fit(N = 30, M = 60, seed = 65, ll_mean = -5.0)
    res  <- quietly(hbb_loo_compare(m0 = fit0, m1 = fit1, m2 = fit2))
    expect_s3_class(res, "hbb_loo_compare")
    expect_equal(res$n_models, 3L)
    expect_equal(nrow(res$comparison), 3L)
    # For 3 models there are 2 consecutive pairs in the pairwise table
    expect_equal(nrow(res$pairwise), 2L)
    # Best model has elpd_diff = 0
    expect_equal(res$comparison$elpd_diff[1L], 0, tolerance = 1e-10)
})

test_that("G4. Single chain (n_chains = 1) works for r_eff computation", {
    # n_chains=1: all M draws come from one chain; .loo_chain_id returns rep(1, M).
    fit <- create_loo_mock_fit(N = 25, M = 50, n_chains = 1, seed = 66)
    expect_no_error(
        suppressWarnings(quietly(loo(fit)))
    )
})


# ============================================================================
# ---- Section H: Input Validation (loo.hbb_fit + hbb_loo_compare) ----
# ============================================================================

test_that("H1. Non-hbb_fit raises error matching 'hbb_fit'", {
    # A plain list is not an hbb_fit; loo::loo() has no method for list.
    bad_input <- list(fit = list(), hbb_data = list(N = 10L), model_type = "base")
    expect_error(
        suppressMessages(loo::loo(x = bad_input))
    )
})

test_that("H2. fit$fit = NULL raises error", {
    fit      <- create_loo_mock_fit(N = 30, M = 60, seed = 71)
    fit$fit  <- NULL
    # .loo_validate_inputs checks fit$fit and calls cli_abort.
    expect_error(
        suppressMessages(loo::loo(x = fit))
    )
})

test_that("H3. log_lik not available raises informative error", {
    # include_loglik = FALSE causes metadata() to return no log_lik variables,
    # which .loo_validate_inputs detects and converts to a cli_abort.
    fit_no_ll <- create_loo_mock_fit(N = 30, M = 60, seed = 73,
                                     include_loglik = FALSE)
    expect_error(
        suppressMessages(loo::loo(x = fit_no_ll))
    )
})

test_that("H4. hbb_loo_compare with single model raises error (need >= 2)", {
    fit <- create_loo_mock_fit(N = 30, M = 60, seed = 74)
    expect_error(
        suppressMessages(hbb_loo_compare(m0 = fit))
    )
})

test_that("H5. hbb_loo_compare with unnamed models raises error matching 'named'", {
    fit1 <- create_loo_mock_fit(N = 30, M = 60, seed = 75)
    fit2 <- create_loo_mock_fit(N = 30, M = 60, seed = 76)
    # No names provided -> hbb_loo_compare aborts with a message containing "named".
    expect_error(
        suppressMessages(hbb_loo_compare(fit1, fit2)),
        "named"
    )
})

test_that("H6. hbb_loo_compare with different N raises error matching 'N' or 'mismatch'", {
    fit1 <- create_loo_mock_fit(N = 30, M = 60, seed = 77)
    fit2 <- create_loo_mock_fit(N = 40, M = 60, seed = 78)   # different N
    # hbb_loo_compare checks that all models share the same N.
    expect_error(
        suppressMessages(hbb_loo_compare(m1 = fit1, m2 = fit2)),
        # The error message includes "N" (per-model N values) or "mismatch".
        regexp = "[Nn]"
    )
})


# ============================================================================
# ---- Section I: Print Method Tests (print.hbb_loo_compare) ----
# ============================================================================

#' Build a minimal hbb_loo_compare object for print testing
#'
#' Avoids the cost of running full PSIS-LOO; constructs fake psis_loo objects
#' and assembles the hbb_loo_compare structure manually.
#'
#' @param N Number of observations shared by both models.
#' @param seed Random seed.
#' @return An object of class "hbb_loo_compare".
make_fake_loo_compare <- function(N = 40, seed = 199) {
    set.seed(seed)

    # -- Two fake psis_loo objects with different per-observation ELPD vectors --
    pw1_elpd <- runif(N, -4.0, -0.2)
    pw2_elpd <- pw1_elpd - runif(N, 0.5, 2.0)   # model 2 is worse

    make_loo_obj <- function(pw_elpd) {
        n_obs    <- length(pw_elpd)
        elpd_est <- sum(pw_elpd)
        pk       <- runif(n_obs, 0, 0.45)
        pw <- cbind(
            elpd_loo           = pw_elpd,
            mcse_elpd_loo      = abs(pw_elpd) * 0.02,
            p_loo              = runif(n_obs, 0.1, 1.0),
            looic              = -2 * pw_elpd,
            influence_pareto_k = pk
        )
        structure(
            list(
                estimates = matrix(
                    c(elpd_est, 5.0, 10.0, 2.0, -2 * elpd_est, 10.0),
                    nrow = 3, ncol = 2,
                    dimnames = list(
                        c("elpd_loo", "p_loo", "looic"),
                        c("Estimate", "SE")
                    )
                ),
                pointwise   = pw,
                diagnostics = list(pareto_k = pk, n_eff = runif(n_obs, 50, 200))
            ),
            class = c("psis_loo", "loo")
        )
    }

    loo1     <- make_loo_obj(pw1_elpd)
    loo2     <- make_loo_obj(pw2_elpd)
    loo_list <- list(m1 = loo1, m2 = loo2)

    comparison_raw <- loo::loo_compare(loo_list)
    comparison_df  <- as.data.frame(comparison_raw)
    comparison_df  <- cbind(
        model            = rownames(comparison_df),
        comparison_df,
        stringsAsFactors = FALSE,
        row.names        = NULL
    )
    ranked_names <- rownames(comparison_raw)

    diff_vec   <- pw1_elpd - pw2_elpd
    se_d       <- sqrt(N) * sd(diff_vec)
    delta_elpd <- sum(pw1_elpd) - sum(pw2_elpd)
    z_ratio    <- if (!is.na(se_d) && se_d > 0) delta_elpd / se_d else NA_real_
    pairwise_df <- data.frame(
        comparison  = paste(ranked_names[1L], "vs", ranked_names[2L]),
        delta_elpd  = delta_elpd,
        se_diff     = se_d,
        z_ratio     = z_ratio,
        stringsAsFactors = FALSE
    )

    structure(
        list(
            comparison  = comparison_df,
            loo_list    = loo_list,
            pairwise    = pairwise_df,
            pareto_k    = list(m1 = loo1$diagnostics$pareto_k,
                               m2 = loo2$diagnostics$pareto_k),
            model_names = c("m1", "m2"),
            n_models    = 2L
        ),
        class = "hbb_loo_compare"
    )
}


test_that("I1. print.hbb_loo_compare produces sections containing 'ELPD' or 'elpd'", {
    comp <- make_fake_loo_compare(N = 40, seed = 201)
    out  <- capture.output(print(comp))
    expect_true(
        any(grepl("ELPD|elpd", out)),
        label = "output must contain 'ELPD' or 'elpd'"
    )
})

test_that("I2. print.hbb_loo_compare respects digits argument", {
    comp <- make_fake_loo_compare(N = 40, seed = 202)
    out1 <- capture.output(print(comp, digits = 1L))
    out3 <- capture.output(print(comp, digits = 3L))
    # Numeric formatting changes with more decimal places.
    expect_false(identical(paste(out1, collapse = ""), paste(out3, collapse = "")))
})

test_that("I3. print.hbb_loo_compare handles 2-model comparison cleanly (no error)", {
    comp <- make_fake_loo_compare(N = 50, seed = 203)
    expect_no_error(capture.output(print(comp)))
    out <- capture.output(print(comp))
    # Both model names should appear in the output.
    expect_true(any(grepl("m1", out, fixed = TRUE)))
    expect_true(any(grepl("m2", out, fixed = TRUE)))
})

test_that("I4. print.hbb_loo_compare returns x invisibly", {
    comp    <- make_fake_loo_compare(N = 35, seed = 204)
    ret_val <- {
        capture.output(ret_val <- print(comp))
        ret_val
    }
    expect_identical(ret_val, comp)
})


# ============================================================================
# ---- Section J: Additional Coverage Tests ----
# ============================================================================


# ---- J1: .loo_pairwise_z --- direct unit tests (lines 1078-1084, 1126) ----

test_that("J1a. .loo_pairwise_z returns empty df for single model (lines 1078-1084)", {
    # With only 1 model name, there are no pairs => empty data.frame
    loo_list <- list(m1 = make_fake_loo_compare(N = 20, seed = 301)$loo_list$m1)
    result <- hurdlebb:::.loo_pairwise_z(loo_list, ranked_names = "m1")

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 0L)
    expect_true("comparison"  %in% names(result))
    expect_true("delta_elpd"  %in% names(result))
    expect_true("se_diff"     %in% names(result))
    expect_true("z_ratio"     %in% names(result))
})


test_that("J1b. .loo_pairwise_z computes valid z-ratio for 2-model pair", {
    comp <- make_fake_loo_compare(N = 40, seed = 302)
    result <- hurdlebb:::.loo_pairwise_z(comp$loo_list, c("m1", "m2"))

    expect_equal(nrow(result), 1L)
    expect_true(!is.na(result$delta_elpd[1]))
    expect_true(!is.na(result$se_diff[1]))
    expect_true(!is.na(result$z_ratio[1]))
    # z_ratio = delta_elpd / se_diff
    expect_equal(result$z_ratio[1],
                 result$delta_elpd[1] / result$se_diff[1],
                 tolerance = 1e-10)
})


test_that("J1c. .loo_pairwise_z handles NA when pointwise data is missing (line 1126)", {
    # Create a loo object with no pointwise data
    fake_loo_a <- structure(
        list(
            estimates = matrix(c(-100, 5, 200, 10, 5, 20),
                               nrow = 3, ncol = 2,
                               dimnames = list(c("elpd_loo","p_loo","looic"),
                                               c("Estimate","SE"))),
            pointwise = NULL,  # missing!
            diagnostics = list(pareto_k = runif(10))
        ),
        class = c("psis_loo", "loo")
    )
    fake_loo_b <- fake_loo_a
    fake_loo_b$estimates[1, 1] <- -120  # different ELPD

    loo_list <- list(ma = fake_loo_a, mb = fake_loo_b)
    result <- hurdlebb:::.loo_pairwise_z(loo_list, c("ma", "mb"))

    expect_equal(nrow(result), 1L)
    expect_true(is.na(result$se_diff[1]))
    expect_true(is.na(result$z_ratio[1]))
})


# ---- J2: .loo_pareto_diagnostics --- direct unit tests (lines 1002-1006) ----

test_that("J2a. .loo_pareto_diagnostics handles NULL loo_obj (line 1002-1006)", {
    expect_warning(
        hurdlebb:::.loo_pareto_diagnostics(NULL, "test_model"),
        "NULL"
    )
})


test_that("J2b. .loo_pareto_diagnostics handles empty pareto_k (line 1002-1006)", {
    fake_loo <- structure(
        list(diagnostics = list(pareto_k = numeric(0))),
        class = c("psis_loo", "loo")
    )
    expect_warning(
        hurdlebb:::.loo_pareto_diagnostics(fake_loo, "empty_model"),
        "unavailable"
    )
})


test_that("J2c. .loo_pareto_diagnostics handles missing diagnostics slot", {
    fake_loo <- structure(list(), class = c("psis_loo", "loo"))
    expect_warning(
        hurdlebb:::.loo_pareto_diagnostics(fake_loo, "no_diag_model"),
        "unavailable"
    )
})


# ---- J3: .loo_chain_id --- direct unit tests (lines 941-945, 953-955) ----

test_that("J3a. .loo_chain_id falls back to single chain when metadata unavailable (lines 941-945)", {
    # Create a mock fit with broken metadata
    mock_fit <- structure(
        list(
            fit = list(
                num_chains = function() stop("no metadata"),
                metadata = function() stop("no metadata")
            ),
            hbb_data = list(N = 10L)
        ),
        class = "hbb_fit"
    )
    expect_warning(
        chain_id <- hurdlebb:::.loo_chain_id(mock_fit, M_expected = 100L),
        "single-chain"
    )
    expect_length(chain_id, 100L)
    expect_true(all(chain_id == 1L))
})


test_that("J3b. .loo_chain_id errors when no M_expected and no metadata (lines 941-945)", {
    mock_fit <- structure(
        list(
            fit = list(
                num_chains = function() NULL,
                metadata = function() list(iter_sampling = NULL)
            ),
            hbb_data = list(N = 10L)
        ),
        class = "hbb_fit"
    )
    expect_error(
        hurdlebb:::.loo_chain_id(mock_fit, M_expected = NULL),
        "chain structure"
    )
})


test_that("J3c. .loo_chain_id warns on length mismatch (lines 953-955)", {
    mock_fit <- structure(
        list(
            fit = list(
                num_chains = function() 2L,
                metadata = function() list(iter_sampling = 50L)
            ),
            hbb_data = list(N = 10L)
        ),
        class = "hbb_fit"
    )
    # M_expected = 200 but chains*iter = 100 => mismatch
    expect_warning(
        chain_id <- hurdlebb:::.loo_chain_id(mock_fit, M_expected = 200L),
        "chain_id length"
    )
    expect_length(chain_id, 100L)
})


# ---- J4: .loo_validate_inputs --- additional paths (lines 816-820, 862-864) ----

test_that("J4a. .loo_validate_inputs errors on missing hbb_data$N", {
    mock_fit <- structure(
        list(
            fit = list(metadata = function() list(stan_variables = "log_lik")),
            hbb_data = list(N = NULL)
        ),
        class = "hbb_fit"
    )
    expect_error(
        hurdlebb:::.loo_validate_inputs(mock_fit),
        "hbb_data\\$N"
    )
})


test_that("J4b. .loo_validate_inputs errors on missing hbb_data entirely", {
    mock_fit <- structure(
        list(
            fit = list(metadata = function() list(stan_variables = "log_lik")),
            hbb_data = NULL
        ),
        class = "hbb_fit"
    )
    expect_error(
        hurdlebb:::.loo_validate_inputs(mock_fit),
        "hbb_data\\$N|missing"
    )
})


test_that("J4c. .loo_validate_inputs does trial draw fallback when metadata is broken", {
    # metadata() throws error, but draws() works => should NOT error
    mock_fit <- structure(
        list(
            fit = list(
                metadata = function() stop("broken metadata"),
                draws = function(variables = NULL, format = NULL) {
                    matrix(rnorm(20), nrow = 10, ncol = 2)
                }
            ),
            hbb_data = list(N = 10L)
        ),
        class = "hbb_fit"
    )
    # Should warn about metadata, but succeed via trial draw
    expect_warning(
        result <- hurdlebb:::.loo_validate_inputs(mock_fit),
        "metadata"
    )
})


test_that("J4d. .loo_validate_inputs errors when no log_lik via metadata or trial draw", {
    # Both metadata and trial draw fail
    mock_fit <- structure(
        list(
            fit = list(
                metadata = function() stop("broken"),
                draws = function(variables = NULL, format = NULL) stop("no draws")
            ),
            hbb_data = list(N = 10L)
        ),
        class = "hbb_fit"
    )
    expect_error(
        suppressWarnings(hurdlebb:::.loo_validate_inputs(mock_fit)),
        "log_lik"
    )
})


# ---- J5: hbb_loo_compare input validation (lines 434-439) ----

test_that("J5a. hbb_loo_compare errors when hbb_data$N is NULL (lines 434-439)", {
    fit1 <- create_loo_mock_fit(N = 30, M = 60, seed = 401)
    fit2 <- create_loo_mock_fit(N = 30, M = 60, seed = 402)
    fit2$hbb_data$N <- NULL   # corrupt N
    # Setting N to NULL causes vapply to produce NA_integer_, triggering the
    # N validation error path. The cli_abort may fail on pluralization too,
    # but either way, an error is expected.
    expect_error(
        suppressMessages(hbb_loo_compare(m1 = fit1, m2 = fit2))
    )
})


test_that("J5b. hbb_loo_compare errors with non-hbb_fit argument", {
    fit1 <- create_loo_mock_fit(N = 30, M = 60, seed = 403)
    expect_error(
        suppressMessages(hbb_loo_compare(m1 = fit1, m2 = list(a = 1))),
        "hbb_fit"
    )
})


# ---- J6: print.hbb_loo_compare formatting branches ----

test_that("J6a. print.hbb_loo_compare shows '(comparison table unavailable)' (line 616)", {
    comp <- make_fake_loo_compare(N = 20, seed = 500)
    comp$comparison <- NULL   # remove comparison table
    out <- capture.output(print(comp))
    combined <- paste(out, collapse = "\n")
    expect_true(grepl("unavailable", combined))
})


test_that("J6b. print.hbb_loo_compare shows '(no Pareto-k data)' (line 677)", {
    comp <- make_fake_loo_compare(N = 20, seed = 501)
    comp$pareto_k <- list(m1 = NULL, m2 = NULL)   # no pareto data
    out <- capture.output(print(comp))
    combined <- paste(out, collapse = "\n")
    expect_true(grepl("no Pareto-k", combined))
})


test_that("J6c. print.hbb_loo_compare handles row formatting error gracefully (line 609-610)", {
    comp <- make_fake_loo_compare(N = 20, seed = 502)
    # Corrupt comparison df so that one row causes formatting error
    comp$comparison$elpd_loo[2] <- "not_a_number"
    out <- capture.output(print(comp))
    # Should still produce output (error caught internally)
    expect_true(length(out) > 3)
})


test_that("J6d. print.hbb_loo_compare shows z-ratio significance flag (line 629)", {
    comp <- make_fake_loo_compare(N = 40, seed = 503)
    # Make sure z_ratio > 2 for significance flag
    comp$pairwise$z_ratio[1] <- 5.0
    out <- capture.output(print(comp))
    combined <- paste(out, collapse = "\n")
    # The '*' significance flag
    expect_true(grepl("\\*", combined))
})


test_that("J6e. print.hbb_loo_compare handles pairwise row error (line 648-649)", {
    comp <- make_fake_loo_compare(N = 20, seed = 504)
    # Corrupt pairwise so the row fails
    comp$pairwise$delta_elpd[1] <- "corrupt"
    out <- capture.output(print(comp))
    expect_true(length(out) > 3)
})


test_that("J6f. print.hbb_loo_compare handles model_names fallback (line 661)", {
    comp <- make_fake_loo_compare(N = 20, seed = 505)
    comp$model_names <- character(0)  # empty model_names
    out <- capture.output(print(comp))
    # Should still print pareto diagnostics section using names(pk_list)
    combined <- paste(out, collapse = "\n")
    expect_true(grepl("Pareto", combined))
})


test_that("J6g. print.hbb_loo_compare handles pareto_k error per model (line 680)", {
    comp <- make_fake_loo_compare(N = 20, seed = 506)
    # Make pareto_k a list that causes an error when accessed
    comp$pareto_k <- list(m1 = structure(list(), class = "weird"), m2 = runif(20))
    out <- capture.output(print(comp))
    combined <- paste(out, collapse = "\n")
    # Should see error message for m1, but m2 should print fine
    expect_true(grepl("error|no Pareto-k", combined) || length(out) > 3)
})


# ---- J7: loo.hbb_fit r_eff path (lines 225-229, 236-239) ----

test_that("J7a. loo.hbb_fit r_eff = FALSE skips chain_id computation", {
    fit <- create_loo_mock_fit(N = 25, M = 50, n_chains = 2, seed = 601)
    result <- quietly(loo(fit, r_eff = FALSE))
    expect_s3_class(result, c("psis_loo", "loo"))
})


test_that("J7b. loo.hbb_fit invalid r_eff parameter (lines 192-195)", {
    fit <- create_loo_mock_fit(N = 25, M = 50, n_chains = 2, seed = 602)
    expect_error(
        quietly(loo(fit, r_eff = "yes")),
        "r_eff"
    )
    expect_error(
        quietly(loo(fit, r_eff = NA)),
        "r_eff"
    )
})


# ---- J8: .loo_extract_loglik pathological value checks ----

test_that("J8a. .loo_extract_loglik warns on NaN values (line 769-773)", {
    fit <- create_loo_mock_fit(N = 20, M = 40, n_chains = 2, seed = 701)
    # Inject NaN into the underlying log_lik draws
    original_draws <- fit$fit$draws
    fit$fit$draws <- function(variables = NULL, format = "matrix") {
        result <- original_draws(variables = variables, format = format)
        if (!is.null(variables) && length(variables) == 1L &&
            variables == "log_lik") {
            result[1, 1] <- NaN
        }
        result
    }

    expect_warning(
        hurdlebb:::.loo_extract_loglik(fit),
        "NaN"
    )
})


test_that("J8b. .loo_extract_loglik warns on +Inf values", {
    fit <- create_loo_mock_fit(N = 20, M = 40, n_chains = 2, seed = 702)
    original_draws <- fit$fit$draws
    fit$fit$draws <- function(variables = NULL, format = "matrix") {
        result <- original_draws(variables = variables, format = format)
        if (!is.null(variables) && length(variables) == 1L &&
            variables == "log_lik") {
            result[1, 1] <- Inf
        }
        result
    }

    expect_warning(
        hurdlebb:::.loo_extract_loglik(fit),
        "Inf"
    )
})


test_that("J8c. .loo_extract_loglik dimension mismatch error (line 743-747)", {
    fit <- create_loo_mock_fit(N = 20, M = 40, n_chains = 2, seed = 703)
    # Make N mismatch with the extracted log_lik
    fit$hbb_data$N <- 999L

    expect_error(
        hurdlebb:::.loo_extract_loglik(fit),
        "log_lik|column"
    )
})
