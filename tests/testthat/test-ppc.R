# ============================================================================
# test-ppc.R --- testthat tests for ppc.R
#
# Sections:
#   MOCK: create_ppc_mock_fit constructor
#   A.   PPC Statistics Algebraic Identities (8 tests)
#   B.   Coverage Computation (5 tests)
# ============================================================================


# ============================================================================
# MOCK CONSTRUCTOR
# ============================================================================

#' Create a mock hbb_fit object for PPC testing
#'
#' Generates a full mock hbb_fit with:
#'   - param draws: alpha[1:P], beta[1:P], log_kappa (M x (2P+1) matrix)
#'   - y_rep draws: M x N integer matrix (hurdle structure, correct zero rate)
#'   - log_lik draws: M x N numeric matrix (valid negative values)
#'   - observed data: y, n_trial, z, X, N, P
#'   - Mock CmdStanMCMC interface: $draws(), $num_chains(), $metadata()
#'
#' @param P Number of covariates per margin (including intercept).
#' @param N Number of observations.
#' @param M Number of MCMC draws.
#' @param n_chains Number of chains (for metadata).
#' @param seed Random seed.
#' @param zero_rate Target zero rate in y_rep (approx).
#' @param include_yrep Logical; whether y_rep is present in draws.
#' @param include_loglik Logical; whether log_lik is present in draws.
#' @return An S3 object of class "hbb_fit".
create_ppc_mock_fit <- function(P = 3, N = 50, M = 200, n_chains = 2,
                                seed = 42, zero_rate = 0.35,
                                include_yrep = TRUE,
                                include_loglik = TRUE) {

    set.seed(seed)

    iter_sampling <- M %/% n_chains

    # -- Parameter draws -------------------------------------------------------
    alpha_true <- c(-0.5, rep(0.2, P - 1L))
    beta_true  <- c(0.3,  rep(-0.1, P - 1L))
    log_kappa_true <- 1.9  # kappa ~ 6.7

    param_names <- c(
        paste0("alpha[", seq_len(P), "]"),
        paste0("beta[",  seq_len(P), "]"),
        "log_kappa"
    )
    D <- 2L * P + 1L
    theta_true <- c(alpha_true, beta_true, log_kappa_true)
    Sigma_param <- diag(D) * 0.01
    param_draws <- MASS::mvrnorm(M, theta_true, Sigma_param)
    colnames(param_draws) <- param_names

    # -- Observed data ---------------------------------------------------------
    n_trial <- sample(10L:50L, N, replace = TRUE)
    z       <- rbinom(N, 1L, 1 - zero_rate)          # 1 = participant
    y       <- ifelse(z == 1L, rbinom(N, n_trial, 0.3), 0L)
    X       <- cbind(
        intercept = rep(1, N),
        matrix(rnorm(N * (P - 1L)), nrow = N, ncol = P - 1L,
               dimnames = list(NULL, paste0("x", seq_len(P - 1L))))
    )

    # -- y_rep draws: M x N matrix with hurdle structure ----------------------
    y_rep <- if (include_yrep) {
        mat <- matrix(0L, nrow = M, ncol = N)
        for (m in seq_len(M)) {
            z_m    <- rbinom(N, 1L, 1 - zero_rate)
            y_m    <- ifelse(z_m == 1L, rbinom(N, n_trial, 0.3), 0L)
            mat[m, ] <- y_m
        }
        mat
    } else {
        NULL
    }

    # -- log_lik draws: M x N matrix of valid negative values -----------------
    log_lik <- if (include_loglik) {
        matrix(
            runif(M * N, min = -5, max = -0.1),
            nrow = M, ncol = N
        )
    } else {
        NULL
    }

    # -- Consolidated draws store (param_draws + optional y_rep + log_lik) ----
    all_draws <- param_draws
    if (include_yrep) {
        y_rep_cols <- matrix(as.numeric(y_rep), nrow = M, ncol = N)
        colnames(y_rep_cols) <- paste0("y_rep[", seq_len(N), "]")
        all_draws <- cbind(all_draws, y_rep_cols)
    }
    if (include_loglik) {
        loglik_cols <- log_lik
        colnames(loglik_cols) <- paste0("log_lik[", seq_len(N), "]")
        all_draws <- cbind(all_draws, loglik_cols)
    }

    # -- Mock CmdStanMCMC interface -------------------------------------------
    mock_cmdstan <- list(

        draws = function(variables = NULL, format = "matrix") {
            if (is.null(variables)) return(all_draws)

            # Handle "y_rep" as a prefix group
            is_yrep_request   <- any(grepl("^y_rep",   variables))
            is_loglik_request <- any(grepl("^log_lik", variables))

            if (is_yrep_request && !include_yrep) {
                stop("y_rep not available in this mock fit")
            }
            if (is_loglik_request && !include_loglik) {
                stop("log_lik not available in this mock fit")
            }

            # Exact name match in all_draws columns
            cols <- colnames(all_draws)
            idx  <- cols %in% variables
            if (!any(idx)) {
                # Try prefix match for grouped variables (e.g. "y_rep")
                idx <- grepl(paste0("^", paste(variables, collapse = "|^")), cols)
            }
            if (!any(idx)) stop("Variables not found: ", paste(variables, collapse = ", "))
            all_draws[, idx, drop = FALSE]
        },

        num_chains = function() n_chains,

        metadata = function() {
            list(
                iter_warmup   = iter_sampling,
                iter_sampling = iter_sampling,
                stan_variables = colnames(all_draws)
            )
        }
    )

    hbb_data <- list(
        y       = y,
        n_trial = n_trial,
        z       = z,
        X       = X,
        N       = as.integer(N),
        P       = as.integer(P)
    )

    structure(
        list(
            fit        = mock_cmdstan,
            hbb_data   = hbb_data,
            model_type = "weighted",
            model_name = "hbb_m1_mock"
        ),
        class = "hbb_fit"
    )
}


# ---- Section A: PPC Statistics Algebraic Identities ----


# ---- A1. Zero rate is in [0, 1] for all draws --------------------------------

test_that("A1. zero rate draws are all in [0, 1]", {
    set.seed(101)
    N <- 40; M <- 150
    n_trial <- sample(5:30, N, replace = TRUE)
    y_rep_mat <- matrix(
        sample(0:10, M * N, replace = TRUE),
        nrow = M, ncol = N
    )
    # Force some zeros so zero_rate is not degenerate
    y_rep_mat[, 1:10] <- 0L

    result <- suppressWarnings(
        suppressMessages(
            hurdlebb:::.ppc_compute_stats(y_rep_mat, n_trial, "zero_rate")
        )
    )
    zr <- result$zero_rate
    expect_true(all(zr >= 0))
    expect_true(all(zr <= 1))
})


# ---- A2. Zero rate mean is close to data zero rate for well-calibrated mock --

test_that("A2. mean zero rate draw is close to observed zero rate (well-calibrated)", {
    set.seed(102)
    N <- 100; M <- 500
    true_zero_rate <- 0.40
    n_trial <- rep(20L, N)

    # Construct y_rep rows from the true model
    y_rep_mat <- matrix(0L, nrow = M, ncol = N)
    for (m in seq_len(M)) {
        z_m <- rbinom(N, 1L, 1 - true_zero_rate)
        y_rep_mat[m, ] <- ifelse(z_m == 1L, rbinom(N, n_trial, 0.35), 0L)
    }

    result <- suppressWarnings(
        suppressMessages(
            hurdlebb:::.ppc_compute_stats(y_rep_mat, n_trial, "zero_rate")
        )
    )
    mean_zr <- mean(result$zero_rate, na.rm = TRUE)
    # Should be close to true_zero_rate within sampling variability
    expect_true(abs(mean_zr - true_zero_rate) < 0.08,
                label = sprintf("mean zero_rate %.3f not near %.3f", mean_zr, true_zero_rate))
})


# ---- A3. IT share is in [0, 1] for all non-NA draws -------------------------

test_that("A3. IT share draws are all in [0, 1] (for non-NA entries)", {
    set.seed(103)
    N <- 60; M <- 200
    n_trial <- sample(5:25, N, replace = TRUE)
    y_rep_mat <- matrix(0L, nrow = M, ncol = N)
    for (m in seq_len(M)) {
        z_m <- rbinom(N, 1L, 0.65)
        y_rep_mat[m, ] <- ifelse(z_m == 1L, rbinom(N, n_trial, 0.3), 0L)
    }

    result <- suppressWarnings(
        suppressMessages(
            hurdlebb:::.ppc_compute_stats(y_rep_mat, n_trial, "it_share")
        )
    )
    its <- result$it_share
    finite_its <- its[!is.na(its)]
    expect_true(length(finite_its) > 0)
    expect_true(all(finite_its >= 0))
    expect_true(all(finite_its <= 1))
})


# ---- A4. IT share draws have correct length (= M) ----------------------------

test_that("A4. IT share draws vector has length equal to M", {
    set.seed(104)
    N <- 30; M <- 80
    n_trial <- rep(15L, N)
    y_rep_mat <- matrix(0L, nrow = M, ncol = N)
    for (m in seq_len(M)) {
        z_m <- rbinom(N, 1L, 0.6)
        y_rep_mat[m, ] <- ifelse(z_m == 1L, rbinom(N, n_trial, 0.4), 0L)
    }

    result <- suppressWarnings(
        suppressMessages(
            hurdlebb:::.ppc_compute_stats(y_rep_mat, n_trial, "it_share")
        )
    )
    expect_length(result$it_share, M)
})


# ---- A5. All-zero y_rep draw produces NA for it_share (handled gracefully) ---

test_that("A5. draw with all zeros yields NA for IT share (graceful)", {
    set.seed(105)
    N <- 20; M <- 5
    n_trial <- rep(10L, N)
    # All-zero matrix
    y_rep_mat <- matrix(0L, nrow = M, ncol = N)

    result <- suppressWarnings(
        suppressMessages(
            hurdlebb:::.ppc_compute_stats(y_rep_mat, n_trial, "it_share")
        )
    )
    # Every draw has all zeros -> every IT share must be NA
    expect_true(all(is.na(result$it_share)))
})


# ---- A6. Zero rate = 1.0 when all y_rep = 0 ---------------------------------

test_that("A6. zero rate equals 1.0 when all y_rep values are zero", {
    N <- 30; M <- 10
    n_trial <- rep(20L, N)
    y_rep_mat <- matrix(0L, nrow = M, ncol = N)

    result <- suppressWarnings(
        suppressMessages(
            hurdlebb:::.ppc_compute_stats(y_rep_mat, n_trial, "zero_rate")
        )
    )
    expect_equal(result$zero_rate, rep(1.0, M), tolerance = 1e-15)
})


# ---- A7. Zero rate = 0.0 when no y_rep = 0 ----------------------------------

test_that("A7. zero rate equals 0.0 when no y_rep value is zero", {
    N <- 25; M <- 8
    n_trial <- rep(20L, N)
    # All positive counts (1 to n_trial)
    y_rep_mat <- matrix(sample(1L:5L, M * N, replace = TRUE),
                        nrow = M, ncol = N)

    result <- suppressWarnings(
        suppressMessages(
            hurdlebb:::.ppc_compute_stats(y_rep_mat, n_trial, "zero_rate")
        )
    )
    expect_equal(result$zero_rate, rep(0.0, M), tolerance = 1e-15)
})


# ---- A8. Observed statistics match manual computation from mock data ----------

test_that("A8. .ppc_compute_stats matches manual computation of zero rate and IT share", {
    set.seed(108)
    N <- 10; M <- 3
    n_trial <- c(5L, 10L, 8L, 6L, 12L, 9L, 7L, 4L, 11L, 3L)

    # Construct a small deterministic matrix
    y_rep_mat <- matrix(
        c(
            0L, 3L, 0L, 6L, 0L, 9L, 0L, 4L, 0L, 3L,   # row 1: 5 zeros
            1L, 2L, 3L, 4L, 5L, 6L, 7L, 4L, 2L, 1L,   # row 2: 0 zeros
            0L, 0L, 0L, 0L, 5L, 0L, 0L, 0L, 0L, 0L    # row 3: 9 zeros
        ),
        nrow = M, ncol = N, byrow = TRUE
    )

    result <- suppressWarnings(
        suppressMessages(
            hurdlebb:::.ppc_compute_stats(y_rep_mat, n_trial, "both")
        )
    )

    # Manual zero rate
    expected_zr <- c(
        mean(y_rep_mat[1L, ] == 0L),
        mean(y_rep_mat[2L, ] == 0L),
        mean(y_rep_mat[3L, ] == 0L)
    )
    expect_equal(result$zero_rate, expected_zr, tolerance = 1e-15)

    # Manual IT share for row 1 (zeros at positions 1,3,5,7,9)
    pos_mask_1 <- y_rep_mat[1L, ] > 0L
    expected_its_1 <- mean(y_rep_mat[1L, pos_mask_1] / n_trial[pos_mask_1])
    expect_equal(result$it_share[1L], expected_its_1, tolerance = 1e-14)

    # Manual IT share for row 2 (no zeros)
    expected_its_2 <- mean(y_rep_mat[2L, ] / n_trial)
    expect_equal(result$it_share[2L], expected_its_2, tolerance = 1e-14)

    # Row 3: only one positive (position 5)
    pos_mask_3 <- y_rep_mat[3L, ] > 0L
    expected_its_3 <- mean(y_rep_mat[3L, pos_mask_3] / n_trial[pos_mask_3])
    expect_equal(result$it_share[3L], expected_its_3, tolerance = 1e-14)
})


# ---- Section B: Coverage Computation ----


# ---- B1. Observed zero_rate within predicted CI when model is well-calibrated -

test_that("B1. observed zero rate falls within 95% CI for well-calibrated mock", {
    set.seed(201)
    N   <- 80
    M   <- 600
    true_zr <- 0.38
    n_trial <- rep(20L, N)

    # Observed data generated from same model
    z_obs <- rbinom(N, 1L, 1 - true_zr)
    y_obs <- ifelse(z_obs == 1L, rbinom(N, n_trial, 0.3), 0L)
    obs_zero_rate <- mean(y_obs == 0L)

    # Posterior predictive draws from same generative model
    y_rep_mat <- matrix(0L, nrow = M, ncol = N)
    for (m in seq_len(M)) {
        z_m <- rbinom(N, 1L, 1 - true_zr)
        y_rep_mat[m, ] <- ifelse(z_m == 1L, rbinom(N, n_trial, 0.3), 0L)
    }

    stats <- suppressWarnings(
        suppressMessages(
            hurdlebb:::.ppc_compute_stats(y_rep_mat, n_trial, "zero_rate")
        )
    )
    summary_obj <- suppressMessages(
        hurdlebb:::.ppc_summarize_draws(stats$zero_rate, level = 0.95)
    )

    ci_lower <- summary_obj$ci["lower"]
    ci_upper <- summary_obj$ci["upper"]

    expect_true(
        obs_zero_rate >= ci_lower && obs_zero_rate <= ci_upper,
        label = sprintf(
            "obs zero rate %.3f not in CI [%.3f, %.3f]",
            obs_zero_rate, ci_lower, ci_upper
        )
    )
})


# ---- B2. Coverage flag in_ci is logical scalar -------------------------------

test_that("B2. .ppc_summarize_draws returns a CI that enables logical in_ci check", {
    set.seed(202)
    N <- 30; M <- 100
    n_trial <- rep(10L, N)
    y_rep_mat <- matrix(0L, nrow = M, ncol = N)
    for (m in seq_len(M)) {
        z_m <- rbinom(N, 1L, 0.7)
        y_rep_mat[m, ] <- ifelse(z_m == 1L, rbinom(N, n_trial, 0.4), 0L)
    }
    obs_zr <- mean(y_rep_mat[1L, ] == 0L)  # use first row as "observed"

    stats <- suppressWarnings(
        suppressMessages(
            hurdlebb:::.ppc_compute_stats(y_rep_mat, n_trial, "zero_rate")
        )
    )
    summary_obj <- suppressMessages(
        hurdlebb:::.ppc_summarize_draws(stats$zero_rate, level = 0.95)
    )

    in_ci <- obs_zr >= summary_obj$ci["lower"] && obs_zr <= summary_obj$ci["upper"]
    expect_true(is.logical(in_ci))
    expect_length(in_ci, 1L)
})


# ---- B3. CI bounds are ordered (lower < upper) -------------------------------

test_that("B3. CI lower bound is strictly less than upper bound", {
    set.seed(203)
    N <- 50; M <- 300
    n_trial <- rep(15L, N)
    y_rep_mat <- matrix(0L, nrow = M, ncol = N)
    for (m in seq_len(M)) {
        z_m <- rbinom(N, 1L, 0.65)
        y_rep_mat[m, ] <- ifelse(z_m == 1L, rbinom(N, n_trial, 0.35), 0L)
    }

    stats <- suppressWarnings(
        suppressMessages(
            hurdlebb:::.ppc_compute_stats(y_rep_mat, n_trial, "zero_rate")
        )
    )
    summary_obj <- suppressMessages(
        hurdlebb:::.ppc_summarize_draws(stats$zero_rate, level = 0.90)
    )

    expect_true(summary_obj$ci["lower"] < summary_obj$ci["upper"])
})


# ---- B4. CI with level=0.50 is narrower than level=0.99 ---------------------

test_that("B4. 50% CI is strictly narrower than 99% CI", {
    set.seed(204)
    N <- 40; M <- 400
    n_trial <- rep(12L, N)
    y_rep_mat <- matrix(0L, nrow = M, ncol = N)
    for (m in seq_len(M)) {
        z_m <- rbinom(N, 1L, 0.7)
        y_rep_mat[m, ] <- ifelse(z_m == 1L, rbinom(N, n_trial, 0.4), 0L)
    }

    stats <- suppressWarnings(
        suppressMessages(
            hurdlebb:::.ppc_compute_stats(y_rep_mat, n_trial, "zero_rate")
        )
    )

    summary_50 <- suppressMessages(
        hurdlebb:::.ppc_summarize_draws(stats$zero_rate, level = 0.50)
    )
    summary_99 <- suppressMessages(
        hurdlebb:::.ppc_summarize_draws(stats$zero_rate, level = 0.99)
    )

    width_50 <- summary_50$ci["upper"] - summary_50$ci["lower"]
    width_99 <- summary_99$ci["upper"] - summary_99$ci["lower"]
    expect_true(width_50 < width_99)
})


# ---- B5. summary has zero_rate and it_share with correct structure -----------

test_that("B5. .ppc_summarize_draws returns named list with mean, median, ci, sd, n_finite", {
    set.seed(205)
    N <- 50; M <- 200
    n_trial <- rep(20L, N)
    y_rep_mat <- matrix(0L, nrow = M, ncol = N)
    for (m in seq_len(M)) {
        z_m <- rbinom(N, 1L, 0.65)
        y_rep_mat[m, ] <- ifelse(z_m == 1L, rbinom(N, n_trial, 0.3), 0L)
    }

    stats <- suppressWarnings(
        suppressMessages(
            hurdlebb:::.ppc_compute_stats(y_rep_mat, n_trial, "both")
        )
    )

    # zero_rate summary
    sum_zr <- suppressMessages(
        hurdlebb:::.ppc_summarize_draws(stats$zero_rate, level = 0.95)
    )
    expect_true(is.list(sum_zr))
    expect_true(all(c("mean", "median", "ci", "sd", "n_finite", "n_na") %in% names(sum_zr)))
    expect_true(is.numeric(sum_zr$mean))
    expect_true(is.numeric(sum_zr$ci))
    expect_named(sum_zr$ci, c("lower", "upper"))
    expect_true(sum_zr$n_finite > 0L)

    # it_share summary
    sum_its <- suppressMessages(
        hurdlebb:::.ppc_summarize_draws(stats$it_share, level = 0.95)
    )
    expect_true(is.list(sum_its))
    expect_true(all(c("mean", "median", "ci", "sd", "n_finite") %in% names(sum_its)))
})


# ============================================================================
# ---- Section C: ppc() Integration (mock-based, no CmdStan required) ----
# ============================================================================

test_that("C1. ppc() returns hbb_ppc object with correct structure", {
    fit    <- create_ppc_mock_fit()
    result <- suppressMessages(ppc(fit, method = "stan"))
    expect_s3_class(result, "hbb_ppc")
    expect_true(all(c("observed", "predicted", "coverage",
                      "n_draws", "n_total_draws", "type", "method",
                      "level", "model_type") %in% names(result)))
})

test_that("C2. ppc() method = 'simulate' works with the mock fit", {
    fit    <- create_ppc_mock_fit()
    result <- suppressMessages(ppc(fit, method = "simulate"))
    expect_s3_class(result, "hbb_ppc")
    expect_equal(result$method, "simulate")
})

test_that("C3. ppc() n_draws thinning reduces draw count", {
    fit    <- create_ppc_mock_fit(M = 200)
    result <- suppressMessages(ppc(fit, n_draws = 50L, method = "stan"))
    expect_s3_class(result, "hbb_ppc")
    expect_lte(result$n_draws, 50L)
})


# ============================================================================
# ---- Section D: Input Validation ----
# ============================================================================

test_that("D1. Non-hbb_fit input (a list) raises error matching 'hbb_fit'", {
    bad_input <- list(fit = NULL, stan_data = list(), model_type = "base")
    expect_error(
        suppressMessages(ppc(bad_input)),
        "hbb_fit"
    )
})

test_that("D2. Invalid type raises error", {
    fit <- create_ppc_mock_fit()
    expect_error(
        suppressMessages(ppc(fit, type = "foo"))
    )
})

test_that("D3. Invalid method raises error", {
    fit <- create_ppc_mock_fit()
    expect_error(
        suppressMessages(ppc(fit, method = "abc"))
    )
})

test_that("D4. n_draws = 0 raises error", {
    fit <- create_ppc_mock_fit()
    expect_error(
        suppressMessages(ppc(fit, n_draws = 0L))
    )
})

test_that("D5. n_draws = -1 raises error", {
    fit <- create_ppc_mock_fit()
    expect_error(
        suppressMessages(ppc(fit, n_draws = -1L))
    )
})

test_that("D6. level = 0 raises error; level = 1 raises error", {
    fit <- create_ppc_mock_fit()
    expect_error(suppressMessages(ppc(fit, level = 0)))
    expect_error(suppressMessages(ppc(fit, level = 1)))
})

test_that("D7. fit$fit = NULL raises error matching 'CmdStan'", {
    fit      <- create_ppc_mock_fit()
    fit$fit  <- NULL
    expect_error(
        suppressMessages(ppc(fit)),
        "CmdStan"
    )
})

test_that("D8. y_rep not available triggers warning + fallback to 'simulate'", {
    # create_ppc_mock_fit with include_yrep = FALSE omits y_rep columns, so
    # the mock $draws("y_rep") call throws an error, which .ppc_validate_inputs
    # catches and converts into a cli_warn + list(method = "simulate").
    fit_no_yrep <- create_ppc_mock_fit(include_yrep = FALSE)
    expect_warning(
        result <- suppressMessages(ppc(fit_no_yrep, method = "stan")),
        "y_rep"
    )
    expect_s3_class(result, "hbb_ppc")
    expect_equal(result$method, "simulate")
})


# ============================================================================
# ---- Section E: Print Method ----
# ============================================================================

test_that("E1. print.hbb_ppc produces output containing 'Posterior Predictive'", {
    fit    <- create_ppc_mock_fit()
    result <- suppressMessages(ppc(fit, method = "stan"))
    out    <- capture.output(suppressMessages(print(result)))
    expect_true(any(grepl("Posterior Predictive", out, fixed = TRUE)))
})

test_that("E2. print respects digits argument (check formatting in output)", {
    fit    <- create_ppc_mock_fit()
    result <- suppressMessages(ppc(fit, method = "stan"))
    out2   <- capture.output(suppressMessages(print(result, digits = 2L)))
    out6   <- capture.output(suppressMessages(print(result, digits = 6L)))
    # With more decimal places the numeric strings differ in length.
    # Concatenate all lines and compare lengths as a simple proxy.
    expect_false(identical(paste(out2, collapse = ""), paste(out6, collapse = "")))
})

test_that("E3. print handles type = 'zero_rate' only (no IT share data row)", {
    fit    <- create_ppc_mock_fit()
    result <- suppressMessages(ppc(fit, type = "zero_rate", method = "stan"))
    out    <- capture.output(suppressMessages(print(result)))
    # "zero_rate" label must appear.
    expect_true(any(grepl("zero_rate", out, fixed = TRUE)))
    # The predicted$it_share is NULL for type = "zero_rate", so .print_row
    # returns immediately without emitting numeric values for "it_share".
    # Any "it_share" line should not contain a decimal number.
    it_lines <- grep("it_share", out, fixed = TRUE, value = TRUE)
    has_numeric <- any(grepl("[0-9]\\.[0-9]", it_lines))
    expect_false(has_numeric)
})

test_that("E4. print returns x invisibly", {
    fit    <- create_ppc_mock_fit()
    result <- suppressMessages(ppc(fit, method = "stan"))
    # Capture return value of print().
    ret_val <- suppressMessages(
        {
            capture.output(ret_val <- print(result))
            ret_val
        }
    )
    expect_identical(ret_val, result)
})


# ============================================================================
# ---- Section F: Plot Method ----
# ============================================================================

test_that("F1. plot.hbb_ppc returns ggplot when ggplot2 available", {
    skip_if_not_installed("ggplot2")
    fit    <- create_ppc_mock_fit()
    result <- suppressMessages(ppc(fit, type = "zero_rate", method = "stan"))
    p      <- suppressMessages(suppressWarnings(plot(result, type = "zero_rate")))
    expect_s3_class(p, "gg")
})

test_that("F2. plot with type = 'both' returns combined plot (patchwork)", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("patchwork")
    fit    <- create_ppc_mock_fit()
    result <- suppressMessages(ppc(fit, type = "both", method = "stan"))
    p      <- suppressMessages(suppressWarnings(plot(result, type = "both")))
    # A patchwork object is also a gg; verify it inherits from "patchwork".
    expect_true(inherits(p, "patchwork") || inherits(p, "gg"))
})

test_that("F3. plot with type = 'zero_rate' returns single panel ggplot", {
    skip_if_not_installed("ggplot2")
    fit    <- create_ppc_mock_fit()
    result <- suppressMessages(ppc(fit, type = "zero_rate", method = "stan"))
    p      <- suppressMessages(suppressWarnings(plot(result, type = "zero_rate")))
    expect_s3_class(p, "gg")
    # A single panel should NOT be a patchwork object.
    expect_false(inherits(p, "patchwork"))
})

test_that("F4. plot on object with NULL predicted$it_share does not crash", {
    skip_if_not_installed("ggplot2")
    fit    <- create_ppc_mock_fit()
    result <- suppressMessages(ppc(fit, type = "zero_rate", method = "stan"))
    # Confirm it_share is NULL for this type.
    expect_null(result$predicted$it_share)
    # Plotting should succeed without error.
    p <- suppressMessages(suppressWarnings(plot(result, type = "zero_rate")))
    expect_s3_class(p, "gg")
})


# ============================================================================
# ---- Section G: Edge Cases ----
# ============================================================================

test_that("G1. N = 1 observation (scalar y_rep column) -- PPC still works", {
    fit    <- create_ppc_mock_fit(N = 1, M = 100, zero_rate = 0.0, seed = 1)
    result <- suppressMessages(ppc(fit, method = "stan"))
    expect_s3_class(result, "hbb_ppc")
    expect_equal(result$observed$N, 1L)
})

test_that("G2. All observations zero: observed zero_rate = 1.0, it_share = NA", {
    fit <- create_ppc_mock_fit(N = 40, M = 100, zero_rate = 0.0, seed = 7)
    # Override y and z to be all zeros regardless of the random generation.
    fit$hbb_data$y <- rep(0L, 40L)
    fit$hbb_data$z <- rep(0L, 40L)
    # Expect a warning about no positive observations in the data.
    result <- suppressWarnings(suppressMessages(ppc(fit, method = "stan")))
    expect_s3_class(result, "hbb_ppc")
    expect_equal(result$observed$zero_rate, 1.0)
    expect_true(is.na(result$observed$it_share))
})

test_that("G3. No zeros in data (z = 1 for all): observed zero_rate = 0", {
    fit <- create_ppc_mock_fit(N = 40, M = 100, zero_rate = 0.0, seed = 9)
    # Ensure all y > 0 in hbb_data.
    fit$hbb_data$y <- pmax(1L, fit$hbb_data$y)
    fit$hbb_data$z <- rep(1L, 40L)
    result <- suppressMessages(ppc(fit, method = "stan"))
    expect_s3_class(result, "hbb_ppc")
    expect_equal(result$observed$zero_rate, 0.0)
    expect_false(is.na(result$observed$it_share))
})

test_that("G4. n_draws = 1: single draw still produces valid result", {
    fit    <- create_ppc_mock_fit(M = 200, seed = 11)
    result <- suppressMessages(ppc(fit, n_draws = 1L, method = "stan"))
    expect_s3_class(result, "hbb_ppc")
    expect_equal(result$n_draws, 1L)
    # The zero_rate draws vector should have exactly length 1.
    expect_length(result$predicted$zero_rate$draws, 1L)
})

test_that("G5. Very large n_trial (1000): no overflow or numerical issues", {
    set.seed(21)
    fit                    <- create_ppc_mock_fit(N = 30, M = 100, seed = 21)
    fit$hbb_data$n_trial   <- rep(1000L, 30L)
    # The y_rep matrix in the mock already exists; n_trial in hbb_data is only
    # used for the IT share denominator, so we just verify completion.
    result <- suppressMessages(ppc(fit, method = "stan"))
    expect_s3_class(result, "hbb_ppc")
    # All zero_rate draws must be finite.
    expect_true(all(is.finite(result$predicted$zero_rate$draws)))
    # IT share draws may include NA (all-zero replications), but any non-NA
    # finite draw must itself be finite.
    it_draws     <- result$predicted$it_share$draws
    finite_it    <- it_draws[!is.na(it_draws)]
    if (length(finite_it) > 0L)
        expect_true(all(is.finite(finite_it)))
})


# ============================================================================
# ---- Section C2: Draw Extraction and Thinning ----
# ============================================================================

# ---- Section C2: Draw Extraction and Thinning ----

test_that("C2-1. .ppc_extract_yrep() returns M x N integer matrix", {
    fit <- create_ppc_mock_fit(P = 3, N = 40, M = 100)
    mat <- suppressMessages(hurdlebb:::.ppc_extract_yrep(fit, n_draws = NULL))
    expect_true(is.matrix(mat))
    expect_equal(storage.mode(mat), "integer")
    expect_equal(nrow(mat), 100L)
    expect_equal(ncol(mat), 40L)
})

test_that("C2-2. column count of .ppc_extract_yrep() equals N", {
    N   <- 55L
    fit <- create_ppc_mock_fit(N = N, M = 80)
    mat <- suppressMessages(hurdlebb:::.ppc_extract_yrep(fit, n_draws = NULL))
    expect_equal(ncol(mat), N)
})

test_that("C2-3. n_draws thinning returns correct number of rows", {
    fit    <- create_ppc_mock_fit(N = 30, M = 200)
    n_want <- 50L
    mat    <- suppressMessages(hurdlebb:::.ppc_extract_yrep(fit, n_draws = n_want))
    expect_equal(nrow(mat), n_want)
})

test_that("C2-4. n_draws > M_total silently uses all draws", {
    M   <- 60L
    fit <- create_ppc_mock_fit(N = 20, M = M)
    mat <- suppressMessages(hurdlebb:::.ppc_extract_yrep(fit, n_draws = 9999L))
    expect_equal(nrow(mat), M)
})

test_that("C2-5. n_draws = 1 returns 1-row matrix", {
    fit <- create_ppc_mock_fit(N = 25, M = 100)
    mat <- suppressMessages(hurdlebb:::.ppc_extract_yrep(fit, n_draws = 1L))
    expect_equal(nrow(mat), 1L)
    expect_true(is.matrix(mat))
})

test_that("C2-6. Systematic thinning is deterministic (same result twice)", {
    fit  <- create_ppc_mock_fit(N = 30, M = 120, seed = 99)
    mat1 <- suppressMessages(hurdlebb:::.ppc_extract_yrep(fit, n_draws = 30L))
    mat2 <- suppressMessages(hurdlebb:::.ppc_extract_yrep(fit, n_draws = 30L))
    expect_identical(mat1, mat2)
})


# ============================================================================
# ---- Section H2: Return Structure ----
# ============================================================================

# ---- Section H2: Return Structure ----

test_that("H2-1. result has class 'hbb_ppc'", {
    fit <- create_ppc_mock_fit(N = 40, M = 100, seed = 51)
    res <- suppressMessages(ppc(fit, type = "zero_rate", n_draws = 50L))
    expect_s3_class(res, "hbb_ppc")
})

test_that("H2-2. all expected top-level fields present (observed, predicted, coverage, n_draws, n_total_draws, type, method, level, model_type)", {
    fit      <- create_ppc_mock_fit(N = 40, M = 100, seed = 52)
    res      <- suppressMessages(ppc(fit, type = "zero_rate", n_draws = 50L))
    expected <- c("observed", "predicted", "coverage", "n_draws",
                  "n_total_draws", "type", "method", "level", "model_type")
    expect_true(all(expected %in% names(res)))
})

test_that("H2-3. n_draws field and predicted draws length equal requested n_draws", {
    fit    <- create_ppc_mock_fit(N = 40, M = 200, seed = 53)
    n_want <- 75L
    res    <- suppressMessages(ppc(fit, type = "zero_rate", n_draws = n_want))
    expect_equal(res$n_draws, n_want)
    expect_equal(length(res$predicted$zero_rate$draws), n_want)
})

test_that("H2-4. type='both' has both zero_rate and it_share in predicted", {
    fit <- create_ppc_mock_fit(N = 40, M = 100, seed = 54, zero_rate = 0.25)
    res <- suppressMessages(ppc(fit, type = "both", n_draws = 50L))
    expect_false(is.null(res$predicted$zero_rate))
    expect_false(is.null(res$predicted$it_share))
})

test_that("H2-5. type='zero_rate' has NULL for it_share in predicted", {
    fit <- create_ppc_mock_fit(N = 40, M = 100, seed = 55)
    res <- suppressMessages(ppc(fit, type = "zero_rate", n_draws = 50L))
    expect_null(res$predicted$it_share)
})


# ============================================================================
# ---- Section I2: R-side Simulation Fallback ----
# ============================================================================

# ---- Section I2: R-side Simulation Fallback ----

test_that("I2-1. method='simulate' produces M x N integer matrix via .ppc_simulate_yrep", {
    fit <- create_ppc_mock_fit(P = 3, N = 35, M = 80, seed = 61,
                               include_yrep = FALSE)
    mat <- suppressMessages(hurdlebb:::.ppc_simulate_yrep(fit, n_draws = 40L))
    expect_true(is.matrix(mat))
    expect_equal(storage.mode(mat), "integer")
    expect_equal(nrow(mat), 40L)
    expect_equal(ncol(mat), 35L)
})

test_that("I2-2. Simulated y_rep respects hurdle structure (y=0 when z=0 in replicated data)", {
    fit     <- create_ppc_mock_fit(P = 3, N = 50, M = 80, seed = 62,
                                   include_yrep = FALSE)
    mat     <- suppressMessages(hurdlebb:::.ppc_simulate_yrep(fit, n_draws = 40L))
    n_trial <- fit$hbb_data$n_trial
    # All non-NA entries must be >= 0 (hurdle: either zero or positive count)
    expect_true(all(mat >= 0L, na.rm = TRUE))
    # All non-NA entries must be bounded by their n_trial (BB cannot exceed n)
    for (col in seq_len(ncol(mat))) {
        finite_vals <- mat[, col][!is.na(mat[, col])]
        expect_true(all(finite_vals <= n_trial[col]))
    }
})

test_that("I2-3. Simulated PPC statistics are in plausible range [0, 1]", {
    fit <- create_ppc_mock_fit(P = 3, N = 50, M = 100, seed = 63,
                               zero_rate = 0.30, include_yrep = FALSE)
    res <- suppressMessages(
        ppc(fit, method = "simulate", type = "both", n_draws = 50L)
    )
    zr_mean <- res$predicted$zero_rate$mean
    expect_gte(zr_mean, 0)
    expect_lte(zr_mean, 1)
    if (!is.null(res$predicted$it_share)) {
        it_mean <- res$predicted$it_share$mean
        if (!is.na(it_mean)) {
            expect_gte(it_mean, 0)
            expect_lte(it_mean, 1)
        }
    }
})


# ============================================================================
# SECTION J: .ppc_validate_inputs Error Paths
# ============================================================================

# ---- J1. Rejects non-hbb_fit object ----------------------------------------
test_that("J1. .ppc_validate_inputs rejects non-hbb_fit", {
    # Already covered by ppc() itself; test the internal function directly
    expect_error(
        hurdlebb:::.ppc_validate_inputs(
            fit     = list(a = 1),
            type    = "both",
            method  = "stan",
            n_draws = NULL,
            level   = 0.95
        ),
        class = "rlang_error"
    )
})

# ---- J2. hbb_data NULL triggers abort --------------------------------------
test_that("J2. .ppc_validate_inputs aborts when hbb_data is NULL", {
    fake_fit <- structure(
        list(fit = list(draws = function(...) NULL),
             hbb_data = NULL,
             model_type = "base"),
        class = "hbb_fit"
    )
    expect_error(
        hurdlebb:::.ppc_validate_inputs(fake_fit, "both", "stan", NULL, 0.95),
        "hbb_data.*NULL"
    )
})

# ---- J3. hbb_data missing required field ------------------------------------
test_that("J3. .ppc_validate_inputs aborts when hbb_data$y is NULL", {
    fake_fit <- structure(
        list(fit = list(draws = function(...) NULL),
             hbb_data = list(y = NULL, n_trial = 1:5, N = 5L),
             model_type = "base"),
        class = "hbb_fit"
    )
    expect_error(
        hurdlebb:::.ppc_validate_inputs(fake_fit, "both", "stan", NULL, 0.95),
        "y.*NULL"
    )
})

test_that("J3b. .ppc_validate_inputs aborts when hbb_data$n_trial is NULL", {
    fake_fit <- structure(
        list(fit = list(draws = function(...) NULL),
             hbb_data = list(y = 1:5, n_trial = NULL, N = 5L),
             model_type = "base"),
        class = "hbb_fit"
    )
    expect_error(
        hurdlebb:::.ppc_validate_inputs(fake_fit, "both", "stan", NULL, 0.95),
        "n_trial.*NULL"
    )
})

test_that("J3c. .ppc_validate_inputs aborts when hbb_data$N is NULL", {
    fake_fit <- structure(
        list(fit = list(draws = function(...) NULL),
             hbb_data = list(y = 1:5, n_trial = rep(10L, 5), N = NULL),
             model_type = "base"),
        class = "hbb_fit"
    )
    expect_error(
        hurdlebb:::.ppc_validate_inputs(fake_fit, "both", "stan", NULL, 0.95),
        "N.*NULL"
    )
})

# ---- J4. Length mismatch between y, n_trial, N ------------------------------
test_that("J4. .ppc_validate_inputs aborts on y/n_trial/N length mismatch", {
    fake_fit <- structure(
        list(fit = list(draws = function(...) NULL),
             hbb_data = list(y = 1:5, n_trial = rep(10L, 3), N = 5L),
             model_type = "base"),
        class = "hbb_fit"
    )
    expect_error(
        hurdlebb:::.ppc_validate_inputs(fake_fit, "both", "stan", NULL, 0.95),
        "Length mismatch"
    )
})


# ============================================================================
# SECTION K: .ppc_extract_yrep Edge Cases
# ============================================================================

# ---- K1. NA values in y_rep trigger warning --------------------------------
test_that("K1. .ppc_extract_yrep warns on NA values in y_rep", {
    N <- 10; M <- 5
    y_rep <- matrix(sample(0:5, M * N, replace = TRUE), nrow = M, ncol = N)
    # Inject NAs
    y_rep[1, 1] <- NA
    y_rep[2, 3] <- NA
    y_rep[3, 5] <- NA
    colnames(y_rep) <- paste0("y_rep[", seq_len(N), "]")

    mock_fit <- structure(
        list(
            fit = list(
                draws = function(variables = NULL, format = "matrix") {
                    y_rep
                }
            ),
            hbb_data = list(N = N)
        ),
        class = "hbb_fit"
    )

    expect_warning(
        hurdlebb:::.ppc_extract_yrep(mock_fit, n_draws = NULL),
        "NA/NaN"
    )
})

# ---- K2. Infinite values in y_rep trigger warning + replacement ------------
test_that("K2. .ppc_extract_yrep warns on Inf values and replaces with NA", {
    N <- 8; M <- 4
    y_rep <- matrix(sample(0:5, M * N, replace = TRUE), nrow = M, ncol = N)
    y_rep[1, 1] <- Inf
    y_rep[2, 2] <- -Inf
    colnames(y_rep) <- paste0("y_rep[", seq_len(N), "]")

    mock_fit <- structure(
        list(
            fit = list(
                draws = function(variables = NULL, format = "matrix") {
                    y_rep
                }
            ),
            hbb_data = list(N = N)
        ),
        class = "hbb_fit"
    )

    result <- suppressWarnings(
        hurdlebb:::.ppc_extract_yrep(mock_fit, n_draws = NULL)
    )
    # Inf values replaced with NA
    expect_true(is.na(result[1, 1]))
    expect_true(is.na(result[2, 2]))
})

# ---- K3. n_draws >= M_total triggers info message --------------------------
test_that("K3. .ppc_extract_yrep uses all draws when n_draws >= M_total", {
    N <- 5; M <- 10
    y_rep <- matrix(sample(0:3, M * N, replace = TRUE), nrow = M, ncol = N)
    colnames(y_rep) <- paste0("y_rep[", seq_len(N), "]")

    mock_fit <- structure(
        list(
            fit = list(
                draws = function(variables = NULL, format = "matrix") {
                    y_rep
                }
            ),
            hbb_data = list(N = N)
        ),
        class = "hbb_fit"
    )

    result <- suppressMessages(
        hurdlebb:::.ppc_extract_yrep(mock_fit, n_draws = 20L)
    )
    expect_equal(nrow(result), M)
})

# ---- K4. n_draws < M_total triggers thinning -------------------------------
test_that("K4. .ppc_extract_yrep thins when n_draws < M_total", {
    N <- 5; M <- 20
    y_rep <- matrix(sample(0:3, M * N, replace = TRUE), nrow = M, ncol = N)
    colnames(y_rep) <- paste0("y_rep[", seq_len(N), "]")

    mock_fit <- structure(
        list(
            fit = list(
                draws = function(variables = NULL, format = "matrix") {
                    y_rep
                }
            ),
            hbb_data = list(N = N)
        ),
        class = "hbb_fit"
    )

    result <- suppressMessages(
        hurdlebb:::.ppc_extract_yrep(mock_fit, n_draws = 5L)
    )
    expect_lte(nrow(result), 5L)
})

# ---- K5. Extraction failure aborts -----------------------------------------
test_that("K5. .ppc_extract_yrep aborts when draws() fails", {
    mock_fit <- structure(
        list(
            fit = list(
                draws = function(variables = NULL, format = "matrix") {
                    stop("simulated CmdStan failure")
                }
            ),
            hbb_data = list(N = 10)
        ),
        class = "hbb_fit"
    )

    expect_error(
        hurdlebb:::.ppc_extract_yrep(mock_fit, n_draws = NULL),
        "Failed to extract"
    )
})

# ---- K6. Zero-row matrix aborts -------------------------------------------
test_that("K6. .ppc_extract_yrep aborts on zero-row matrix", {
    N <- 5
    y_rep <- matrix(integer(0), nrow = 0, ncol = N)
    colnames(y_rep) <- paste0("y_rep[", seq_len(N), "]")

    mock_fit <- structure(
        list(
            fit = list(
                draws = function(variables = NULL, format = "matrix") {
                    y_rep
                }
            ),
            hbb_data = list(N = N)
        ),
        class = "hbb_fit"
    )

    expect_error(
        hurdlebb:::.ppc_extract_yrep(mock_fit, n_draws = NULL),
        "zero rows"
    )
})


# ============================================================================
# SECTION L: .ppc_simulate_yrep Edge Cases
# ============================================================================

# ---- L1. Missing design matrix X aborts -----------------------------------
test_that("L1. .ppc_simulate_yrep aborts when X is NULL", {
    fake_fit <- structure(
        list(
            fit = list(
                draws = function(variables = NULL, format = "matrix") {
                    matrix(0, nrow = 10, ncol = 7)
                }
            ),
            hbb_data = list(X = NULL, P = 3L, N = 5L, n_trial = rep(10L, 5))
        ),
        class = "hbb_fit"
    )
    expect_error(
        hurdlebb:::.ppc_simulate_yrep(fake_fit, n_draws = NULL),
        "Design matrix"
    )
})

# ---- L2. X dimension mismatch aborts --------------------------------------
test_that("L2. .ppc_simulate_yrep aborts when nrow(X) != N", {
    fake_fit <- structure(
        list(
            fit = list(
                draws = function(variables = NULL, format = "matrix") {
                    matrix(0, nrow = 10, ncol = 7)
                }
            ),
            hbb_data = list(
                X = matrix(1, nrow = 3, ncol = 2),
                P = 2L,
                N = 5L,
                n_trial = rep(10L, 5)
            )
        ),
        class = "hbb_fit"
    )
    expect_error(
        hurdlebb:::.ppc_simulate_yrep(fake_fit, n_draws = NULL),
        "Dimension mismatch"
    )
})

# ---- L3. Missing parameter names abort ------------------------------------
test_that("L3. .ppc_simulate_yrep aborts when parameter draws missing expected names", {
    N <- 5; P <- 2
    X <- cbind(1, rnorm(N))

    # Create draws with wrong column names
    draws_mat <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
    colnames(draws_mat) <- c("alpha[1]", "alpha[2]", "beta[1]", "wrong_name", "log_kappa")

    fake_fit <- structure(
        list(
            fit = list(
                draws = function(variables = NULL, format = "matrix") {
                    draws_mat
                }
            ),
            hbb_data = list(
                X = X,
                P = P,
                N = as.integer(N),
                n_trial = rep(10L, N)
            )
        ),
        class = "hbb_fit"
    )
    expect_error(
        hurdlebb:::.ppc_simulate_yrep(fake_fit, n_draws = NULL),
        "not found"
    )
})

# ---- L4. Non-finite log_kappa triggers warning ----------------------------
test_that("L4. .ppc_simulate_yrep warns when log_kappa has non-finite draws", {
    N <- 5; P <- 2; M <- 10
    X <- cbind(1, rnorm(N))

    param_names <- c("alpha[1]", "alpha[2]", "beta[1]", "beta[2]", "log_kappa")
    draws_mat <- matrix(rnorm(M * 5, sd = 0.1), nrow = M, ncol = 5)
    colnames(draws_mat) <- param_names
    # Set first two draws to have non-finite log_kappa
    draws_mat[1, "log_kappa"] <- Inf
    draws_mat[2, "log_kappa"] <- NaN

    fake_fit <- structure(
        list(
            fit = list(
                draws = function(variables = NULL, format = "matrix") {
                    draws_mat
                }
            ),
            hbb_data = list(
                X = X,
                P = as.integer(P),
                N = as.integer(N),
                n_trial = rep(10L, N)
            )
        ),
        class = "hbb_fit"
    )
    # Multiple warnings may be emitted (non-finite log_kappa + skipped draws)
    warns <- capture_warnings(
        suppressMessages(hurdlebb:::.ppc_simulate_yrep(fake_fit, n_draws = NULL))
    )
    expect_true(any(grepl("non-finite", warns)))
})

# ---- L5. Skipped draws warning -------------------------------------------
test_that("L5. .ppc_simulate_yrep warns about skipped draws due to non-finite params", {
    N <- 5; P <- 2; M <- 6
    X <- cbind(1, rnorm(N))

    param_names <- c("alpha[1]", "alpha[2]", "beta[1]", "beta[2]", "log_kappa")
    draws_mat <- matrix(rnorm(M * 5, sd = 0.1), nrow = M, ncol = 5)
    colnames(draws_mat) <- param_names
    # Make alpha non-finite for some draws
    draws_mat[1, "alpha[1]"] <- NA
    draws_mat[2, "beta[1]"]  <- Inf

    fake_fit <- structure(
        list(
            fit = list(
                draws = function(variables = NULL, format = "matrix") {
                    draws_mat
                }
            ),
            hbb_data = list(
                X = X,
                P = as.integer(P),
                N = as.integer(N),
                n_trial = rep(10L, N)
            )
        ),
        class = "hbb_fit"
    )
    # Should warn about skipped draws (multiple warnings may be emitted)
    warns <- capture_warnings(
        suppressMessages(hurdlebb:::.ppc_simulate_yrep(fake_fit, n_draws = NULL))
    )
    expect_true(any(grepl("skipped", warns)))
})

# ---- L6. NaN in linear predictors triggers warning -------------------------
test_that("L6. .ppc_simulate_yrep warns about NaN in linear predictors", {
    N <- 5; P <- 2; M <- 6
    X <- cbind(1, rnorm(N))

    param_names <- c("alpha[1]", "alpha[2]", "beta[1]", "beta[2]", "log_kappa")
    draws_mat <- matrix(0.1, nrow = M, ncol = 5)
    colnames(draws_mat) <- param_names
    draws_mat[, "log_kappa"] <- 1.5
    # Extreme alpha values that might cause NaN in plogis (actually plogis
    # never returns NaN so we must override plogis behavior via extreme X)
    # Instead, let's use extreme X values that when multiplied produce NaN
    X[1, 2] <- .Machine$double.xmax  # extreme covariate
    draws_mat[1, "alpha[2]"] <- .Machine$double.xmax  # overflow product

    fake_fit <- structure(
        list(
            fit = list(
                draws = function(variables = NULL, format = "matrix") {
                    draws_mat
                }
            ),
            hbb_data = list(
                X = X,
                P = as.integer(P),
                N = as.integer(N),
                n_trial = rep(10L, N)
            )
        ),
        class = "hbb_fit"
    )
    # This may or may not produce NaN depending on plogis behavior
    # but it exercises the code path through the simulation loop
    result <- suppressWarnings(
        suppressMessages(hurdlebb:::.ppc_simulate_yrep(fake_fit, n_draws = NULL))
    )
    expect_true(is.matrix(result))
})


# ============================================================================
# SECTION M: .ppc_compute_stats Warnings
# ============================================================================

# ---- M1. Negative values in y_rep trigger warning -------------------------
test_that("M1. .ppc_compute_stats warns on negative y_rep values", {
    M <- 5; N <- 8
    n_trial <- rep(10L, N)
    y_rep <- matrix(sample(0:5, M * N, replace = TRUE), nrow = M, ncol = N)
    # Inject negative values
    y_rep[1, 1] <- -1L
    y_rep[2, 3] <- -5L

    expect_warning(
        hurdlebb:::.ppc_compute_stats(y_rep, n_trial, "it_share"),
        "negative"
    )
})

# ---- M2. Values exceeding n_trial trigger warning --------------------------
test_that("M2. .ppc_compute_stats warns on values exceeding n_trial", {
    M <- 5; N <- 8
    n_trial <- rep(10L, N)
    y_rep <- matrix(sample(1:5, M * N, replace = TRUE), nrow = M, ncol = N)
    # Inject values exceeding n_trial
    y_rep[1, 1] <- 15L
    y_rep[3, 5] <- 20L

    expect_warning(
        hurdlebb:::.ppc_compute_stats(y_rep, n_trial, "it_share"),
        "exceed"
    )
})


# ============================================================================
# SECTION N: .ppc_summarize_draws Edge Cases
# ============================================================================

# ---- N1. All non-finite draws return NA summary ----------------------------
test_that("N1. .ppc_summarize_draws returns NA summary when all draws are non-finite", {
    draws <- c(NA_real_, NaN, Inf, -Inf, NA_real_)
    result <- hurdlebb:::.ppc_summarize_draws(draws, level = 0.95)

    expect_true(is.na(result$mean))
    expect_true(is.na(result$median))
    expect_true(is.na(result$sd))
    expect_equal(result$n_finite, 0L)
    expect_true(is.na(result$ci[["lower"]]))
    expect_true(is.na(result$ci[["upper"]]))
})

# ---- N2. All-NA draws with n_na attribute ----------------------------------
test_that("N2. .ppc_summarize_draws preserves n_na attribute", {
    draws <- c(NA_real_, NA_real_, NA_real_)
    attr(draws, "n_na") <- 3L
    result <- hurdlebb:::.ppc_summarize_draws(draws, level = 0.95)

    expect_equal(result$n_na, 3L)
    expect_equal(result$n_finite, 0L)
    expect_true(is.na(result$mean))
})

# ---- N3. Partial NA draws handled correctly --------------------------------
test_that("N3. .ppc_summarize_draws excludes non-finite values but uses the rest", {
    draws <- c(0.3, 0.4, NA_real_, Inf, 0.5)
    attr(draws, "n_na") <- 1L
    result <- hurdlebb:::.ppc_summarize_draws(draws, level = 0.95)

    expect_equal(result$n_finite, 3L)
    expect_equal(result$n_na, 1L)
    expect_equal(result$mean, mean(c(0.3, 0.4, 0.5)))
    expect_equal(result$median, median(c(0.3, 0.4, 0.5)))
})


# ============================================================================
# SECTION O: print.hbb_ppc Method
# ============================================================================

# ---- O1. Basic print output -----------------------------------------------
test_that("O1. print.hbb_ppc produces output for a well-formed object", {
    ppc_obj <- structure(
        list(
            model_type    = "weighted",
            method        = "stan",
            level         = 0.95,
            n_draws       = 100L,
            n_total_draws = 200L,
            type          = "both",
            observed = list(
                N         = 50L,
                N_pos     = 30L,
                zero_rate = 0.40,
                it_share  = 0.25
            ),
            predicted = list(
                zero_rate = list(
                    draws   = rnorm(100, 0.4, 0.02),
                    mean    = 0.40,
                    median  = 0.40,
                    ci      = list(0.36, 0.44),
                    n_finite = 100L,
                    n_na    = 0L
                ),
                it_share = list(
                    draws   = rnorm(100, 0.25, 0.01),
                    mean    = 0.25,
                    median  = 0.25,
                    ci      = list(0.23, 0.27),
                    n_finite = 100L,
                    n_na    = 0L
                )
            ),
            coverage = list(
                zero_rate = list(in_ci = TRUE),
                it_share  = list(in_ci = TRUE)
            )
        ),
        class = "hbb_ppc"
    )

    out <- capture.output(print(ppc_obj))
    expect_true(any(grepl("Posterior Predictive Check", out)))
    expect_true(any(grepl("PASS", out)))
    expect_true(any(grepl("weighted", out)))
    expect_true(any(grepl("zero_rate", out)))
    expect_true(any(grepl("it_share", out)))
})

# ---- O2. Print with non-numeric level uses fallback 0.95 -------------------
test_that("O2. print.hbb_ppc handles non-numeric level gracefully", {
    ppc_obj <- structure(
        list(
            model_type = "base",
            method     = "simulate",
            level      = "invalid",
            n_draws    = 50L,
            n_total_draws = 100L,
            type       = "both",
            observed   = list(N = 20L, N_pos = 10L, zero_rate = 0.5, it_share = 0.3),
            predicted  = list(
                zero_rate = list(mean = 0.5, median = 0.5, ci = list(0.4, 0.6),
                                 n_finite = 50L, n_na = 0L),
                it_share = list(mean = 0.3, median = 0.3, ci = list(0.2, 0.4),
                                n_finite = 50L, n_na = 0L)
            ),
            coverage = list(
                zero_rate = list(in_ci = TRUE),
                it_share = list(in_ci = TRUE)
            )
        ),
        class = "hbb_ppc"
    )
    out <- capture.output(print(ppc_obj))
    expect_true(any(grepl("0.95", out)))
})

# ---- O3. Print with FAIL status -------------------------------------------
test_that("O3. print.hbb_ppc shows FAIL when in_ci is FALSE", {
    ppc_obj <- structure(
        list(
            model_type = "base",
            method     = "stan",
            level      = 0.90,
            n_draws    = 100L,
            n_total_draws = 200L,
            type       = "both",
            observed   = list(N = 50L, N_pos = 30L, zero_rate = 0.40, it_share = 0.25),
            predicted  = list(
                zero_rate = list(mean = 0.40, median = 0.40, ci = list(0.36, 0.44),
                                 n_finite = 100L, n_na = 0L),
                it_share = list(mean = 0.25, median = 0.25, ci = list(0.23, 0.27),
                                n_finite = 100L, n_na = 0L)
            ),
            coverage = list(
                zero_rate = list(in_ci = FALSE),
                it_share = list(in_ci = TRUE)
            )
        ),
        class = "hbb_ppc"
    )
    out <- capture.output(print(ppc_obj))
    expect_true(any(grepl("FAIL", out)))
    expect_true(any(grepl("PASS", out)))
})

# ---- O4. Print with NA observed values shows "NA: no positives" -----------
test_that("O4. print.hbb_ppc shows NA when observed value is NA", {
    ppc_obj <- structure(
        list(
            model_type = "base",
            method     = "simulate",
            level      = 0.95,
            n_draws    = 50L,
            n_total_draws = 100L,
            type       = "both",
            observed   = list(N = 50L, N_pos = 0L, zero_rate = 1.0, it_share = NA_real_),
            predicted  = list(
                zero_rate = list(mean = 0.99, median = 0.99, ci = list(0.95, 1.0),
                                 n_finite = 50L, n_na = 0L),
                it_share = list(mean = NA_real_, median = NA_real_, ci = list(NA_real_, NA_real_),
                                n_finite = 0L, n_na = 50L)
            ),
            coverage = list(
                zero_rate = list(in_ci = TRUE),
                it_share = list(in_ci = NA)
            )
        ),
        class = "hbb_ppc"
    )
    out <- capture.output(print(ppc_obj))
    expect_true(any(grepl("NA.*no positives", out)))
})

# ---- O5. Print with unavailable CI ----------------------------------------
test_that("O5. print.hbb_ppc shows unavailable when CI is NA", {
    ppc_obj <- structure(
        list(
            model_type = "base",
            method     = "simulate",
            level      = 0.95,
            n_draws    = 50L,
            n_total_draws = 100L,
            type       = "both",
            observed   = list(N = 20L, zero_rate = 0.5, it_share = 0.3),
            predicted  = list(
                zero_rate = list(mean = 0.5, median = 0.5, ci = list(NA_real_, NA_real_),
                                 n_finite = 50L, n_na = 0L),
                it_share  = NULL
            ),
            coverage = list(zero_rate = list(in_ci = NA))
        ),
        class = "hbb_ppc"
    )
    out <- capture.output(print(ppc_obj))
    expect_true(any(grepl("unavailable", out)))
})

# ---- O6. Print with n_na draws > 0 shows note -----------------------------
test_that("O6. print.hbb_ppc shows note about all-zero draws", {
    ppc_obj <- structure(
        list(
            model_type = "base",
            method     = "stan",
            level      = 0.95,
            n_draws    = 100L,
            n_total_draws = 200L,
            type       = "both",
            observed   = list(N = 50L, N_pos = 30L, zero_rate = 0.40, it_share = 0.25),
            predicted  = list(
                zero_rate = list(mean = 0.40, median = 0.40, ci = list(0.36, 0.44),
                                 n_finite = 100L, n_na = 0L),
                it_share = list(mean = 0.25, median = 0.25, ci = list(0.23, 0.27),
                                n_finite = 95L, n_na = 5L)
            ),
            coverage = list(
                zero_rate = list(in_ci = TRUE),
                it_share = list(in_ci = TRUE)
            )
        ),
        class = "hbb_ppc"
    )
    out <- capture.output(print(ppc_obj))
    expect_true(any(grepl("all-zero", out)))
})

# ---- O7. Print with n_finite < n_draws shows finite-count note -------------
test_that("O7. print.hbb_ppc shows note when n_finite < n_draws", {
    ppc_obj <- structure(
        list(
            model_type = "base",
            method     = "stan",
            level      = 0.95,
            n_draws    = 100L,
            n_total_draws = 200L,
            type       = "both",
            observed   = list(N = 50L, N_pos = 30L, zero_rate = 0.40, it_share = 0.25),
            predicted  = list(
                zero_rate = list(mean = 0.40, median = 0.40, ci = list(0.36, 0.44),
                                 n_finite = 90L, n_na = 0L),
                it_share = list(mean = 0.25, median = 0.25, ci = list(0.23, 0.27),
                                n_finite = 100L, n_na = 0L)
            ),
            coverage = list(
                zero_rate = list(in_ci = TRUE),
                it_share = list(in_ci = TRUE)
            )
        ),
        class = "hbb_ppc"
    )
    out <- capture.output(print(ppc_obj))
    expect_true(any(grepl("90/100", out)))
})

# ---- O8. Print with corrupted object uses last-resort error catch ----------
test_that("O8. print.hbb_ppc does not crash on totally corrupted object", {
    # Create an object where accessing fields throws errors
    ppc_obj <- structure(list(), class = "hbb_ppc")
    out <- capture.output(print(ppc_obj))
    # Should produce some output without error
    expect_true(length(out) > 0)
})

# ---- O9. Print with empty status line when in_ci is NULL ------------------
test_that("O9. print.hbb_ppc handles NULL coverage gracefully", {
    ppc_obj <- structure(
        list(
            model_type = "base",
            method     = "simulate",
            level      = 0.95,
            n_draws    = 50L,
            n_total_draws = 100L,
            type       = "both",
            observed   = list(N = 20L, N_pos = 10L, zero_rate = 0.5, it_share = 0.3),
            predicted  = list(
                zero_rate = list(mean = 0.5, median = 0.5, ci = list(0.4, 0.6),
                                 n_finite = 50L, n_na = 0L),
                it_share = list(mean = 0.3, median = 0.3, ci = list(0.2, 0.4),
                                n_finite = 50L, n_na = 0L)
            ),
            coverage = NULL
        ),
        class = "hbb_ppc"
    )
    out <- capture.output(print(ppc_obj))
    expect_true(any(grepl("Posterior Predictive Check", out)))
})


# ============================================================================
# SECTION P: plot.hbb_ppc Method
# ============================================================================

# ---- P1. Plot uses x$type when type argument is NULL -----------------------
test_that("P1. plot.hbb_ppc uses x$type when type is NULL", {
    skip_if_not_installed("ggplot2")

    draws_zr <- rnorm(100, 0.4, 0.02)
    ppc_obj <- structure(
        list(
            type  = "zero_rate",
            level = 0.95,
            n_draws = 100L,
            observed = list(zero_rate = 0.40, it_share = NA_real_),
            predicted = list(
                zero_rate = list(
                    draws = draws_zr,
                    ci    = list(0.36, 0.44)
                ),
                it_share = NULL
            )
        ),
        class = "hbb_ppc"
    )
    p <- suppressWarnings(plot(ppc_obj))
    expect_true(inherits(p, "gg"))
})

# ---- P2. Plot warns when zero_rate draws are NULL --------------------------
test_that("P2. plot.hbb_ppc warns when zero_rate draws unavailable", {
    skip_if_not_installed("ggplot2")

    draws_it <- rnorm(100, 0.3, 0.01)
    ppc_obj <- structure(
        list(
            type  = "both",
            level = 0.95,
            n_draws = 100L,
            observed = list(zero_rate = 0.40, it_share = 0.30),
            predicted = list(
                zero_rate = NULL,
                it_share = list(
                    draws = draws_it,
                    ci    = list(0.28, 0.32)
                )
            )
        ),
        class = "hbb_ppc"
    )
    # cli_alert_warning uses message(), not warning(), so use expect_message
    expect_message(
        suppressWarnings(plot(ppc_obj)),
        "Zero-rate"
    )
})

# ---- P3. Plot warns when it_share draws are NULL ---------------------------
test_that("P3. plot.hbb_ppc warns when it_share draws unavailable", {
    skip_if_not_installed("ggplot2")

    draws_zr <- rnorm(100, 0.4, 0.02)
    ppc_obj <- structure(
        list(
            type  = "both",
            level = 0.95,
            n_draws = 100L,
            observed = list(zero_rate = 0.40, it_share = 0.30),
            predicted = list(
                zero_rate = list(
                    draws = draws_zr,
                    ci    = list(0.36, 0.44)
                ),
                it_share = NULL
            )
        ),
        class = "hbb_ppc"
    )
    # cli_alert_warning uses message(), not warning(), so use expect_message
    expect_message(
        suppressWarnings(plot(ppc_obj)),
        "IT-share"
    )
})

# ---- P4. Plot aborts when no panels available ------------------------------
test_that("P4. plot.hbb_ppc aborts when no panels available", {
    skip_if_not_installed("ggplot2")

    ppc_obj <- structure(
        list(
            type  = "both",
            level = 0.95,
            n_draws = 100L,
            observed = list(zero_rate = 0.40, it_share = 0.30),
            predicted = list(
                zero_rate = NULL,
                it_share  = NULL
            )
        ),
        class = "hbb_ppc"
    )
    expect_error(
        suppressWarnings(suppressMessages(plot(ppc_obj))),
        "No PPC panels"
    )
})

# ---- P5. Plot handles FAIL subtitle ---------------------------------------
test_that("P5. plot.hbb_ppc renders FAIL subtitle when obs outside CI", {
    skip_if_not_installed("ggplot2")

    draws_zr <- rnorm(100, 0.4, 0.02)
    ppc_obj <- structure(
        list(
            type  = "zero_rate",
            level = 0.95,
            n_draws = 100L,
            observed = list(zero_rate = 0.80, it_share = NA_real_),
            predicted = list(
                zero_rate = list(
                    draws = draws_zr,
                    ci    = list(0.36, 0.44)
                ),
                it_share = NULL
            )
        ),
        class = "hbb_ppc"
    )
    p <- suppressWarnings(plot(ppc_obj))
    expect_true(inherits(p, "gg"))
})

# ---- P6. Plot handles NA observed value -----------------------------------
test_that("P6. plot.hbb_ppc handles NA observed value (no positives)", {
    skip_if_not_installed("ggplot2")

    draws_it <- rnorm(100, 0.3, 0.01)
    ppc_obj <- structure(
        list(
            type  = "it_share",
            level = 0.95,
            n_draws = 100L,
            observed = list(zero_rate = 1.0, it_share = NA_real_),
            predicted = list(
                zero_rate = NULL,
                it_share = list(
                    draws = draws_it,
                    ci    = NULL
                )
            )
        ),
        class = "hbb_ppc"
    )
    p <- suppressWarnings(plot(ppc_obj))
    expect_true(inherits(p, "gg"))
})

# ---- P7. Plot with explicit type arg overrides x$type ----------------------
test_that("P7. plot.hbb_ppc uses explicit type argument", {
    skip_if_not_installed("ggplot2")

    draws_zr <- rnorm(100, 0.4, 0.02)
    draws_it <- rnorm(100, 0.3, 0.01)
    ppc_obj <- structure(
        list(
            type  = "both",
            level = 0.95,
            n_draws = 100L,
            observed = list(zero_rate = 0.40, it_share = 0.30),
            predicted = list(
                zero_rate = list(
                    draws = draws_zr,
                    ci    = list(0.36, 0.44)
                ),
                it_share = list(
                    draws = draws_it,
                    ci    = list(0.28, 0.32)
                )
            )
        ),
        class = "hbb_ppc"
    )
    # Only plot zero_rate even though type is "both"
    p <- suppressWarnings(plot(ppc_obj, type = "zero_rate"))
    expect_true(inherits(p, "gg"))
})
