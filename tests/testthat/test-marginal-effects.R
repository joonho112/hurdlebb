# ============================================================================
# test-marginal-effects.R --- testthat tests for marginal-effects.R
#
# Sections:
#   A. Algebraic identities     (total = ext + int, product rule)
#   B. Decomposition table      (13 cols, shares, sign patterns, CIs)
#   C. Wald delta method        (gradient, SE, CI, zero-V, symmetry)
#   D. Input validation         (classes, levels, n_draws, missing fields)
#   E. Print method             (output sections, digits, corrupted objects)
#   F. Edge cases               (P=1, n_draws=1, large P, N=1, all fields)
#   G. Reversal probability     (known reversal, known reinforcement)
#   H. Subsampling              (n_draws, systematic thinning, dimensions)
#   I. Point estimates          (named P-vectors, total = ext + int)
#   J. Mean_q and mean_mu       (bounded, correct length, plogis)
# ============================================================================


# ---- Mock object constructors -----------------------------------------------

#' Create a mock hbb_fit object for AME testing
#'
#' Generates draws from MVN(theta_true, Sigma_true) and wraps in
#' a CmdStanR-like fit interface with a design matrix X.
#'
#' @param P Number of covariates per margin (including intercept).
#' @param N Number of observations.
#' @param M Number of MCMC draws.
#' @param seed Random seed.
#' @param alpha_true Optional true alpha vector of length P.
#' @param beta_true Optional true beta vector of length P.
#' @param X Optional N x P design matrix.
#' @return An S3 object of class "hbb_fit".
create_ame_mock_fit <- function(P = 3, N = 100, M = 2000, seed = 42,
                                alpha_true = NULL, beta_true = NULL,
                                X = NULL) {
    D <- 2L * P + 1L

    param_names <- c(
        paste0("alpha[", seq_len(P), "]"),
        paste0("beta[", seq_len(P), "]"),
        "log_kappa"
    )

    set.seed(seed)

    if (is.null(alpha_true)) {
        alpha_true <- c(-0.5, seq(-0.3, 0.3, length.out = max(P - 1L, 1)))
    }
    if (is.null(beta_true)) {
        beta_true <- c(0.2, seq(-0.2, 0.2, length.out = max(P - 1L, 1)))
    }
    if (length(alpha_true) < P) {
        alpha_true <- rep_len(alpha_true, P)
    } else if (length(alpha_true) > P) {
        alpha_true <- alpha_true[seq_len(P)]
    }
    if (length(beta_true) < P) {
        beta_true <- rep_len(beta_true, P)
    } else if (length(beta_true) > P) {
        beta_true <- beta_true[seq_len(P)]
    }

    theta_true <- c(alpha_true, beta_true, 1.5)
    Sigma_true <- diag(D) * 0.01

    draws <- MASS::mvrnorm(M, theta_true, Sigma_true)
    colnames(draws) <- param_names

    # Design matrix
    if (is.null(X)) {
        X <- cbind(
            intercept = rep(1, N),
            matrix(rnorm(N * (P - 1L)), nrow = N, ncol = P - 1L,
                   dimnames = list(NULL, paste0("x", seq_len(P - 1L))))
        )
    }

    mock_cmdstan <- list(
        draws = function(variables = NULL, format = "matrix") {
            if (!is.null(variables)) draws[, variables, drop = FALSE]
            else draws
        }
    )

    z       <- rbinom(N, 1, 0.6)
    n_trial <- sample(10:50, N, replace = TRUE)
    y       <- ifelse(z == 1, rbinom(N, n_trial, 0.3), 0L)

    hbb_data <- list(
        X       = X,
        N       = as.integer(N),
        P       = as.integer(P),
        z       = z,
        n_trial = n_trial,
        y       = y
    )

    structure(
        list(
            fit        = mock_cmdstan,
            hbb_data   = hbb_data,
            model_type = "weighted"
        ),
        class = "hbb_fit"
    )
}


#' Create a mock hbb_cholesky object for AME testing
#'
#' @param fit A mock hbb_fit object.
#' @param theta_hat Optional point estimate override.
#' @return An S3 object of class "hbb_cholesky".
create_ame_mock_cholesky <- function(fit, theta_hat = NULL) {
    P <- fit$hbb_data$P
    D <- 2L * P + 1L
    param_names <- c(
        paste0("alpha[", seq_len(P), "]"),
        paste0("beta[", seq_len(P), "]"),
        "log_kappa"
    )

    theta_corrected <- fit$fit$draws(variables = param_names, format = "matrix")
    if (is.null(theta_hat)) theta_hat <- colMeans(theta_corrected)

    structure(
        list(
            theta_corrected = theta_corrected,
            theta_hat       = theta_hat,
            D               = D,
            M               = nrow(theta_corrected),
            P               = P
        ),
        class = "hbb_cholesky"
    )
}


#' Create a mock hbb_sandwich object for AME testing
#'
#' @param P Number of covariates per margin.
#' @param V_sand Optional sandwich variance matrix (D x D).
#' @return An S3 object of class "hbb_sandwich".
create_ame_mock_sandwich <- function(P = 3, V_sand = NULL) {
    D <- 2L * P + 1L

    if (is.null(V_sand)) {
        V_sand <- diag(D) * 0.01 + 0.001
    }

    plabs <- c(
        paste0("alpha_", c("intercept", paste0("x", seq_len(P - 1L)))),
        paste0("beta_",  c("intercept", paste0("x", seq_len(P - 1L)))),
        "log_kappa"
    )

    structure(
        list(
            V_sand       = V_sand,
            H_obs_inv    = diag(D) * 0.005,
            Sigma_MCMC   = diag(D) * 0.5,
            DER          = diag(V_sand) / (diag(D) * 0.005),
            param_labels = plabs,
            D            = D,
            P            = P,
            N            = 100L
        ),
        class = "hbb_sandwich"
    )
}


#' Suppress cli messages during testing
quietly <- function(expr) {
    suppressMessages(expr)
}


# ============================================================================
# Section A: Algebraic identities (total = ext + int, product rule)
# ============================================================================

test_that("A1. total AME equals ext + int for every draw and covariate", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 500, seed = 1)
    res <- quietly(ame(fit))

    diff_mat <- res$total_ame_draws - (res$ext_ame_draws + res$int_ame_draws)
    expect_equal(max(abs(diff_mat)), 0, tolerance = 1e-14)
})

test_that("A2. total AME point estimate equals ext + int at theta_hat", {
    fit <- create_ame_mock_fit(P = 4, N = 80, M = 300, seed = 2)
    res <- quietly(ame(fit))

    expect_equal(
        unname(res$ame_total_hat),
        unname(res$ame_ext_hat + res$ame_int_hat),
        tolerance = 1e-14
    )
})

test_that("A3. posterior mean of total draws = mean(ext) + mean(int)", {
    fit <- create_ame_mock_fit(P = 3, N = 100, M = 1000, seed = 3)
    res <- quietly(ame(fit))

    expect_equal(
        res$total_summary$post_mean,
        res$ext_summary$post_mean + res$int_summary$post_mean,
        tolerance = 1e-12
    )
})

test_that("A4. total = ext + int with cholesky-corrected draws", {
    fit <- create_ame_mock_fit(P = 4, N = 80, M = 300, seed = 4)
    chol <- create_ame_mock_cholesky(fit)
    res <- quietly(ame(fit, cholesky = chol))

    diff_mat <- res$total_ame_draws - (res$ext_ame_draws + res$int_ame_draws)
    expect_equal(max(abs(diff_mat)), 0, tolerance = 1e-14)
})

test_that("A5. manual product-rule verification at theta_hat", {
    set.seed(100)
    P <- 2L
    N <- 30L
    X <- cbind(intercept = rep(1, N), x1 = rnorm(N))

    alpha <- c(-0.5, 0.3)
    beta  <- c(0.2, -0.4)

    q  <- plogis(as.numeric(X %*% alpha))
    mu <- plogis(as.numeric(X %*% beta))

    # Manual AME via product rule
    ext_manual <- alpha * mean(q * (1 - q) * mu)
    int_manual <- beta  * mean(mu * (1 - mu) * q)

    # Use .ame_at_theta
    theta <- c(alpha, beta, 1.5)
    result <- .ame_at_theta(theta, X, P)

    expect_equal(result$ext,   ext_manual, tolerance = 1e-14)
    expect_equal(result$int,   int_manual, tolerance = 1e-14)
    expect_equal(result$total, ext_manual + int_manual, tolerance = 1e-14)
})

test_that("A6. AME at zero coefficients is zero", {
    P <- 3L
    N <- 50L
    set.seed(200)
    X <- cbind(1, matrix(rnorm(N * (P - 1L)), N, P - 1L))
    colnames(X) <- c("intercept", "x1", "x2")

    theta_zero <- c(rep(0, 2L * P), 1.5)
    result <- .ame_at_theta(theta_zero, X, P)

    expect_equal(result$ext,   rep(0, P))
    expect_equal(result$int,   rep(0, P))
    expect_equal(result$total, rep(0, P))
    expect_equal(result$mean_q,  0.5)
    expect_equal(result$mean_mu, 0.5)
})

test_that("A7. ext AME sign follows alpha sign", {
    set.seed(301)
    P <- 3L
    N <- 100L
    X <- cbind(1, matrix(rnorm(N * 2), N, 2))
    colnames(X) <- c("intercept", "x1", "x2")

    alpha <- c(0, 0.5, -0.5)
    beta  <- c(0, 0, 0)
    theta <- c(alpha, beta, 1.5)

    result <- .ame_at_theta(theta, X, P)
    expect_true(result$ext[2] > 0)
    expect_true(result$ext[3] < 0)
})

test_that("A8. int AME sign follows beta sign", {
    set.seed(302)
    P <- 3L
    N <- 100L
    X <- cbind(1, matrix(rnorm(N * 2), N, 2))
    colnames(X) <- c("intercept", "x1", "x2")

    alpha <- c(0, 0, 0)
    beta  <- c(0, 0.8, -0.3)
    theta <- c(alpha, beta, 1.5)

    result <- .ame_at_theta(theta, X, P)
    expect_true(result$int[2] > 0)
    expect_true(result$int[3] < 0)
})

test_that("A9. analytical AME matches for constant X (intercept only)", {
    P <- 1L
    N <- 50L
    D <- 3L
    M <- 10L

    alpha_true <- 1.0
    beta_true  <- -0.5
    theta_true <- c(alpha_true, beta_true, 1.0)

    draws <- matrix(rep(theta_true, each = M), nrow = M, ncol = D)
    param_names <- c("alpha[1]", "beta[1]", "log_kappa")
    colnames(draws) <- param_names

    X <- matrix(1, nrow = N, ncol = 1)
    colnames(X) <- "intercept"

    mock_cmdstan <- list(
        draws = function(variables = NULL, format = "matrix") {
            if (!is.null(variables)) draws[, variables, drop = FALSE]
            else draws
        }
    )

    fit <- structure(
        list(
            fit = mock_cmdstan,
            model_type = "weighted",
            hbb_data = list(X = X, N = N, P = 1L, z = rep(1L, N),
                            n_trial = rep(10L, N), y = rep(5L, N))
        ),
        class = "hbb_fit"
    )

    res <- quietly(ame(fit))

    q  <- plogis(alpha_true)
    mu <- plogis(beta_true)
    expected_ext   <- alpha_true * q * (1 - q) * mu
    expected_int   <- beta_true  * mu * (1 - mu) * q
    expected_total <- expected_ext + expected_int

    expect_equal(unname(res$ame_ext_hat[1]),   expected_ext,   tolerance = 1e-12)
    expect_equal(unname(res$ame_int_hat[1]),   expected_int,   tolerance = 1e-12)
    expect_equal(unname(res$ame_total_hat[1]), expected_total, tolerance = 1e-12)

    # All draws identical => posterior mean equals any single draw
    expect_equal(
        unname(res$ext_ame_draws[1, 1]), expected_ext, tolerance = 1e-12
    )
})


# ============================================================================
# Section B: Decomposition table (shares, sign patterns, structure)
# ============================================================================

test_that("B1. decomposition table has exactly 13 columns", {
    fit <- create_ame_mock_fit(P = 4, N = 50, M = 300, seed = 10)
    res <- quietly(ame(fit))
    expect_equal(ncol(res$decomp_table), 13L)
})

test_that("B2. decomposition table has correct column names", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 200, seed = 11)
    res <- quietly(ame(fit))

    expected_cols <- c(
        "covariate", "ext_ame", "ext_ci_lo", "ext_ci_hi",
        "int_ame", "int_ci_lo", "int_ci_hi",
        "total_ame", "total_ci_lo", "total_ci_hi",
        "ext_share", "int_share", "sign_pattern"
    )
    expect_equal(names(res$decomp_table), expected_cols)
})

test_that("B3. decomposition table has P-1 rows (non-intercept only)", {
    fit <- create_ame_mock_fit(P = 4, N = 50, M = 200, seed = 12)
    res <- quietly(ame(fit))

    expect_equal(nrow(res$decomp_table), res$P - 1L)
    expect_false("intercept" %in% res$decomp_table$covariate)
})

test_that("B4. ext_share + int_share = 100 for each covariate", {
    fit <- create_ame_mock_fit(P = 5, N = 100, M = 500, seed = 13)
    res <- quietly(ame(fit))

    dt <- res$decomp_table
    share_sum <- dt$ext_share + dt$int_share
    expect_equal(share_sum, rep(100, nrow(dt)), tolerance = 1e-10)
})

test_that("B5. shares are non-negative", {
    fit <- create_ame_mock_fit(P = 5, N = 50, M = 200, seed = 14)
    res <- quietly(ame(fit))

    dt <- res$decomp_table
    expect_true(all(dt$ext_share >= 0))
    expect_true(all(dt$int_share >= 0))
    expect_true(all(dt$ext_share <= 100))
    expect_true(all(dt$int_share <= 100))
})

test_that("B6. sign_pattern is 'opposing' when ext and int have different signs", {
    P <- 2L
    N <- 100L
    set.seed(501)
    X <- cbind(intercept = 1, poverty = rnorm(N))

    fit <- create_ame_mock_fit(
        P = P, N = N, M = 1000, seed = 502,
        alpha_true = c(-0.5, -0.8),
        beta_true  = c(0.2, 0.6),
        X = X
    )

    res <- quietly(ame(fit))
    expect_equal(res$decomp_table$sign_pattern[1], "opposing")
})

test_that("B7. sign_pattern is 'reinforcing' when ext and int have same sign", {
    P <- 2L
    N <- 100L
    set.seed(601)
    X <- cbind(intercept = 1, x1 = rnorm(N))

    fit <- create_ame_mock_fit(
        P = P, N = N, M = 1000, seed = 602,
        alpha_true = c(0, 0.5),
        beta_true  = c(0, 0.5),
        X = X
    )
    res <- quietly(ame(fit))
    expect_equal(res$decomp_table$sign_pattern[1], "reinforcing")
})

test_that("B8. total_ame = ext_ame + int_ame in decomp_table", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 300, seed = 15)
    res <- quietly(ame(fit))

    dt <- res$decomp_table
    expect_equal(dt$total_ame, dt$ext_ame + dt$int_ame, tolerance = 1e-12)
})

test_that("B9. CI ordering: ci_lo <= point <= ci_hi in decomp_table", {
    fit <- create_ame_mock_fit(P = 4, N = 100, M = 500, seed = 16)
    res <- quietly(ame(fit))

    dt <- res$decomp_table
    expect_true(all(dt$ext_ci_lo   <= dt$ext_ci_hi))
    expect_true(all(dt$int_ci_lo   <= dt$int_ci_hi))
    expect_true(all(dt$total_ci_lo <= dt$total_ci_hi))
})

test_that("B10. ame_decomposition extracts same table as decomp_table field", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 200, seed = 17)
    res <- quietly(ame(fit))

    dt <- ame_decomposition(res)
    expect_identical(dt, res$decomp_table)
})

test_that("B11. P=1 produces 0-row decomposition table", {
    P <- 1L
    N <- 30L
    X <- matrix(1, nrow = N, ncol = 1, dimnames = list(NULL, "intercept"))

    fit <- create_ame_mock_fit(
        P = P, N = N, M = 100, seed = 18,
        alpha_true = c(-0.5), beta_true = c(0.3), X = X
    )
    res <- quietly(ame(fit))

    expect_equal(nrow(res$decomp_table), 0L)
    expect_equal(ncol(res$decomp_table), 13L)
})

test_that("B12. sign_pattern is only 'reinforcing' or 'opposing'", {
    fit <- create_ame_mock_fit(P = 5, N = 100, M = 200, seed = 19)
    res <- quietly(ame(fit))

    expect_true(all(res$decomp_table$sign_pattern %in%
                        c("reinforcing", "opposing")))
})

test_that("B13. decomp_table matches total_summary for non-intercept rows", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 200, seed = 20)
    res <- quietly(ame(fit))

    dt <- res$decomp_table
    ts <- res$total_summary

    for (i in seq_len(nrow(dt))) {
        cov_name <- dt$covariate[i]
        ts_row <- ts[ts$covariate == cov_name, ]
        expect_equal(dt$total_ame[i], ts_row$post_mean, tolerance = 1e-14)
    }
})


# ============================================================================
# Section C: Wald delta method (gradient, SE, CI consistency)
# ============================================================================

test_that("C1. Wald AME is NULL when no sandwich provided", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 300, seed = 30)
    res <- quietly(ame(fit))
    expect_null(res$wald_summary)
})

test_that("C2. Wald AME is computed when sandwich provided", {
    fit  <- create_ame_mock_fit(P = 3, N = 50, M = 300, seed = 31)
    sand <- create_ame_mock_sandwich(P = 3)
    res  <- quietly(ame(fit, sandwich = sand))

    expect_false(is.null(res$wald_summary))
    expect_s3_class(res$wald_summary, "data.frame")
})

test_that("C3. Wald summary has correct 13 columns", {
    fit  <- create_ame_mock_fit(P = 3, N = 50, M = 300, seed = 32)
    sand <- create_ame_mock_sandwich(P = 3)
    res  <- quietly(ame(fit, sandwich = sand))

    ws <- res$wald_summary
    required_cols <- c(
        "covariate", "ext_point", "ext_se", "ext_lo", "ext_hi",
        "int_point", "int_se", "int_lo", "int_hi",
        "total_point", "total_se", "total_lo", "total_hi"
    )
    for (col in required_cols) {
        expect_true(col %in% names(ws), info = paste("Missing:", col))
    }
})

test_that("C4. Wald summary has P rows", {
    P <- 4L
    fit  <- create_ame_mock_fit(P = P, N = 80, M = 300, seed = 33)
    sand <- create_ame_mock_sandwich(P = P)
    res  <- quietly(ame(fit, sandwich = sand))

    expect_equal(nrow(res$wald_summary), P)
})

test_that("C5. Wald total_point = ext_point + int_point", {
    fit  <- create_ame_mock_fit(P = 4, N = 80, M = 300, seed = 34)
    sand <- create_ame_mock_sandwich(P = 4)
    res  <- quietly(ame(fit, sandwich = sand))

    ws <- res$wald_summary
    expect_equal(ws$total_point, ws$ext_point + ws$int_point, tolerance = 1e-12)
})

test_that("C6. Wald SEs are non-negative", {
    fit  <- create_ame_mock_fit(P = 3, N = 50, M = 300, seed = 35)
    sand <- create_ame_mock_sandwich(P = 3)
    res  <- quietly(ame(fit, sandwich = sand))

    ws <- res$wald_summary
    expect_true(all(ws$ext_se >= 0))
    expect_true(all(ws$int_se >= 0))
    expect_true(all(ws$total_se >= 0))
})

test_that("C7. Wald CIs are symmetric around point estimate", {
    fit  <- create_ame_mock_fit(P = 3, N = 50, M = 300, seed = 36)
    sand <- create_ame_mock_sandwich(P = 3)
    res  <- quietly(ame(fit, sandwich = sand))

    ws <- res$wald_summary

    mid_ext   <- (ws$ext_lo   + ws$ext_hi)   / 2
    mid_int   <- (ws$int_lo   + ws$int_hi)   / 2
    mid_total <- (ws$total_lo + ws$total_hi) / 2

    expect_equal(mid_ext,   ws$ext_point,   tolerance = 1e-12)
    expect_equal(mid_int,   ws$int_point,   tolerance = 1e-12)
    expect_equal(mid_total, ws$total_point, tolerance = 1e-12)
})

test_that("C8. Wald SE = 0 when V_sand is zero matrix", {
    P <- 2L
    D <- 2L * P + 1L
    fit  <- create_ame_mock_fit(P = P, N = 30, M = 200, seed = 37)
    sand <- create_ame_mock_sandwich(P = P, V_sand = matrix(0, D, D))

    res <- quietly(ame(fit, sandwich = sand))
    ws  <- res$wald_summary

    expect_true(all(ws$ext_se == 0))
    expect_true(all(ws$int_se == 0))
    expect_true(all(ws$total_se == 0))
})

test_that("C9. Wald point estimates match ame_*_hat", {
    fit  <- create_ame_mock_fit(P = 3, N = 50, M = 300, seed = 38)
    sand <- create_ame_mock_sandwich(P = 3)
    res  <- quietly(ame(fit, sandwich = sand))

    expect_equal(
        unname(res$wald_summary$total_point),
        unname(res$ame_total_hat),
        tolerance = 1e-10
    )
    expect_equal(
        unname(res$wald_summary$ext_point),
        unname(res$ame_ext_hat),
        tolerance = 1e-10
    )
    expect_equal(
        unname(res$wald_summary$int_point),
        unname(res$ame_int_hat),
        tolerance = 1e-10
    )
})

test_that("C10. Wald SE positive for non-trivial V_sand", {
    P <- 2L
    D <- 2L * P + 1L
    set.seed(390)
    X <- cbind(1, rnorm(50))
    colnames(X) <- c("intercept", "x1")
    fit  <- create_ame_mock_fit(P = P, N = 50, M = 200, seed = 39, X = X)
    sand <- create_ame_mock_sandwich(P = P, V_sand = diag(D) * 0.05)

    res <- quietly(ame(fit, sandwich = sand))
    ws  <- res$wald_summary

    expect_true(all(ws$total_se > 0))
    expect_true(all(ws$ext_se > 0))
})

test_that("C11. Wald CI width increases with level", {
    fit  <- create_ame_mock_fit(P = 2, N = 50, M = 200, seed = 310)
    sand <- create_ame_mock_sandwich(P = 2)

    r90 <- quietly(ame(fit, sandwich = sand, level = 0.90))
    r95 <- quietly(ame(fit, sandwich = sand, level = 0.95))
    r99 <- quietly(ame(fit, sandwich = sand, level = 0.99))

    w90 <- r90$wald_summary$total_hi - r90$wald_summary$total_lo
    w95 <- r95$wald_summary$total_hi - r95$wald_summary$total_lo
    w99 <- r99$wald_summary$total_hi - r99$wald_summary$total_lo

    expect_true(all(w90 <= w95 + 1e-12))
    expect_true(all(w95 <= w99 + 1e-12))
})


# ============================================================================
# Section D: Input validation (classes, levels, n_draws)
# ============================================================================

test_that("D1. rejects non-hbb_fit object", {
    expect_error(ame(list(a = 1)), "hbb_fit")
    expect_error(ame("not_a_fit"), "hbb_fit")
    expect_error(ame(42), "hbb_fit")
})

test_that("D2. rejects non-hbb_cholesky object for cholesky arg", {
    fit <- create_ame_mock_fit()
    expect_error(ame(fit, cholesky = list(a = 1)),  "hbb_cholesky")
    expect_error(ame(fit, cholesky = "not_chol"),   "hbb_cholesky")
})

test_that("D3. rejects non-hbb_sandwich object for sandwich arg", {
    fit <- create_ame_mock_fit()
    expect_error(
        ame(fit, sandwich = list(V_sand = diag(5))),
        "hbb_sandwich"
    )
})

test_that("D4. rejects level = 0", {
    fit <- create_ame_mock_fit()
    expect_error(ame(fit, level = 0), "strictly")
})

test_that("D5. rejects level = 1", {
    fit <- create_ame_mock_fit()
    expect_error(ame(fit, level = 1), "strictly")
})

test_that("D6. rejects negative level", {
    fit <- create_ame_mock_fit()
    expect_error(ame(fit, level = -0.5), "strictly")
})

test_that("D7. rejects string level", {
    fit <- create_ame_mock_fit()
    expect_error(ame(fit, level = "a"), "single numeric")
})

test_that("D8. rejects vector level", {
    fit <- create_ame_mock_fit()
    expect_error(ame(fit, level = c(0.9, 0.95)), "single numeric")
})

test_that("D9. rejects NA level", {
    fit <- create_ame_mock_fit()
    expect_error(ame(fit, level = NA), "single numeric")
})

test_that("D10. rejects n_draws = 0", {
    fit <- create_ame_mock_fit()
    expect_error(ame(fit, n_draws = 0), "positive integer")
})

test_that("D11. rejects negative n_draws", {
    fit <- create_ame_mock_fit()
    expect_error(ame(fit, n_draws = -1), "positive integer")
})

test_that("D12. rejects fractional n_draws", {
    fit <- create_ame_mock_fit()
    expect_error(ame(fit, n_draws = 1.5), "positive integer")
})

test_that("D13. rejects string n_draws", {
    fit <- create_ame_mock_fit()
    expect_error(ame(fit, n_draws = "a"), "positive integer")
})

test_that("D14. rejects NA n_draws", {
    fit <- create_ame_mock_fit()
    expect_error(ame(fit, n_draws = NA), "positive integer")
})

test_that("D15. rejects fit with missing hbb_data$X", {
    fit <- create_ame_mock_fit()
    fit$hbb_data$X <- NULL
    expect_error(quietly(ame(fit)), "X|matrix")
})

test_that("D16. rejects fit with X not a matrix", {
    fit <- create_ame_mock_fit()
    fit$hbb_data$X <- "not_a_matrix"
    expect_error(quietly(ame(fit)), "matrix")
})

test_that("D17. rejects fit$fit NULL when no cholesky", {
    fit <- create_ame_mock_fit()
    fit$fit <- NULL
    expect_error(ame(fit), "NULL")
})

test_that("D18. rejects fit with NULL hbb_data", {
    bad_fit <- structure(
        list(fit = list(), model_type = "weighted", hbb_data = NULL),
        class = "hbb_fit"
    )
    expect_error(ame(bad_fit), "NULL")
})

test_that("D19. rejects cholesky with wrong number of columns", {
    fit <- create_ame_mock_fit(P = 3)
    chol <- create_ame_mock_cholesky(fit)
    chol$theta_corrected <- cbind(chol$theta_corrected, 0)  # add extra col
    expect_error(quietly(ame(fit, cholesky = chol)), "columns")
})

test_that("D20. ame_decomposition rejects non-hbb_ame input", {
    expect_error(ame_decomposition(list(a = 1)), "hbb_ame")
    expect_error(ame_decomposition("not_ame"),   "hbb_ame")
    expect_error(ame_decomposition(42),          "hbb_ame")
})

test_that("D21. rejects X with wrong dimensions", {
    fit <- create_ame_mock_fit(P = 2, N = 50, M = 100)
    fit$hbb_data$X <- matrix(0, 30, 2)  # wrong N
    expect_error(ame(fit), "dimensions")
})

test_that("D22. accepts valid level values", {
    fit <- create_ame_mock_fit(P = 2, N = 20, M = 100, seed = 400)
    for (lvl in c(0.50, 0.90, 0.95, 0.99)) {
        res <- quietly(ame(fit, level = lvl))
        expect_equal(res$level, lvl)
    }
})


# ============================================================================
# Section E: Print method (output sections, digits, corrupted objects)
# ============================================================================

test_that("E1. print runs without error and produces output", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 200, seed = 50)
    res <- quietly(ame(fit))

    out <- capture.output(print(res))
    expect_true(length(out) > 0)
})

test_that("E2. print returns object invisibly", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 200, seed = 51)
    res <- quietly(ame(fit))

    out <- capture.output(ret <- print(res))
    expect_identical(ret, res)
})

test_that("E3. print output contains key section headers", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 200, seed = 52)
    res <- quietly(ame(fit))

    out  <- capture.output(print(res))
    full <- paste(out, collapse = "\n")
    expect_true(grepl("Average Marginal Effects", full))
    expect_true(grepl("Decomposition", full))
    expect_true(grepl("Observations", full))
})

test_that("E4. print shows Wald section when sandwich provided", {
    fit  <- create_ame_mock_fit(P = 3, N = 50, M = 200, seed = 53)
    sand <- create_ame_mock_sandwich(P = 3)
    res  <- quietly(ame(fit, sandwich = sand))

    out  <- capture.output(print(res))
    full <- paste(out, collapse = "\n")
    expect_true(grepl("Wald AME", full))
    expect_true(grepl("agreement", full, ignore.case = TRUE))
})

test_that("E5. print does not show Wald section when no sandwich", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 200, seed = 54)
    res <- quietly(ame(fit))

    out  <- capture.output(print(res))
    full <- paste(out, collapse = "\n")
    expect_false(grepl("Wald AME", full))
})

test_that("E6. print respects digits argument", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 200, seed = 55)
    res <- quietly(ame(fit))

    out2 <- capture.output(print(res, digits = 2))
    out6 <- capture.output(print(res, digits = 6))
    expect_false(identical(out2, out6))
})

test_that("E7. print handles corrupted hbb_ame object gracefully", {
    corrupted <- structure(
        list(
            N = 100, P = 3, D = 7,
            M_use = 500, M_total = 500,
            level = 0.95,
            mean_q = runif(10), mean_mu = runif(10),
            decomp_table = NULL,
            reversal_probs = NULL,
            wald_summary = NULL
        ),
        class = "hbb_ame"
    )
    # tryCatch inside print protects from error
    expect_no_error(capture.output(print(corrupted)))
})

test_that("E8. print shows sign patterns in decomposition", {
    P <- 3L
    N <- 100L
    set.seed(560)
    X <- cbind(1, rnorm(N), rnorm(N))
    colnames(X) <- c("intercept", "poverty", "urban")

    fit <- create_ame_mock_fit(
        P = P, N = N, M = 500, seed = 56,
        alpha_true = c(0, -0.8, 0.5),
        beta_true  = c(0,  0.6, 0.5),
        X = X
    )
    res <- quietly(ame(fit))

    out  <- capture.output(print(res))
    full <- paste(out, collapse = "\n")
    expect_true(grepl("opposing",    full))
    expect_true(grepl("reinforcing", full))
})

test_that("E9. print shows Reversal probabilities", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 200, seed = 57)
    res <- quietly(ame(fit))

    out  <- capture.output(print(res))
    full <- paste(out, collapse = "\n")
    expect_true(grepl("Reversal", full, ignore.case = TRUE))
})

test_that("E10. print shows dimensions", {
    fit <- create_ame_mock_fit(P = 3, N = 80, M = 200, seed = 58)
    res <- quietly(ame(fit))

    out  <- capture.output(print(res))
    full <- paste(out, collapse = "\n")
    expect_true(grepl("80", full))   # N
    expect_true(grepl("200", full))  # M_use
})


# ============================================================================
# Section F: Edge cases (P=1, n_draws=1, large P, all required fields)
# ============================================================================

test_that("F1. minimal P=1 (intercept only) works correctly", {
    P <- 1L
    N <- 30L
    X <- matrix(1, nrow = N, ncol = 1, dimnames = list(NULL, "intercept"))

    fit <- create_ame_mock_fit(
        P = P, N = N, M = 500, seed = 70,
        alpha_true = c(-0.5), beta_true = c(0.3), X = X
    )
    res <- quietly(ame(fit))

    expect_s3_class(res, "hbb_ame")
    expect_equal(res$P, 1L)
    expect_equal(res$D, 3L)
    expect_equal(nrow(res$decomp_table), 0L)
    expect_equal(nrow(res$ext_summary), 1L)
    expect_equal(nrow(res$int_summary), 1L)
    expect_equal(nrow(res$total_summary), 1L)
    expect_equal(length(res$reversal_probs), 0L)
})

test_that("F2. n_draws = 1 works without error", {
    fit <- create_ame_mock_fit(P = 2, N = 30, M = 100, seed = 71)
    res <- quietly(ame(fit, n_draws = 1L))

    expect_equal(res$M_use, 1L)
    expect_equal(nrow(res$ext_ame_draws), 1L)
    expect_equal(nrow(res$int_ame_draws), 1L)
    expect_equal(nrow(res$total_ame_draws), 1L)

    # Identity still holds
    diff <- res$total_ame_draws - (res$ext_ame_draws + res$int_ame_draws)
    expect_equal(max(abs(diff)), 0, tolerance = 1e-14)
})

test_that("F3. n_draws > M_total uses all draws", {
    M <- 200L
    fit <- create_ame_mock_fit(P = 2, N = 30, M = M, seed = 72)
    res <- quietly(ame(fit, n_draws = 5000L))

    expect_equal(res$M_use, M)
    expect_equal(res$M_total, M)
})

test_that("F4. large P=8 works correctly", {
    P <- 8L
    N <- 100L
    M <- 500L
    set.seed(770)
    X <- cbind(1, matrix(rnorm(N * (P - 1L)), N, P - 1L))
    colnames(X) <- c("intercept", paste0("x", 1:(P - 1L)))

    fit <- create_ame_mock_fit(P = P, N = N, M = M, seed = 77, X = X)
    res <- quietly(ame(fit))

    expect_equal(res$P, P)
    expect_equal(res$D, 2L * P + 1L)
    expect_equal(nrow(res$decomp_table), P - 1L)
    expect_equal(ncol(res$ext_ame_draws), P)
    expect_equal(length(res$reversal_probs), P - 1L)

    # Identity still holds
    diff_mat <- res$total_ame_draws - (res$ext_ame_draws + res$int_ame_draws)
    expect_equal(max(abs(diff_mat)), 0, tolerance = 1e-14)
})

test_that("F5. hbb_ame return has all required fields", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 200, seed = 74)
    res <- quietly(ame(fit))

    required_fields <- c(
        "theta_draws_used", "theta_hat",
        "ext_ame_draws", "int_ame_draws", "total_ame_draws",
        "ext_summary", "int_summary", "total_summary",
        "decomp_table", "reversal_probs", "wald_summary",
        "ame_ext_hat", "ame_int_hat", "ame_total_hat",
        "mean_q", "mean_mu",
        "N", "P", "D", "M_total", "M_use", "level", "cov_labels"
    )
    for (field in required_fields) {
        expect_true(field %in% names(res),
                    info = paste("Missing field:", field))
    }
})

test_that("F6. summary data frames have correct dimensions (P rows x 7 cols)", {
    P <- 4L
    fit <- create_ame_mock_fit(P = P, N = 50, M = 200, seed = 75)
    res <- quietly(ame(fit))

    expected_cols <- c(
        "covariate", "post_mean", "post_median",
        "ci_lo", "ci_hi", "post_sd", "pr_positive"
    )

    for (summ in list(res$ext_summary, res$int_summary, res$total_summary)) {
        expect_equal(nrow(summ), P)
        expect_equal(ncol(summ), 7L)
        expect_true(all(expected_cols %in% names(summ)))
    }
})

test_that("F7. pr_positive is in [0, 1]", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 500, seed = 76)
    res <- quietly(ame(fit))

    for (summ in list(res$ext_summary, res$int_summary, res$total_summary)) {
        expect_true(all(summ$pr_positive >= 0 & summ$pr_positive <= 1))
    }
})

test_that("F8. cholesky draws are used when provided", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 500, seed = 78)
    chol <- create_ame_mock_cholesky(fit)

    known_hat <- colMeans(chol$theta_corrected) + 0.01
    chol$theta_hat <- known_hat

    res <- quietly(ame(fit, cholesky = chol))

    expect_equal(unname(res$theta_hat), unname(known_hat), tolerance = 1e-14)
})

test_that("F9. deterministic output for same inputs", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 200, seed = 79)

    r1 <- quietly(ame(fit))
    r2 <- quietly(ame(fit))

    expect_identical(r1$ext_ame_draws, r2$ext_ame_draws)
    expect_identical(r1$int_ame_draws, r2$int_ame_draws)
    expect_identical(r1$decomp_table, r2$decomp_table)
})

test_that("F10. extreme positive coefficients don't cause NaN", {
    P <- 2L
    fit <- create_ame_mock_fit(
        P = P, N = 50, M = 100, seed = 710,
        alpha_true = c(0, 5.0), beta_true = c(0, 5.0)
    )
    res <- quietly(ame(fit))

    expect_false(any(is.nan(res$ext_ame_draws)))
    expect_false(any(is.nan(res$int_ame_draws)))
    expect_false(any(is.nan(res$total_ame_draws)))
})

test_that("F11. extreme negative coefficients don't cause NaN", {
    P <- 2L
    fit <- create_ame_mock_fit(
        P = P, N = 50, M = 100, seed = 711,
        alpha_true = c(0, -5.0), beta_true = c(0, -5.0)
    )
    res <- quietly(ame(fit))

    expect_false(any(is.nan(res$ext_ame_draws)))
    expect_false(any(is.nan(res$int_ame_draws)))
})

test_that("F12. very small N = 5 works", {
    fit <- create_ame_mock_fit(P = 2, N = 5, M = 50, seed = 712)
    res <- quietly(ame(fit))

    expect_equal(res$N, 5L)
    expect_equal(nrow(res$ext_ame_draws), 50L)

    diff_mat <- res$total_ame_draws - (res$ext_ame_draws + res$int_ame_draws)
    expect_equal(max(abs(diff_mat)), 0, tolerance = 1e-14)
})

test_that("F13. N = 1 does not crash", {
    P <- 2L
    N <- 1L
    M <- 50L
    D <- 5L

    X <- matrix(c(1, 0.5), nrow = 1, ncol = 2)
    colnames(X) <- c("intercept", "x1")

    set.seed(713)
    draws <- MASS::mvrnorm(M, rep(0, D), diag(D) * 0.01)
    param_names <- c("alpha[1]", "alpha[2]", "beta[1]", "beta[2]", "log_kappa")
    colnames(draws) <- param_names

    mock <- list(
        draws = function(variables = NULL, format = "matrix") {
            if (!is.null(variables)) draws[, variables, drop = FALSE]
            else draws
        }
    )

    fit <- structure(
        list(
            fit = mock, model_type = "weighted",
            hbb_data = list(X = X, N = N, P = P, z = 1L,
                            n_trial = 10L, y = 3L)
        ),
        class = "hbb_fit"
    )

    res <- quietly(ame(fit))
    expect_s3_class(res, "hbb_ame")
    expect_equal(res$N, 1L)

    diff_mat <- res$total_ame_draws - (res$ext_ame_draws + res$int_ame_draws)
    expect_equal(max(abs(diff_mat)), 0, tolerance = 1e-14)
})

test_that("F14. M = 2 (minimum practical draws) works", {
    fit <- create_ame_mock_fit(P = 2, N = 30, M = 2, seed = 714)
    res <- quietly(ame(fit))

    expect_equal(res$M_use, 2L)
    expect_equal(nrow(res$ext_ame_draws), 2L)
})


# ============================================================================
# Section G: Reversal probability (known reversal, known reinforcement)
# ============================================================================

test_that("G1. reversal probability near 1 for known opposing signs", {
    P <- 2L
    N <- 50L
    set.seed(801)
    X <- cbind(1, rnorm(N))
    colnames(X) <- c("intercept", "poverty")

    fit <- create_ame_mock_fit(
        P = P, N = N, M = 2000, seed = 80,
        alpha_true = c(-0.5, -1.0),
        beta_true  = c(0.2,  0.8),
        X = X
    )

    res <- quietly(ame(fit))
    expect_true(res$reversal_probs[["poverty"]] > 0.95)
})

test_that("G2. reversal probability near 0 for known same signs", {
    P <- 2L
    N <- 50L
    set.seed(811)
    X <- cbind(1, rnorm(N))
    colnames(X) <- c("intercept", "urban")

    fit <- create_ame_mock_fit(
        P = P, N = N, M = 2000, seed = 81,
        alpha_true = c(0.5, 0.8),
        beta_true  = c(0.2, 0.6),
        X = X
    )

    res <- quietly(ame(fit))
    expect_true(res$reversal_probs[["urban"]] < 0.05)
})

test_that("G3. reversal probs are named with covariate labels", {
    P <- 4L
    N <- 50L
    set.seed(820)
    X <- cbind(1, matrix(rnorm(N * 3), N, 3))
    colnames(X) <- c("intercept", "poverty", "urban", "black")

    fit <- create_ame_mock_fit(P = P, N = N, M = 200, seed = 82, X = X)
    res <- quietly(ame(fit))

    expect_equal(sort(names(res$reversal_probs)),
                 sort(c("poverty", "urban", "black")))
})

test_that("G4. reversal probability is in [0, 1]", {
    fit <- create_ame_mock_fit(P = 5, N = 50, M = 500, seed = 83)
    res <- quietly(ame(fit))

    for (nm in names(res$reversal_probs)) {
        pr <- res$reversal_probs[[nm]]
        expect_true(pr >= 0 && pr <= 1,
                    info = paste("Reversal prob for", nm, "=", pr))
    }
})

test_that("G5. number of reversal probs equals P - 1", {
    for (P in c(2, 3, 5)) {
        fit <- create_ame_mock_fit(P = P, N = 30, M = 100, seed = 84 + P)
        res <- quietly(ame(fit))
        expect_equal(length(res$reversal_probs), P - 1L)
    }
})

test_that("G6. reversal probs are empty for P = 1", {
    P <- 1L
    N <- 30L
    X <- matrix(1, nrow = N, ncol = 1, dimnames = list(NULL, "intercept"))

    fit <- create_ame_mock_fit(
        P = P, N = N, M = 100, seed = 85,
        alpha_true = c(-0.5), beta_true = c(0.3), X = X
    )
    res <- quietly(ame(fit))

    expect_equal(length(res$reversal_probs), 0L)
    expect_type(res$reversal_probs, "double")
})

test_that("G7. reversal probability near 1 for poverty reversal (strong signal)", {
    # Mimics the paper's finding: alpha < 0, beta > 0 with very tight draws
    P <- 2L
    N <- 200L
    D <- 5L
    M <- 500L

    set.seed(860)
    X <- cbind(1, rnorm(N))
    colnames(X) <- c("intercept", "poverty")

    # Deterministic draws (no noise) at opposing truth
    draws <- matrix(0, nrow = M, ncol = D)
    draws[, 1] <- -0.5
    draws[, 2] <- -0.8 + rnorm(M, 0, 0.01)  # alpha_poverty strongly negative
    draws[, 3] <- 0.2
    draws[, 4] <- 0.6 + rnorm(M, 0, 0.01)   # beta_poverty strongly positive
    draws[, 5] <- 1.5
    param_names <- c("alpha[1]", "alpha[2]", "beta[1]", "beta[2]", "log_kappa")
    colnames(draws) <- param_names

    mock <- list(
        draws = function(variables = NULL, format = "matrix") {
            if (!is.null(variables)) draws[, variables, drop = FALSE]
            else draws
        }
    )

    fit <- structure(
        list(
            fit = mock, model_type = "weighted",
            hbb_data = list(X = X, N = N, P = P, z = rep(1L, N),
                            n_trial = rep(10L, N), y = rep(3L, N))
        ),
        class = "hbb_fit"
    )

    res <- quietly(ame(fit))
    expect_true(res$reversal_probs[["poverty"]] > 0.99)
})


# ============================================================================
# Section H: Subsampling (n_draws, systematic thinning, dimensions)
# ============================================================================

test_that("H1. n_draws controls number of draws used", {
    M_total <- 1000L
    fit <- create_ame_mock_fit(P = 3, N = 50, M = M_total, seed = 90)

    res_all <- quietly(ame(fit))
    res_500 <- quietly(ame(fit, n_draws = 500L))
    res_100 <- quietly(ame(fit, n_draws = 100L))

    expect_equal(res_all$M_use, M_total)
    expect_equal(res_all$M_total, M_total)

    expect_true(res_500$M_use <= 500L)
    expect_true(res_500$M_use > 0L)

    expect_true(res_100$M_use <= 100L)
    expect_true(res_100$M_use > 0L)
})

test_that("H2. systematic thinning is deterministic (not random)", {
    M <- 100L
    fit <- create_ame_mock_fit(P = 2, N = 30, M = M, seed = 91)

    res1 <- quietly(ame(fit, n_draws = 25L))
    res2 <- quietly(ame(fit, n_draws = 25L))

    expect_identical(res1$ext_ame_draws, res2$ext_ame_draws)
    expect_identical(res1$theta_draws_used, res2$theta_draws_used)
})

test_that("H3. n_draws = NULL uses all draws", {
    M <- 300L
    fit <- create_ame_mock_fit(P = 2, N = 30, M = M, seed = 92)
    res <- quietly(ame(fit, n_draws = NULL))
    expect_equal(res$M_use, M)
})

test_that("H4. theta_draws_used dimensions match M_use x D", {
    M <- 500L
    P <- 3L
    fit <- create_ame_mock_fit(P = P, N = 50, M = M, seed = 93)

    res_all <- quietly(ame(fit))
    expect_equal(nrow(res_all$theta_draws_used), M)
    expect_equal(ncol(res_all$theta_draws_used), 2L * P + 1L)

    res_100 <- quietly(ame(fit, n_draws = 100L))
    expect_equal(nrow(res_100$theta_draws_used), res_100$M_use)
    expect_equal(ncol(res_100$theta_draws_used), 2L * P + 1L)
})

test_that("H5. AME draw matrices have M_use rows and P cols", {
    M <- 400L
    P <- 4L
    n_draws <- 50L
    fit <- create_ame_mock_fit(P = P, N = 50, M = M, seed = 94)
    res <- quietly(ame(fit, n_draws = n_draws))

    expect_equal(nrow(res$ext_ame_draws),   res$M_use)
    expect_equal(ncol(res$ext_ame_draws),   P)
    expect_equal(nrow(res$int_ame_draws),   res$M_use)
    expect_equal(ncol(res$int_ame_draws),   P)
    expect_equal(nrow(res$total_ame_draws), res$M_use)
    expect_equal(ncol(res$total_ame_draws), P)
})


# ============================================================================
# Section I: Point estimates (named P-vectors, total = ext + int)
# ============================================================================

test_that("I1. point estimates are named P-vectors", {
    P <- 4L
    fit <- create_ame_mock_fit(P = P, N = 50, M = 200, seed = 110)
    res <- quietly(ame(fit))

    expect_length(res$ame_ext_hat,   P)
    expect_length(res$ame_int_hat,   P)
    expect_length(res$ame_total_hat, P)

    expect_false(is.null(names(res$ame_ext_hat)))
    expect_false(is.null(names(res$ame_int_hat)))
    expect_false(is.null(names(res$ame_total_hat)))
})

test_that("I2. point estimate total = ext + int", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 200, seed = 111)
    res <- quietly(ame(fit))

    expect_equal(
        unname(res$ame_total_hat),
        unname(res$ame_ext_hat + res$ame_int_hat),
        tolerance = 1e-14
    )
})

test_that("I3. point estimates are close to posterior means", {
    fit <- create_ame_mock_fit(P = 3, N = 100, M = 500, seed = 112)
    res <- quietly(ame(fit))

    max_diff_ext <- max(abs(unname(res$ame_ext_hat) -
                                res$ext_summary$post_mean))
    max_diff_int <- max(abs(unname(res$ame_int_hat) -
                                res$int_summary$post_mean))

    # With small MCMC noise (Sigma = 0.01*I), point estimate and
    # posterior mean should be very close
    expect_true(max_diff_ext < 0.05,
                info = sprintf("Max ext diff: %.6f", max_diff_ext))
    expect_true(max_diff_int < 0.05,
                info = sprintf("Max int diff: %.6f", max_diff_int))
})

test_that("I4. point estimate names match cov_labels", {
    P <- 3L
    N <- 50L
    set.seed(1130)
    X <- cbind(intercept = 1, poverty = rnorm(N), urban = rnorm(N))

    fit <- create_ame_mock_fit(P = P, N = N, M = 200, seed = 113, X = X)
    res <- quietly(ame(fit))

    expect_equal(names(res$ame_ext_hat),   res$cov_labels)
    expect_equal(names(res$ame_int_hat),   res$cov_labels)
    expect_equal(names(res$ame_total_hat), res$cov_labels)
})


# ============================================================================
# Section J: Mean_q and mean_mu diagnostics (bounded, length, plogis)
# ============================================================================

test_that("J1. mean_q values are in (0, 1)", {
    fit <- create_ame_mock_fit(P = 3, N = 100, M = 200, seed = 120)
    res <- quietly(ame(fit))

    expect_true(all(res$mean_q > 0))
    expect_true(all(res$mean_q < 1))
})

test_that("J2. mean_mu values are in (0, 1)", {
    fit <- create_ame_mock_fit(P = 3, N = 100, M = 200, seed = 121)
    res <- quietly(ame(fit))

    expect_true(all(res$mean_mu > 0))
    expect_true(all(res$mean_mu < 1))
})

test_that("J3. mean_q and mean_mu have correct length M_use", {
    M <- 300L
    n_draws <- 75L
    fit <- create_ame_mock_fit(P = 3, N = 50, M = M, seed = 122)
    res <- quietly(ame(fit, n_draws = n_draws))

    expect_equal(length(res$mean_q),  res$M_use)
    expect_equal(length(res$mean_mu), res$M_use)
})

test_that("J4. mean_q consistent with plogis for intercept-only model", {
    P <- 1L
    N <- 100L
    alpha_true <- c(0.5)
    beta_true  <- c(-0.3)
    X <- matrix(1, nrow = N, ncol = 1)
    colnames(X) <- "intercept"

    fit <- create_ame_mock_fit(
        P = P, N = N, M = 50, seed = 123,
        alpha_true = alpha_true, beta_true = beta_true, X = X
    )
    res <- quietly(ame(fit))

    # Intercept-only => all obs have identical q = plogis(alpha)
    expected_q <- plogis(alpha_true[1])
    expect_equal(mean(res$mean_q), expected_q, tolerance = 0.05)

    expected_mu <- plogis(beta_true[1])
    expect_equal(mean(res$mean_mu), expected_mu, tolerance = 0.05)
})

test_that("J5. cov_labels are stored and used consistently", {
    P <- 3L
    N <- 50L
    set.seed(1240)
    X <- cbind(intercept = 1, poverty = rnorm(N), urban = rnorm(N))

    fit <- create_ame_mock_fit(P = P, N = N, M = 200, seed = 124, X = X)
    res <- quietly(ame(fit))

    expect_equal(res$cov_labels, c("intercept", "poverty", "urban"))
    expect_equal(colnames(res$ext_ame_draws), res$cov_labels)
    expect_equal(colnames(res$int_ame_draws), res$cov_labels)
    expect_equal(colnames(res$total_ame_draws), res$cov_labels)
    expect_equal(res$ext_summary$covariate, res$cov_labels)
    expect_equal(res$decomp_table$covariate, c("poverty", "urban"))
})

test_that("J6. level parameter affects CI width", {
    fit <- create_ame_mock_fit(P = 3, N = 50, M = 500, seed = 125)

    r90 <- quietly(ame(fit, level = 0.90))
    r95 <- quietly(ame(fit, level = 0.95))
    r99 <- quietly(ame(fit, level = 0.99))

    for (i in seq_len(r90$P)) {
        w90 <- r90$total_summary$ci_hi[i] - r90$total_summary$ci_lo[i]
        w95 <- r95$total_summary$ci_hi[i] - r95$total_summary$ci_lo[i]
        w99 <- r99$total_summary$ci_hi[i] - r99$total_summary$ci_lo[i]

        expect_true(w90 <= w95 + 1e-10,
                    info = paste("Width 90 > 95 for cov", i))
        expect_true(w95 <= w99 + 1e-10,
                    info = paste("Width 95 > 99 for cov", i))
    }
})

test_that("J7. class attribute is hbb_ame", {
    fit <- create_ame_mock_fit(P = 2, N = 30, M = 100, seed = 126)
    res <- quietly(ame(fit))
    expect_s3_class(res, "hbb_ame")
})

test_that("J8. CI_lo <= post_mean <= CI_hi in summaries", {
    fit <- create_ame_mock_fit(P = 3, N = 100, M = 2000, seed = 127)
    res <- quietly(ame(fit))

    for (summ in list(res$ext_summary, res$int_summary, res$total_summary)) {
        expect_true(all(summ$ci_lo <= summ$post_mean + 1e-10))
        expect_true(all(summ$ci_hi >= summ$post_mean - 1e-10))
    }
})
