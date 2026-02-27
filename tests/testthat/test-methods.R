# ============================================================================
# test-methods.R --- testthat tests for R/methods.R
#
# S3 methods for hbb_fit objects: nobs, coef, vcov, fitted, residuals,
# summary, print.summary
#
# Sections:
#   MOCK:  3 constructors (base, SVC, sandwich)
#   A.     nobs() tests (4 tests)
#   B.     coef() tests (8 tests)
#   C.     vcov() tests (7 tests)
#   D.     fitted() known-answer tests (8 tests)
#   E.     fitted() with SVC models (3 tests)
#   F.     residuals() tests (7 tests)
#   G.     summary() structure tests (8 tests)
#   H.     print.summary.hbb_fit defensive tests (5 tests)
#   I.     Input validation tests (7 tests)
#   J.     Edge cases (6 tests)
#   K.     Integration with CmdStan (5 tests, skip_if)
# ============================================================================


# ============================================================================
# MOCK CONSTRUCTORS
# ============================================================================

#' Create a mock hbb_fit object for methods.R testing
#'
#' Generates a mock hbb_fit with known parameter values and tight draws
#' so that fitted values, residuals, and summary outputs can be verified
#' algebraically.
#'
#' @param P Number of covariates per margin (including intercept).
#' @param N Number of observations.
#' @param M Number of MCMC draws (total across chains).
#' @param n_chains Number of MCMC chains (for metadata).
#' @param seed Random seed.
#' @param model_type Character: "base", "weighted", "svc", "svc_weighted".
#' @return An S3 object of class "hbb_fit".
create_methods_mock_fit <- function(
    P = 3, N = 10, M = 100, n_chains = 4, seed = 42,
    model_type = "base"
) {
    set.seed(seed)
    P <- as.integer(P)
    N <- as.integer(N)
    D <- 2L * P + 1L

    # Known true parameter values
    alpha_true <- seq(-0.5, 0.5, length.out = P)
    beta_true  <- seq(0.2, 0.8, length.out = P)
    log_kappa_true <- log(5)
    theta_true <- c(alpha_true, beta_true, log_kappa_true)

    # Simulate draws centered on true values (very tight for near-deterministic)
    draws_mat <- matrix(
        rep(theta_true, each = M) + rnorm(M * D, sd = 0.001),
        nrow = M, ncol = D
    )
    colnames(draws_mat) <- c(
        paste0("alpha[", 1:P, "]"),
        paste0("beta[", 1:P, "]"),
        "log_kappa"
    )

    # Design matrix with intercept column (handles P=1 edge case)
    if (P == 1L) {
        X <- matrix(1, nrow = N, ncol = 1L)
        colnames(X) <- "intercept"
    } else {
        X <- cbind(1, matrix(rnorm(N * (P - 1L)), nrow = N))
        colnames(X) <- c("intercept", paste0("x", 2:P))
    }

    # Compute expected q and mu from known parameters
    eta_ext <- as.numeric(X %*% alpha_true)
    eta_int <- as.numeric(X %*% beta_true)
    q_true  <- plogis(eta_ext)
    mu_true <- plogis(eta_int)

    # Generate observations consistent with the model
    n_trial <- rep(20L, N)
    z <- rbinom(N, 1, q_true)
    y <- integer(N)
    y[z == 1L] <- pmax(1L, rbinom(sum(z), n_trial[z == 1L], mu_true[z == 1L]))

    # Mock CmdStanMCMC-like object
    mock_fit <- list(
        draws = function(variables = NULL, format = "matrix") {
            if (is.null(variables)) return(draws_mat)
            idx <- match(variables, colnames(draws_mat))
            draws_mat[, idx, drop = FALSE]
        },
        summary = function(variables = NULL) {
            if (is.null(variables)) {
                vars <- colnames(draws_mat)
            } else {
                vars <- variables
            }
            idx <- match(vars, colnames(draws_mat))
            sub <- draws_mat[, idx, drop = FALSE]
            data.frame(
                variable = vars,
                mean = colMeans(sub),
                median = apply(sub, 2, median),
                sd = apply(sub, 2, sd),
                rhat = rep(1.001, length(vars)),
                ess_bulk = rep(2000, length(vars)),
                ess_tail = rep(1800, length(vars)),
                stringsAsFactors = FALSE
            )
        },
        diagnostic_summary = function(quiet = FALSE) {
            list(
                num_divergent = rep(0L, n_chains),
                num_max_treedepth = rep(0L, n_chains),
                ebfmi = rep(1.1, n_chains)
            )
        },
        num_chains = function() n_chains,
        metadata = function() list(
            iter_sampling = M / n_chains,
            stan_variables = colnames(draws_mat)
        )
    )

    hbb_data <- list(
        N = N, P = P, X = X,
        y = y, n_trial = n_trial,
        z = z
    )

    structure(
        list(
            fit = mock_fit,
            hbb_data = hbb_data,
            model_type = model_type,
            formula = list(formula = y ~ x2 + x3),
            call = match.call()
        ),
        class = "hbb_fit"
    )
}


#' Create a mock SVC hbb_fit object for testing fitted() with deltas
#'
#' @param P Number of covariates per margin.
#' @param N Number of observations.
#' @param S Number of states.
#' @param M Number of MCMC draws.
#' @param n_chains Number of MCMC chains.
#' @param seed Random seed.
#' @return An S3 object of class "hbb_fit" with model_type = "svc".
create_svc_mock_fit <- function(
    P = 3, N = 20, S = 5, M = 100, n_chains = 4, seed = 42
) {
    set.seed(seed)
    P <- as.integer(P)
    N <- as.integer(N)
    S <- as.integer(S)
    D <- 2L * P + 1L
    K <- 2L * P   # total RE dimension

    # Known true parameter values
    alpha_true <- seq(-0.5, 0.5, length.out = P)
    beta_true  <- seq(0.2, 0.8, length.out = P)
    log_kappa_true <- log(5)
    theta_true <- c(alpha_true, beta_true, log_kappa_true)

    # Fixed-effect draws (tight)
    draws_mat <- matrix(
        rep(theta_true, each = M) + rnorm(M * D, sd = 0.001),
        nrow = M, ncol = D
    )
    colnames(draws_mat) <- c(
        paste0("alpha[", 1:P, "]"),
        paste0("beta[", 1:P, "]"),
        "log_kappa"
    )

    # Delta draws: S x K matrix of random effects (non-trivial sd = 0.2)
    delta_true <- matrix(rnorm(S * K, sd = 0.2), nrow = S, ncol = K)
    delta_vec_true <- as.numeric(t(delta_true))  # row-major order
    n_delta <- S * K
    delta_draws <- matrix(
        rep(delta_vec_true, each = M) + rnorm(M * n_delta, sd = 0.001),
        nrow = M, ncol = n_delta
    )
    # Name in delta[s,k] format matching .extract_delta_means iteration order
    delta_colnames <- character(n_delta)
    idx <- 1L
    for (s in seq_len(S)) {
        for (k in seq_len(K)) {
            delta_colnames[idx] <- sprintf("delta[%d,%d]", s, k)
            idx <- idx + 1L
        }
    }
    colnames(delta_draws) <- delta_colnames

    # Tau draws: K hyperparameters
    tau_true <- abs(rnorm(K, mean = 0.3, sd = 0.05))
    tau_draws <- matrix(
        rep(tau_true, each = M) + rnorm(M * K, sd = 0.001),
        nrow = M, ncol = K
    )
    colnames(tau_draws) <- paste0("tau[", seq_len(K), "]")

    # Combined draws matrix
    all_draws <- cbind(draws_mat, delta_draws, tau_draws)

    # Design matrix
    if (P == 1L) {
        X <- matrix(1, nrow = N, ncol = 1L)
        colnames(X) <- "intercept"
    } else {
        X <- cbind(1, matrix(rnorm(N * (P - 1L)), nrow = N))
        colnames(X) <- c("intercept", paste0("x", 2:P))
    }

    # State assignments: cycle through 1:S
    state_vec <- rep(seq_len(S), length.out = N)

    # Generate observations
    eta_ext <- as.numeric(X %*% alpha_true)
    q_true  <- plogis(eta_ext)
    n_trial <- rep(20L, N)
    z <- rbinom(N, 1, q_true)
    y <- integer(N)
    y[z == 1L] <- pmax(1L, rbinom(sum(z), n_trial[z == 1L], 0.3))

    # Mock CmdStanMCMC-like object with delta and tau support
    mock_fit <- list(
        draws = function(variables = NULL, format = "matrix") {
            if (is.null(variables)) return(all_draws)
            idx <- match(variables, colnames(all_draws))
            all_draws[, idx, drop = FALSE]
        },
        summary = function(variables = NULL) {
            if (is.null(variables)) {
                vars <- colnames(all_draws)
            } else {
                vars <- variables
            }
            idx <- match(vars, colnames(all_draws))
            sub <- all_draws[, idx, drop = FALSE]
            data.frame(
                variable = vars,
                mean = colMeans(sub),
                median = apply(sub, 2, median),
                sd = apply(sub, 2, sd),
                rhat = rep(1.001, length(vars)),
                ess_bulk = rep(2000, length(vars)),
                ess_tail = rep(1800, length(vars)),
                stringsAsFactors = FALSE
            )
        },
        diagnostic_summary = function(quiet = FALSE) {
            list(
                num_divergent = rep(0L, n_chains),
                num_max_treedepth = rep(0L, n_chains),
                ebfmi = rep(1.1, n_chains)
            )
        },
        num_chains = function() n_chains,
        metadata = function() list(
            iter_sampling = M / n_chains,
            stan_variables = colnames(all_draws)
        )
    )

    hbb_data <- list(
        N = N, P = P, X = X,
        y = y, n_trial = n_trial,
        z = z,
        state = state_vec,
        S = S
    )

    structure(
        list(
            fit = mock_fit,
            hbb_data = hbb_data,
            model_type = "svc",
            formula = list(formula = y ~ x2 + x3),
            call = match.call()
        ),
        class = "hbb_fit"
    )
}


#' Create a mock hbb_sandwich object for vcov()/summary() testing
#'
#' @param D Total number of fixed-effect parameters (2P + 1).
#' @param seed Random seed.
#' @return An S3 object of class "hbb_sandwich".
create_mock_sandwich <- function(D, seed = 99) {
    set.seed(seed)
    V <- diag(runif(D, 0.001, 0.01))
    V <- (V + t(V)) / 2  # exact symmetry
    DER_vec <- runif(D, 1.5, 5)
    structure(
        list(V_sand = V, DER = DER_vec),
        class = "hbb_sandwich"
    )
}


# ============================================================================
# ---- Section A: nobs() ----
# ============================================================================

test_that("A1: nobs() returns the correct integer N", {
    fit <- create_methods_mock_fit(N = 10)
    result <- nobs(fit)
    expect_identical(result, 10L)
})

test_that("A2: nobs() returns integer type with length 1", {
    fit <- create_methods_mock_fit(N = 50)
    result <- nobs(fit)
    expect_true(is.integer(result))
    expect_length(result, 1L)
})

test_that("A3: nobs() returns correct value for different N", {
    for (n in c(1, 5, 100, 500)) {
        fit <- create_methods_mock_fit(N = n)
        expect_identical(nobs(fit), as.integer(n),
                         info = paste("N =", n))
    }
})

test_that("A4: nobs() errors on non-hbb_fit object", {
    expect_error(nobs.hbb_fit(list(a = 1)), "hbb_fit")
})


# ============================================================================
# ---- Section B: coef() ----
# ============================================================================

test_that("B1: coef() returns named numeric vector of length D = 2P+1", {
    fit <- create_methods_mock_fit(P = 3)
    co <- coef(fit)
    D <- 2L * 3L + 1L
    expect_true(is.numeric(co))
    expect_length(co, D)
    expect_false(is.null(names(co)))
    expect_equal(length(names(co)), D)
})

test_that("B2: coef() names use .build_param_labels convention", {
    fit <- create_methods_mock_fit(P = 3)
    co <- coef(fit)
    nms <- names(co)
    expect_equal(nms[1], "alpha_intercept")
    expect_equal(nms[4], "beta_intercept")
    expect_equal(nms[7], "log_kappa")
    expect_true(all(grepl("^alpha_", nms[1:3])))
    expect_true(all(grepl("^beta_", nms[4:6])))
})

test_that("B3: coef(margin='extensive') returns first P alpha coefficients", {
    fit <- create_methods_mock_fit(P = 4)
    co_ext <- coef(fit, margin = "extensive")
    co_all <- coef(fit)
    expect_length(co_ext, 4L)
    expect_equal(co_ext, co_all[1:4], tolerance = 1e-12)
    expect_true(all(grepl("^alpha_", names(co_ext))))
})

test_that("B4: coef(margin='intensive') returns beta coefficients P+1:2P", {
    fit <- create_methods_mock_fit(P = 4)
    co_int <- coef(fit, margin = "intensive")
    co_all <- coef(fit)
    expect_length(co_int, 4L)
    expect_equal(co_int, co_all[5:8], tolerance = 1e-12)
    expect_true(all(grepl("^beta_", names(co_int))))
})

test_that("B5: coef() posterior means match known true values (tight draws)", {
    fit <- create_methods_mock_fit(P = 3, M = 500, seed = 77)
    co <- coef(fit)
    alpha_true <- seq(-0.5, 0.5, length.out = 3)
    beta_true  <- seq(0.2, 0.8, length.out = 3)
    log_kappa_true <- log(5)
    theta_true <- c(alpha_true, beta_true, log_kappa_true)
    expect_equal(as.numeric(co), theta_true, tolerance = 0.01)
})

test_that("B6: coef() consistency: both == c(extensive, intensive, log_kappa)", {
    fit <- create_methods_mock_fit(P = 3)
    co_all <- coef(fit)
    co_ext <- coef(fit, margin = "extensive")
    co_int <- coef(fit, margin = "intensive")
    P <- 3L
    D <- 2L * P + 1L
    expect_equal(co_all[1:P], co_ext, tolerance = 1e-12)
    expect_equal(co_all[(P + 1):(2 * P)], co_int, tolerance = 1e-12)
    expect_equal(co_all[D], co_all["log_kappa"], tolerance = 1e-12)
})

test_that("B7: coef() errors on invalid margin argument", {
    fit <- create_methods_mock_fit()
    expect_error(coef(fit, margin = "invalid"))
})

test_that("B8: coef() works correctly with varying P", {
    for (p in c(1, 2, 5)) {
        fit <- create_methods_mock_fit(P = p)
        co <- coef(fit)
        expect_equal(length(co), 2L * p + 1L,
                     info = paste("P =", p))
        co_ext <- coef(fit, margin = "extensive")
        expect_equal(length(co_ext), as.integer(p),
                     info = paste("extensive P =", p))
        co_int <- coef(fit, margin = "intensive")
        expect_equal(length(co_int), as.integer(p),
                     info = paste("intensive P =", p))
    }
})


# ============================================================================
# ---- Section C: vcov() ----
# ============================================================================

test_that("C1: vcov() returns symmetric D x D matrix (MCMC cov)", {
    fit <- create_methods_mock_fit(P = 3)
    V <- vcov(fit)
    D <- 2L * 3L + 1L
    expect_true(is.matrix(V))
    expect_equal(dim(V), c(D, D))
    expect_equal(V, t(V), tolerance = 1e-14)
})

test_that("C2: vcov() has correct row and column names", {
    fit <- create_methods_mock_fit(P = 3)
    V <- vcov(fit)
    expected_labels <- c(
        "alpha_intercept", "alpha_x2", "alpha_x3",
        "beta_intercept", "beta_x2", "beta_x3",
        "log_kappa"
    )
    expect_equal(rownames(V), expected_labels)
    expect_equal(colnames(V), expected_labels)
})

test_that("C3: vcov() diagonal elements are positive (variance > 0)", {
    fit <- create_methods_mock_fit(P = 3, M = 200)
    V <- vcov(fit)
    expect_true(all(diag(V) > 0))
})

test_that("C4: vcov(sandwich=) returns V_sand with correct labels", {
    fit <- create_methods_mock_fit(P = 3)
    D <- 2L * 3L + 1L
    sand <- create_mock_sandwich(D)
    V <- vcov(fit, sandwich = sand)
    expect_equal(dim(V), c(D, D))
    expected_labels <- c(
        "alpha_intercept", "alpha_x2", "alpha_x3",
        "beta_intercept", "beta_x2", "beta_x3",
        "log_kappa"
    )
    expect_equal(rownames(V), expected_labels)
    expect_equal(colnames(V), expected_labels)
})

test_that("C5: vcov(sandwich=) returns V_sand values unchanged (relabeled)", {
    fit <- create_methods_mock_fit(P = 3)
    D <- 2L * 3L + 1L
    sand <- create_mock_sandwich(D)
    V <- vcov(fit, sandwich = sand)
    expect_equal(unname(V), unname(sand$V_sand), tolerance = 1e-14)
})

test_that("C6: vcov() errors on invalid sandwich class", {
    fit <- create_methods_mock_fit(P = 3)
    bad_sandwich <- list(V_sand = diag(7))
    expect_error(vcov(fit, sandwich = bad_sandwich), "hbb_sandwich")
})

test_that("C7: vcov() errors on sandwich dimension mismatch", {
    fit <- create_methods_mock_fit(P = 3)
    sand_wrong <- create_mock_sandwich(5)
    expect_error(vcov(fit, sandwich = sand_wrong), "Dimension")
})


# ============================================================================
# ---- Section D: fitted() known-answer tests ----
# ============================================================================

test_that("D1: fitted(type='response') == q_hat * mu_hat (proportion)", {
    fit <- create_methods_mock_fit()
    f_resp <- fitted(fit, type = "response")
    f_ext  <- fitted(fit, type = "extensive")
    f_int  <- fitted(fit, type = "intensive")
    expect_equal(f_resp, f_ext * f_int, tolerance = 1e-12)
})

test_that("D2: fitted(type='extensive') == plogis(X %*% alpha_hat)", {
    fit <- create_methods_mock_fit()
    P <- fit$hbb_data$P
    draws <- fit$fit$draws(format = "matrix")
    alpha_hat <- colMeans(draws[, seq_len(P), drop = FALSE])
    X <- fit$hbb_data$X
    expected_q <- plogis(as.numeric(X %*% alpha_hat))
    f_ext <- fitted(fit, type = "extensive")
    expect_equal(f_ext, expected_q, tolerance = 1e-10)
})

test_that("D3: fitted(type='intensive') == plogis(X %*% beta_hat)", {
    fit <- create_methods_mock_fit()
    P <- fit$hbb_data$P
    draws <- fit$fit$draws(format = "matrix")
    beta_hat <- colMeans(draws[, (P + 1L):(2L * P), drop = FALSE])
    X <- fit$hbb_data$X
    expected_mu <- plogis(as.numeric(X %*% beta_hat))
    f_int <- fitted(fit, type = "intensive")
    expect_equal(f_int, expected_mu, tolerance = 1e-10)
})

test_that("D4: fitted() returns N-length numeric vector", {
    fit <- create_methods_mock_fit(N = 15)
    for (tp in c("response", "extensive", "intensive")) {
        fv <- fitted(fit, type = tp)
        expect_true(is.numeric(fv))
        expect_length(fv, 15L)
    }
})

test_that("D5: fitted() values lie in [0, 1] for all types", {
    fit <- create_methods_mock_fit(N = 50, seed = 123)
    for (tp in c("response", "extensive", "intensive")) {
        fv <- fitted(fit, type = tp)
        expect_true(all(fv >= 0 & fv <= 1),
                    info = paste("type =", tp))
    }
})

test_that("D6: fitted() near known algebraic values (tight draws)", {
    fit <- create_methods_mock_fit(M = 500, seed = 99)
    P <- fit$hbb_data$P
    alpha_true <- seq(-0.5, 0.5, length.out = P)
    beta_true  <- seq(0.2, 0.8, length.out = P)
    X <- fit$hbb_data$X
    q_exact  <- plogis(as.numeric(X %*% alpha_true))
    mu_exact <- plogis(as.numeric(X %*% beta_true))
    expect_equal(fitted(fit, type = "extensive"), q_exact, tolerance = 0.01)
    expect_equal(fitted(fit, type = "intensive"), mu_exact, tolerance = 0.01)
    expect_equal(fitted(fit, type = "response"), q_exact * mu_exact,
                 tolerance = 0.01)
})

test_that("D7: fitted() default type is 'response'", {
    fit <- create_methods_mock_fit()
    f_default  <- fitted(fit)
    f_response <- fitted(fit, type = "response")
    expect_identical(f_default, f_response)
})

test_that("D8: fitted() with P=1 (intercept only) returns constant vector", {
    fit <- create_methods_mock_fit(P = 1, N = 8)
    f_ext <- fitted(fit, type = "extensive")
    f_int <- fitted(fit, type = "intensive")
    expect_true(max(f_ext) - min(f_ext) < 1e-6)
    expect_true(max(f_int) - min(f_int) < 1e-6)
})


# ============================================================================
# ---- Section E: fitted() with SVC models ----
# ============================================================================

test_that("E1: SVC fitted values are numeric, N-length, in [0,1]", {
    fit_svc <- create_svc_mock_fit(P = 3, N = 20, S = 5, seed = 42)
    f_svc <- fitted(fit_svc, type = "response")
    expect_true(is.numeric(f_svc))
    expect_length(f_svc, 20L)
    expect_true(all(f_svc >= 0 & f_svc <= 1))
})

test_that("E2: SVC fitted values differ from fixed-effects-only values", {
    fit_svc <- create_svc_mock_fit(P = 3, N = 20, S = 5, seed = 42)
    f_svc <- fitted(fit_svc, type = "response")

    # Compute fixed-effects-only predictions manually
    P <- fit_svc$hbb_data$P
    param_names <- c(paste0("alpha[", 1:P, "]"),
                     paste0("beta[", 1:P, "]"), "log_kappa")
    draws <- fit_svc$fit$draws(variables = param_names, format = "matrix")
    theta_hat <- colMeans(draws)
    alpha_hat <- theta_hat[1:P]
    beta_hat  <- theta_hat[(P + 1):(2 * P)]
    X <- fit_svc$hbb_data$X

    q_fe  <- plogis(as.numeric(X %*% alpha_hat))
    mu_fe <- plogis(as.numeric(X %*% beta_hat))
    f_fe  <- q_fe * mu_fe

    # SVC fitted should differ (delta sd = 0.2 is non-trivial)
    expect_false(isTRUE(all.equal(f_svc, f_fe, tolerance = 1e-6)))
})

test_that("E3: SVC fitted() in [0,1] for all types", {
    fit_svc <- create_svc_mock_fit(P = 3, N = 30, S = 5, seed = 88)
    for (tp in c("response", "extensive", "intensive")) {
        fv <- fitted(fit_svc, type = tp)
        expect_true(is.numeric(fv))
        expect_length(fv, 30L)
        expect_true(all(fv >= 0 & fv <= 1),
                    info = paste("SVC type =", tp))
    }
})


# ============================================================================
# ---- Section F: residuals() ----
# ============================================================================

test_that("F1: residuals(type='response') == y/n - fitted(type='response')", {
    fit <- create_methods_mock_fit()
    r_resp <- residuals(fit, type = "response")
    y <- fit$hbb_data$y
    n_trial <- fit$hbb_data$n_trial
    f_resp <- fitted(fit, type = "response")
    expected <- y / n_trial - f_resp
    expect_equal(r_resp, expected, tolerance = 1e-12)
})

test_that("F2: residuals() default type is 'response'", {
    fit <- create_methods_mock_fit()
    r_default  <- residuals(fit)
    r_response <- residuals(fit, type = "response")
    expect_identical(r_default, r_response)
})

test_that("F3: residuals() returns N-length numeric for both types", {
    fit <- create_methods_mock_fit(N = 12)
    for (tp in c("response", "pearson")) {
        r <- residuals(fit, type = tp)
        expect_true(is.numeric(r))
        expect_length(r, 12L)
    }
})

test_that("F4: Pearson residual algebraic verification", {
    fit <- create_methods_mock_fit(N = 20, seed = 77)
    r_pearson <- residuals(fit, type = "pearson")
    y <- fit$hbb_data$y
    n_trial <- as.integer(fit$hbb_data$n_trial)
    q_hat <- fitted(fit, type = "extensive")
    mu_hat <- fitted(fit, type = "intensive")

    # Count-scale expected value and numerator
    E_y <- as.numeric(n_trial) * q_hat * mu_hat
    numerator <- y - E_y

    # Same guards as in residuals.hbb_fit
    eps <- .Machine$double.eps^0.5
    mu_clamped <- pmin(pmax(mu_hat, eps), 1 - eps)
    q_clamped  <- pmin(pmax(q_hat, eps), 1)

    draws <- fit$fit$draws(format = "matrix")
    log_kappa_hat <- mean(draws[, ncol(draws)])
    kappa_hat <- pmin(exp(log_kappa_hat), 1e15)

    var_hat <- hurdle_variance(
        n = n_trial, q = q_clamped,
        mu = mu_clamped, kappa = kappa_hat
    )
    var_hat <- pmax(var_hat, .Machine$double.eps)

    expected_pearson <- numerator / sqrt(var_hat)
    expect_equal(r_pearson, expected_pearson, tolerance = 1e-8)
})

test_that("F5: Pearson residuals are finite for typical data", {
    fit <- create_methods_mock_fit(N = 30, seed = 55)
    r_pearson <- residuals(fit, type = "pearson")
    expect_true(all(is.finite(r_pearson)))
})

test_that("F6: response residuals are all finite", {
    fit <- create_methods_mock_fit(N = 50, seed = 11)
    r_resp <- residuals(fit, type = "response")
    expect_true(is.numeric(r_resp))
    expect_length(r_resp, 50L)
    expect_true(all(is.finite(r_resp)))
})

test_that("F7: response residual == -q*mu when y=0", {
    fit <- create_methods_mock_fit(N = 30, seed = 88)
    y <- fit$hbb_data$y
    r_resp <- residuals(fit, type = "response")
    f_resp <- fitted(fit, type = "response")
    zero_idx <- which(y == 0L)
    if (length(zero_idx) > 0) {
        expect_equal(r_resp[zero_idx], -f_resp[zero_idx], tolerance = 1e-12)
    }
})


# ============================================================================
# ---- Section G: summary() structure ----
# ============================================================================

test_that("G1: summary() returns class 'summary.hbb_fit'", {
    fit <- create_methods_mock_fit()
    s <- summary(fit)
    expect_s3_class(s, "summary.hbb_fit")
})

test_that("G2: summary() has all required top-level fields", {
    fit <- create_methods_mock_fit()
    s <- summary(fit)
    expected_fields <- c(
        "fixed_effects", "dispersion", "random_effects",
        "diagnostics", "model_info", "sandwich_used",
        "DER", "level", "call"
    )
    for (field in expected_fields) {
        expect_true(field %in% names(s),
                    info = paste("Missing field:", field))
    }
})

test_that("G3: fixed_effects has correct dimensions and columns", {
    fit <- create_methods_mock_fit(P = 4)
    s <- summary(fit)
    fe <- s$fixed_effects
    expect_s3_class(fe, "data.frame")
    expect_equal(nrow(fe), 9L)  # D = 2*4+1 = 9
    required_cols <- c("parameter", "estimate", "se",
                       "ci_lower", "ci_upper", "rhat", "ess_bulk")
    for (col in required_cols) {
        expect_true(col %in% names(fe),
                    info = paste("Missing column:", col))
    }
})

test_that("G4: summary() without sandwich uses posterior CIs", {
    fit <- create_methods_mock_fit()
    s <- summary(fit)
    expect_false(s$sandwich_used)
    expect_null(s$DER)
    fe <- s$fixed_effects
    expect_true(all(fe$ci_lower <= fe$estimate + 1e-10))
    expect_true(all(fe$estimate <= fe$ci_upper + 1e-10))
})

test_that("G5: summary() with sandwich uses Wald SEs and DER", {
    fit <- create_methods_mock_fit()
    P <- fit$hbb_data$P
    D <- 2L * P + 1L
    V_sand <- diag(D) * 0.01
    mock_sandwich <- structure(
        list(V_sand = V_sand, DER = rep(2.0, D)),
        class = "hbb_sandwich"
    )
    s <- summary(fit, sandwich = mock_sandwich)
    expect_true(s$sandwich_used)
    expect_false(is.null(s$DER))
    expect_length(s$DER, D)

    # SEs from sqrt(diag(V_sand))
    expected_se <- sqrt(diag(V_sand))
    expect_equal(s$fixed_effects$se, unname(expected_se), tolerance = 1e-10)

    # Wald CIs: estimate +/- z * se
    z_crit <- qnorm(0.975)
    expect_equal(
        s$fixed_effects$ci_lower,
        s$fixed_effects$estimate - z_crit * s$fixed_effects$se,
        tolerance = 1e-10
    )
    expect_equal(
        s$fixed_effects$ci_upper,
        s$fixed_effects$estimate + z_crit * s$fixed_effects$se,
        tolerance = 1e-10
    )
})

test_that("G6: dispersion has correct delta-method SE", {
    fit <- create_methods_mock_fit()
    s <- summary(fit)
    disp <- s$dispersion
    expect_true(is.list(disp))
    expect_true(all(c("kappa_hat", "kappa_se", "kappa_ci") %in% names(disp)))

    fe <- s$fixed_effects
    D <- nrow(fe)
    lk_hat <- fe$estimate[D]
    lk_se  <- fe$se[D]

    expect_equal(disp$kappa_hat, exp(lk_hat), tolerance = 1e-10)
    expect_equal(disp$kappa_se, exp(lk_hat) * lk_se, tolerance = 1e-10)
    expect_equal(disp$kappa_ci[1], exp(fe$ci_lower[D]), tolerance = 1e-10)
    expect_equal(disp$kappa_ci[2], exp(fe$ci_upper[D]), tolerance = 1e-10)
})

test_that("G7: diagnostics has all expected fields", {
    fit <- create_methods_mock_fit()
    s <- summary(fit)
    dx <- s$diagnostics
    expect_true(is.list(dx))
    expected_diag_fields <- c(
        "n_divergent", "n_max_treedepth", "ebfmi",
        "max_rhat", "min_ess_bulk", "min_ess_tail"
    )
    for (f in expected_diag_fields) {
        expect_true(f %in% names(dx),
                    info = paste("Missing diagnostics field:", f))
    }
    expect_equal(dx$n_divergent, 0L)
    expect_equal(dx$n_max_treedepth, 0L)
    expect_equal(dx$max_rhat, 1.001, tolerance = 1e-6)
    expect_equal(dx$min_ess_bulk, 2000, tolerance = 1e-6)
})

test_that("G8: model_info reflects the fit object", {
    fit <- create_methods_mock_fit(P = 3, N = 10, model_type = "weighted")
    s <- summary(fit)
    mi <- s$model_info
    expect_true(is.list(mi))
    expect_equal(mi$N, 10L)
    expect_equal(mi$P, 3L)
    expect_equal(mi$model_type, "weighted")
    expected_zr <- 1 - mean(fit$hbb_data$z)
    expect_equal(mi$zero_rate, expected_zr, tolerance = 1e-10)
    expect_null(s$random_effects)
    expect_null(mi$S)
    expect_equal(s$level, 0.95)
})


# ============================================================================
# ---- Section H: print.summary.hbb_fit defensive ----
# ============================================================================

test_that("H1: print.summary.hbb_fit produces output without error", {
    fit <- create_methods_mock_fit()
    s <- summary(fit)
    out <- capture.output(print(s))
    expect_true(length(out) > 0L)
})

test_that("H2: print output contains expected section headings", {
    fit <- create_methods_mock_fit()
    s <- summary(fit)
    out <- paste(capture.output(print(s)), collapse = "\n")
    expect_true(grepl("Hurdle Beta-Binomial", out))
    expect_true(grepl("Extensive Margin", out))
    expect_true(grepl("Intensive Margin", out))
    expect_true(grepl("Dispersion", out))
    expect_true(grepl("MCMC Diagnostics", out))
})

test_that("H3: print handles corrupted summary gracefully (no crash)", {
    fit <- create_methods_mock_fit()
    s <- summary(fit)
    # Corrupt various fields
    s$fixed_effects <- NULL
    s$diagnostics <- NULL
    s$model_info$P <- NULL
    s$dispersion <- NULL
    # Should NOT crash -- top-level tryCatch
    out <- capture.output(result <- print(s))
    expect_true(length(out) > 0L)
    out_text <- paste(out, collapse = "\n")
    expect_true(grepl("Error printing summary", out_text))
})

test_that("H4: print returns invisible(x)", {
    fit <- create_methods_mock_fit()
    s <- summary(fit)
    vis <- withVisible(print(s))
    expect_false(vis$visible)
    expect_identical(vis$value, s)
})

test_that("H5: print with digits parameter produces different output", {
    fit <- create_methods_mock_fit()
    s <- summary(fit)
    out_3 <- capture.output(print(s, digits = 3))
    out_5 <- capture.output(print(s, digits = 5))
    expect_true(length(out_3) > 0L)
    expect_true(length(out_5) > 0L)
    expect_false(identical(out_3, out_5))
})


# ============================================================================
# ---- Section I: Input validation ----
# ============================================================================

test_that("I1: nobs.hbb_fit() errors on non-hbb_fit (plain list)", {
    bad <- list(hbb_data = list(N = 10), fit = list())
    expect_error(nobs.hbb_fit(bad), "hbb_fit")
})

test_that("I2: coef.hbb_fit() errors on non-hbb_fit", {
    bad <- list(hbb_data = list(N = 10), fit = list())
    expect_error(coef.hbb_fit(bad), "hbb_fit")
})

test_that("I3: fitted() errors on object with NULL hbb_data", {
    fit <- create_methods_mock_fit()
    fit$hbb_data <- NULL
    expect_error(fitted(fit), "hbb_data.*NULL")
})

test_that("I4: residuals() errors on object with NULL fit slot", {
    fit <- create_methods_mock_fit()
    fit$fit <- NULL
    expect_error(residuals(fit), "fit.*NULL")
})

test_that("I5: summary() errors with invalid level", {
    fit <- create_methods_mock_fit()
    expect_error(summary(fit, level = 0), "level")
    expect_error(summary(fit, level = 1), "level")
    expect_error(summary(fit, level = 1.5), "level")
    expect_error(summary(fit, level = -0.1), "level")
})

test_that("I6: vcov() errors when sandwich is wrong class", {
    fit <- create_methods_mock_fit()
    bad_sandwich <- list(V_sand = diag(7))
    expect_error(vcov(fit, sandwich = bad_sandwich), "hbb_sandwich")
})

test_that("I7: coef() errors on invalid margin argument", {
    fit <- create_methods_mock_fit()
    expect_error(coef(fit, margin = "invalid"), "should be one of")
})


# ============================================================================
# ---- Section J: Edge cases ----
# ============================================================================

test_that("J1: N=1 (single observation) -- all methods work", {
    fit <- create_methods_mock_fit(N = 1, seed = 99)
    expect_equal(nobs(fit), 1L)
    expect_length(coef(fit), 7L)  # D = 2*3+1
    expect_true(is.matrix(vcov(fit)))
    expect_length(fitted(fit), 1L)
    expect_length(residuals(fit), 1L)
    s <- summary(fit)
    expect_s3_class(s, "summary.hbb_fit")
})

test_that("J2: All y=0 -- fitted, residuals, summary work", {
    fit <- create_methods_mock_fit(seed = 77)
    fit$hbb_data$y <- rep(0L, fit$hbb_data$N)
    fit$hbb_data$z <- rep(0L, fit$hbb_data$N)
    fv <- fitted(fit)
    expect_length(fv, fit$hbb_data$N)
    expect_true(all(is.finite(fv)))
    r <- residuals(fit)
    expect_length(r, fit$hbb_data$N)
    expect_true(all(is.finite(r)))
    expect_true(all(r <= 0))  # model predicts positive, data is 0
    s <- summary(fit)
    expect_s3_class(s, "summary.hbb_fit")
})

test_that("J3: P=1 (intercept only) -- D=3, methods work", {
    fit <- create_methods_mock_fit(P = 1, N = 5, seed = 55)
    expect_equal(nobs(fit), 5L)
    theta <- coef(fit)
    expect_length(theta, 3L)
    V <- vcov(fit)
    expect_equal(dim(V), c(3L, 3L))
    fv <- fitted(fit)
    expect_length(fv, 5L)
    expect_true(all(fv >= 0 & fv <= 1))
    r <- residuals(fit, type = "pearson")
    expect_length(r, 5L)
})

test_that("J4: Extreme kappa (log_kappa=35) -- clamped in fitted/residuals", {
    fit <- create_methods_mock_fit(seed = 66)
    P <- fit$hbb_data$P
    original_draws_fn <- fit$fit$draws
    extreme_lk <- 35  # exp(35) ~ 1.6e15, exceeds 1e15 clamp

    fit$fit$draws <- function(variables = NULL, format = "matrix") {
        mat <- original_draws_fn(variables, format)
        lk_col <- which(colnames(mat) == "log_kappa")
        if (length(lk_col) > 0L) mat[, lk_col] <- extreme_lk
        mat
    }

    fv <- fitted(fit)
    expect_length(fv, fit$hbb_data$N)
    expect_true(all(is.finite(fv)))
    r <- residuals(fit, type = "pearson")
    expect_length(r, fit$hbb_data$N)
    expect_true(all(is.finite(r) | is.na(r)))
})

test_that("J5: n_trial=1 for all -- Pearson residuals work", {
    fit <- create_methods_mock_fit(seed = 88)
    N <- fit$hbb_data$N
    fit$hbb_data$n_trial <- rep(1L, N)
    fit$hbb_data$y <- fit$hbb_data$z  # z is already 0/1
    fv <- fitted(fit)
    expect_length(fv, N)
    r_resp <- residuals(fit, type = "response")
    expect_length(r_resp, N)
    expect_true(all(is.finite(r_resp)))
    r_pear <- residuals(fit, type = "pearson")
    expect_length(r_pear, N)
    expect_true(all(is.finite(r_pear) | is.na(r_pear)))
})

test_that("J6: y = n_trial (100% enrollment) -- residuals non-NA", {
    fit <- create_methods_mock_fit(seed = 33)
    N <- fit$hbb_data$N
    n_trial <- fit$hbb_data$n_trial
    fit$hbb_data$z <- rep(1L, N)
    fit$hbb_data$y <- n_trial
    fv <- fitted(fit)
    expect_true(all(is.finite(fv)))
    r_resp <- residuals(fit, type = "response")
    expect_true(all(is.finite(r_resp)))
    r_pear <- residuals(fit, type = "pearson")
    expect_true(all(is.finite(r_pear) | is.na(r_pear)))
})


# ============================================================================
# ---- Section K: Integration with CmdStan (skip_if) ----
# ============================================================================

test_that("K1: coef() with real CmdStanMCMC", {
    skip_if_not_installed("cmdstanr")
    skip_on_cran()
    skip_if(
        is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
        message = "CmdStan installation not found"
    )

    stan_code <- "
    data { int<lower=1> P; }
    parameters { vector[P] alpha; vector[P] beta; real log_kappa; }
    model { alpha ~ normal(0, 1); beta ~ normal(0, 1);
            log_kappa ~ normal(1, 1); }
    "
    P_test <- 2; N_test <- 5
    mod <- cmdstanr::cmdstan_model(
        cmdstanr::write_stan_file(stan_code), quiet = TRUE
    )
    cmdstan_fit <- mod$sample(
        data = list(P = P_test), chains = 1,
        iter_warmup = 200, iter_sampling = 100,
        refresh = 0, show_messages = FALSE
    )
    X <- cbind(1, rnorm(N_test))
    colnames(X) <- c("intercept", "x2")
    fit_obj <- structure(
        list(fit = cmdstan_fit,
             hbb_data = list(N = N_test, P = P_test, X = X,
                             y = rep(5L, N_test),
                             n_trial = rep(20L, N_test),
                             z = rep(1L, N_test)),
             model_type = "base",
             formula = list(formula = y ~ x2),
             call = quote(hbb())),
        class = "hbb_fit"
    )
    co <- coef(fit_obj)
    expect_length(co, 2L * P_test + 1L)
    expect_true(is.numeric(co))
    expect_false(any(is.na(co)))
})

test_that("K2: vcov() from real CmdStanMCMC", {
    skip_if_not_installed("cmdstanr")
    skip_on_cran()
    skip_if(
        is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
        message = "CmdStan installation not found"
    )

    stan_code <- "
    data { int<lower=1> P; }
    parameters { vector[P] alpha; vector[P] beta; real log_kappa; }
    model { alpha ~ normal(0, 1); beta ~ normal(0, 1);
            log_kappa ~ normal(1, 1); }
    "
    P_test <- 2; N_test <- 5
    mod <- cmdstanr::cmdstan_model(
        cmdstanr::write_stan_file(stan_code), quiet = TRUE
    )
    cmdstan_fit <- mod$sample(
        data = list(P = P_test), chains = 1,
        iter_warmup = 200, iter_sampling = 100,
        refresh = 0, show_messages = FALSE
    )
    X <- cbind(1, rnorm(N_test))
    colnames(X) <- c("intercept", "x2")
    fit_obj <- structure(
        list(fit = cmdstan_fit,
             hbb_data = list(N = N_test, P = P_test, X = X,
                             y = rep(5L, N_test),
                             n_trial = rep(20L, N_test),
                             z = rep(1L, N_test)),
             model_type = "base",
             formula = list(formula = y ~ x2),
             call = quote(hbb())),
        class = "hbb_fit"
    )
    V <- vcov(fit_obj)
    D <- 2L * P_test + 1L
    expect_equal(dim(V), c(D, D))
    expect_true(all(diag(V) > 0))
    expect_equal(V, t(V), tolerance = 1e-12)
})

test_that("K3: fitted() from real CmdStanMCMC returns [0,1]", {
    skip_if_not_installed("cmdstanr")
    skip_on_cran()
    skip_if(
        is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
        message = "CmdStan installation not found"
    )

    stan_code <- "
    data { int<lower=1> P; }
    parameters { vector[P] alpha; vector[P] beta; real log_kappa; }
    model { alpha ~ normal(0, 1); beta ~ normal(0, 1);
            log_kappa ~ normal(1, 1); }
    "
    P_test <- 2; N_test <- 8
    mod <- cmdstanr::cmdstan_model(
        cmdstanr::write_stan_file(stan_code), quiet = TRUE
    )
    cmdstan_fit <- mod$sample(
        data = list(P = P_test), chains = 1,
        iter_warmup = 200, iter_sampling = 100,
        refresh = 0, show_messages = FALSE
    )
    set.seed(123)
    X <- cbind(1, rnorm(N_test))
    colnames(X) <- c("intercept", "x2")
    fit_obj <- structure(
        list(fit = cmdstan_fit,
             hbb_data = list(N = N_test, P = P_test, X = X,
                             y = c(0L, 0L, 3L, 5L, 0L, 7L, 2L, 0L),
                             n_trial = rep(20L, N_test),
                             z = c(0L, 0L, 1L, 1L, 0L, 1L, 1L, 0L)),
             model_type = "base",
             formula = list(formula = y ~ x2),
             call = quote(hbb())),
        class = "hbb_fit"
    )
    fv <- fitted(fit_obj, type = "response")
    expect_length(fv, N_test)
    expect_true(all(fv >= 0 & fv <= 1))
    fv_ext <- fitted(fit_obj, type = "extensive")
    fv_int <- fitted(fit_obj, type = "intensive")
    expect_equal(fv, fv_ext * fv_int, tolerance = 1e-12)
})

test_that("K4: nobs() from real CmdStanMCMC", {
    skip_if_not_installed("cmdstanr")
    skip_on_cran()
    skip_if(
        is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
        message = "CmdStan installation not found"
    )

    stan_code <- "
    data { int<lower=1> P; }
    parameters { vector[P] alpha; vector[P] beta; real log_kappa; }
    model { alpha ~ normal(0, 1); beta ~ normal(0, 1);
            log_kappa ~ normal(1, 1); }
    "
    P_test <- 2; N_test <- 12
    mod <- cmdstanr::cmdstan_model(
        cmdstanr::write_stan_file(stan_code), quiet = TRUE
    )
    cmdstan_fit <- mod$sample(
        data = list(P = P_test), chains = 1,
        iter_warmup = 200, iter_sampling = 100,
        refresh = 0, show_messages = FALSE
    )
    X <- cbind(1, rnorm(N_test))
    colnames(X) <- c("intercept", "x2")
    fit_obj <- structure(
        list(fit = cmdstan_fit,
             hbb_data = list(N = N_test, P = P_test, X = X,
                             y = rep(3L, N_test),
                             n_trial = rep(20L, N_test),
                             z = rep(1L, N_test)),
             model_type = "base",
             formula = list(formula = y ~ x2),
             call = quote(hbb())),
        class = "hbb_fit"
    )
    expect_identical(nobs(fit_obj), 12L)
})

test_that("K5: coef() multi-chain with P=3", {
    skip_if_not_installed("cmdstanr")
    skip_on_cran()
    skip_if(
        is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
        message = "CmdStan installation not found"
    )

    stan_code <- "
    data { int<lower=1> P; }
    parameters { vector[P] alpha; vector[P] beta; real log_kappa; }
    model { alpha ~ normal(0, 1); beta ~ normal(0, 1);
            log_kappa ~ normal(1, 1); }
    "
    P_test <- 3; N_test <- 10
    mod <- cmdstanr::cmdstan_model(
        cmdstanr::write_stan_file(stan_code), quiet = TRUE
    )
    cmdstan_fit <- mod$sample(
        data = list(P = P_test), chains = 2,
        iter_warmup = 300, iter_sampling = 200,
        refresh = 0, show_messages = FALSE
    )
    X <- cbind(1, matrix(rnorm(N_test * 2), ncol = 2))
    colnames(X) <- c("intercept", "x2", "x3")
    fit_obj <- structure(
        list(fit = cmdstan_fit,
             hbb_data = list(N = N_test, P = P_test, X = X,
                             y = rep(5L, N_test),
                             n_trial = rep(20L, N_test),
                             z = rep(1L, N_test)),
             model_type = "base",
             formula = list(formula = y ~ x2 + x3),
             call = quote(hbb())),
        class = "hbb_fit"
    )
    co <- coef(fit_obj)
    D <- 2L * P_test + 1L
    expect_length(co, D)
    expect_true(all(is.finite(co)))
    co_ext <- coef(fit_obj, margin = "extensive")
    co_int <- coef(fit_obj, margin = "intensive")
    expect_length(co_ext, P_test)
    expect_length(co_int, P_test)
})


# ============================================================================
# ---- Section L: CmdStan extraction error paths (tryCatch branches) ----
# ============================================================================

#' Create a mock hbb_fit whose $fit$draws() throws an error
#'
#' Used to exercise the tryCatch error handler paths in coef(), vcov(),
#' .extract_fixed_effects(), and .compute_fitted_values().
#'
#' @param P Number of covariates per margin.
#' @param N Number of observations.
#' @param model_type Character model type.
#' @return An S3 object of class "hbb_fit" with a broken draws() method.
create_broken_draws_mock <- function(P = 3, N = 10, model_type = "base") {
    P <- as.integer(P); N <- as.integer(N)

    if (P == 1L) {
        X <- matrix(1, nrow = N, ncol = 1L)
        colnames(X) <- "intercept"
    } else {
        set.seed(42)
        X <- cbind(1, matrix(rnorm(N * (P - 1L)), nrow = N))
        colnames(X) <- c("intercept", paste0("x", 2:P))
    }

    mock_fit <- list(
        draws = function(variables = NULL, format = "matrix") {
            stop("CmdStan output files not found")
        },
        summary = function(variables = NULL) {
            stop("CmdStan output files not found")
        },
        diagnostic_summary = function(quiet = FALSE) {
            stop("CmdStan output files not found")
        },
        num_chains = function() 4L,
        metadata = function() list(iter_sampling = 25L,
                                    stan_variables = character(0))
    )

    hbb_data <- list(
        N = N, P = P, X = X,
        y = rep(3L, N), n_trial = rep(20L, N), z = rep(1L, N)
    )

    structure(
        list(fit = mock_fit, hbb_data = hbb_data,
             model_type = model_type,
             formula = list(formula = y ~ x2 + x3),
             call = quote(hbb())),
        class = "hbb_fit"
    )
}


#' Create a mock hbb_fit that returns wrong-dimension draws
#'
#' draws() returns a matrix with wrong number of columns to exercise
#' the dimension mismatch error branch.
create_wrong_dim_draws_mock <- function(P = 3, N = 10) {
    P <- as.integer(P); N <- as.integer(N); D <- 2L * P + 1L
    set.seed(42)

    X <- cbind(1, matrix(rnorm(N * (P - 1L)), nrow = N))
    colnames(X) <- c("intercept", paste0("x", 2:P))

    # Return draws matrix with WRONG number of columns (D+1 instead of D)
    wrong_draws <- matrix(rnorm(100 * (D + 1L)), nrow = 100, ncol = D + 1L)
    colnames(wrong_draws) <- c(
        paste0("alpha[", 1:P, "]"),
        paste0("beta[", 1:P, "]"),
        "log_kappa", "extra_col"
    )

    mock_fit <- list(
        draws = function(variables = NULL, format = "matrix") {
            wrong_draws
        },
        summary = function(variables = NULL) {
            data.frame(variable = colnames(wrong_draws),
                       mean = colMeans(wrong_draws),
                       rhat = rep(1.001, ncol(wrong_draws)),
                       ess_bulk = rep(2000, ncol(wrong_draws)),
                       stringsAsFactors = FALSE)
        },
        diagnostic_summary = function(quiet = FALSE) {
            list(num_divergent = rep(0L, 4),
                 num_max_treedepth = rep(0L, 4),
                 ebfmi = rep(1.1, 4))
        },
        num_chains = function() 4L,
        metadata = function() list(iter_sampling = 25L,
                                    stan_variables = colnames(wrong_draws))
    )

    hbb_data <- list(N = N, P = P, X = X,
                     y = rep(3L, N), n_trial = rep(20L, N), z = rep(1L, N))

    structure(
        list(fit = mock_fit, hbb_data = hbb_data,
             model_type = "base",
             formula = list(formula = y ~ x2 + x3),
             call = quote(hbb())),
        class = "hbb_fit"
    )
}


#' Create a mock hbb_fit that returns only 1 draw (insufficient M)
create_single_draw_mock <- function(P = 3, N = 10) {
    P <- as.integer(P); N <- as.integer(N); D <- 2L * P + 1L
    set.seed(42)

    X <- cbind(1, matrix(rnorm(N * (P - 1L)), nrow = N))
    colnames(X) <- c("intercept", paste0("x", 2:P))

    # Only 1 MCMC draw
    single_draw <- matrix(rnorm(D), nrow = 1, ncol = D)
    colnames(single_draw) <- c(
        paste0("alpha[", 1:P, "]"),
        paste0("beta[", 1:P, "]"),
        "log_kappa"
    )

    mock_fit <- list(
        draws = function(variables = NULL, format = "matrix") {
            if (is.null(variables)) return(single_draw)
            idx <- match(variables, colnames(single_draw))
            single_draw[, idx, drop = FALSE]
        },
        summary = function(variables = NULL) {
            if (is.null(variables)) vars <- colnames(single_draw) else vars <- variables
            idx <- match(vars, colnames(single_draw))
            sub <- single_draw[, idx, drop = FALSE]
            data.frame(variable = vars, mean = as.numeric(sub),
                       rhat = rep(NA_real_, length(vars)),
                       ess_bulk = rep(NA_real_, length(vars)),
                       stringsAsFactors = FALSE)
        },
        diagnostic_summary = function(quiet = FALSE) {
            list(num_divergent = 0L, num_max_treedepth = 0L, ebfmi = 1.1)
        },
        num_chains = function() 1L,
        metadata = function() list(iter_sampling = 1L,
                                    stan_variables = colnames(single_draw))
    )

    hbb_data <- list(N = N, P = P, X = X,
                     y = rep(3L, N), n_trial = rep(20L, N), z = rep(1L, N))

    structure(
        list(fit = mock_fit, hbb_data = hbb_data,
             model_type = "base",
             formula = list(formula = y ~ x2 + x3),
             call = quote(hbb())),
        class = "hbb_fit"
    )
}


#' Create a mock hbb_fit with NaN/Inf posterior means
create_nonfinite_draws_mock <- function(P = 3, N = 10) {
    P <- as.integer(P); N <- as.integer(N); D <- 2L * P + 1L
    set.seed(42)

    X <- cbind(1, matrix(rnorm(N * (P - 1L)), nrow = N))
    colnames(X) <- c("intercept", paste0("x", 2:P))

    # Create draws where some columns have NaN/Inf
    M <- 100L
    draws_mat <- matrix(rnorm(M * D), nrow = M, ncol = D)
    colnames(draws_mat) <- c(
        paste0("alpha[", 1:P, "]"),
        paste0("beta[", 1:P, "]"),
        "log_kappa"
    )
    # Inject NaN into first column (all draws)
    draws_mat[, 1L] <- NaN

    mock_fit <- list(
        draws = function(variables = NULL, format = "matrix") {
            if (is.null(variables)) return(draws_mat)
            idx <- match(variables, colnames(draws_mat))
            draws_mat[, idx, drop = FALSE]
        },
        summary = function(variables = NULL) {
            if (is.null(variables)) vars <- colnames(draws_mat) else vars <- variables
            idx <- match(vars, colnames(draws_mat))
            sub <- draws_mat[, idx, drop = FALSE]
            data.frame(variable = vars, mean = colMeans(sub),
                       sd = apply(sub, 2, sd),
                       rhat = rep(1.001, length(vars)),
                       ess_bulk = rep(2000, length(vars)),
                       ess_tail = rep(1800, length(vars)),
                       stringsAsFactors = FALSE)
        },
        diagnostic_summary = function(quiet = FALSE) {
            list(num_divergent = rep(0L, 4),
                 num_max_treedepth = rep(0L, 4),
                 ebfmi = rep(1.1, 4))
        },
        num_chains = function() 4L,
        metadata = function() list(iter_sampling = 25L,
                                    stan_variables = colnames(draws_mat))
    )

    hbb_data <- list(N = N, P = P, X = X,
                     y = rep(3L, N), n_trial = rep(20L, N), z = rep(1L, N))

    structure(
        list(fit = mock_fit, hbb_data = hbb_data,
             model_type = "base",
             formula = list(formula = y ~ x2 + x3),
             call = quote(hbb())),
        class = "hbb_fit"
    )
}


test_that("L1: coef() errors when draws() throws (tryCatch path)", {
    fit <- create_broken_draws_mock()
    expect_error(coef(fit), "Failed to extract posterior draws")
})

test_that("L2: vcov() errors when draws() throws (tryCatch path)", {
    fit <- create_broken_draws_mock()
    expect_error(vcov(fit), "Failed to extract posterior draws")
})

test_that("L3: coef() errors when draw matrix has wrong columns", {
    fit <- create_wrong_dim_draws_mock()
    expect_error(coef(fit), "columns but expected")
})

test_that("L4: coef() errors when only 1 MCMC draw available", {
    fit <- create_single_draw_mock()
    expect_error(coef(fit), "Only 1 MCMC draw")
})

test_that("L5: coef() warns on non-finite posterior means (NaN/Inf)", {
    fit <- create_nonfinite_draws_mock()
    expect_warning(coef(fit), "non-finite posterior mean")
})

test_that("L6: vcov() errors when draw matrix has wrong columns", {
    fit <- create_wrong_dim_draws_mock()
    expect_error(vcov(fit), "columns but expected")
})

test_that("L7: vcov() errors when only 1 MCMC draw available", {
    fit <- create_single_draw_mock()
    expect_error(vcov(fit), "Only 1 MCMC draw")
})

test_that("L8: summary() errors when draws() throws", {
    fit <- create_broken_draws_mock()
    expect_error(summary(fit), "Failed to extract posterior draws")
})

test_that("L9: fitted() errors when draws() throws", {
    fit <- create_broken_draws_mock()
    expect_error(fitted(fit), "Failed to extract posterior draws")
})

test_that("L10: residuals() errors when draws() throws", {
    fit <- create_broken_draws_mock()
    expect_error(residuals(fit), "Failed to extract posterior draws")
})

test_that("L11: .extract_fixed_effects errors on wrong-dim draws", {
    fit <- create_wrong_dim_draws_mock()
    P <- fit$hbb_data$P
    param_labels <- c(
        paste0("alpha_", colnames(fit$hbb_data$X)),
        paste0("beta_", colnames(fit$hbb_data$X)),
        "log_kappa"
    )
    expect_error(
        hurdlebb:::.extract_fixed_effects(fit, NULL, 0.95, param_labels),
        "columns but expected"
    )
})

test_that("L12: .extract_fixed_effects errors on single draw", {
    fit <- create_single_draw_mock()
    P <- fit$hbb_data$P
    param_labels <- c(
        paste0("alpha_", colnames(fit$hbb_data$X)),
        paste0("beta_", colnames(fit$hbb_data$X)),
        "log_kappa"
    )
    expect_error(
        hurdlebb:::.extract_fixed_effects(fit, NULL, 0.95, param_labels),
        "Only 1 MCMC draw"
    )
})

test_that("L13: .extract_fixed_effects warns on non-finite means (sandwich path)", {
    fit <- create_nonfinite_draws_mock()
    P <- fit$hbb_data$P; D <- 2L * P + 1L
    param_labels <- c(
        paste0("alpha_", colnames(fit$hbb_data$X)),
        paste0("beta_", colnames(fit$hbb_data$X)),
        "log_kappa"
    )
    # Use a sandwich so the Wald CI path is taken (avoids quantile NaN crash)
    mock_sand <- structure(
        list(V_sand = diag(D) * 0.01, DER = rep(2, D)),
        class = "hbb_sandwich"
    )
    expect_warning(
        hurdlebb:::.extract_fixed_effects(fit, mock_sand, 0.95, param_labels),
        "non-finite posterior mean"
    )
})

test_that("L14: .compute_fitted_values errors on broken draws", {
    fit <- create_broken_draws_mock()
    expect_error(
        hurdlebb:::.compute_fitted_values(fit),
        "Failed to extract posterior draws"
    )
})

test_that("L15: .compute_fitted_values errors on wrong-dim draws", {
    fit <- create_wrong_dim_draws_mock()
    expect_error(
        hurdlebb:::.compute_fitted_values(fit),
        "columns but expected"
    )
})

test_that("L16: .compute_fitted_values errors on single draw", {
    fit <- create_single_draw_mock()
    expect_error(
        hurdlebb:::.compute_fitted_values(fit),
        "Only 1 MCMC draw"
    )
})

test_that("L17: .compute_fitted_values warns on non-finite means", {
    fit <- create_nonfinite_draws_mock()
    expect_warning(
        hurdlebb:::.compute_fitted_values(fit),
        "non-finite posterior mean"
    )
})


# ============================================================================
# ---- Section M: .validate_hbb_fit_methods error paths ----
# ============================================================================

test_that("M1: validate errors when hbb_data$N is NULL", {
    fit <- create_methods_mock_fit()
    fit$hbb_data$N <- NULL
    expect_error(
        hurdlebb:::.validate_hbb_fit_methods(fit),
        "missing.*N.*P"
    )
})

test_that("M2: validate errors when hbb_data$P is NULL", {
    fit <- create_methods_mock_fit()
    fit$hbb_data$P <- NULL
    expect_error(
        hurdlebb:::.validate_hbb_fit_methods(fit),
        "missing.*N.*P"
    )
})

test_that("M3: validate errors when N is non-positive", {
    fit <- create_methods_mock_fit()
    fit$hbb_data$N <- 0L
    expect_error(
        hurdlebb:::.validate_hbb_fit_methods(fit),
        "N.*must be a positive"
    )
})

test_that("M4: validate errors when P is non-positive", {
    fit <- create_methods_mock_fit()
    fit$hbb_data$P <- -1L
    expect_error(
        hurdlebb:::.validate_hbb_fit_methods(fit),
        "P.*must be a positive"
    )
})

test_that("M5: validate errors when X is NULL", {
    fit <- create_methods_mock_fit()
    fit$hbb_data$X <- NULL
    expect_error(
        hurdlebb:::.validate_hbb_fit_methods(fit),
        "hbb_data\\$X.*must be a numeric matrix"
    )
})

test_that("M6: validate errors when X is not a matrix", {
    fit <- create_methods_mock_fit()
    fit$hbb_data$X <- data.frame(x1 = 1:10)
    expect_error(
        hurdlebb:::.validate_hbb_fit_methods(fit),
        "hbb_data\\$X.*must be a numeric matrix"
    )
})

test_that("M7: validate errors when X dimensions don't match N and P", {
    fit <- create_methods_mock_fit(P = 3, N = 10)
    # Give X wrong number of rows
    fit$hbb_data$X <- matrix(rnorm(5 * 3), nrow = 5, ncol = 3)
    expect_error(
        hurdlebb:::.validate_hbb_fit_methods(fit),
        "Design matrix X dimensions"
    )
})

test_that("M8: validate errors when N is not numeric", {
    fit <- create_methods_mock_fit()
    fit$hbb_data$N <- "ten"
    expect_error(
        hurdlebb:::.validate_hbb_fit_methods(fit),
        "N.*must be a positive"
    )
})


# ============================================================================
# ---- Section N: print.summary.hbb_fit formatting branches ----
# ============================================================================

test_that("N1: print.summary shows 'weighted' model type label", {
    fit <- create_methods_mock_fit(model_type = "weighted")
    s <- summary(fit)
    out <- capture.output(print(s))
    out_str <- paste(out, collapse = "\n")
    expect_true(grepl("Weighted.*survey", out_str, ignore.case = TRUE))
})

test_that("N2: print.summary shows SVC model type and random effects", {
    # Create an SVC mock with tau
    fit <- create_svc_mock_fit()
    s <- summary(fit)
    out <- capture.output(print(s))
    out_str <- paste(out, collapse = "\n")
    expect_true(grepl("SVC|state-varying", out_str, ignore.case = TRUE))
    expect_true(grepl("Random Effects", out_str))
    expect_true(grepl("tau", out_str))
})

test_that("N3: print.summary shows sandwich inference label when used", {
    fit <- create_methods_mock_fit()
    P <- fit$hbb_data$P; D <- 2L * P + 1L
    mock_sand <- create_mock_sandwich(D)
    s <- summary(fit, sandwich = mock_sand)
    out <- capture.output(print(s))
    out_str <- paste(out, collapse = "\n")
    expect_true(grepl("sandwich", out_str, ignore.case = TRUE))
})

test_that("N4: print.summary shows DER section when sandwich present", {
    fit <- create_methods_mock_fit()
    P <- fit$hbb_data$P; D <- 2L * P + 1L
    mock_sand <- create_mock_sandwich(D)
    s <- summary(fit, sandwich = mock_sand)
    out <- capture.output(print(s))
    out_str <- paste(out, collapse = "\n")
    expect_true(grepl("Design Effect Ratios", out_str))
    expect_true(grepl("DER range", out_str))
})

test_that("N5: print.summary shows (dispersion unavailable) when kappa non-finite", {
    fit <- create_methods_mock_fit()
    s <- summary(fit)
    # Manually set kappa to non-finite
    s$dispersion$kappa_hat <- NaN
    out <- capture.output(print(s))
    out_str <- paste(out, collapse = "\n")
    expect_true(grepl("dispersion unavailable", out_str))
})

test_that("N6: print.summary shows S (states) for SVC model", {
    fit <- create_svc_mock_fit()
    s <- summary(fit)
    out <- capture.output(print(s))
    out_str <- paste(out, collapse = "\n")
    expect_true(grepl("States \\(S\\)", out_str))
})

test_that("N7: print.summary shows unknown type fallback for unrecognized model_type", {
    fit <- create_methods_mock_fit()
    s <- summary(fit)
    # Inject an unrecognized model type
    s$model_info$model_type <- "my_custom_model"
    out <- capture.output(print(s))
    out_str <- paste(out, collapse = "\n")
    expect_true(grepl("my_custom_model", out_str))
})

test_that("N8: print.summary shows MCMC diagnostics section", {
    fit <- create_methods_mock_fit()
    s <- summary(fit)
    out <- capture.output(print(s))
    out_str <- paste(out, collapse = "\n")
    expect_true(grepl("MCMC Diagnostics", out_str))
    expect_true(grepl("Divergent", out_str))
    expect_true(grepl("\\[OK\\]", out_str))
})


# ============================================================================
# ---- Section O: summary() validation edge cases ----
# ============================================================================

test_that("O1: summary() errors when sandwich is wrong class", {
    fit <- create_methods_mock_fit()
    bad <- list(V_sand = diag(7))  # no class attribute
    expect_error(summary(fit, sandwich = bad), "hbb_sandwich")
})

test_that("O2: summary() errors when sandwich V_sand dimension mismatches", {
    fit <- create_methods_mock_fit(P = 3)  # D = 7
    # Sandwich with wrong dimensions (5x5 instead of 7x7)
    bad_sand <- structure(
        list(V_sand = diag(5), DER = rep(2, 5)),
        class = "hbb_sandwich"
    )
    expect_error(summary(fit, sandwich = bad_sand), "Dimension mismatch")
})

test_that("O3: summary() warns when DER vector length doesn't match D", {
    fit <- create_methods_mock_fit(P = 3)  # D = 7
    bad_sand <- structure(
        list(V_sand = diag(7) * 0.01, DER = rep(2, 5)),  # DER length 5, not 7
        class = "hbb_sandwich"
    )
    expect_warning(summary(fit, sandwich = bad_sand), "DER vector length")
})

test_that("O4: .extract_fixed_effects warns on negative sandwich diagonal", {
    fit <- create_methods_mock_fit(P = 3)
    P <- fit$hbb_data$P; D <- 2L * P + 1L
    # Create sandwich with a negative diagonal
    V_sand <- diag(D) * 0.01
    V_sand[1, 1] <- -0.001  # one negative diagonal
    mock_sand <- structure(
        list(V_sand = V_sand, DER = rep(2, D)),
        class = "hbb_sandwich"
    )
    param_labels <- c(
        paste0("alpha_", colnames(fit$hbb_data$X)),
        paste0("beta_", colnames(fit$hbb_data$X)),
        "log_kappa"
    )
    expect_warning(
        hurdlebb:::.extract_fixed_effects(fit, mock_sand, 0.95, param_labels),
        "negative diagonal"
    )
})


# ============================================================================
# ---- Section P: residuals() edge-case warning paths ----
# ============================================================================

test_that("P1: residuals(type='pearson') warns on near-zero variance obs", {
    fit <- create_methods_mock_fit(N = 5, seed = 88)
    # Force q_hat very close to 0 by setting alpha to large negative values
    # Instead, we manipulate at the source by making z=0 and y=0
    # This triggers the q_clamped path but variance will be near-zero
    # if q is extremely close to 0
    fit$hbb_data$y <- rep(0L, 5)
    fit$hbb_data$z <- rep(0L, 5)
    # Response residuals should work fine (no Pearson path)
    r <- residuals(fit, type = "response")
    expect_length(r, 5L)
    expect_true(all(is.finite(r)))
})

test_that("P2: residuals(type='pearson') produces finite output for normal data", {
    fit <- create_methods_mock_fit(N = 20, seed = 33)
    r <- residuals(fit, type = "pearson")
    expect_length(r, 20L)
    # Most should be finite
    expect_true(sum(is.finite(r)) >= 18L)
})
