# ============================================================================
# test-predict.R --- testthat tests for R/predict.R
#
# Sections:
#   MOCK:  3 constructors (base with center/scale, SVC, nontrivial center/scale)
#   A.     In-sample point predictions (7 tests)
#   B.     Credible intervals (5 tests)
#   C.     Out-of-sample newdata (6 tests)
#   D.     SVC state-specific prediction (5 tests)
#   E.     Type argument (5 tests)
#   F.     Input validation (8 tests)
#   G.     Edge cases (5 tests)
#   H.     Standardization centering/scaling (4 tests)
#
# Total: 45 tests across 8 sections.
# ============================================================================


# ============================================================================
# MOCK CONSTRUCTORS
# ============================================================================

#' Create a mock hbb_fit object for predict.R testing
#'
#' Canonical version combining x_center/x_scale fields (identity by
#' default: center = 0, scale = 1) with generalised formula
#' construction and draws(format = "array") support.
#'
#' @param P Number of covariates per margin (including intercept).
#' @param N Number of observations.
#' @param M Number of MCMC draws (total across chains).
#' @param n_chains Number of MCMC chains.
#' @param seed Random seed.
#' @param model_type Character: "base" or "weighted".
#' @return An S3 object of class "hbb_fit".
create_predict_mock_fit <- function(P = 3, N = 10, M = 100, n_chains = 4,
                                     seed = 42, model_type = "base") {
    set.seed(seed)
    P <- as.integer(P); N <- as.integer(N); D <- 2L * P + 1L

    alpha_true <- seq(-0.5, 0.5, length.out = P)
    beta_true <- seq(0.2, 0.8, length.out = P)
    log_kappa_true <- log(5)
    theta_true <- c(alpha_true, beta_true, log_kappa_true)

    draws_mat <- matrix(rep(theta_true, each = M) + rnorm(M * D, sd = 0.001),
                         nrow = M, ncol = D)
    colnames(draws_mat) <- c(paste0("alpha[", 1:P, "]"),
                              paste0("beta[", 1:P, "]"),
                              "log_kappa")

    iter_sampling <- M / n_chains

    # Create 3-D array: [iter_sampling, n_chains, D]
    draws_array <- array(NA_real_, dim = c(iter_sampling, n_chains, D))
    for (ch in seq_len(n_chains)) {
        row_start <- (ch - 1L) * iter_sampling + 1L
        row_end   <- ch * iter_sampling
        draws_array[, ch, ] <- draws_mat[row_start:row_end, ]
    }
    dimnames(draws_array) <- list(
        NULL,
        paste0("chain:", seq_len(n_chains)),
        colnames(draws_mat)
    )

    fixed_names <- if (P > 1L) paste0("x", 2:P) else character(0)

    if (P == 1L) {
        X <- matrix(1, nrow = N, ncol = 1L); colnames(X) <- "intercept"
    } else {
        X <- cbind(1, matrix(rnorm(N * (P - 1L)), nrow = N))
        colnames(X) <- c("intercept", fixed_names)
    }

    eta_ext <- as.numeric(X %*% alpha_true)
    eta_int <- as.numeric(X %*% beta_true)
    q_true <- plogis(eta_ext); mu_true <- plogis(eta_int)
    n_trial <- rep(20L, N); z <- rbinom(N, 1, q_true); y <- integer(N)
    y[z == 1L] <- pmax(1L, rbinom(sum(z), n_trial[z == 1L], mu_true[z == 1L]))

    mock_fit <- list(
        draws = function(variables = NULL, format = "matrix") {
            if (format == "array") {
                if (is.null(variables)) return(draws_array)
                idx <- match(variables, colnames(draws_mat))
                return(draws_array[, , idx, drop = FALSE])
            }
            if (is.null(variables)) return(draws_mat)
            idx <- match(variables, colnames(draws_mat))
            draws_mat[, idx, drop = FALSE]
        },
        summary = function(variables = NULL) {
            if (is.null(variables)) vars <- colnames(draws_mat) else vars <- variables
            idx <- match(vars, colnames(draws_mat))
            sub <- draws_mat[, idx, drop = FALSE]
            data.frame(variable = vars,
                       mean = colMeans(sub),
                       median = apply(sub, 2, median),
                       sd = apply(sub, 2, sd),
                       rhat = rep(1.001, length(vars)),
                       ess_bulk = rep(2000, length(vars)),
                       ess_tail = rep(1800, length(vars)),
                       q5 = apply(sub, 2, quantile, 0.05),
                       q95 = apply(sub, 2, quantile, 0.95),
                       stringsAsFactors = FALSE)
        },
        diagnostic_summary = function(quiet = FALSE) {
            list(num_divergent = rep(0L, n_chains),
                 num_max_treedepth = rep(0L, n_chains),
                 ebfmi = rep(1.1, n_chains))
        },
        num_chains = function() n_chains,
        metadata = function() list(iter_sampling = iter_sampling,
                                    stan_variables = colnames(draws_mat))
    )

    # Identity centering and scaling: center = 0, scale = 1
    # This means newdata is used as-is (no transformation)
    x_center <- setNames(rep(0, length(fixed_names)), fixed_names)
    x_scale <- setNames(rep(1, length(fixed_names)), fixed_names)

    hbb_data <- list(N = N, P = P, X = X, y = y, n_trial = n_trial, z = z,
                     x_center = x_center, x_scale = x_scale)

    structure(list(fit = mock_fit, hbb_data = hbb_data,
                   model_type = model_type,
                   formula = list(
                       fixed = fixed_names,
                       formula = if (P > 1) {
                           as.formula(paste("y ~", paste0("x", 2:P, collapse = " + ")))
                       } else {
                           y ~ 1
                       }
                   ),
                   call = quote(hbb())),
              class = "hbb_fit")
}


#' Create a mock SVC hbb_fit object for predict.R testing
#'
#' Includes delta random effects, state structure, and centering/scaling.
#'
#' @param P Number of covariates per margin.
#' @param N Number of observations.
#' @param S Number of states.
#' @param M Number of MCMC draws.
#' @param n_chains Number of MCMC chains.
#' @param seed Random seed.
#' @return An S3 object of class "hbb_fit" with model_type = "svc".
create_svc_predict_mock_fit <- function(P = 3, N = 20, S = 5, M = 100,
                                         n_chains = 4, seed = 42) {
    set.seed(seed)
    P <- as.integer(P); N <- as.integer(N); S <- as.integer(S)
    D <- 2L * P + 1L; K <- 2L * P

    alpha_true <- seq(-0.5, 0.5, length.out = P)
    beta_true <- seq(0.2, 0.8, length.out = P)
    log_kappa_true <- log(5)
    theta_true <- c(alpha_true, beta_true, log_kappa_true)

    draws_mat <- matrix(rep(theta_true, each = M) + rnorm(M * D, sd = 0.001),
                         nrow = M, ncol = D)
    colnames(draws_mat) <- c(paste0("alpha[", 1:P, "]"),
                              paste0("beta[", 1:P, "]"),
                              "log_kappa")

    delta_true <- matrix(rnorm(S * K, sd = 0.2), nrow = S, ncol = K)
    delta_vec_true <- as.numeric(t(delta_true))
    n_delta <- S * K
    delta_draws <- matrix(rep(delta_vec_true, each = M) +
                              rnorm(M * n_delta, sd = 0.001),
                           nrow = M, ncol = n_delta)
    delta_colnames <- character(n_delta)
    idx <- 1L
    for (s in seq_len(S)) {
        for (k in seq_len(K)) {
            delta_colnames[idx] <- sprintf("delta[%d,%d]", s, k)
            idx <- idx + 1L
        }
    }
    colnames(delta_draws) <- delta_colnames

    all_draws <- cbind(draws_mat, delta_draws)
    iter_sampling <- M / n_chains

    fixed_names <- if (P > 1L) paste0("x", 2:P) else character(0)

    if (P == 1L) {
        X <- matrix(1, nrow = N, ncol = 1L); colnames(X) <- "intercept"
    } else {
        X <- cbind(1, matrix(rnorm(N * (P - 1L)), nrow = N))
        colnames(X) <- c("intercept", fixed_names)
    }

    state_vec <- rep(seq_len(S), length.out = N)
    state_labels <- paste0("State_", LETTERS[seq_len(S)])
    n_trial <- rep(20L, N); z <- rbinom(N, 1, 0.7); y <- integer(N)
    y[z == 1L] <- pmax(1L, rbinom(sum(z), n_trial[z == 1L], 0.3))

    mock_fit <- list(
        draws = function(variables = NULL, format = "matrix") {
            if (is.null(variables)) return(all_draws)
            idx <- match(variables, colnames(all_draws))
            all_draws[, idx, drop = FALSE]
        },
        summary = function(variables = NULL) {
            if (is.null(variables)) vars <- colnames(all_draws) else vars <- variables
            idx <- match(vars, colnames(all_draws))
            sub <- all_draws[, idx, drop = FALSE]
            data.frame(variable = vars,
                       mean = colMeans(sub),
                       median = apply(sub, 2, median),
                       sd = apply(sub, 2, sd),
                       rhat = rep(1.001, length(vars)),
                       ess_bulk = rep(2000, length(vars)),
                       ess_tail = rep(1800, length(vars)),
                       q5 = apply(sub, 2, quantile, 0.05),
                       q95 = apply(sub, 2, quantile, 0.95),
                       stringsAsFactors = FALSE)
        },
        diagnostic_summary = function(quiet = FALSE) {
            list(num_divergent = rep(0L, n_chains),
                 num_max_treedepth = rep(0L, n_chains),
                 ebfmi = rep(1.1, n_chains))
        },
        num_chains = function() n_chains,
        metadata = function() list(iter_sampling = iter_sampling,
                                    stan_variables = colnames(all_draws))
    )

    x_center <- setNames(rep(0, length(fixed_names)), fixed_names)
    x_scale <- setNames(rep(1, length(fixed_names)), fixed_names)

    hbb_data <- list(N = N, P = P, S = S, K = K, X = X,
                     y = y, n_trial = n_trial, z = z,
                     state = state_vec,
                     group_levels = state_labels,
                     x_center = x_center, x_scale = x_scale)

    structure(list(fit = mock_fit, hbb_data = hbb_data,
                   model_type = "svc",
                   formula = list(
                       fixed = fixed_names,
                       formula = if (P == 3) y ~ x2 + x3 else y ~ 1,
                       group = "state_id",
                       svc = TRUE
                   ),
                   call = quote(hbb())),
              class = "hbb_fit")
}


#' Create a mock hbb_fit with non-trivial centering and scaling
#'
#' The training data was centered at specific means and scaled by specific SDs.
#' This allows testing that predict() correctly applies the transformation.
#'
#' @param seed Random seed.
#' @return An S3 object of class "hbb_fit" with non-trivial x_center and x_scale.
create_centered_predict_mock_fit <- function(seed = 42) {
    set.seed(seed)
    P <- 3L; N <- 20L; M <- 100L; n_chains <- 4L; D <- 2L * P + 1L

    # Define non-trivial centering/scaling from hypothetical training data
    fixed_names <- c("x2", "x3")
    train_means <- c(x2 = 2.5, x3 = -1.0)
    train_sds <- c(x2 = 1.5, x3 = 0.8)

    alpha_true <- seq(-0.5, 0.5, length.out = P)
    beta_true <- seq(0.2, 0.8, length.out = P)
    log_kappa_true <- log(5)
    theta_true <- c(alpha_true, beta_true, log_kappa_true)

    draws_mat <- matrix(rep(theta_true, each = M) + rnorm(M * D, sd = 0.001),
                         nrow = M, ncol = D)
    colnames(draws_mat) <- c(paste0("alpha[", 1:P, "]"),
                              paste0("beta[", 1:P, "]"),
                              "log_kappa")

    iter_sampling <- M / n_chains

    # Build training X with non-trivial means and SDs
    # Generate raw, then center and scale
    x2_raw <- rnorm(N, mean = 2.5, sd = 1.5)
    x3_raw <- rnorm(N, mean = -1.0, sd = 0.8)
    x2_std <- (x2_raw - train_means["x2"]) / train_sds["x2"]
    x3_std <- (x3_raw - train_means["x3"]) / train_sds["x3"]
    X <- cbind(1, x2_std, x3_std)
    colnames(X) <- c("intercept", "x2", "x3")

    eta_ext <- as.numeric(X %*% alpha_true)
    eta_int <- as.numeric(X %*% beta_true)
    q_true <- plogis(eta_ext); mu_true <- plogis(eta_int)
    n_trial <- rep(20L, N); z <- rbinom(N, 1, q_true); y <- integer(N)
    y[z == 1L] <- pmax(1L, rbinom(sum(z), n_trial[z == 1L], mu_true[z == 1L]))

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
            data.frame(variable = vars,
                       mean = colMeans(sub),
                       median = apply(sub, 2, median),
                       sd = apply(sub, 2, sd),
                       rhat = rep(1.001, length(vars)),
                       ess_bulk = rep(2000, length(vars)),
                       ess_tail = rep(1800, length(vars)),
                       q5 = apply(sub, 2, quantile, 0.05),
                       q95 = apply(sub, 2, quantile, 0.95),
                       stringsAsFactors = FALSE)
        },
        diagnostic_summary = function(quiet = FALSE) {
            list(num_divergent = rep(0L, n_chains),
                 num_max_treedepth = rep(0L, n_chains),
                 ebfmi = rep(1.1, n_chains))
        },
        num_chains = function() n_chains,
        metadata = function() list(iter_sampling = iter_sampling,
                                    stan_variables = colnames(draws_mat))
    )

    hbb_data <- list(N = N, P = P, X = X, y = y, n_trial = n_trial, z = z,
                     x_center = train_means,
                     x_scale = train_sds)

    # Store raw training data for verification in tests
    train_df <- data.frame(x2 = x2_raw, x3 = x3_raw)

    obj <- structure(list(fit = mock_fit, hbb_data = hbb_data,
                          model_type = "base",
                          formula = list(
                              fixed = fixed_names,
                              formula = y ~ x2 + x3
                          ),
                          call = quote(hbb())),
                     class = "hbb_fit")

    # Return both the fit object and the raw training data for test verification
    list(fit = obj, train_df = train_df)
}


# ============================================================================
# Section A: In-sample point predictions (7 tests)
# ============================================================================

test_that("A1: predict(fit) returns data.frame with column 'fit'", {
    fit  <- create_predict_mock_fit(P = 3)
    pred <- predict(fit)

    expect_s3_class(pred, "data.frame")
    expect_true("fit" %in% names(pred))
})


test_that("A2: predictions have length N", {
    N   <- 15L
    fit <- create_predict_mock_fit(P = 3, N = N)
    pred <- predict(fit)

    expect_equal(nrow(pred), N)
})


test_that("A3: response predictions match fitted(fit, type='response') within tol", {
    fit  <- create_predict_mock_fit(P = 3, seed = 123)
    pred <- predict(fit, type = "response")
    fv   <- fitted(fit, type = "response")

    # predict() and fitted() should be nearly identical for in-sample
    # (both use posterior means of fixed effects applied to training X)
    expect_equal(pred$fit, fv, tolerance = 1e-6)
})


test_that("A4: extensive predictions match fitted(fit, type='extensive')", {
    fit  <- create_predict_mock_fit(P = 3, seed = 123)
    pred <- predict(fit, type = "extensive")
    fv   <- fitted(fit, type = "extensive")

    expect_equal(pred$fit, fv, tolerance = 1e-6)
})


test_that("A5: intensive predictions match fitted(fit, type='intensive')", {
    fit  <- create_predict_mock_fit(P = 3, seed = 123)
    pred <- predict(fit, type = "intensive")
    fv   <- fitted(fit, type = "intensive")

    expect_equal(pred$fit, fv, tolerance = 1e-6)
})


test_that("A6: all response predictions are in [0, 1]", {
    fit  <- create_predict_mock_fit(P = 3, N = 50, seed = 77)
    pred <- predict(fit, type = "response")

    expect_true(all(pred$fit >= 0))
    expect_true(all(pred$fit <= 1))
})


test_that("A7: P=1 (intercept-only) works correctly", {
    fit  <- create_predict_mock_fit(P = 1, N = 10)
    pred <- predict(fit, type = "response")

    expect_s3_class(pred, "data.frame")
    expect_equal(nrow(pred), 10L)
    expect_true(all(pred$fit >= 0 & pred$fit <= 1))

    # For intercept-only, all predictions should be the same
    # (since all rows of X are just the intercept)
    expect_equal(length(unique(round(pred$fit, 8))), 1L)
})


# ============================================================================
# Section B: Credible intervals (5 tests)
# ============================================================================

test_that("B1: interval='credible' returns 3 columns: fit, lwr, upr", {
    fit  <- create_predict_mock_fit(P = 3, M = 200)
    pred <- predict(fit, interval = "credible", level = 0.95)

    expect_s3_class(pred, "data.frame")
    expect_true(all(c("fit", "lwr", "upr") %in% names(pred)))
    expect_equal(ncol(pred), 3L)
})


test_that("B2: lwr <= fit <= upr for all rows", {
    fit  <- create_predict_mock_fit(P = 3, M = 200, seed = 55)
    pred <- predict(fit, interval = "credible", level = 0.95)

    expect_true(all(pred$lwr <= pred$fit + 1e-10))
    expect_true(all(pred$fit <= pred$upr + 1e-10))
})


test_that("B3: level=0.50 produces narrower CIs than level=0.99", {
    fit <- create_predict_mock_fit(P = 3, M = 400, seed = 99)

    pred50 <- predict(fit, interval = "credible", level = 0.50)
    pred99 <- predict(fit, interval = "credible", level = 0.99)

    width50 <- pred50$upr - pred50$lwr
    width99 <- pred99$upr - pred99$lwr

    # 99% interval must be at least as wide as 50% for every observation
    expect_true(all(width99 >= width50 - 1e-10))

    # At least some should be strictly wider (draws have non-zero variance)
    expect_true(any(width99 > width50 + 1e-10))
})


test_that("B4: ndraws=10 works and returns valid data.frame", {
    fit  <- create_predict_mock_fit(P = 3, M = 200)
    pred <- predict(fit, interval = "credible", level = 0.95, ndraws = 10)

    expect_s3_class(pred, "data.frame")
    expect_true(all(c("fit", "lwr", "upr") %in% names(pred)))
    expect_equal(nrow(pred), 10L)  # N = 10 from default
    expect_true(all(pred$lwr <= pred$fit + 1e-10))
    expect_true(all(pred$fit <= pred$upr + 1e-10))
})


test_that("B5: ndraws=NULL (use all draws) works", {
    fit <- create_predict_mock_fit(P = 3, M = 100)

    # ndraws = NULL is the default: should use all M draws
    pred <- predict(fit, interval = "credible", level = 0.95, ndraws = NULL)

    expect_s3_class(pred, "data.frame")
    expect_true(all(c("fit", "lwr", "upr") %in% names(pred)))
    expect_true(all(pred$lwr <= pred$fit + 1e-10))
    expect_true(all(pred$fit <= pred$upr + 1e-10))
})


# ============================================================================
# Section C: Out-of-sample newdata (6 tests)
# ============================================================================

test_that("C1: predict(fit, newdata = train_df) returns data.frame", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)

    # Build newdata matching the training covariates (same values as design matrix)
    train_df <- data.frame(
        x2 = fit$hbb_data$X[, "x2"],
        x3 = fit$hbb_data$X[, "x3"]
    )

    result <- predict(fit, newdata = train_df)

    expect_s3_class(result, "data.frame")
    expect_true("fit" %in% names(result),
                info = "Prediction data.frame should have a 'fit' column")
})


test_that("C2: newdata with novel values works and returns [0,1] predictions", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)

    # Novel covariate values not in training data
    novel_df <- data.frame(
        x2 = c(-3.0, 0.0, 3.0, 5.0, -5.0),
        x3 = c(2.0, -2.0, 0.0, 1.5, -1.5)
    )

    result <- predict(fit, newdata = novel_df)

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), nrow(novel_df))
    # All predictions should be in [0, 1]
    expect_true(all(result$fit >= 0 & result$fit <= 1),
                info = "All predictions should be clamped to [0, 1]")
})


test_that("C3: Missing column in newdata raises error", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)

    # Provide only x2 but model needs x2 and x3
    incomplete_df <- data.frame(x2 = c(0.1, 0.2, 0.3))

    expect_error(
        predict(fit, newdata = incomplete_df),
        regexp = "missing|x3|required",
        info = "Missing covariate column should raise an informative error"
    )
})


test_that("C4: Extra columns in newdata are silently ignored", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)

    # Include all required columns plus extras
    extra_df <- data.frame(
        x2 = c(0.1, 0.2),
        x3 = c(-0.1, 0.3),
        x_extra = c(100, 200),
        another_col = c("a", "b")
    )

    # Should work without error or warning about extra columns
    result <- predict(fit, newdata = extra_df)

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 2L)
})


test_that("C5: newdata with 1 row works", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)

    single_row <- data.frame(x2 = 0.5, x3 = -0.3)

    result <- predict(fit, newdata = single_row)

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 1L,
                 info = "Single-row newdata should produce single-row prediction")
    expect_true(result$fit >= 0 && result$fit <= 1)
})


test_that("C6: newdata predictions have correct number of rows", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)

    for (n_rows in c(1L, 5L, 50L)) {
        nd <- data.frame(
            x2 = rnorm(n_rows),
            x3 = rnorm(n_rows)
        )
        result <- predict(fit, newdata = nd)
        expect_equal(nrow(result), n_rows,
                     info = paste("newdata with", n_rows,
                                  "rows should yield", n_rows, "predictions"))
    }
})


# ============================================================================
# Section D: SVC state-specific prediction (5 tests)
# ============================================================================

test_that("D1: SVC model with state = 'State_A' returns state-specific predictions", {
    svc_fit <- create_svc_predict_mock_fit(P = 3, N = 20, S = 5, M = 100)

    nd <- data.frame(
        x2 = c(0.1, 0.2, 0.3),
        x3 = c(-0.1, 0.0, 0.1)
    )

    # Provide state as a vector of length N_new
    result <- predict(svc_fit, newdata = nd,
                      state = rep("State_A", 3))

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 3L)
    expect_true(all(result$fit >= 0 & result$fit <= 1))
})


test_that("D2: SVC model without state gives population-average (with warning)", {
    svc_fit <- create_svc_predict_mock_fit(P = 3, N = 20, S = 5, M = 100)

    nd <- data.frame(
        x2 = c(0.1, 0.2),
        x3 = c(-0.1, 0.0)
    )

    # state = NULL on SVC model should warn and use fixed effects only
    expect_warning(
        result <- predict(svc_fit, newdata = nd),
        regexp = "state.*NULL|population|fixed effects",
        info = "SVC model without state should warn about population-average"
    )

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 2L)
})


test_that("D3: Population-average differs from state-specific predictions", {
    svc_fit <- create_svc_predict_mock_fit(P = 3, N = 20, S = 5, M = 100)

    nd <- data.frame(
        x2 = c(0.5, -0.5, 1.0),
        x3 = c(0.0, 0.5, -0.5)
    )

    # Population average (no state, expect warning)
    pa_result <- suppressWarnings(predict(svc_fit, newdata = nd))

    # State-specific (State_A)
    state_result <- predict(svc_fit, newdata = nd,
                            state = rep("State_A", 3))

    # Since delta random effects are non-zero (sd = 0.2), the state-specific
    # predictions should differ from population-average (which ignores delta)
    expect_false(isTRUE(all.equal(pa_result$fit, state_result$fit)),
                 info = paste("Population-average and state-specific predictions",
                              "should differ when random effects are non-zero"))
})


test_that("D4: Unknown state label raises error", {
    svc_fit <- create_svc_predict_mock_fit(P = 3, N = 20, S = 5, M = 100)

    nd <- data.frame(x2 = 0.1, x3 = -0.1)

    expect_error(
        predict(svc_fit, newdata = nd, state = "Nonexistent_State"),
        regexp = "unknown|Nonexistent_State|state",
        info = "Unknown state label should produce an informative error"
    )
})


test_that("D5: state argument on base model raises error", {
    base_fit <- create_predict_mock_fit(P = 3, N = 10, M = 100,
                                         model_type = "base")

    nd <- data.frame(x2 = 0.1, x3 = -0.1)

    expect_error(
        predict(base_fit, newdata = nd, state = "State_A"),
        regexp = "state.*SVC|model type|base",
        info = "Using state argument on base model should raise an error"
    )
})


# ============================================================================
# Section E: Type argument (5 tests)
# ============================================================================

test_that("E1: type = 'response' returns q * mu (known-answer check)", {
    fit <- create_predict_mock_fit(P = 3, N = 15, M = 100)
    pred_resp <- predict(fit, type = "response")
    pred_ext  <- predict(fit, type = "extensive")
    pred_int  <- predict(fit, type = "intensive")

    # response should be q * mu (element-wise product of ext and int)
    expected <- pred_ext$fit * pred_int$fit
    expect_equal(pred_resp$fit, expected, tolerance = 1e-10)
})

test_that("E2: type = 'extensive' returns q only (all in (0,1))", {
    fit <- create_predict_mock_fit(P = 3, N = 20, M = 100)
    pred <- predict(fit, type = "extensive")

    expect_true(is.data.frame(pred))
    expect_true("fit" %in% names(pred))
    # q = logistic(X %*% alpha) should be strictly in (0, 1)
    expect_true(all(pred$fit > 0))
    expect_true(all(pred$fit < 1))
})

test_that("E3: type = 'intensive' returns mu only (all in (0,1))", {
    fit <- create_predict_mock_fit(P = 3, N = 20, M = 100)
    pred <- predict(fit, type = "intensive")

    expect_true(is.data.frame(pred))
    expect_true("fit" %in% names(pred))
    # mu = logistic(X %*% beta) should be strictly in (0, 1)
    expect_true(all(pred$fit > 0))
    expect_true(all(pred$fit < 1))
})

test_that("E4: invalid type raises error", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)
    expect_error(
        predict(fit, type = "invalid_type"),
        regexp = "type"
    )
})

test_that("E5: all three types produce correct lengths equal to N", {
    N <- 25
    fit <- create_predict_mock_fit(P = 3, N = N, M = 100)

    pred_resp <- predict(fit, type = "response")
    pred_ext  <- predict(fit, type = "extensive")
    pred_int  <- predict(fit, type = "intensive")

    expect_equal(nrow(pred_resp), N)
    expect_equal(nrow(pred_ext),  N)
    expect_equal(nrow(pred_int),  N)
})


# ============================================================================
# Section F: Input validation (8 tests)
# ============================================================================

test_that("F1: non-hbb_fit input raises error", {
    expect_error(
        predict.hbb_fit(list(a = 1)),
        regexp = "hbb_fit"
    )
    expect_error(
        predict.hbb_fit(42),
        regexp = "hbb_fit"
    )
})

test_that("F2: invalid interval raises error", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)
    expect_error(
        predict(fit, interval = "confidence"),
        regexp = "interval"
    )
    expect_error(
        predict(fit, interval = "prediction"),
        regexp = "interval"
    )
})

test_that("F3: level = 0 or level = 1 raises error (via .validate_level)", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)

    # level = 0 should fail
    expect_error(
        predict(fit, interval = "credible", level = 0),
        regexp = "level"
    )

    # level = 1 should fail
    expect_error(
        predict(fit, interval = "credible", level = 1),
        regexp = "level"
    )
})

test_that("F4: negative ndraws raises error", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)
    expect_error(
        predict(fit, interval = "credible", ndraws = -5),
        regexp = "ndraws"
    )
})

test_that("F5: ndraws = 0 raises error", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)
    expect_error(
        predict(fit, interval = "credible", ndraws = 0),
        regexp = "ndraws"
    )
})

test_that("F6: newdata must be data.frame (give a matrix, expect error)", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)
    # Provide a matrix instead of data.frame
    mat_newdata <- matrix(rnorm(15), nrow = 5, ncol = 3)
    colnames(mat_newdata) <- c("intercept", "x2", "x3")
    expect_error(
        predict(fit, newdata = mat_newdata),
        regexp = "newdata.*data\\.frame|data\\.frame"
    )
})

test_that("F7: state argument on base model raises error", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100, model_type = "base")
    expect_error(
        predict(fit, state = c("AL", "CA", "NY", "TX", "FL",
                               "OH", "PA", "GA", "NC", "MI")),
        regexp = "state.*SVC|SVC.*model"
    )
})

test_that("F8: newdata with 0 rows raises error", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)
    empty_df <- data.frame(x2 = numeric(0), x3 = numeric(0))
    expect_error(
        predict(fit, newdata = empty_df),
        regexp = "zero rows|0 rows"
    )
})


# ============================================================================
# Section G: Edge cases (5 tests)
# ============================================================================

test_that("G1: N=1 observation: predict returns data.frame with 1 row", {
    fit <- create_predict_mock_fit(P = 3, N = 1, M = 100)
    pred <- predict(fit)

    expect_true(is.data.frame(pred))
    expect_equal(nrow(pred), 1L)
    expect_true("fit" %in% names(pred))
    # The single prediction should be a valid probability
    expect_true(pred$fit >= 0 && pred$fit <= 1)
})

test_that("G2: P=1 (intercept only): predict works", {
    fit <- create_predict_mock_fit(P = 1, N = 15, M = 100)
    pred <- predict(fit)

    expect_true(is.data.frame(pred))
    expect_equal(nrow(pred), 15L)
    # All predictions should be the same since X is just an intercept
    expect_true(all(pred$fit > 0))
    expect_true(all(pred$fit < 1))
    # With intercept only, all predictions should be nearly identical
    expect_true(sd(pred$fit) < 1e-6)
})

test_that("G3: very large kappa (log_kappa = 30) does not crash", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)
    # Override the draws to have very large log_kappa
    P <- fit$hbb_data$P
    D <- 2L * P + 1L
    original_draws <- fit$fit$draws(format = "matrix")
    # Set log_kappa column to 30 (kappa ~ 1e13)
    modified_draws <- original_draws
    modified_draws[, D] <- 30

    # Rebuild the mock fit with modified draws
    fit$fit$draws <- function(variables = NULL, format = "matrix") {
        if (format == "array") {
            # Not needed for predict, but handle gracefully
            return(NULL)
        }
        if (is.null(variables)) return(modified_draws)
        idx <- match(variables, colnames(modified_draws))
        modified_draws[, idx, drop = FALSE]
    }

    # predict should still work -- kappa only affects the variance model,
    # not the point predictions (which are q * mu)
    pred <- predict(fit, type = "response")
    expect_true(is.data.frame(pred))
    expect_equal(nrow(pred), 10L)
    # Predictions should still be valid probabilities
    expect_true(all(pred$fit >= 0 & pred$fit <= 1, na.rm = TRUE))
})

test_that("G4: M=4 draws (minimum) with interval='credible' works", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 4, n_chains = 4)
    pred <- predict(fit, interval = "credible", level = 0.90)

    expect_true(is.data.frame(pred))
    expect_equal(nrow(pred), 10L)
    expect_true("fit" %in% names(pred))
    expect_true("lwr" %in% names(pred))
    expect_true("upr" %in% names(pred))
    # lwr <= fit <= upr
    expect_true(all(pred$lwr <= pred$fit + 1e-10))
    expect_true(all(pred$fit <= pred$upr + 1e-10))
})

test_that("G5: predict(fit, interval='none') is deterministic (same on repeated calls)", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)
    pred1 <- predict(fit, interval = "none")
    pred2 <- predict(fit, interval = "none")

    expect_equal(pred1$fit, pred2$fit, tolerance = 0)
})


# ============================================================================
# Section H: Standardization centering/scaling (4 tests)
# ============================================================================

test_that("H1: predict(fit, newdata = training_data) matches predict(fit) when center/scale stored", {
    mock_info <- create_centered_predict_mock_fit(seed = 42)
    fit <- mock_info$fit
    train_df <- mock_info$train_df

    # In-sample predictions (uses stored X directly)
    insample <- predict(fit)

    # Out-of-sample using the raw training data (should apply center/scale
    # to reconstruct the same standardised design matrix)
    oos <- predict(fit, newdata = train_df)

    expect_equal(nrow(insample), nrow(oos),
                 info = "In-sample and newdata predictions should have same number of rows")

    # The predictions should be very close (same data, same model)
    expect_equal(insample$fit, oos$fit,
                 tolerance = 1e-4,
                 info = paste("predict(fit) and predict(fit, newdata=train_df)",
                              "should produce nearly identical predictions"))
})


test_that("H2: Without x_center/x_scale, newdata predictions still work", {
    # Create a mock where centering/scaling info is absent
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)
    # Remove centering/scaling
    fit$hbb_data$x_center <- NULL
    fit$hbb_data$x_scale <- NULL

    nd <- data.frame(x2 = c(0.1, 0.2), x3 = c(-0.1, 0.3))

    # Should still produce valid predictions (no transformation applied)
    result <- predict(fit, newdata = nd)

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 2L)
    expect_true(all(result$fit >= 0 & result$fit <= 1))
})


test_that("H3: Centering changes prediction values vs uncentered", {
    # Create two mock fits: one with non-trivial centering, one without
    set.seed(123)
    P <- 3L; N <- 10L; M <- 100L; n_chains <- 4L; D <- 2L * P + 1L
    fixed_names <- c("x2", "x3")

    alpha_true <- seq(-0.5, 0.5, length.out = P)
    beta_true <- seq(0.2, 0.8, length.out = P)
    log_kappa_true <- log(5)
    theta_true <- c(alpha_true, beta_true, log_kappa_true)

    draws_mat <- matrix(rep(theta_true, each = M) + rnorm(M * D, sd = 0.001),
                         nrow = M, ncol = D)
    colnames(draws_mat) <- c(paste0("alpha[", 1:P, "]"),
                              paste0("beta[", 1:P, "]"),
                              "log_kappa")

    # Shared design matrix (already standardised for the "uncentered" fit)
    X <- cbind(1, matrix(rnorm(N * 2), nrow = N, ncol = 2))
    colnames(X) <- c("intercept", "x2", "x3")

    eta_ext <- as.numeric(X %*% alpha_true)
    eta_int <- as.numeric(X %*% beta_true)
    q_true <- plogis(eta_ext); mu_true <- plogis(eta_int)
    n_trial <- rep(20L, N); z <- rbinom(N, 1, q_true); y <- integer(N)
    y[z == 1L] <- pmax(1L, rbinom(sum(z), n_trial[z == 1L], mu_true[z == 1L]))

    make_mock_fit_obj <- function(x_center, x_scale) {
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
                           median = apply(sub, 2, median),
                           sd = apply(sub, 2, sd),
                           rhat = rep(1.001, length(vars)),
                           ess_bulk = rep(2000, length(vars)),
                           ess_tail = rep(1800, length(vars)),
                           q5 = apply(sub, 2, quantile, 0.05),
                           q95 = apply(sub, 2, quantile, 0.95),
                           stringsAsFactors = FALSE)
            },
            diagnostic_summary = function(quiet = FALSE) {
                list(num_divergent = rep(0L, n_chains),
                     num_max_treedepth = rep(0L, n_chains),
                     ebfmi = rep(1.1, n_chains))
            },
            num_chains = function() n_chains,
            metadata = function() list(iter_sampling = M / n_chains,
                                        stan_variables = colnames(draws_mat))
        )
        hbb_data <- list(N = N, P = P, X = X, y = y,
                         n_trial = n_trial, z = z,
                         x_center = x_center, x_scale = x_scale)
        structure(list(fit = mock_fit, hbb_data = hbb_data,
                       model_type = "base",
                       formula = list(fixed = fixed_names,
                                      formula = y ~ x2 + x3),
                       call = quote(hbb())),
                  class = "hbb_fit")
    }

    # Fit A: no centering (center = 0, scale = 1)
    fit_nocenter <- make_mock_fit_obj(
        x_center = setNames(c(0, 0), fixed_names),
        x_scale = setNames(c(1, 1), fixed_names)
    )

    # Fit B: non-trivial centering
    fit_centered <- make_mock_fit_obj(
        x_center = setNames(c(5.0, -3.0), fixed_names),
        x_scale = setNames(c(2.0, 0.5), fixed_names)
    )

    # Same raw newdata
    nd <- data.frame(x2 = c(1.0, 3.0), x3 = c(0.5, -0.5))

    pred_nocenter <- predict(fit_nocenter, newdata = nd)
    pred_centered <- predict(fit_centered, newdata = nd)

    # Predictions should differ because the design matrices differ
    # (centering/scaling transforms the covariates differently)
    expect_false(isTRUE(all.equal(pred_nocenter$fit, pred_centered$fit)),
                 info = paste("Non-trivial centering should change predictions",
                              "compared to no centering"))
})


test_that("H4: newdata with column names matching formula works correctly", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)
    expected_names <- fit$formula$fixed  # c("x2", "x3")

    # Provide newdata with exactly matching column names
    nd <- data.frame(x2 = c(0.1, 0.5), x3 = c(-0.2, 0.8))

    result <- predict(fit, newdata = nd)

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 2L)
    expect_true(all(result$fit >= 0 & result$fit <= 1))

    # Verify that formula$fixed matches what we provided
    expect_true(all(expected_names %in% names(nd)),
                info = "All formula$fixed names should be present in newdata")

    # Also verify column order does not matter: swap x2 and x3
    nd_swapped <- data.frame(x3 = c(-0.2, 0.8), x2 = c(0.1, 0.5))
    result_swapped <- predict(fit, newdata = nd_swapped)

    expect_equal(result$fit, result_swapped$fit,
                 tolerance = 1e-10,
                 info = "Column order in newdata should not affect predictions")
})


# ============================================================================
# ---- Section I: CmdStan extraction error paths in predict ----
# ============================================================================

#' Create a mock hbb_fit whose $fit$draws() throws an error
#' (predict.R variant with formula$fixed, x_center, x_scale)
create_broken_draws_predict_mock <- function(P = 3, N = 10, model_type = "base") {
    P <- as.integer(P); N <- as.integer(N)

    fixed_names <- if (P > 1L) paste0("x", 2:P) else character(0)

    if (P == 1L) {
        X <- matrix(1, nrow = N, ncol = 1L)
        colnames(X) <- "intercept"
    } else {
        set.seed(42)
        X <- cbind(1, matrix(rnorm(N * (P - 1L)), nrow = N))
        colnames(X) <- c("intercept", fixed_names)
    }

    x_center <- setNames(rep(0, length(fixed_names)), fixed_names)
    x_scale <- setNames(rep(1, length(fixed_names)), fixed_names)

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

    n_trial <- rep(20L, N); z <- rep(1L, N); y <- rep(5L, N)
    hbb_data <- list(N = N, P = P, X = X, y = y, n_trial = n_trial, z = z,
                     x_center = x_center, x_scale = x_scale)

    structure(
        list(fit = mock_fit, hbb_data = hbb_data,
             model_type = model_type,
             formula = list(
                 fixed = fixed_names,
                 formula = if (P > 1) {
                     as.formula(paste("y ~", paste0("x", 2:P, collapse = " + ")))
                 } else {
                     y ~ 1
                 }
             ),
             call = quote(hbb())),
        class = "hbb_fit"
    )
}


#' Create a mock hbb_fit that returns only 1 draw (predict.R variant)
create_single_draw_predict_mock <- function(P = 3, N = 10) {
    P <- as.integer(P); N <- as.integer(N); D <- 2L * P + 1L
    set.seed(42)

    fixed_names <- if (P > 1L) paste0("x", 2:P) else character(0)
    X <- cbind(1, matrix(rnorm(N * (P - 1L)), nrow = N))
    colnames(X) <- c("intercept", fixed_names)

    single_draw <- matrix(rnorm(D, sd = 0.3), nrow = 1, ncol = D)
    colnames(single_draw) <- c(
        paste0("alpha[", 1:P, "]"),
        paste0("beta[", 1:P, "]"),
        "log_kappa"
    )

    x_center <- setNames(rep(0, length(fixed_names)), fixed_names)
    x_scale <- setNames(rep(1, length(fixed_names)), fixed_names)

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

    n_trial <- rep(20L, N); z <- rep(1L, N); y <- rep(5L, N)
    hbb_data <- list(N = N, P = P, X = X, y = y, n_trial = n_trial, z = z,
                     x_center = x_center, x_scale = x_scale)

    structure(
        list(fit = mock_fit, hbb_data = hbb_data,
             model_type = "base",
             formula = list(
                 fixed = fixed_names,
                 formula = as.formula(paste("y ~", paste0("x", 2:P, collapse = " + ")))
             ),
             call = quote(hbb())),
        class = "hbb_fit"
    )
}


#' Create a mock with NaN draws for predict testing
create_nonfinite_draws_predict_mock <- function(P = 3, N = 10) {
    P <- as.integer(P); N <- as.integer(N); D <- 2L * P + 1L
    set.seed(42)

    fixed_names <- if (P > 1L) paste0("x", 2:P) else character(0)
    X <- cbind(1, matrix(rnorm(N * (P - 1L)), nrow = N))
    colnames(X) <- c("intercept", fixed_names)

    M <- 100L
    draws_mat <- matrix(rnorm(M * D, sd = 0.3), nrow = M, ncol = D)
    colnames(draws_mat) <- c(
        paste0("alpha[", 1:P, "]"),
        paste0("beta[", 1:P, "]"),
        "log_kappa"
    )
    # Inject NaN into first column
    draws_mat[, 1L] <- NaN

    x_center <- setNames(rep(0, length(fixed_names)), fixed_names)
    x_scale <- setNames(rep(1, length(fixed_names)), fixed_names)

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

    n_trial <- rep(20L, N); z <- rep(1L, N); y <- rep(5L, N)
    hbb_data <- list(N = N, P = P, X = X, y = y, n_trial = n_trial, z = z,
                     x_center = x_center, x_scale = x_scale)

    structure(
        list(fit = mock_fit, hbb_data = hbb_data,
             model_type = "base",
             formula = list(
                 fixed = fixed_names,
                 formula = as.formula(paste("y ~", paste0("x", 2:P, collapse = " + ")))
             ),
             call = quote(hbb())),
        class = "hbb_fit"
    )
}


test_that("I1: predict() errors when draws() throws", {
    fit <- create_broken_draws_predict_mock()
    expect_error(predict(fit), "Failed to extract posterior draws")
})

test_that("I2: predict() with interval='credible' errors when draws() throws", {
    fit <- create_broken_draws_predict_mock()
    expect_error(
        predict(fit, interval = "credible"),
        "Failed to extract posterior draws"
    )
})

test_that("I3: .predict_compute_point errors on broken draws", {
    fit <- create_broken_draws_predict_mock()
    X_new <- fit$hbb_data$X
    expect_error(
        hurdlebb:::.predict_compute_point(fit, X_new, fit$hbb_data$P, NULL, "response"),
        "Failed to extract posterior draws"
    )
})

test_that("I4: .predict_compute_point warns on non-finite posterior means", {
    fit <- create_nonfinite_draws_predict_mock()
    X_new <- fit$hbb_data$X
    expect_warning(
        hurdlebb:::.predict_compute_point(fit, X_new, fit$hbb_data$P, NULL, "response"),
        "non-finite posterior mean"
    )
})

test_that("I5: .predict_compute_interval errors on broken draws", {
    fit <- create_broken_draws_predict_mock()
    X_new <- fit$hbb_data$X
    expect_error(
        hurdlebb:::.predict_compute_interval(
            fit, X_new, fit$hbb_data$P, NULL, "response", 0.95, NULL
        ),
        "Failed to extract posterior draws"
    )
})

test_that("I6: .predict_compute_interval errors when only 1 draw", {
    fit <- create_single_draw_predict_mock()
    X_new <- fit$hbb_data$X
    expect_error(
        hurdlebb:::.predict_compute_interval(
            fit, X_new, fit$hbb_data$P, NULL, "response", 0.95, NULL
        ),
        "Only 1 MCMC draw"
    )
})


# ============================================================================
# ---- Section J: .predict_build_X newdata validation paths ----
# ============================================================================

test_that("J1: .predict_build_X errors when newdata missing required columns", {
    fit <- create_predict_mock_fit(P = 3)
    bad_nd <- data.frame(x2 = c(1.0, 2.0))  # missing x3
    expect_error(
        hurdlebb:::.predict_build_X(fit, bad_nd),
        "missing required covariate"
    )
})

test_that("J2: .predict_build_X errors when newdata column is non-numeric", {
    fit <- create_predict_mock_fit(P = 3)
    bad_nd <- data.frame(x2 = c("a", "b"), x3 = c(0.5, -0.5),
                          stringsAsFactors = FALSE)
    expect_error(
        hurdlebb:::.predict_build_X(fit, bad_nd),
        "must be numeric"
    )
})

test_that("J3: .predict_build_X warns on non-finite values in newdata", {
    fit <- create_predict_mock_fit(P = 3)
    bad_nd <- data.frame(x2 = c(1.0, NaN), x3 = c(0.5, Inf))
    expect_warning(
        hurdlebb:::.predict_build_X(fit, bad_nd),
        "non-finite value"
    )
})

test_that("J4: .predict_build_X handles intercept-only model (P=1)", {
    fit <- create_predict_mock_fit(P = 1, N = 5)
    nd <- data.frame(dummy = c(1, 2, 3))
    X_new <- hurdlebb:::.predict_build_X(fit, nd)
    expect_equal(ncol(X_new), 1L)
    expect_equal(nrow(X_new), 3L)
    expect_true(all(X_new[, 1] == 1))
})

test_that("J5: .predict_build_X warns when centering info is missing for a covariate", {
    fit <- create_predict_mock_fit(P = 3)
    # Manipulate x_center so it is missing one covariate
    fit$hbb_data$x_center <- c(x2 = 0)  # missing x3
    nd <- data.frame(x2 = c(1.0, 2.0), x3 = c(0.5, -0.5))
    expect_warning(
        hurdlebb:::.predict_build_X(fit, nd),
        "Centering info missing"
    )
})

test_that("J6: .predict_build_X warns when scaling info is missing for a covariate", {
    fit <- create_predict_mock_fit(P = 3)
    # Manipulate x_scale so it is missing one covariate
    fit$hbb_data$x_scale <- c(x2 = 1)  # missing x3
    nd <- data.frame(x2 = c(1.0, 2.0), x3 = c(0.5, -0.5))
    expect_warning(
        hurdlebb:::.predict_build_X(fit, nd),
        "Scaling info missing"
    )
})

test_that("J7: .predict_build_X warns when training SD is zero for a covariate", {
    fit <- create_predict_mock_fit(P = 3)
    # Set scale for x2 to zero
    fit$hbb_data$x_scale <- c(x2 = 0, x3 = 1)
    nd <- data.frame(x2 = c(1.0, 2.0), x3 = c(0.5, -0.5))
    expect_warning(
        hurdlebb:::.predict_build_X(fit, nd),
        "Training-set SD is zero"
    )
})

test_that("J8: .predict_build_X returns correct dimensions with proper newdata", {
    fit <- create_predict_mock_fit(P = 3)
    nd <- data.frame(x2 = c(1.0, 2.0, 3.0), x3 = c(0.5, -0.5, 0.1))
    X_new <- hurdlebb:::.predict_build_X(fit, nd)
    expect_equal(ncol(X_new), fit$hbb_data$P)
    expect_equal(nrow(X_new), 3L)
    expect_equal(colnames(X_new)[1], "(Intercept)")
})

test_that("J9: .predict_build_X returns training X when newdata is NULL", {
    fit <- create_predict_mock_fit(P = 3, N = 10)
    X_new <- hurdlebb:::.predict_build_X(fit, NULL)
    expect_identical(X_new, fit$hbb_data$X)
})

test_that("J10: .predict_build_X errors when final ncol doesn't match P", {
    # Create a fit where P=4 but only 2 fixed covariates in formula
    # so that intercept + 2 = 3 != P=4
    fit <- create_predict_mock_fit(P = 3)
    fit$hbb_data$P <- 4L  # artificially set P too high
    nd <- data.frame(x2 = c(1.0, 2.0), x3 = c(0.5, -0.5))
    expect_error(
        hurdlebb:::.predict_build_X(fit, nd),
        "column mismatch"
    )
})


# ============================================================================
# ---- Section K: .predict_resolve_delta paths ----
# ============================================================================

test_that("K1: .predict_resolve_delta returns NULL for non-SVC model", {
    fit <- create_predict_mock_fit(P = 3, model_type = "base")
    result <- hurdlebb:::.predict_resolve_delta(fit, NULL, 10L)
    expect_null(result)
})

test_that("K2: .predict_resolve_delta warns when state given for non-SVC model", {
    fit <- create_predict_mock_fit(P = 3, model_type = "base")
    expect_warning(
        hurdlebb:::.predict_resolve_delta(fit, "AL", 10L),
        "state.*ignored.*non-SVC"
    )
})

test_that("K3: .predict_resolve_delta warns when SVC but state is NULL", {
    fit <- create_svc_predict_mock_fit(P = 3, N = 20, S = 5)
    expect_warning(
        hurdlebb:::.predict_resolve_delta(fit, NULL, 20L),
        "SVC model.*state.*NULL"
    )
})

test_that("K4: .predict_resolve_delta resolves character state vector", {
    fit <- create_svc_predict_mock_fit(P = 3, N = 20, S = 5)
    state_labels <- fit$hbb_data$group_levels
    # Provide a vector of state labels matching N_new
    state_vec <- rep(state_labels, length.out = 20L)
    result <- hurdlebb:::.predict_resolve_delta(fit, state_vec, 20L)
    expect_true(is.list(result))
    expect_true("delta_hat" %in% names(result))
    expect_true("state_idx" %in% names(result))
    expect_length(result$state_idx, 20L)
})

test_that("K5: .predict_resolve_delta uses column name from newdata", {
    fit <- create_svc_predict_mock_fit(P = 3, N = 20, S = 5)
    state_labels <- fit$hbb_data$group_levels
    nd <- data.frame(
        x2 = rnorm(5), x3 = rnorm(5),
        my_state = state_labels
    )
    result <- hurdlebb:::.predict_resolve_delta(fit, "my_state", 5L, newdata = nd)
    expect_true(is.list(result))
    expect_length(result$state_idx, 5L)
})

test_that("K6: .predict_resolve_delta errors on unknown state labels", {
    fit <- create_svc_predict_mock_fit(P = 3, N = 20, S = 5)
    bad_states <- rep("Unknown_State_ZZ", 20L)
    expect_error(
        hurdlebb:::.predict_resolve_delta(fit, bad_states, 20L),
        "unknown state label"
    )
})

test_that("K7: .predict_resolve_delta errors when state column not in newdata", {
    fit <- create_svc_predict_mock_fit(P = 3, N = 20, S = 5)
    nd <- data.frame(x2 = rnorm(5), x3 = rnorm(5))
    expect_error(
        hurdlebb:::.predict_resolve_delta(fit, "nonexistent_col", 5L, newdata = nd),
        "not found in"
    )
})

test_that("K8: .predict_resolve_delta errors on wrong-length state argument", {
    fit <- create_svc_predict_mock_fit(P = 3, N = 20, S = 5)
    # Provide a state vector of length 3 when N_new is 20
    expect_error(
        hurdlebb:::.predict_resolve_delta(fit, c("A", "B", "C"), 20L),
        "length 1.*or.*character vector"
    )
})

test_that("K9: .predict_resolve_delta single label repeated for in-sample", {
    fit <- create_svc_predict_mock_fit(P = 3, N = 20, S = 5)
    state_labels <- fit$hbb_data$group_levels
    # Single state label, no newdata -> repeat for N_new
    result <- hurdlebb:::.predict_resolve_delta(
        fit, state_labels[1], 20L, newdata = NULL
    )
    expect_true(is.list(result))
    expect_length(result$state_idx, 20L)
    expect_true(all(result$state_idx == 1L))
})

test_that("K10: .predict_resolve_delta errors when group_levels is NULL", {
    fit <- create_svc_predict_mock_fit(P = 3, N = 20, S = 5)
    fit$hbb_data$group_levels <- NULL
    state_labels <- paste0("State_", LETTERS[1:5])
    expect_error(
        hurdlebb:::.predict_resolve_delta(
            fit, rep(state_labels[1], 20L), 20L
        ),
        "group_levels.*NULL"
    )
})


# ============================================================================
# ---- Section L: predict() newdata validation (top-level) ----
# ============================================================================

test_that("L1: predict() warns on NA rows in newdata", {
    fit <- create_predict_mock_fit(P = 3, N = 10)
    nd <- data.frame(x2 = c(1.0, NA), x3 = c(0.5, -0.5))
    expect_warning(
        predict(fit, newdata = nd),
        "contain NA values"
    )
})

test_that("L2: predict() errors when state used with non-SVC model", {
    fit <- create_predict_mock_fit(P = 3, model_type = "base")
    expect_error(
        predict(fit, state = "AL"),
        "state.*can only be used with SVC"
    )
})

test_that("L3: predict() errors on zero-row newdata", {
    fit <- create_predict_mock_fit(P = 3)
    nd <- data.frame(x2 = numeric(0), x3 = numeric(0))
    expect_error(
        predict(fit, newdata = nd),
        "zero rows"
    )
})

test_that("L4: predict() errors when newdata is not a data.frame", {
    fit <- create_predict_mock_fit(P = 3)
    expect_error(
        predict(fit, newdata = matrix(1:4, nrow = 2)),
        "must be a data.frame"
    )
})

test_that("L5: predict() with ndraws > M warns and uses all draws", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)
    # ndraws = 200 exceeds M = 100
    expect_warning(
        predict(fit, interval = "credible", ndraws = 200L),
        "exceeds available draws"
    )
})

test_that("L6: predict() with credible interval returns lwr and upr columns", {
    fit <- create_predict_mock_fit(P = 3, N = 10, M = 100)
    result <- predict(fit, interval = "credible", level = 0.90)
    expect_true("lwr" %in% names(result))
    expect_true("upr" %in% names(result))
    expect_true(all(result$lwr <= result$fit + 1e-10))
    expect_true(all(result$fit <= result$upr + 1e-10))
})

test_that("L7: predict() with type='extensive' returns probabilities in [0,1]", {
    fit <- create_predict_mock_fit(P = 3, N = 10)
    result <- predict(fit, type = "extensive")
    expect_true(all(result$fit >= 0 & result$fit <= 1))
})

test_that("L8: predict() with type='intensive' returns probabilities in [0,1]", {
    fit <- create_predict_mock_fit(P = 3, N = 10)
    result <- predict(fit, type = "intensive")
    expect_true(all(result$fit >= 0 & result$fit <= 1))
})
