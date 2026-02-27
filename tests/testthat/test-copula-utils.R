# ============================================================================
# test-copula-utils.R — Tests for internal copula utility functions
#
# Covers:
#   A. .convert_to_z       — rank-based z-score transformation
#   B. .make_inv_ecdf      — inverse ECDF builder (linear / constant)
#   C. .copula_generate    — Gaussian copula data generation
#   D. .extract_copula_params — parameter extraction from real data
#   E. Integration tests   — full pipeline roundtrip
#
# All functions are internal (not exported), so use hurdlebb::: prefix.
# ============================================================================


# ============================================================================
# A. .convert_to_z — rank-based z-score transformation
# ============================================================================

test_that(".convert_to_z: output length matches input length", {
    set.seed(100)
    y <- rnorm(50)
    z <- hurdlebb:::.convert_to_z(y)
    expect_length(z, length(y))
})

test_that(".convert_to_z: output is approximately standard normal for large n", {
    set.seed(101)
    y <- rexp(5000, rate = 2)  # non-normal input
    z <- hurdlebb:::.convert_to_z(y)
    # Mean should be approximately 0, sd approximately 1
    expect_equal(mean(z), 0, tolerance = 0.05)
    expect_equal(sd(z), 1, tolerance = 0.05)
})

test_that(".convert_to_z: monotone — sorted input gives (approximately) sorted output", {
    set.seed(102)
    y <- sort(runif(200))
    z <- hurdlebb:::.convert_to_z(y)
    # Allow small deviations from monotonicity due to tie-breaking,
    # but the overall trend should be increasing (rank correlation ~ 1)
    expect_gt(cor(y, z, method = "spearman"), 0.99)
})

test_that(".convert_to_z: handles NA — NAs propagated, non-NA values valid", {
    set.seed(103)
    y <- c(1, 2, NA, 4, 5, NA, 7)
    z <- hurdlebb:::.convert_to_z(y)
    expect_length(z, length(y))
    expect_true(is.na(z[3]))
    expect_true(is.na(z[6]))
    # Non-NA values should be finite
    non_na <- !is.na(z)
    expect_true(all(is.finite(z[non_na])))
})

test_that(".convert_to_z: constant vector returns all zeros (with warning)", {
    y <- rep(5, 20)
    expect_warning(
        z <- hurdlebb:::.convert_to_z(y),
        regexp = NULL  # any warning is acceptable
    )
    expect_true(all(z == 0))
})

test_that(".convert_to_z: single value returns 0 (with warning)", {
    expect_warning(
        z <- hurdlebb:::.convert_to_z(42),
        regexp = NULL
    )
    expect_equal(z, 0)
})

test_that(".convert_to_z: range is bounded by qnorm(1/(N+1)) and qnorm(N/(N+1))", {
    set.seed(104)
    n <- 100
    y <- rnorm(n)
    z <- hurdlebb:::.convert_to_z(y)
    # z-scores should be in approximately [qnorm(1/(n+1)), qnorm(n/(n+1))]
    lower_bound <- qnorm(0.5 / n)   # slightly looser than 1/(n+1)
    upper_bound <- qnorm(1 - 0.5 / n)
    expect_true(all(z >= lower_bound - 0.5))
    expect_true(all(z <= upper_bound + 0.5))
})

test_that(".convert_to_z: ties produce valid z-scores (no Inf, no NaN)", {
    y <- c(1, 1, 1, 2, 2, 3, 3, 3, 3, 4)
    z <- hurdlebb:::.convert_to_z(y)
    expect_true(all(is.finite(z)))
    expect_false(any(is.nan(z)))
    expect_false(any(is.infinite(z)))
})

test_that(".convert_to_z: preserves rank order for distinct values", {
    y <- c(10, 20, 30, 40, 50)
    z <- hurdlebb:::.convert_to_z(y)
    # Strictly distinct inputs → strictly ordered z
    expect_true(all(diff(z) > 0))
})

test_that(".convert_to_z: large n (10000) completes quickly", {
    set.seed(105)
    y <- rnorm(10000)
    t0 <- proc.time()["elapsed"]
    z <- hurdlebb:::.convert_to_z(y)
    elapsed <- proc.time()["elapsed"] - t0
    expect_lt(elapsed, 5)  # should finish well within 5 seconds
    expect_length(z, 10000)
})

test_that(".convert_to_z: integer input works", {
    y <- c(1L, 3L, 5L, 7L, 9L)
    z <- hurdlebb:::.convert_to_z(y)
    expect_length(z, 5)
    expect_true(all(is.finite(z)))
})

test_that(".convert_to_z: negative values handled correctly", {
    y <- c(-10, -5, 0, 5, 10)
    z <- hurdlebb:::.convert_to_z(y)
    # Rank order preserved
    expect_true(all(diff(z) > 0))
})

test_that(".convert_to_z: output is numeric (not integer)", {
    y <- 1:20
    z <- hurdlebb:::.convert_to_z(y)
    expect_true(is.numeric(z))
    expect_true(is.double(z))
})

test_that(".convert_to_z: two-element vector works", {
    y <- c(3, 7)
    z <- hurdlebb:::.convert_to_z(y)
    expect_length(z, 2)
    expect_true(z[1] < z[2])
})

test_that(".convert_to_z: all-NA input returns all NA", {
    y <- c(NA_real_, NA_real_, NA_real_)
    z <- suppressWarnings(hurdlebb:::.convert_to_z(y))
    expect_true(all(is.na(z)))
})

test_that(".convert_to_z: reproducibility with set.seed (ties use random breaking)", {
    y <- c(1, 1, 2, 2, 3, 3)
    set.seed(200)
    z1 <- hurdlebb:::.convert_to_z(y)
    set.seed(200)
    z2 <- hurdlebb:::.convert_to_z(y)
    expect_equal(z1, z2)
})


# ============================================================================
# B. .make_inv_ecdf — inverse ECDF builder
# ============================================================================

test_that(".make_inv_ecdf: roundtrip for continuous data", {
    set.seed(201)
    x <- rnorm(500)
    inv_f <- hurdlebb:::.make_inv_ecdf(x, method = "linear")
    # ecdf(x) maps x -> [0,1]; inv_ecdf should approximately invert
    F_x <- ecdf(x)
    u <- F_x(x)
    x_recovered <- inv_f(u)
    # Correlation should be very high (not exact due to downsampling)
    expect_gt(cor(x, x_recovered), 0.95)
})

test_that(".make_inv_ecdf: range preservation — output within original range", {
    set.seed(202)
    x <- runif(1000, min = 2, max = 8)
    inv_f <- hurdlebb:::.make_inv_ecdf(x, method = "linear")
    # At boundaries
    expect_gte(inv_f(0), min(x) - 0.1)  # allow tiny numerical slack
    expect_lte(inv_f(1), max(x) + 0.1)
    # Interior points
    u_grid <- seq(0.01, 0.99, by = 0.01)
    x_grid <- inv_f(u_grid)
    expect_true(all(x_grid >= min(x) - 0.1))
    expect_true(all(x_grid <= max(x) + 0.1))
})

test_that(".make_inv_ecdf: monotone — non-decreasing for increasing input", {
    set.seed(203)
    x <- rnorm(500)
    inv_f <- hurdlebb:::.make_inv_ecdf(x, method = "linear")
    u_grid <- seq(0, 1, length.out = 200)
    x_grid <- inv_f(u_grid)
    expect_true(all(diff(x_grid) >= -1e-10))
})

test_that(".make_inv_ecdf: 'linear' method produces intermediate values", {
    x <- c(0, 1, 2, 3, 4, 5)
    inv_f <- hurdlebb:::.make_inv_ecdf(x, method = "linear")
    # At u = 0.5, linear interpolation should give something between
    # the observed quantiles, potentially non-integer
    val <- inv_f(0.5)
    expect_true(val >= 0 && val <= 5)
    # Check that some interior u values produce non-integer output
    u_test <- seq(0.1, 0.9, by = 0.1)
    x_test <- inv_f(u_test)
    has_non_integer <- any(abs(x_test - round(x_test)) > 0.001)
    expect_true(has_non_integer)
})

test_that(".make_inv_ecdf: 'constant' method produces only observed values", {
    x <- c(0, 1, 5, 10, 20)
    inv_f <- hurdlebb:::.make_inv_ecdf(x, method = "constant")
    u_grid <- seq(0, 1, length.out = 100)
    x_grid <- inv_f(u_grid)
    # All output values should be within the set of observed values
    for (val in x_grid) {
        dists <- abs(val - x)
        expect_true(min(dists) < 1e-10,
                    label = sprintf("Value %.4f not in observed set", val))
    }
})

test_that(".make_inv_ecdf: binary variable with 'constant' maps to {0, 1}", {
    x <- c(rep(0, 300), rep(1, 700))
    inv_f <- hurdlebb:::.make_inv_ecdf(x, method = "constant")
    u_grid <- seq(0, 1, length.out = 100)
    x_grid <- inv_f(u_grid)
    # Every value should be 0 or 1
    expect_true(all(x_grid == 0 | x_grid == 1))
})

test_that(".make_inv_ecdf: integer data with 'constant' preserves integer values", {
    set.seed(204)
    x <- rpois(500, lambda = 3)
    inv_f <- hurdlebb:::.make_inv_ecdf(x, method = "constant")
    u_grid <- seq(0.01, 0.99, length.out = 200)
    x_grid <- inv_f(u_grid)
    # All outputs should be integers (within tolerance)
    expect_true(all(abs(x_grid - round(x_grid)) < 1e-10))
})

test_that(".make_inv_ecdf: single unique value returns constant function", {
    x <- rep(3.14, 50)
    inv_f <- hurdlebb:::.make_inv_ecdf(x, method = "linear")
    expect_equal(inv_f(0), 3.14, tolerance = 1e-10)
    expect_equal(inv_f(0.5), 3.14, tolerance = 1e-10)
    expect_equal(inv_f(1), 3.14, tolerance = 1e-10)
})

test_that(".make_inv_ecdf: returned function works on vectors", {
    set.seed(205)
    x <- rnorm(200)
    inv_f <- hurdlebb:::.make_inv_ecdf(x, method = "linear")
    u <- c(0.1, 0.25, 0.5, 0.75, 0.9)
    result <- inv_f(u)
    expect_length(result, 5)
    expect_true(all(is.finite(result)))
})

test_that(".make_inv_ecdf: 300-point downsampling works with large input", {
    set.seed(206)
    # 10000 unique values — should be downsampled to ~300
    x <- rnorm(10000)
    inv_f <- hurdlebb:::.make_inv_ecdf(x, method = "linear")
    # Function should still produce reasonable output
    u_grid <- seq(0.01, 0.99, length.out = 50)
    x_grid <- inv_f(u_grid)
    expect_true(all(is.finite(x_grid)))
    # Should be approximately monotone
    expect_true(all(diff(x_grid) >= -0.5))  # allow small non-monotonicity from downsampling
    # Mean and range should be roughly preserved
    expect_equal(mean(x_grid), mean(x), tolerance = 0.5)
})

test_that(".make_inv_ecdf: 'linear' returns a function object", {
    x <- 1:10
    inv_f <- hurdlebb:::.make_inv_ecdf(x, method = "linear")
    expect_true(is.function(inv_f))
})

test_that(".make_inv_ecdf: 'constant' returns a function object", {
    x <- 1:10
    inv_f <- hurdlebb:::.make_inv_ecdf(x, method = "constant")
    expect_true(is.function(inv_f))
})

test_that(".make_inv_ecdf: extreme quantiles preserved", {
    set.seed(207)
    x <- rnorm(1000)
    inv_f <- hurdlebb:::.make_inv_ecdf(x, method = "linear")
    # Very low quantile ~ min(x), very high ~ max(x)
    expect_equal(inv_f(0), min(x), tolerance = 0.5)
    expect_equal(inv_f(1), max(x), tolerance = 0.5)
})

test_that(".make_inv_ecdf: monotone for 'constant' method", {
    set.seed(208)
    x <- sample(1:20, 500, replace = TRUE)
    inv_f <- hurdlebb:::.make_inv_ecdf(x, method = "constant")
    u_grid <- seq(0, 1, length.out = 200)
    x_grid <- inv_f(u_grid)
    expect_true(all(diff(x_grid) >= -1e-10))
})

test_that(".make_inv_ecdf: two distinct values produce step function with 'constant'", {
    x <- c(rep(10, 400), rep(20, 600))
    inv_f <- hurdlebb:::.make_inv_ecdf(x, method = "constant")
    # Low quantiles should map to 10, high quantiles to 20
    expect_equal(inv_f(0.1), 10)
    expect_equal(inv_f(0.9), 20)
})


# ============================================================================
# C. .copula_generate — Gaussian copula data generation
# ============================================================================

test_that(".copula_generate: output dimensions are n rows x p columns", {
    set.seed(301)
    inv_ecdfs <- list(
        x1 = hurdlebb:::.make_inv_ecdf(rnorm(500), "linear"),
        x2 = hurdlebb:::.make_inv_ecdf(rnorm(500), "linear"),
        x3 = hurdlebb:::.make_inv_ecdf(rnorm(500), "linear")
    )
    cor_mat <- diag(3)
    result <- hurdlebb:::.copula_generate(n = 200, cor_matrix = cor_mat,
                                          inv_ecdfs = inv_ecdfs, seed = 42)
    expect_equal(nrow(result), 200)
    expect_equal(ncol(result), 3)
})

test_that(".copula_generate: column names match names(inv_ecdfs)", {
    set.seed(302)
    inv_ecdfs <- list(
        var_a = hurdlebb:::.make_inv_ecdf(rnorm(200), "linear"),
        var_b = hurdlebb:::.make_inv_ecdf(rnorm(200), "linear")
    )
    cor_mat <- diag(2)
    result <- hurdlebb:::.copula_generate(n = 50, cor_matrix = cor_mat,
                                          inv_ecdfs = inv_ecdfs, seed = 1)
    expect_equal(names(result), c("var_a", "var_b"))
})

test_that(".copula_generate: returns data.frame", {
    set.seed(303)
    inv_ecdfs <- list(
        x1 = hurdlebb:::.make_inv_ecdf(rnorm(200), "linear")
    )
    cor_mat <- matrix(1, 1, 1)
    result <- hurdlebb:::.copula_generate(n = 50, cor_matrix = cor_mat,
                                          inv_ecdfs = inv_ecdfs, seed = 1)
    expect_s3_class(result, "data.frame")
})

test_that(".copula_generate: reproducibility — same seed gives same output", {
    set.seed(304)
    x_vals <- rnorm(500)
    inv_ecdfs <- list(
        x1 = hurdlebb:::.make_inv_ecdf(x_vals, "linear"),
        x2 = hurdlebb:::.make_inv_ecdf(runif(500), "linear")
    )
    cor_mat <- matrix(c(1, 0.5, 0.5, 1), 2, 2)

    result1 <- hurdlebb:::.copula_generate(n = 300, cor_matrix = cor_mat,
                                           inv_ecdfs = inv_ecdfs, seed = 999)
    result2 <- hurdlebb:::.copula_generate(n = 300, cor_matrix = cor_mat,
                                           inv_ecdfs = inv_ecdfs, seed = 999)
    expect_identical(result1, result2)
})

test_that(".copula_generate: different seeds give different output", {
    set.seed(305)
    inv_ecdfs <- list(
        x1 = hurdlebb:::.make_inv_ecdf(rnorm(300), "linear"),
        x2 = hurdlebb:::.make_inv_ecdf(rnorm(300), "linear")
    )
    cor_mat <- matrix(c(1, 0.3, 0.3, 1), 2, 2)

    result1 <- hurdlebb:::.copula_generate(n = 200, cor_matrix = cor_mat,
                                           inv_ecdfs = inv_ecdfs, seed = 1)
    result2 <- hurdlebb:::.copula_generate(n = 200, cor_matrix = cor_mat,
                                           inv_ecdfs = inv_ecdfs, seed = 2)
    expect_false(identical(result1, result2))
})

test_that(".copula_generate: correlation preservation — Spearman within tolerance", {
    set.seed(306)
    x1_raw <- rnorm(2000)
    x2_raw <- rnorm(2000)
    inv_ecdfs <- list(
        x1 = hurdlebb:::.make_inv_ecdf(x1_raw, "linear"),
        x2 = hurdlebb:::.make_inv_ecdf(x2_raw, "linear")
    )
    target_cor <- 0.6
    cor_mat <- matrix(c(1, target_cor, target_cor, 1), 2, 2)

    result <- hurdlebb:::.copula_generate(n = 5000, cor_matrix = cor_mat,
                                          inv_ecdfs = inv_ecdfs, seed = 42)
    observed_cor <- cor(result$x1, result$x2, method = "spearman")
    expect_equal(observed_cor, target_cor, tolerance = 0.1)
})

test_that(".copula_generate: marginal preservation via KS test", {
    set.seed(307)
    x_raw <- rexp(2000, rate = 1.5)
    inv_ecdfs <- list(
        x1 = hurdlebb:::.make_inv_ecdf(x_raw, "linear")
    )
    cor_mat <- matrix(1, 1, 1)

    result <- hurdlebb:::.copula_generate(n = 3000, cor_matrix = cor_mat,
                                          inv_ecdfs = inv_ecdfs, seed = 42)
    # KS test: generated vs original distribution
    ks_result <- ks.test(result$x1, x_raw)
    expect_gt(ks_result$p.value, 0.01)
})

test_that(".copula_generate: p = 1 case works for single variable", {
    set.seed(308)
    inv_ecdfs <- list(
        only_var = hurdlebb:::.make_inv_ecdf(rnorm(300), "linear")
    )
    cor_mat <- matrix(1, 1, 1)

    result <- hurdlebb:::.copula_generate(n = 100, cor_matrix = cor_mat,
                                          inv_ecdfs = inv_ecdfs, seed = 42)
    expect_equal(nrow(result), 100)
    expect_equal(ncol(result), 1)
    expect_equal(names(result), "only_var")
})

test_that(".copula_generate: validates dimension mismatch", {
    set.seed(309)
    inv_ecdfs <- list(
        x1 = hurdlebb:::.make_inv_ecdf(rnorm(200), "linear"),
        x2 = hurdlebb:::.make_inv_ecdf(rnorm(200), "linear")
    )
    # 3x3 cor_mat but only 2 inv_ecdfs
    cor_mat <- diag(3)
    expect_error(
        hurdlebb:::.copula_generate(n = 50, cor_matrix = cor_mat,
                                    inv_ecdfs = inv_ecdfs, seed = 1)
    )
})

test_that(".copula_generate: no NaN, NA, or Inf in output for valid inputs", {
    set.seed(310)
    inv_ecdfs <- list(
        x1 = hurdlebb:::.make_inv_ecdf(rnorm(500), "linear"),
        x2 = hurdlebb:::.make_inv_ecdf(rpois(500, 3), "constant"),
        x3 = hurdlebb:::.make_inv_ecdf(rbeta(500, 2, 5), "linear")
    )
    cor_mat <- matrix(c(1, 0.3, -0.2,
                        0.3, 1, 0.1,
                        -0.2, 0.1, 1), 3, 3)

    result <- hurdlebb:::.copula_generate(n = 1000, cor_matrix = cor_mat,
                                          inv_ecdfs = inv_ecdfs, seed = 42)
    for (col_name in names(result)) {
        expect_false(any(is.na(result[[col_name]])),
                     label = paste0("NA in column ", col_name))
        expect_false(any(is.nan(result[[col_name]])),
                     label = paste0("NaN in column ", col_name))
        expect_false(any(is.infinite(result[[col_name]])),
                     label = paste0("Inf in column ", col_name))
    }
})

test_that(".copula_generate: identity cor_matrix gives uncorrelated output", {
    set.seed(311)
    inv_ecdfs <- list(
        x1 = hurdlebb:::.make_inv_ecdf(rnorm(1000), "linear"),
        x2 = hurdlebb:::.make_inv_ecdf(rnorm(1000), "linear")
    )
    cor_mat <- diag(2)

    result <- hurdlebb:::.copula_generate(n = 5000, cor_matrix = cor_mat,
                                          inv_ecdfs = inv_ecdfs, seed = 42)
    observed_cor <- cor(result$x1, result$x2, method = "spearman")
    expect_equal(observed_cor, 0, tolerance = 0.05)
})

test_that(".copula_generate: negative correlation preserved", {
    set.seed(312)
    inv_ecdfs <- list(
        x1 = hurdlebb:::.make_inv_ecdf(rnorm(1000), "linear"),
        x2 = hurdlebb:::.make_inv_ecdf(rnorm(1000), "linear")
    )
    target_cor <- -0.5
    cor_mat <- matrix(c(1, target_cor, target_cor, 1), 2, 2)

    result <- hurdlebb:::.copula_generate(n = 5000, cor_matrix = cor_mat,
                                          inv_ecdfs = inv_ecdfs, seed = 42)
    observed_cor <- cor(result$x1, result$x2, method = "spearman")
    expect_equal(observed_cor, target_cor, tolerance = 0.1)
})

test_that(".copula_generate: non-PD matrix produces warning and still works", {
    set.seed(313)
    inv_ecdfs <- list(
        x1 = hurdlebb:::.make_inv_ecdf(rnorm(300), "linear"),
        x2 = hurdlebb:::.make_inv_ecdf(rnorm(300), "linear"),
        x3 = hurdlebb:::.make_inv_ecdf(rnorm(300), "linear")
    )
    # A non-PD correlation matrix (contradictory correlations)
    cor_mat <- matrix(c(1,   0.9, 0.9,
                        0.9, 1,  -0.9,
                        0.9, -0.9, 1), 3, 3)

    expect_warning(
        result <- hurdlebb:::.copula_generate(n = 100, cor_matrix = cor_mat,
                                              inv_ecdfs = inv_ecdfs, seed = 42),
        regexp = NULL
    )
    # Should still return a data.frame with correct dimensions
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 100)
    expect_equal(ncol(result), 3)
})

test_that(".copula_generate: high-dimensional case (p = 6)", {
    set.seed(314)
    p <- 6
    inv_ecdfs <- setNames(
        lapply(1:p, function(i) {
            hurdlebb:::.make_inv_ecdf(rnorm(300), "linear")
        }),
        paste0("x", 1:p)
    )
    cor_mat <- diag(p)
    # Add some modest correlations
    cor_mat[1, 2] <- cor_mat[2, 1] <- 0.3
    cor_mat[3, 4] <- cor_mat[4, 3] <- -0.2

    result <- hurdlebb:::.copula_generate(n = 500, cor_matrix = cor_mat,
                                          inv_ecdfs = inv_ecdfs, seed = 42)
    expect_equal(nrow(result), 500)
    expect_equal(ncol(result), 6)
    expect_equal(names(result), paste0("x", 1:6))
})

test_that(".copula_generate: mixed continuous and discrete marginals", {
    set.seed(315)
    inv_ecdfs <- list(
        continuous = hurdlebb:::.make_inv_ecdf(rnorm(500), "linear"),
        discrete   = hurdlebb:::.make_inv_ecdf(rpois(500, 5), "constant")
    )
    cor_mat <- matrix(c(1, 0.4, 0.4, 1), 2, 2)

    result <- hurdlebb:::.copula_generate(n = 1000, cor_matrix = cor_mat,
                                          inv_ecdfs = inv_ecdfs, seed = 42)
    # Discrete column should be integer-valued
    expect_true(all(abs(result$discrete - round(result$discrete)) < 1e-10))
    # Continuous column should have non-integer values
    expect_true(any(abs(result$continuous - round(result$continuous)) > 0.01))
})


# ============================================================================
# D. .extract_copula_params — parameter extraction from real data
# ============================================================================

test_that(".extract_copula_params: returns correct structure", {
    set.seed(401)
    df <- data.frame(
        a = rnorm(100),
        b = rpois(100, 3),
        c = rbeta(100, 2, 5)
    )
    params <- hurdlebb:::.extract_copula_params(
        data = df,
        var_names = c("a", "b", "c"),
        types = c("continuous", "discrete", "continuous")
    )
    expect_true(is.list(params))
    expect_true("cor_matrix" %in% names(params))
    expect_true("inv_ecdfs" %in% names(params))
    expect_true("n_complete" %in% names(params))
    expect_true("var_names" %in% names(params))
})

test_that(".extract_copula_params: cor_matrix is symmetric and PD", {
    set.seed(402)
    df <- data.frame(
        x1 = rnorm(200),
        x2 = rnorm(200),
        x3 = rnorm(200)
    )
    params <- hurdlebb:::.extract_copula_params(
        data = df,
        var_names = c("x1", "x2", "x3"),
        types = c("continuous", "continuous", "continuous")
    )
    R <- params$cor_matrix
    # Symmetric
    expect_equal(R, t(R), tolerance = 1e-12)
    # PD: all eigenvalues > 0
    eigs <- eigen(R, symmetric = TRUE)$values
    expect_true(all(eigs > 0))
    # Correct dimensions
    expect_equal(nrow(R), 3)
    expect_equal(ncol(R), 3)
    # Diagonal = 1
    expect_equal(unname(diag(R)), rep(1, 3), tolerance = 1e-10)
})

test_that(".extract_copula_params: inv_ecdfs is a list of functions", {
    set.seed(403)
    df <- data.frame(
        v1 = rnorm(150),
        v2 = rpois(150, 4)
    )
    params <- hurdlebb:::.extract_copula_params(
        data = df,
        var_names = c("v1", "v2"),
        types = c("continuous", "discrete")
    )
    expect_true(is.list(params$inv_ecdfs))
    expect_length(params$inv_ecdfs, 2)
    expect_true(is.function(params$inv_ecdfs[[1]]))
    expect_true(is.function(params$inv_ecdfs[[2]]))
})

test_that(".extract_copula_params: handles NAs — reduces to complete cases", {
    set.seed(404)
    n <- 100
    df <- data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n)
    )
    # Introduce 10 NAs in x1 and 5 different NAs in x2
    df$x1[1:10] <- NA
    df$x2[15:19] <- NA

    params <- hurdlebb:::.extract_copula_params(
        data = df,
        var_names = c("x1", "x2"),
        types = c("continuous", "continuous")
    )
    # n_complete should be less than 100
    expect_lt(params$n_complete, n)
    expect_equal(params$n_complete, sum(complete.cases(df[, c("x1", "x2")])))
})

test_that(".extract_copula_params: types mapping — continuous -> linear, discrete -> constant", {
    set.seed(405)
    df <- data.frame(
        cont_var = rnorm(200),
        disc_var = rpois(200, 5)
    )
    params <- hurdlebb:::.extract_copula_params(
        data = df,
        var_names = c("cont_var", "disc_var"),
        types = c("continuous", "discrete")
    )
    # Test that continuous inv_ecdf produces non-integer values
    u_test <- seq(0.1, 0.9, by = 0.1)
    cont_vals <- params$inv_ecdfs[["cont_var"]](u_test)
    has_non_int <- any(abs(cont_vals - round(cont_vals)) > 0.001)
    expect_true(has_non_int)

    # Test that discrete inv_ecdf produces only integer values
    disc_vals <- params$inv_ecdfs[["disc_var"]](u_test)
    expect_true(all(abs(disc_vals - round(disc_vals)) < 1e-10))
})

test_that(".extract_copula_params: error on invalid var_names", {
    df <- data.frame(x1 = 1:10, x2 = 11:20)
    expect_error(
        hurdlebb:::.extract_copula_params(
            data = df,
            var_names = c("x1", "nonexistent"),
            types = c("continuous", "continuous")
        )
    )
})

test_that(".extract_copula_params: error on mismatched types/var_names length", {
    df <- data.frame(x1 = 1:10, x2 = 11:20)
    expect_error(
        hurdlebb:::.extract_copula_params(
            data = df,
            var_names = c("x1", "x2"),
            types = c("continuous")  # length 1 vs 2
        )
    )
})

test_that(".extract_copula_params: var_names preserved in output", {
    set.seed(406)
    df <- data.frame(
        poverty = runif(100),
        urban = rbinom(100, 1, 0.6)
    )
    params <- hurdlebb:::.extract_copula_params(
        data = df,
        var_names = c("poverty", "urban"),
        types = c("continuous", "discrete")
    )
    expect_equal(params$var_names, c("poverty", "urban"))
    expect_equal(names(params$inv_ecdfs), c("poverty", "urban"))
})

test_that(".extract_copula_params: n_complete correct with no NAs", {
    set.seed(407)
    n <- 250
    df <- data.frame(
        x1 = rnorm(n),
        x2 = runif(n)
    )
    params <- hurdlebb:::.extract_copula_params(
        data = df,
        var_names = c("x1", "x2"),
        types = c("continuous", "continuous")
    )
    expect_equal(params$n_complete, n)
})

test_that(".extract_copula_params: single variable works", {
    set.seed(408)
    df <- data.frame(x1 = rnorm(100))
    params <- hurdlebb:::.extract_copula_params(
        data = df,
        var_names = "x1",
        types = "continuous"
    )
    expect_equal(unname(params$cor_matrix), matrix(1, 1, 1))
    expect_length(params$inv_ecdfs, 1)
})

test_that(".extract_copula_params: cor_matrix rownames and colnames match var_names", {
    set.seed(409)
    df <- data.frame(alpha = rnorm(100), beta = rnorm(100))
    params <- hurdlebb:::.extract_copula_params(
        data = df,
        var_names = c("alpha", "beta"),
        types = c("continuous", "continuous")
    )
    expect_equal(rownames(params$cor_matrix), c("alpha", "beta"))
    expect_equal(colnames(params$cor_matrix), c("alpha", "beta"))
})


# ============================================================================
# E. Integration tests — full pipeline roundtrip
# ============================================================================

test_that("integration: extract -> generate -> compare rank correlations", {
    set.seed(501)
    n_orig <- 2000
    # Create correlated data via Cholesky
    L <- matrix(c(1, 0, 0,
                  0.5, sqrt(1 - 0.25), 0,
                  -0.3, 0.2, sqrt(1 - 0.09 - 0.04)), 3, 3, byrow = TRUE)
    z <- matrix(rnorm(n_orig * 3), n_orig, 3) %*% t(L)
    orig_data <- data.frame(
        x1 = z[, 1],
        x2 = z[, 2],
        x3 = pnorm(z[, 3]) * 10  # skewed marginal
    )

    params <- hurdlebb:::.extract_copula_params(
        data = orig_data,
        var_names = c("x1", "x2", "x3"),
        types = c("continuous", "continuous", "continuous")
    )

    synth_data <- hurdlebb:::.copula_generate(
        n = 5000,
        cor_matrix = params$cor_matrix,
        inv_ecdfs = params$inv_ecdfs,
        seed = 42
    )

    # Check Spearman correlations within +/- 0.1
    orig_spearman <- cor(orig_data, method = "spearman")
    synth_spearman <- cor(synth_data, method = "spearman")

    for (i in 1:3) {
        for (j in 1:3) {
            if (i != j) {
                expect_equal(synth_spearman[i, j], orig_spearman[i, j],
                             tolerance = 0.1,
                             label = sprintf("cor(%d,%d)", i, j))
            }
        }
    }
})

test_that("integration: marginal moments preserved (mean within 1 SD)", {
    set.seed(502)
    n_orig <- 1000
    orig_data <- data.frame(
        income = rlnorm(n_orig, meanlog = 3, sdlog = 0.5),
        count  = rpois(n_orig, lambda = 7)
    )

    params <- hurdlebb:::.extract_copula_params(
        data = orig_data,
        var_names = c("income", "count"),
        types = c("continuous", "discrete")
    )

    synth_data <- hurdlebb:::.copula_generate(
        n = 5000,
        cor_matrix = params$cor_matrix,
        inv_ecdfs = params$inv_ecdfs,
        seed = 42
    )

    # Mean within 1 SD of original
    for (vn in c("income", "count")) {
        orig_mean <- mean(orig_data[[vn]])
        orig_sd   <- sd(orig_data[[vn]])
        synth_mean <- mean(synth_data[[vn]])
        expect_true(
            abs(synth_mean - orig_mean) < orig_sd,
            label = sprintf("%s: synth mean %.2f vs orig mean %.2f +/- %.2f",
                            vn, synth_mean, orig_mean, orig_sd)
        )
    }
})

test_that("integration: marginal variance preserved (within factor of 2)", {
    set.seed(503)
    n_orig <- 1000
    orig_data <- data.frame(
        x1 = rnorm(n_orig, mean = 5, sd = 2),
        x2 = rgamma(n_orig, shape = 2, rate = 1)
    )

    params <- hurdlebb:::.extract_copula_params(
        data = orig_data,
        var_names = c("x1", "x2"),
        types = c("continuous", "continuous")
    )

    synth_data <- hurdlebb:::.copula_generate(
        n = 5000,
        cor_matrix = params$cor_matrix,
        inv_ecdfs = params$inv_ecdfs,
        seed = 42
    )

    for (vn in c("x1", "x2")) {
        ratio <- var(synth_data[[vn]]) / var(orig_data[[vn]])
        expect_true(ratio > 0.5 && ratio < 2.0,
                    label = sprintf("%s: variance ratio = %.2f", vn, ratio))
    }
})

test_that("integration: zero-one bounded variables stay in [0, 1]", {
    set.seed(504)
    n_orig <- 500
    orig_data <- data.frame(
        proportion = rbeta(n_orig, 2, 5)
    )

    params <- hurdlebb:::.extract_copula_params(
        data = orig_data,
        var_names = "proportion",
        types = "continuous"
    )

    synth_data <- hurdlebb:::.copula_generate(
        n = 2000,
        cor_matrix = params$cor_matrix,
        inv_ecdfs = params$inv_ecdfs,
        seed = 42
    )

    expect_true(all(synth_data$proportion >= 0 - 1e-10))
    expect_true(all(synth_data$proportion <= 1 + 1e-10))
})

test_that("integration: integer variables remain integer-ish after 'constant'", {
    set.seed(505)
    n_orig <- 500
    orig_data <- data.frame(
        int_var = rpois(n_orig, lambda = 10)
    )

    params <- hurdlebb:::.extract_copula_params(
        data = orig_data,
        var_names = "int_var",
        types = "discrete"
    )

    synth_data <- hurdlebb:::.copula_generate(
        n = 1000,
        cor_matrix = params$cor_matrix,
        inv_ecdfs = params$inv_ecdfs,
        seed = 42
    )

    # All values should be integers (within tolerance)
    expect_true(
        all(abs(synth_data$int_var - round(synth_data$int_var)) < 1e-10),
        label = "integer preservation after constant method"
    )
})

test_that("integration: binary variable preserved with correct proportions", {
    set.seed(506)
    n_orig <- 1000
    true_prob <- 0.65
    orig_data <- data.frame(
        binary = rbinom(n_orig, 1, true_prob)
    )

    params <- hurdlebb:::.extract_copula_params(
        data = orig_data,
        var_names = "binary",
        types = "discrete"
    )

    synth_data <- hurdlebb:::.copula_generate(
        n = 5000,
        cor_matrix = params$cor_matrix,
        inv_ecdfs = params$inv_ecdfs,
        seed = 42
    )

    # All values 0 or 1
    expect_true(all(synth_data$binary == 0 | synth_data$binary == 1))
    # Proportion close to original
    synth_prob <- mean(synth_data$binary)
    orig_prob  <- mean(orig_data$binary)
    expect_equal(synth_prob, orig_prob, tolerance = 0.05)
})

test_that("integration: multivariate pipeline with 5 variables", {
    set.seed(507)
    n_orig <- 800
    orig_data <- data.frame(
        normal   = rnorm(n_orig, 0, 1),
        uniform  = runif(n_orig, 0, 100),
        poisson  = rpois(n_orig, 5),
        binary   = rbinom(n_orig, 1, 0.3),
        lognorm  = rlnorm(n_orig, 1, 0.5)
    )

    params <- hurdlebb:::.extract_copula_params(
        data = orig_data,
        var_names = c("normal", "uniform", "poisson", "binary", "lognorm"),
        types = c("continuous", "continuous", "discrete", "discrete", "continuous")
    )

    synth_data <- hurdlebb:::.copula_generate(
        n = 2000,
        cor_matrix = params$cor_matrix,
        inv_ecdfs = params$inv_ecdfs,
        seed = 42
    )

    expect_equal(ncol(synth_data), 5)
    expect_equal(nrow(synth_data), 2000)
    expect_equal(names(synth_data),
                 c("normal", "uniform", "poisson", "binary", "lognorm"))

    # No missing values in output
    expect_false(any(is.na(synth_data)))

    # Discrete columns should be integer-valued
    expect_true(all(abs(synth_data$poisson - round(synth_data$poisson)) < 1e-10))
    expect_true(all(synth_data$binary == 0 | synth_data$binary == 1))
})

test_that("integration: extract_copula_params -> copula_generate roundtrip preserves shape", {
    # Generate data from a known skewed distribution,
    # verify the synthetic data has similar skewness
    set.seed(508)
    n_orig <- 2000
    orig_data <- data.frame(
        skewed = rexp(n_orig, rate = 0.5)
    )

    params <- hurdlebb:::.extract_copula_params(
        data = orig_data,
        var_names = "skewed",
        types = "continuous"
    )

    synth_data <- hurdlebb:::.copula_generate(
        n = 5000,
        cor_matrix = params$cor_matrix,
        inv_ecdfs = params$inv_ecdfs,
        seed = 42
    )

    # Skewness should be positive for both (exponential is right-skewed)
    orig_skew  <- mean((orig_data$skewed  - mean(orig_data$skewed))^3) /
                  sd(orig_data$skewed)^3
    synth_skew <- mean((synth_data$skewed - mean(synth_data$skewed))^3) /
                  sd(synth_data$skewed)^3

    expect_true(orig_skew > 0)
    expect_true(synth_skew > 0)
    # Skewness should be in the same ballpark
    expect_equal(synth_skew, orig_skew, tolerance = 0.5)
})
