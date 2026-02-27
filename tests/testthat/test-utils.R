# ============================================================================
# test-utils.R — Tests for internal utility functions
#
# All functions are internal (not exported), so use hurdlebb::: prefix.
# ============================================================================


# ============================================================================
# A. inv_logit
# ============================================================================

test_that("inv_logit: special values", {
    expect_equal(hurdlebb:::inv_logit(0), 0.5)
    expect_equal(hurdlebb:::inv_logit(-Inf), 0)
    expect_equal(hurdlebb:::inv_logit(Inf), 1)
})

test_that("inv_logit: matches plogis", {
    x <- seq(-10, 10, length.out = 100)
    expect_identical(hurdlebb:::inv_logit(x), plogis(x))
})

test_that("inv_logit: monotonically increasing", {
    x <- seq(-20, 20, length.out = 1000)
    y <- hurdlebb:::inv_logit(x)
    expect_true(all(diff(y) >= 0))
})

test_that("inv_logit: output in [0, 1] for finite input", {
    x <- c(-100, -10, -1, 0, 1, 10, 100)
    y <- hurdlebb:::inv_logit(x)
    expect_true(all(y >= 0 & y <= 1))
    # Moderate inputs are strictly in (0, 1)
    x_mod <- c(-10, -1, 0, 1, 10)
    y_mod <- hurdlebb:::inv_logit(x_mod)
    expect_true(all(y_mod > 0 & y_mod < 1))
})

test_that("inv_logit: vectorised", {
    x <- c(-2, 0, 2)
    expect_length(hurdlebb:::inv_logit(x), 3)
})


# ============================================================================
# B. logit
# ============================================================================

test_that("logit: special values", {
    expect_equal(hurdlebb:::logit(0.5), 0)
    expect_equal(hurdlebb:::logit(0), -Inf)
    expect_equal(hurdlebb:::logit(1), Inf)
})

test_that("logit: matches qlogis", {
    p <- seq(0.01, 0.99, length.out = 100)
    expect_identical(hurdlebb:::logit(p), qlogis(p))
})

test_that("logit: roundtrip with inv_logit", {
    x <- seq(-10, 10, length.out = 50)
    expect_equal(hurdlebb:::logit(hurdlebb:::inv_logit(x)), x, tolerance = 1e-12)

    p <- seq(0.01, 0.99, length.out = 50)
    expect_equal(hurdlebb:::inv_logit(hurdlebb:::logit(p)), p, tolerance = 1e-12)
})

test_that("logit: monotonically increasing", {
    p <- seq(0.01, 0.99, length.out = 1000)
    y <- hurdlebb:::logit(p)
    expect_true(all(diff(y) > 0))
})

test_that("logit: vectorised", {
    p <- c(0.1, 0.5, 0.9)
    expect_length(hurdlebb:::logit(p), 3)
})


# ============================================================================
# C. normalize_weights
# ============================================================================

test_that("normalize_weights: sum equals length", {
    w <- c(1, 2, 3, 4, 5)
    wn <- hurdlebb:::normalize_weights(w)
    expect_equal(sum(wn), length(w), tolerance = 1e-12)
})

test_that("normalize_weights: preserves ratios", {
    w <- c(2, 4, 8)
    wn <- hurdlebb:::normalize_weights(w)
    expect_equal(wn[2] / wn[1], 2, tolerance = 1e-12)
    expect_equal(wn[3] / wn[1], 4, tolerance = 1e-12)
})

test_that("normalize_weights: uniform weights unchanged", {
    w <- rep(3, 10)
    wn <- hurdlebb:::normalize_weights(w)
    expect_equal(wn, rep(1, 10), tolerance = 1e-12)
})

test_that("normalize_weights: single weight → 1", {
    expect_equal(hurdlebb:::normalize_weights(7), 1, tolerance = 1e-12)
})

test_that("normalize_weights: large vector", {
    set.seed(1)
    w <- runif(10000, min = 0.1, max = 100)
    wn <- hurdlebb:::normalize_weights(w)
    expect_equal(sum(wn), length(w), tolerance = 1e-8)
    expect_true(all(wn > 0))
})

test_that("normalize_weights: errors on zero weights", {
    expect_error(hurdlebb:::normalize_weights(c(1, 0, 3)))
})

test_that("normalize_weights: errors on negative weights", {
    expect_error(hurdlebb:::normalize_weights(c(1, -2, 3)))
})

test_that("normalize_weights: errors on NA", {
    expect_error(hurdlebb:::normalize_weights(c(1, NA, 3)))
})

test_that("normalize_weights: errors on empty vector", {
    expect_error(hurdlebb:::normalize_weights(numeric(0)))
})


# ============================================================================
# D. log1mexp
# ============================================================================

test_that("log1mexp: correctness for known values", {
    # log(1 - exp(-1)) = log(1 - 0.3679...) = log(0.6321...)
    expected <- log(1 - exp(-1))
    expect_equal(hurdlebb:::log1mexp(1), expected, tolerance = 1e-14)

    # log(1 - exp(-5)) ≈ log(1 - 0.00674)
    expected <- log(1 - exp(-5))
    expect_equal(hurdlebb:::log1mexp(5), expected, tolerance = 1e-14)
})

test_that("log1mexp: x = 0 returns -Inf", {
    expect_equal(hurdlebb:::log1mexp(0), -Inf)
})

test_that("log1mexp: large x returns approximately -exp(-x)", {
    # For large x: log(1 - exp(-x)) ≈ -exp(-x) to high relative accuracy
    x <- 20
    result <- hurdlebb:::log1mexp(x)
    approx <- -exp(-x)
    # Both are ~-2.06e-9; relative error is O(exp(-x)) ≈ 2e-9
    expect_equal(result, approx, tolerance = 1e-6)
})

test_that("log1mexp: vectorised", {
    x <- c(0.1, 0.5, 1, 2, 5, 10)
    result <- hurdlebb:::log1mexp(x)
    expected <- log(1 - exp(-x))
    expect_equal(result, expected, tolerance = 1e-12)
})

test_that("log1mexp: small x near zero (high relative error regime)", {
    # This is exactly where the two-branch matters: x <= log(2)
    x <- c(1e-10, 1e-8, 1e-5, 1e-3, 0.1, 0.5, log(2))
    result <- hurdlebb:::log1mexp(x)
    # Compare with high-precision reference: log(1 - exp(-x))
    # For very small x: 1 - exp(-x) ≈ x, so log1mexp(x) ≈ log(x)
    for (i in seq_along(x)) {
        ref <- log(-expm1(-x[i]))  # the accurate branch
        expect_equal(result[i], ref, tolerance = 1e-14,
                     label = sprintf("x = %g", x[i]))
    }
})

test_that("log1mexp: branch boundary at log(2)", {
    # Both branches should give the same answer at the boundary
    x <- log(2)
    # Small branch: log(-expm1(-x))
    small_ans <- log(-expm1(-x))
    # Large branch: log1p(-exp(-x))
    large_ans <- log1p(-exp(-x))
    result <- hurdlebb:::log1mexp(x)
    # All three should match closely
    expect_equal(result, small_ans, tolerance = 1e-14)
    expect_equal(small_ans, large_ans, tolerance = 1e-12)
})

test_that("log1mexp: monotonically increasing", {
    x <- seq(0.01, 20, length.out = 1000)
    y <- hurdlebb:::log1mexp(x)
    expect_true(all(diff(y) >= 0))
})

test_that("log1mexp: NA propagation", {
    result <- hurdlebb:::log1mexp(c(1, NA, 3))
    expect_equal(result[1], log(1 - exp(-1)), tolerance = 1e-14)
    expect_true(is.na(result[2]))
    expect_equal(result[3], log(1 - exp(-3)), tolerance = 1e-14)
})

test_that("log1mexp: output always <= 0", {
    x <- c(0.001, 0.01, 0.1, 0.5, 1, 2, 5, 10, 50)
    result <- hurdlebb:::log1mexp(x)
    expect_true(all(result <= 0))
})
