# ============================================================================
# test-distributions.R — Tests for dbetabinom, pbetabinom, rbetabinom,
#                         dztbetabinom, rztbetabinom, compute_p0,
#                         compute_ztbb_mean
# ============================================================================


# ============================================================================
# A. dbetabinom — PMF of Beta-Binomial
# ============================================================================

test_that("dbetabinom: PMF sums to 1 over grid of parameters", {
    ns     <- c(1, 5, 10, 20, 50)
    mus    <- c(0.1, 0.3, 0.5, 0.7, 0.9)
    kappas <- c(0.5, 2, 10, 50, 200)
    for (ni in ns) {
        for (mui in mus) {
            for (ki in kappas) {
                pmf <- dbetabinom(0:ni, n = ni, mu = mui, kappa = ki)
                expect_equal(sum(pmf), 1, tolerance = 1e-10,
                             label = sprintf("n=%d, mu=%.1f, kappa=%.0f", ni, mui, ki))
            }
        }
    }
})

test_that("dbetabinom: hand-computed for n = 1 (Bernoulli case)", {
    # BetaBin(y|1, mu, kappa) = Bern(mu) regardless of kappa
    # P(Y=0) = (1-mu), P(Y=1) = mu
    for (ki in c(1, 5, 100)) {
        p0 <- dbetabinom(0, n = 1, mu = 0.3, kappa = ki)
        p1 <- dbetabinom(1, n = 1, mu = 0.3, kappa = ki)
        expect_equal(p0, 0.7, tolerance = 1e-10)
        expect_equal(p1, 0.3, tolerance = 1e-10)
    }
})

test_that("dbetabinom: hand-computed for n = 2", {
    # a = 0.3 * 10 = 3, b = 0.7 * 10 = 7
    # P(Y=0) = B(7+2, 3)/B(7, 3) = B(9,3)/B(7,3)
    a <- 3
    b <- 7
    n <- 2
    p0_expected <- exp(lbeta(b + n, a) - lbeta(b, a))  # choose(2,0) = 1
    p1_expected <- choose(2, 1) * exp(lbeta(1 + a, 1 + b) - lbeta(a, b))
    p2_expected <- exp(lbeta(2 + a, b) - lbeta(a, b))  # choose(2,2) = 1

    expect_equal(dbetabinom(0, 2, 0.3, 10), p0_expected, tolerance = 1e-12)
    expect_equal(dbetabinom(1, 2, 0.3, 10), p1_expected, tolerance = 1e-12)
    expect_equal(dbetabinom(2, 2, 0.3, 10), p2_expected, tolerance = 1e-12)
})

test_that("dbetabinom: log/linear consistency", {
    y <- 0:10
    pmf     <- dbetabinom(y, n = 10, mu = 0.4, kappa = 5)
    log_pmf <- dbetabinom(y, n = 10, mu = 0.4, kappa = 5, log = TRUE)
    expect_equal(log(pmf), log_pmf, tolerance = 1e-14)
})

test_that("dbetabinom: out-of-support returns 0 / -Inf", {
    expect_equal(dbetabinom(-1, n = 5, mu = 0.3, kappa = 10), 0)
    expect_equal(dbetabinom(6,  n = 5, mu = 0.3, kappa = 10), 0)
    expect_equal(dbetabinom(-1, n = 5, mu = 0.3, kappa = 10, log = TRUE), -Inf)
    expect_equal(dbetabinom(6,  n = 5, mu = 0.3, kappa = 10, log = TRUE), -Inf)
})

test_that("dbetabinom: mu = 0 gives point mass at y = 0", {
    pmf <- dbetabinom(0:5, n = 5, mu = 0, kappa = 10)
    expect_equal(pmf[1], 1)
    expect_true(all(pmf[-1] == 0))
})

test_that("dbetabinom: mu = 1 gives point mass at y = n", {
    pmf <- dbetabinom(0:5, n = 5, mu = 1, kappa = 10)
    expect_equal(pmf[6], 1)    # y = 5 = n
    expect_true(all(pmf[1:5] == 0))
})

test_that("dbetabinom: vectorisation works correctly", {
    # Vectorised over mu
    p <- dbetabinom(c(0, 1), n = 2, mu = c(0.3, 0.7), kappa = 5)
    expect_length(p, 2)
    expect_equal(p[1], dbetabinom(0, 2, 0.3, 5))
    expect_equal(p[2], dbetabinom(1, 2, 0.7, 5))
})

test_that("dbetabinom: symmetry P(Y=y|mu) = P(Y=n-y|1-mu)", {
    n <- 8
    mu <- 0.35
    kappa <- 12
    for (y in 0:n) {
        p1 <- dbetabinom(y,     n, mu,     kappa)
        p2 <- dbetabinom(n - y, n, 1 - mu, kappa)
        expect_equal(p1, p2, tolerance = 1e-12,
                     label = sprintf("y=%d", y))
    }
})

test_that("dbetabinom: validates inputs", {
    expect_error(dbetabinom(0, n = 5, mu = -0.1, kappa = 10))
    expect_error(dbetabinom(0, n = 5, mu = 1.1,  kappa = 10))
    expect_error(dbetabinom(0, n = 5, mu = 0.3,  kappa = -1))
    expect_error(dbetabinom(0, n = -1, mu = 0.3, kappa = 10))
})


# ============================================================================
# B. pbetabinom — CDF of Beta-Binomial
# ============================================================================

test_that("pbetabinom: CDF at boundaries", {
    expect_equal(pbetabinom(-1,  n = 5, mu = 0.3, kappa = 10), 0)
    expect_equal(pbetabinom(5,   n = 5, mu = 0.3, kappa = 10), 1)
    expect_equal(pbetabinom(100, n = 5, mu = 0.3, kappa = 10), 1)
})

test_that("pbetabinom: monotonically non-decreasing", {
    cdf <- pbetabinom(0:10, n = 10, mu = 0.5, kappa = 5)
    expect_true(all(diff(cdf) >= -1e-15))
})

test_that("pbetabinom: consistent with cumsum of PMF", {
    n <- 8
    mu <- 0.4
    kappa <- 7
    pmf <- dbetabinom(0:n, n = n, mu = mu, kappa = kappa)
    cdf_expected <- cumsum(pmf)
    cdf_actual   <- pbetabinom(0:n, n = n, mu = mu, kappa = kappa)
    expect_equal(cdf_actual, cdf_expected, tolerance = 1e-10)
})

test_that("pbetabinom: upper tail = 1 - lower tail", {
    p_lower <- pbetabinom(3, n = 10, mu = 0.5, kappa = 5, lower.tail = TRUE)
    p_upper <- pbetabinom(3, n = 10, mu = 0.5, kappa = 5, lower.tail = FALSE)
    expect_equal(p_lower + p_upper, 1, tolerance = 1e-12)
})


# ============================================================================
# C. rbetabinom — Random generation
# ============================================================================

test_that("rbetabinom: all draws in [0, n]", {
    set.seed(123)
    x <- rbetabinom(10000, n = 15, mu = 0.4, kappa = 3)
    expect_true(all(x >= 0))
    expect_true(all(x <= 15))
})

test_that("rbetabinom: LLN for mean", {
    set.seed(42)
    n <- 20
    mu <- 0.3
    kappa <- 10
    x <- rbetabinom(50000, n = n, mu = mu, kappa = kappa)
    expect_equal(mean(x), n * mu, tolerance = 0.1)
})

test_that("rbetabinom: LLN for variance", {
    set.seed(42)
    n <- 20
    mu <- 0.3
    kappa <- 10
    x <- rbetabinom(50000, n = n, mu = mu, kappa = kappa)
    # Var_BB = n*mu*(1-mu)*(n+kappa)/(1+kappa)
    var_expected <- n * mu * (1 - mu) * (n + kappa) / (1 + kappa)
    expect_equal(var(x), var_expected, tolerance = 0.3)
})

test_that("rbetabinom: empirical PMF matches dbetabinom for small n", {
    set.seed(1)
    n <- 5
    mu <- 0.4
    kappa <- 8
    x <- rbetabinom(100000, n = n, mu = mu, kappa = kappa)
    emp_pmf <- tabulate(x + 1L, nbins = n + 1L) / length(x)
    the_pmf <- dbetabinom(0:n, n = n, mu = mu, kappa = kappa)
    expect_equal(emp_pmf, the_pmf, tolerance = 0.01)
})

test_that("rbetabinom: nn = 0 returns integer(0)", {
    expect_identical(rbetabinom(0, n = 5, mu = 0.3, kappa = 10), integer(0L))
})

test_that("rbetabinom: vectorised parameters", {
    set.seed(99)
    x <- rbetabinom(4, n = c(5, 10), mu = c(0.2, 0.8), kappa = c(3, 20))
    expect_length(x, 4)
    expect_true(all(x[c(1, 3)] <= 5))
    expect_true(all(x[c(2, 4)] <= 10))
})


# ============================================================================
# D. dztbetabinom — Zero-Truncated Beta-Binomial PMF
# ============================================================================

test_that("dztbetabinom: PMF sums to 1 over grid", {
    ns     <- c(1, 5, 10, 20)
    mus    <- c(0.1, 0.3, 0.5, 0.7, 0.9)
    kappas <- c(0.5, 2, 10, 50)
    for (ni in ns) {
        for (mui in mus) {
            for (ki in kappas) {
                pmf <- dztbetabinom(1:ni, n = ni, mu = mui, kappa = ki)
                expect_equal(sum(pmf), 1, tolerance = 1e-10,
                             label = sprintf("n=%d, mu=%.1f, kappa=%.0f", ni, mui, ki))
            }
        }
    }
})

test_that("dztbetabinom: y = 0 returns 0 / -Inf", {
    expect_equal(dztbetabinom(0, n = 5, mu = 0.3, kappa = 10), 0)
    expect_equal(dztbetabinom(0, n = 5, mu = 0.3, kappa = 10, log = TRUE), -Inf)
})

test_that("dztbetabinom: equals BB/(1 - p0) for y >= 1", {
    n <- 8
    mu <- 0.4
    kappa <- 12
    p0 <- compute_p0(n, mu, kappa)
    for (y in 1:n) {
        zt_pmf <- dztbetabinom(y, n, mu, kappa)
        bb_pmf <- dbetabinom(y, n, mu, kappa)
        expect_equal(zt_pmf, bb_pmf / (1 - p0), tolerance = 1e-12,
                     label = sprintf("y=%d", y))
    }
})

test_that("dztbetabinom: log/linear consistency", {
    y <- 1:10
    pmf     <- dztbetabinom(y, n = 10, mu = 0.4, kappa = 5)
    log_pmf <- dztbetabinom(y, n = 10, mu = 0.4, kappa = 5, log = TRUE)
    expect_equal(log(pmf), log_pmf, tolerance = 1e-14)
})

test_that("dztbetabinom: rejects mu = 0 (checkmate)", {
    expect_error(dztbetabinom(1, n = 5, mu = 0, kappa = 10))
})

test_that("dztbetabinom: out-of-support returns 0", {
    expect_equal(dztbetabinom(0,  n = 5, mu = 0.3, kappa = 10), 0)
    expect_equal(dztbetabinom(6,  n = 5, mu = 0.3, kappa = 10), 0)
    expect_equal(dztbetabinom(-1, n = 5, mu = 0.3, kappa = 10), 0)
})


# ============================================================================
# E. rztbetabinom — Random generation from ZT-Beta-Binomial
# ============================================================================

test_that("rztbetabinom: all draws in [1, n]", {
    set.seed(123)
    x <- rztbetabinom(10000, n = 15, mu = 0.4, kappa = 3)
    expect_true(all(x >= 1))
    expect_true(all(x <= 15))
})

test_that("rztbetabinom: no zeros", {
    set.seed(42)
    x <- rztbetabinom(10000, n = 10, mu = 0.1, kappa = 2)
    expect_true(all(x > 0))
})

test_that("rztbetabinom: LLN for mean", {
    set.seed(42)
    n <- 20
    mu <- 0.3
    kappa <- 10
    x <- rztbetabinom(50000, n = n, mu = mu, kappa = kappa)
    analytical_mean <- compute_ztbb_mean(n, mu, kappa)
    expect_equal(mean(x), analytical_mean, tolerance = 0.1)
})

test_that("rztbetabinom: nn = 0 returns integer(0)", {
    expect_identical(rztbetabinom(0, n = 5, mu = 0.3, kappa = 10), integer(0L))
})

test_that("rztbetabinom: vectorised parameters", {
    set.seed(99)
    x <- rztbetabinom(6, n = c(5, 10, 20), mu = c(0.2, 0.5, 0.8),
                      kappa = c(3, 10, 50))
    expect_length(x, 6)
    expect_true(all(x[c(1, 4)] >= 1 & x[c(1, 4)] <= 5))
    expect_true(all(x[c(2, 5)] >= 1 & x[c(2, 5)] <= 10))
    expect_true(all(x[c(3, 6)] >= 1 & x[c(3, 6)] <= 20))
})


# ============================================================================
# F. compute_p0
# ============================================================================

test_that("compute_p0: matches lgamma formula", {
    n <- 10
    mu <- 0.3
    kappa <- 8
    b <- (1 - mu) * kappa
    expected <- exp(lgamma(b + n) + lgamma(kappa) - lgamma(b) - lgamma(kappa + n))
    expect_equal(compute_p0(n, mu, kappa), expected, tolerance = 1e-14)
})

test_that("compute_p0: matches dbetabinom(0, ...)", {
    ns     <- c(1, 5, 10, 30)
    mus    <- c(0.1, 0.5, 0.9)
    kappas <- c(1, 10, 100)
    for (ni in ns) {
        for (mui in mus) {
            for (ki in kappas) {
                p0_direct <- compute_p0(ni, mui, ki)
                p0_pmf    <- dbetabinom(0, n = ni, mu = mui, kappa = ki)
                expect_equal(p0_direct, p0_pmf, tolerance = 1e-12,
                             label = sprintf("n=%d, mu=%.1f, kappa=%.0f", ni, mui, ki))
            }
        }
    }
})

test_that("compute_p0: Binomial limit (kappa -> Inf)", {
    # As kappa -> Inf, BetaBin -> Binom, p0 -> (1-mu)^n
    n <- 10
    mu <- 0.4
    kappa <- 1e8
    expected <- (1 - mu)^n
    expect_equal(compute_p0(n, mu, kappa), expected, tolerance = 1e-4)
})

test_that("compute_p0: result in [0, 1]", {
    ns     <- c(1, 10, 100, 1000)
    mus    <- c(0.01, 0.5, 0.99)
    kappas <- c(0.1, 10, 1000)
    for (ni in ns) {
        for (mui in mus) {
            for (ki in kappas) {
                p0 <- compute_p0(ni, mui, ki)
                expect_true(p0 >= 0 && p0 <= 1,
                            label = sprintf("n=%d, mu=%.2f, kappa=%.1f: p0=%.6f",
                                            ni, mui, ki, p0))
            }
        }
    }
})

test_that("compute_p0: boundary mu = 0 returns 1", {
    expect_equal(compute_p0(5, mu = 0, kappa = 10), 1)
    expect_equal(compute_p0(100, mu = 0, kappa = 1), 1)
})

test_that("compute_p0: boundary mu = 1 returns 0", {
    expect_equal(compute_p0(5, mu = 1, kappa = 10), 0)
    expect_equal(compute_p0(100, mu = 1, kappa = 1), 0)
})

test_that("compute_p0: boundary n = 0 returns 1", {
    expect_equal(compute_p0(0, mu = 0.3, kappa = 10), 1)
    expect_equal(compute_p0(0, mu = 0.9, kappa = 50), 1)
})


# ============================================================================
# G. compute_ztbb_mean
# ============================================================================

test_that("compute_ztbb_mean: formula = n*mu/(1-p0)", {
    n <- 15
    mu <- 0.4
    kappa <- 7
    p0 <- compute_p0(n, mu, kappa)
    expected <- n * mu / (1 - p0)
    expect_equal(compute_ztbb_mean(n, mu, kappa), expected, tolerance = 1e-12)
})

test_that("compute_ztbb_mean: inflated relative to n*mu", {
    n <- 10
    mu <- 0.3
    kappa <- 5
    ztbb_mean <- compute_ztbb_mean(n, mu, kappa)
    expect_true(ztbb_mean > n * mu)
})

test_that("compute_ztbb_mean: LLN with rztbetabinom", {
    set.seed(42)
    n <- 20
    mu <- 0.3
    kappa <- 10
    x <- rztbetabinom(100000, n = n, mu = mu, kappa = kappa)
    analytical <- compute_ztbb_mean(n, mu, kappa)
    expect_equal(mean(x), analytical, tolerance = 0.05)
})


# ============================================================================
# H. Edge cases
# ============================================================================

test_that("edge: n = 0 for dbetabinom", {
    # Only y = 0 is in support
    expect_equal(dbetabinom(0, n = 0, mu = 0.5, kappa = 10), 1)
    expect_equal(dbetabinom(1, n = 0, mu = 0.5, kappa = 10), 0)
})

test_that("edge: n = 1 for dztbetabinom (only y = 1)", {
    pmf <- dztbetabinom(1, n = 1, mu = 0.5, kappa = 10)
    expect_equal(pmf, 1, tolerance = 1e-12)
})

test_that("edge: very large kappa approaches Binomial", {
    n <- 10
    mu <- 0.4
    kappa <- 1e7
    bb_pmf   <- dbetabinom(0:n, n = n, mu = mu, kappa = kappa)
    binom_pmf <- dbinom(0:n, size = n, prob = mu)
    expect_equal(bb_pmf, binom_pmf, tolerance = 1e-4)
})

test_that("edge: mu near 0 — most mass at y = 0", {
    mu <- 0.001
    n <- 10
    kappa <- 5
    p0 <- dbetabinom(0, n = n, mu = mu, kappa = kappa)
    expect_true(p0 > 0.95)
})

test_that("edge: mu near 1 — most mass at y = n", {
    mu <- 0.999
    n <- 10
    kappa <- 5
    pn <- dbetabinom(n, n = n, mu = mu, kappa = kappa)
    expect_true(pn > 0.95)
})

test_that("edge: n = 1000 (large n), PMF sums to 1", {
    n <- 1000
    mu <- 0.5
    kappa <- 50
    # Only check boundary + a few points; full sum via pbetabinom
    cdf_n <- pbetabinom(n, n = n, mu = mu, kappa = kappa)
    expect_equal(cdf_n, 1, tolerance = 1e-8)
})

test_that("edge: very small kappa (high overdispersion)", {
    n <- 10
    mu <- 0.5
    kappa <- 0.01  # extreme overdispersion
    pmf <- dbetabinom(0:n, n = n, mu = mu, kappa = kappa)
    expect_equal(sum(pmf), 1, tolerance = 1e-10)
    expect_true(all(pmf >= 0))
})
