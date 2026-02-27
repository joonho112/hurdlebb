# ============================================================================
# test-hurdle.R — Tests for dhurdle_betabinom, rhurdle_betabinom,
#                 hurdle_mean, hurdle_variance
# ============================================================================


# ============================================================================
# A. dhurdle_betabinom — PMF
# ============================================================================

test_that("dhurdle_betabinom: P(Y=0) = 1 - q", {
    qs <- c(0, 0.1, 0.5, 0.7, 1)
    for (qi in qs) {
        p0 <- dhurdle_betabinom(0, n = 10, q = qi, mu = 0.3, kappa = 5)
        expect_equal(p0, 1 - qi, tolerance = 1e-14,
                     label = sprintf("q = %.1f", qi))
    }
})

test_that("dhurdle_betabinom: PMF sums to 1 over grid", {
    ns     <- c(1, 5, 10, 20)
    qs     <- c(0.1, 0.5, 0.9)
    mus    <- c(0.1, 0.3, 0.5, 0.7, 0.9)
    kappas <- c(0.5, 5, 50)
    for (ni in ns) {
        for (qi in qs) {
            for (mui in mus) {
                for (ki in kappas) {
                    pmf <- dhurdle_betabinom(0:ni, n = ni, q = qi,
                                             mu = mui, kappa = ki)
                    expect_equal(
                        sum(pmf), 1, tolerance = 1e-10,
                        label = sprintf("n=%d, q=%.1f, mu=%.1f, kappa=%.1f",
                                        ni, qi, mui, ki)
                    )
                }
            }
        }
    }
})

test_that("dhurdle_betabinom: log/linear consistency", {
    y <- 0:10
    pmf     <- dhurdle_betabinom(y, n = 10, q = 0.7, mu = 0.4, kappa = 8)
    log_pmf <- dhurdle_betabinom(y, n = 10, q = 0.7, mu = 0.4, kappa = 8,
                                 log = TRUE)
    # For y = 0, pmf = 0.3, log_pmf = log(0.3)
    # For y > 0, check log consistency
    positive <- pmf > 0
    expect_equal(log(pmf[positive]), log_pmf[positive], tolerance = 1e-14)
    # For y where pmf = 0, log_pmf should be -Inf
    expect_true(all(log_pmf[!positive] == -Inf))
})

test_that("dhurdle_betabinom: out-of-support returns 0", {
    expect_equal(dhurdle_betabinom(-1, n = 5, q = 0.7, mu = 0.3, kappa = 10), 0)
    expect_equal(dhurdle_betabinom(6,  n = 5, q = 0.7, mu = 0.3, kappa = 10), 0)
})

test_that("dhurdle_betabinom: q = 0 gives all mass at y = 0", {
    pmf <- dhurdle_betabinom(0:5, n = 5, q = 0, mu = 0.3, kappa = 10)
    expect_equal(pmf[1], 1)      # P(Y=0) = 1
    expect_true(all(pmf[-1] == 0))
})

test_that("dhurdle_betabinom: q = 1 gives zero mass at y = 0", {
    p0 <- dhurdle_betabinom(0, n = 5, q = 1, mu = 0.3, kappa = 10)
    expect_equal(p0, 0)
    # Positive part should sum to 1
    pmf_pos <- dhurdle_betabinom(1:5, n = 5, q = 1, mu = 0.3, kappa = 10)
    expect_equal(sum(pmf_pos), 1, tolerance = 1e-10)
})

test_that("dhurdle_betabinom: q = 1 matches ZT-BB exactly", {
    n <- 8
    mu <- 0.4
    kappa <- 12
    for (y in 1:n) {
        hurdle_p <- dhurdle_betabinom(y, n = n, q = 1, mu = mu, kappa = kappa)
        ztbb_p   <- dztbetabinom(y, n = n, mu = mu, kappa = kappa)
        expect_equal(hurdle_p, ztbb_p, tolerance = 1e-14,
                     label = sprintf("y=%d", y))
    }
})

test_that("dhurdle_betabinom: vectorisation over q", {
    p <- dhurdle_betabinom(c(0, 0), n = 5, q = c(0.3, 0.8),
                           mu = 0.4, kappa = 10)
    expect_equal(p[1], 1 - 0.3, tolerance = 1e-14)
    expect_equal(p[2], 1 - 0.8, tolerance = 1e-14)
})


# ============================================================================
# B. rhurdle_betabinom — Random generation
# ============================================================================

test_that("rhurdle_betabinom: zero rate approximately 1 - q", {
    set.seed(42)
    q <- 0.6
    x <- rhurdle_betabinom(50000, n = 10, q = q, mu = 0.4, kappa = 5)
    zero_rate <- mean(x == 0)
    expect_equal(zero_rate, 1 - q, tolerance = 0.02)
})

test_that("rhurdle_betabinom: all draws in [0, n]", {
    set.seed(123)
    x <- rhurdle_betabinom(10000, n = 15, q = 0.7, mu = 0.4, kappa = 3)
    expect_true(all(x >= 0))
    expect_true(all(x <= 15))
})

test_that("rhurdle_betabinom: LLN for mean", {
    set.seed(42)
    n <- 20
    q <- 0.7
    mu <- 0.3
    kappa <- 10
    x <- rhurdle_betabinom(100000, n = n, q = q, mu = mu, kappa = kappa)
    analytical <- hurdle_mean(n, q, mu, kappa)
    expect_equal(mean(x), analytical, tolerance = 0.05)
})

test_that("rhurdle_betabinom: LLN for variance", {
    set.seed(42)
    n <- 15
    q <- 0.6
    mu <- 0.4
    kappa <- 8
    x <- rhurdle_betabinom(100000, n = n, q = q, mu = mu, kappa = kappa)
    analytical <- hurdle_variance(n, q, mu, kappa)
    expect_equal(var(x), analytical, tolerance = 0.3)
})

test_that("rhurdle_betabinom: q = 0 gives all zeros", {
    set.seed(1)
    x <- rhurdle_betabinom(1000, n = 10, q = 0, mu = 0.3, kappa = 5)
    expect_true(all(x == 0))
})

test_that("rhurdle_betabinom: q = 1 gives no zeros", {
    set.seed(1)
    x <- rhurdle_betabinom(10000, n = 10, q = 1, mu = 0.3, kappa = 5)
    expect_true(all(x >= 1))
})

test_that("rhurdle_betabinom: nn = 0 returns integer(0)", {
    expect_identical(
        rhurdle_betabinom(0, n = 5, q = 0.5, mu = 0.3, kappa = 10),
        integer(0L)
    )
})

test_that("rhurdle_betabinom: vectorised parameters", {
    set.seed(99)
    x <- rhurdle_betabinom(6, n = c(5, 10, 20), q = c(0.5, 0.8, 0.9),
                           mu = c(0.2, 0.5, 0.8), kappa = c(3, 10, 50))
    expect_length(x, 6)
    # Elements with n=5 (positions 1, 4) must be in [0, 5]
    expect_true(all(x[c(1, 4)] >= 0 & x[c(1, 4)] <= 5))
    expect_true(all(x[c(2, 5)] >= 0 & x[c(2, 5)] <= 10))
    expect_true(all(x[c(3, 6)] >= 0 & x[c(3, 6)] <= 20))
})

test_that("rhurdle_betabinom: empirical PMF for small n", {
    set.seed(1)
    n <- 5
    q <- 0.6
    mu <- 0.4
    kappa <- 8
    x <- rhurdle_betabinom(200000, n = n, q = q, mu = mu, kappa = kappa)
    emp_pmf <- tabulate(x + 1L, nbins = n + 1L) / length(x)
    the_pmf <- dhurdle_betabinom(0:n, n = n, q = q, mu = mu, kappa = kappa)
    expect_equal(emp_pmf, the_pmf, tolerance = 0.01)
})


# ============================================================================
# C. hurdle_mean
# ============================================================================

test_that("hurdle_mean: formula = q * n * mu / (1 - p0)", {
    n <- 10
    q <- 0.7
    mu <- 0.3
    kappa <- 5
    p0 <- compute_p0(n, mu, kappa)
    expected <- q * n * mu / (1 - p0)
    expect_equal(hurdle_mean(n, q, mu, kappa), expected, tolerance = 1e-12)
})

test_that("hurdle_mean: matches PMF sum", {
    n <- 8
    q <- 0.65
    mu <- 0.4
    kappa <- 12
    pmf <- dhurdle_betabinom(0:n, n = n, q = q, mu = mu, kappa = kappa)
    pmf_mean <- sum((0:n) * pmf)
    analytical <- hurdle_mean(n, q, mu, kappa)
    expect_equal(analytical, pmf_mean, tolerance = 1e-10)
})

test_that("hurdle_mean: LLN with rhurdle_betabinom", {
    set.seed(42)
    n <- 15
    q <- 0.8
    mu <- 0.5
    kappa <- 10
    x <- rhurdle_betabinom(100000, n = n, q = q, mu = mu, kappa = kappa)
    expect_equal(mean(x), hurdle_mean(n, q, mu, kappa), tolerance = 0.05)
})

test_that("hurdle_mean: q = 0 returns 0", {
    expect_equal(hurdle_mean(10, q = 0, mu = 0.3, kappa = 5), 0)
})

test_that("hurdle_mean: q = 1 returns ZT-BB mean", {
    n <- 10
    mu <- 0.3
    kappa <- 5
    expect_equal(
        hurdle_mean(n, q = 1, mu, kappa),
        compute_ztbb_mean(n, mu, kappa),
        tolerance = 1e-12
    )
})

test_that("hurdle_mean: vectorised", {
    result <- hurdle_mean(
        n = c(5, 10),
        q = c(0.5, 0.8),
        mu = c(0.3, 0.6),
        kappa = c(5, 10)
    )
    expect_length(result, 2)
    expect_equal(result[1], hurdle_mean(5, 0.5, 0.3, 5), tolerance = 1e-12)
    expect_equal(result[2], hurdle_mean(10, 0.8, 0.6, 10), tolerance = 1e-12)
})

test_that("hurdle_mean: increases with q", {
    n <- 10
    mu <- 0.4
    kappa <- 8
    qs <- seq(0, 1, by = 0.1)
    means <- hurdle_mean(n, q = qs, mu = mu, kappa = kappa)
    expect_true(all(diff(means) >= -1e-15))
})


# ============================================================================
# D. hurdle_variance
# ============================================================================

test_that("hurdle_variance: matches PMF formula", {
    n <- 8
    q <- 0.65
    mu <- 0.4
    kappa <- 12
    pmf <- dhurdle_betabinom(0:n, n = n, q = q, mu = mu, kappa = kappa)
    pmf_mean <- sum((0:n) * pmf)
    pmf_var  <- sum(((0:n) - pmf_mean)^2 * pmf)
    analytical <- hurdle_variance(n, q, mu, kappa)
    expect_equal(analytical, pmf_var, tolerance = 1e-10)
})

test_that("hurdle_variance: LLN with rhurdle_betabinom", {
    set.seed(42)
    n <- 15
    q <- 0.6
    mu <- 0.4
    kappa <- 8
    x <- rhurdle_betabinom(200000, n = n, q = q, mu = mu, kappa = kappa)
    analytical <- hurdle_variance(n, q, mu, kappa)
    expect_equal(var(x), analytical, tolerance = 0.3)
})

test_that("hurdle_variance: non-negative", {
    ns     <- c(1, 5, 10)
    qs     <- c(0, 0.3, 0.7, 1)
    mus    <- c(0.1, 0.5, 0.9)
    kappas <- c(0.5, 5, 50)
    for (ni in ns) {
        for (qi in qs) {
            for (mui in mus) {
                for (ki in kappas) {
                    v <- hurdle_variance(ni, qi, mui, ki)
                    expect_true(v >= 0,
                                label = sprintf("n=%d, q=%.1f, mu=%.1f, kappa=%.1f",
                                                ni, qi, mui, ki))
                }
            }
        }
    }
})

test_that("hurdle_variance: q = 0 returns 0", {
    expect_equal(hurdle_variance(10, q = 0, mu = 0.3, kappa = 5), 0)
})

test_that("hurdle_variance: vectorised", {
    result <- hurdle_variance(
        n = c(5, 10),
        q = c(0.5, 0.8),
        mu = c(0.3, 0.6),
        kappa = c(5, 10)
    )
    expect_length(result, 2)
    expect_equal(result[1], hurdle_variance(5, 0.5, 0.3, 5), tolerance = 1e-12)
    expect_equal(result[2], hurdle_variance(10, 0.8, 0.6, 10), tolerance = 1e-12)
})

test_that("hurdle_variance: law of total variance decomposition", {
    # Verify: Var[Y] = E[Var[Y|Z]] + Var[E[Y|Z]]
    # where Z ~ Bern(q), Y|Z=0 = 0, Y|Z=1 ~ ZT-BB
    n <- 10
    q <- 0.7
    mu <- 0.4
    kappa <- 8

    p0 <- compute_p0(n, mu, kappa)
    E_zt <- n * mu / (1 - p0)

    # E_BB[Y^2] = Var_BB + (E_BB)^2
    V_bb  <- n * mu * (1 - mu) * (n + kappa) / (1 + kappa)
    E_bb  <- n * mu
    E2_bb <- V_bb + E_bb^2

    V_zt <- E2_bb / (1 - p0) - E_zt^2

    # E[Var[Y|Z]] = q * V_zt + (1-q) * 0 = q * V_zt
    # Var[E[Y|Z]] = E[E[Y|Z]^2] - (E[E[Y|Z]])^2
    #             = q * E_zt^2 - (q * E_zt)^2 = q*(1-q)*E_zt^2
    total_var <- q * V_zt + q * (1 - q) * E_zt^2

    expect_equal(hurdle_variance(n, q, mu, kappa), total_var,
                 tolerance = 1e-12)
})


# ============================================================================
# E. Consistency checks
# ============================================================================

test_that("consistency: q = 1 hurdle reduces to ZT-BB for all moments", {
    n <- 12
    mu <- 0.35
    kappa <- 15

    # Mean
    expect_equal(
        hurdle_mean(n, q = 1, mu, kappa),
        compute_ztbb_mean(n, mu, kappa),
        tolerance = 1e-12
    )

    # PMF
    for (y in 1:n) {
        expect_equal(
            dhurdle_betabinom(y, n, q = 1, mu, kappa),
            dztbetabinom(y, n, mu, kappa),
            tolerance = 1e-14,
            label = sprintf("y=%d", y)
        )
    }
})

test_that("consistency: hurdle_variance q=1 matches ZT-BB variance from PMF", {
    n <- 8
    mu <- 0.4
    kappa <- 10
    # Compute ZT-BB variance from PMF directly
    ztbb_pmf <- dztbetabinom(1:n, n = n, mu = mu, kappa = kappa)
    ztbb_mean <- sum((1:n) * ztbb_pmf)
    ztbb_var  <- sum(((1:n) - ztbb_mean)^2 * ztbb_pmf)
    # hurdle_variance with q=1 should match
    expect_equal(hurdle_variance(n, q = 1, mu, kappa), ztbb_var,
                 tolerance = 1e-10)
})

test_that("consistency: hurdle mean/variance agree across grid", {
    # For a grid of parameters, verify that:
    # 1. hurdle_mean matches E[Y] from PMF sum
    # 2. hurdle_variance matches Var[Y] from PMF sum
    params <- expand.grid(
        n = c(3, 8),
        q = c(0.3, 0.7),
        mu = c(0.2, 0.5, 0.8),
        kappa = c(2, 20)
    )
    for (i in seq_len(nrow(params))) {
        ni <- params$n[i]
        qi <- params$q[i]
        mui <- params$mu[i]
        ki <- params$kappa[i]

        pmf <- dhurdle_betabinom(0:ni, n = ni, q = qi, mu = mui, kappa = ki)
        pmf_mean <- sum((0:ni) * pmf)
        pmf_var  <- sum(((0:ni) - pmf_mean)^2 * pmf)

        expect_equal(hurdle_mean(ni, qi, mui, ki), pmf_mean,
                     tolerance = 1e-10,
                     label = sprintf("mean: n=%d,q=%.1f,mu=%.1f,k=%.0f",
                                     ni, qi, mui, ki))
        expect_equal(hurdle_variance(ni, qi, mui, ki), pmf_var,
                     tolerance = 1e-10,
                     label = sprintf("var: n=%d,q=%.1f,mu=%.1f,k=%.0f",
                                     ni, qi, mui, ki))
    }
})
