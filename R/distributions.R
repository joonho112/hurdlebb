# ============================================================================
# distributions.R — Beta-Binomial and Zero-Truncated Beta-Binomial
#
# Core probability functions for the hurdlebb package.
#
# Contents:
#   1. dbetabinom       — PMF of Beta-Binomial
#   2. pbetabinom       — CDF of Beta-Binomial
#   3. rbetabinom       — Random generation from Beta-Binomial
#   4. dztbetabinom     — PMF of Zero-Truncated Beta-Binomial
#   5. rztbetabinom     — Random generation from ZT-Beta-Binomial
#   6. compute_p0       — P(Y = 0) under Beta-Binomial
#   7. compute_ztbb_mean — E[Y | Y > 0] under Beta-Binomial
# ============================================================================


# -- 1. dbetabinom -----------------------------------------------------------

#' Probability mass function of the Beta-Binomial distribution
#'
#' Computes the PMF of the Beta-Binomial distribution parameterised by
#' mean \eqn{\mu \in [0, 1]} and concentration \eqn{\kappa > 0}. The
#' shape parameters are \eqn{a = \mu \kappa} and \eqn{b = (1 - \mu) \kappa}.
#'
#' @details
#' The PMF is
#' \deqn{
#'   f(y \mid n, \mu, \kappa)
#'   = \binom{n}{y} \frac{B(y + a,\; n - y + b)}{B(a, b)},
#'   \quad y \in \{0, 1, \ldots, n\},
#' }
#' where \eqn{a = \mu \kappa} and \eqn{b = (1 - \mu) \kappa}.
#'
#' **Boundary handling.** When \eqn{\mu = 0} the Beta-Binomial degenerates
#' to a point mass at \eqn{y = 0} (all probability on zero). When
#' \eqn{\mu = 1} it degenerates to a point mass at \eqn{y = n}. These
#' boundaries are handled with explicit branches because the \code{lbeta}
#' formula produces \code{NaN} when \eqn{a = 0} or \eqn{b = 0}.
#'
#' All arguments are recycled to common length via \code{rep_len}.
#'
#' @param y Integer vector of observed counts.
#' @param n Integer vector of trial sizes (\eqn{n \ge 0}).
#' @param mu Numeric vector of means, each in \eqn{[0, 1]}.
#' @param kappa Numeric vector of concentrations, each \eqn{> 0}.
#' @param log Logical; if \code{TRUE}, return log-probabilities.
#'
#' @return A numeric vector of (log-)probabilities, same length as the
#'   recycled inputs.
#'
#' @references
#' Ghosal, S., Ghosh, S., and Moores, M. (2020).
#' \dQuote{Hierarchical beta-binomial models for batch effects in
#' cytometry data.} \emph{Journal of the Royal Statistical Society:
#' Series A}, \strong{183}(4), 1579--1601.
#'
#' @examples
#' # Standard usage
#' dbetabinom(0:5, n = 5, mu = 0.3, kappa = 10)
#' dbetabinom(0:5, n = 5, mu = 0.3, kappa = 10, log = TRUE)
#'
#' # Boundary: mu = 0 gives point mass at y = 0
#' dbetabinom(0:3, n = 3, mu = 0, kappa = 5)
#'
#' # Boundary: mu = 1 gives point mass at y = n
#' dbetabinom(0:3, n = 3, mu = 1, kappa = 5)
#'
#' @family distributions
#' @export
dbetabinom <- function(y, n, mu, kappa, log = FALSE) {

    # -- Input validation ----------------------------------------------------
    assert_flag(log)
    assert_integerish(y)
    assert_integerish(n, lower = 0)
    assert_numeric(mu, lower = 0, upper = 1)
    assert_numeric(kappa, lower = .Machine$double.eps)

    # -- Recycle to common length --------------------------------------------
    len <- max(length(y), length(n), length(mu), length(kappa))
    y     <- rep_len(as.numeric(y), len)
    n     <- rep_len(as.numeric(n), len)
    mu    <- rep_len(mu, len)
    kappa <- rep_len(kappa, len)

    # -- Allocate result -----------------------------------------------------
    log_f <- rep(NA_real_, len)

    # -- Mark invalid support entries ----------------------------------------
    invalid <- (y < 0) | (y > n) | (y != round(y))

    # -- Boundary: mu = 0 → point mass at y = 0 -----------------------------
    at_mu0 <- (mu == 0) & !invalid
    log_f[at_mu0 & (y == 0)] <- 0        # log(1) = 0
    log_f[at_mu0 & (y != 0)] <- -Inf     # log(0) = -Inf

    # -- Boundary: mu = 1 → point mass at y = n -----------------------------
    at_mu1 <- (mu == 1) & !invalid
    log_f[at_mu1 & (y == n)] <- 0
    log_f[at_mu1 & (y != n)] <- -Inf

    # -- Main formula for interior mu ∈ (0, 1) ------------------------------
    interior <- !invalid & (mu > 0) & (mu < 1) & is.na(log_f)
    if (any(interior)) {
        a_i <- mu[interior] * kappa[interior]
        b_i <- (1 - mu[interior]) * kappa[interior]
        log_f[interior] <- lchoose(n[interior], y[interior]) +
            lbeta(y[interior] + a_i, n[interior] - y[interior] + b_i) -
            lbeta(a_i, b_i)
    }

    # -- Invalid entries -----------------------------------------------------
    log_f[invalid] <- -Inf

    if (log) log_f else exp(log_f)
}


# -- 2. pbetabinom -----------------------------------------------------------

#' Cumulative distribution function of the Beta-Binomial distribution
#'
#' Computes the CDF \eqn{P(Y \le q)} by summing the PMF from 0 to
#' \code{floor(q)}, element-wise.
#'
#' @details
#' Numerical stability is maintained by performing the summation in the
#' log domain using the log-sum-exp trick.
#'
#' @param q Numeric vector of quantiles.
#' @param n Integer vector of trial sizes (\eqn{n \ge 0}).
#' @param mu Numeric vector of means, each in \eqn{[0, 1]}.
#' @param kappa Numeric vector of concentrations, each \eqn{> 0}.
#' @param lower.tail Logical; if \code{TRUE} (default), returns
#'   \eqn{P(Y \le q)}; otherwise \eqn{P(Y > q)}.
#'
#' @return A numeric vector of probabilities.
#'
#' @references
#' Ghosal, S., Ghosh, S., and Moores, M. (2020).
#' \dQuote{Hierarchical beta-binomial models for batch effects in
#' cytometry data.} \emph{Journal of the Royal Statistical Society:
#' Series A}, \strong{183}(4), 1579--1601.
#'
#' @examples
#' # CDF at each value
#' pbetabinom(0:5, n = 5, mu = 0.3, kappa = 10)
#'
#' # Upper tail
#' pbetabinom(2, n = 5, mu = 0.3, kappa = 10, lower.tail = FALSE)
#'
#' @family distributions
#' @export
pbetabinom <- function(q, n, mu, kappa, lower.tail = TRUE) {

    assert_flag(lower.tail)
    assert_numeric(mu, lower = 0, upper = 1)
    assert_numeric(kappa, lower = .Machine$double.eps)
    assert_integerish(n, lower = 0)

    # Recycle to common length
    len <- max(length(q), length(n), length(mu), length(kappa))
    q     <- rep_len(as.numeric(q), len)
    n     <- rep_len(as.numeric(n), len)
    mu    <- rep_len(mu, len)
    kappa <- rep_len(kappa, len)

    # Element-wise CDF computation
    cdf <- vapply(seq_len(len), function(i) {
        qi <- floor(q[i])
        if (qi < 0) return(0)
        if (qi >= n[i]) return(1)

        # Sum PMF from 0 to qi using log-sum-exp
        yy <- 0:qi
        log_pmfs <- dbetabinom(yy, n = n[i], mu = mu[i],
                               kappa = kappa[i], log = TRUE)

        # Log-sum-exp for numerical stability
        max_lp <- max(log_pmfs)
        if (is.infinite(max_lp) && max_lp < 0) return(0)
        exp(max_lp + log(sum(exp(log_pmfs - max_lp))))
    }, numeric(1L))

    if (lower.tail) cdf else 1 - cdf
}


# -- 3. rbetabinom -----------------------------------------------------------

#' Random generation from the Beta-Binomial distribution
#'
#' Generates random draws via the hierarchical representation:
#' \eqn{p \sim \text{Beta}(a, b)}, \eqn{Y \mid p \sim \text{Binomial}(n, p)}.
#'
#' @details
#' All parameter vectors are recycled to length \code{nn}. The generation
#' is fully vectorised (no loop).
#'
#' @param nn Number of draws to generate (positive integer).
#' @param n Integer vector of trial sizes (\eqn{n \ge 0}).
#' @param mu Numeric vector of means, each in \eqn{[0, 1]}.
#' @param kappa Numeric vector of concentrations, each \eqn{> 0}.
#'
#' @return An integer vector of length \code{nn}.
#'
#' @references
#' Ghosal, S., Ghosh, S., and Moores, M. (2020).
#' \dQuote{Hierarchical beta-binomial models for batch effects in
#' cytometry data.} \emph{Journal of the Royal Statistical Society:
#' Series A}, \strong{183}(4), 1579--1601.
#'
#' @examples
#' set.seed(42)
#' x <- rbetabinom(1000, n = 10, mu = 0.3, kappa = 5)
#' hist(x, breaks = -0.5:10.5)
#' mean(x)  # approximately 3
#'
#' # Vectorised parameters
#' rbetabinom(4, n = c(5, 10), mu = c(0.2, 0.8), kappa = c(3, 20))
#'
#' @family distributions
#' @export
rbetabinom <- function(nn, n, mu, kappa) {

    assert_count(nn)
    if (nn == 0L) return(integer(0L))

    assert_integerish(n, lower = 0)
    assert_numeric(mu, lower = 0, upper = 1)
    assert_numeric(kappa, lower = .Machine$double.eps)

    n     <- rep_len(as.integer(n), nn)
    mu    <- rep_len(mu, nn)
    kappa <- rep_len(kappa, nn)

    a <- mu * kappa
    b <- (1 - mu) * kappa

    # Hierarchical: p ~ Beta(a, b), Y | p ~ Binom(n, p)
    p <- rbeta(nn, shape1 = a, shape2 = b)
    rbinom(nn, size = n, prob = p)
}


# -- 4. dztbetabinom ---------------------------------------------------------

#' Probability mass function of the Zero-Truncated Beta-Binomial
#'
#' Computes the PMF of the zero-truncated Beta-Binomial distribution
#' for \eqn{y \in \{1, 2, \ldots, n\}}.
#'
#' @details
#' The zero-truncated PMF is
#' \deqn{
#'   f_{\text{ZT}}(y \mid n, \mu, \kappa)
#'   = \frac{f_{\text{BB}}(y \mid n, \mu, \kappa)}{1 - p_0},
#'   \quad y \in \{1, \ldots, n\},
#' }
#' where \eqn{p_0 = P(Y = 0)} under the untruncated Beta-Binomial.
#'
#' Numerical stability is ensured by computing in the log domain and
#' using \code{log1mexp()} for \eqn{\log(1 - p_0)}.
#'
#' @param y Integer vector of observed counts (\eqn{y \ge 1}).
#' @param n Integer vector of trial sizes (\eqn{n \ge 1}).
#' @param mu Numeric vector of means, each in
#'   \eqn{(\varepsilon, 1 - \varepsilon)} where \eqn{\varepsilon} is
#'   machine epsilon.
#' @param kappa Numeric vector of concentrations, each \eqn{> 0}.
#' @param log Logical; if \code{TRUE}, return log-probabilities.
#'
#' @return A numeric vector of (log-)probabilities.
#'
#' @references
#' Ghosal, S., Ghosh, S., and Moores, M. (2020).
#' \dQuote{Hierarchical beta-binomial models for batch effects in
#' cytometry data.} \emph{Journal of the Royal Statistical Society:
#' Series A}, \strong{183}(4), 1579--1601.
#'
#' @examples
#' # ZT-BB PMF (support starts at 1)
#' dztbetabinom(1:5, n = 5, mu = 0.3, kappa = 10)
#'
#' # Sums to 1 over the support
#' sum(dztbetabinom(1:5, n = 5, mu = 0.3, kappa = 10))
#'
#' # Log scale
#' dztbetabinom(1:5, n = 5, mu = 0.3, kappa = 10, log = TRUE)
#'
#' @family distributions
#' @export
dztbetabinom <- function(y, n, mu, kappa, log = FALSE) {

    assert_flag(log)
    assert_integerish(y)
    assert_integerish(n, lower = 1)
    eps <- .Machine$double.eps
    assert_numeric(mu, lower = eps, upper = 1 - eps)
    assert_numeric(kappa, lower = eps)

    # Recycle
    len <- max(length(y), length(n), length(mu), length(kappa))
    y     <- rep_len(as.numeric(y), len)
    n     <- rep_len(as.numeric(n), len)
    mu    <- rep_len(mu, len)
    kappa <- rep_len(kappa, len)

    # log P(Y = 0) via lgamma
    b <- (1 - mu) * kappa
    log_p0 <- lgamma(b + n) + lgamma(kappa) - lgamma(b) - lgamma(kappa + n)

    # log(1 - p0) via log1mexp (numerically stable)
    log_1mp0 <- log1mexp(-log_p0)

    # BB log-PMF (handles out-of-support internally)
    log_f_bb <- dbetabinom(y, n, mu, kappa, log = TRUE)

    # ZT-BB log-PMF
    log_f_zt <- log_f_bb - log_1mp0

    # Enforce support: y in {1, ..., n}
    invalid <- (y < 1) | (y > n) | (y != round(y))
    log_f_zt[invalid] <- -Inf

    if (log) log_f_zt else exp(log_f_zt)
}


# -- 5. rztbetabinom ----------------------------------------------------------

#' Random generation from the Zero-Truncated Beta-Binomial
#'
#' Generates random draws from the zero-truncated Beta-Binomial via
#' rejection sampling on the standard Beta-Binomial.
#'
#' @details
#' The function is vectorised over \code{n}, \code{mu}, and \code{kappa}
#' using a parameter-grouping strategy. Unique parameter combinations are
#' identified, and rejection sampling is performed once per group with
#' smart batch sizing based on the analytical \eqn{p_0}.
#'
#' @param nn Number of draws to generate (positive integer).
#' @param n Integer vector of trial sizes (\eqn{n \ge 1}), recycled to
#'   length \code{nn}.
#' @param mu Numeric vector of means in \eqn{(\varepsilon, 1 - \varepsilon)},
#'   recycled to length \code{nn}.
#' @param kappa Numeric vector of concentrations (\eqn{> 0}), recycled to
#'   length \code{nn}.
#'
#' @return An integer vector of length \code{nn}, each element in
#'   \eqn{\{1, \ldots, n_i\}}.
#'
#' @references
#' Ghosal, S., Ghosh, S., and Moores, M. (2020).
#' \dQuote{Hierarchical beta-binomial models for batch effects in
#' cytometry data.} \emph{Journal of the Royal Statistical Society:
#' Series A}, \strong{183}(4), 1579--1601.
#'
#' @examples
#' set.seed(42)
#' x <- rztbetabinom(1000, n = 10, mu = 0.3, kappa = 5)
#' min(x)  # always >= 1
#' mean(x)
#'
#' # Vectorised parameters
#' rztbetabinom(6, n = c(5, 10, 20), mu = c(0.2, 0.5, 0.8), kappa = c(3, 10, 50))
#'
#' @family distributions
#' @export
rztbetabinom <- function(nn, n, mu, kappa) {

    assert_count(nn)
    if (nn == 0L) return(integer(0L))

    assert_integerish(n, lower = 1)
    eps <- .Machine$double.eps
    assert_numeric(mu, lower = eps, upper = 1 - eps)
    assert_numeric(kappa, lower = eps)

    n     <- rep_len(as.integer(n), nn)
    mu    <- rep_len(mu, nn)
    kappa <- rep_len(kappa, nn)

    # -- Parameter-grouping strategy for vectorised rejection sampling -------
    # Create group keys from unique (n, mu, kappa) combinations
    group_key <- paste(n, mu, kappa, sep = "|")
    unique_keys <- unique(group_key)

    result <- integer(nn)

    for (key in unique_keys) {
        idx <- which(group_key == key)
        need <- length(idx)
        ni <- n[idx[1]]
        mui <- mu[idx[1]]
        ki <- kappa[idx[1]]

        ai <- mui * ki
        bi <- (1 - mui) * ki

        # Compute p0 for smart batch sizing
        p0 <- exp(lgamma(bi + ni) + lgamma(ki) -
                      lgamma(bi) - lgamma(ki + ni))
        accept_rate <- max(1 - p0, 0.01)

        # Rejection sampling
        collected <- integer(0L)
        max_iter <- 1000L
        iter <- 0L

        while (length(collected) < need && iter < max_iter) {
            iter <- iter + 1L
            remaining <- need - length(collected)
            batch_size <- as.integer(ceiling(remaining / accept_rate)) + 10L

            # Hierarchical draw: p ~ Beta, Y | p ~ Binom
            p_draws <- rbeta(batch_size, shape1 = ai, shape2 = bi)
            y_draws <- rbinom(batch_size, size = ni, prob = p_draws)

            # Keep non-zero draws
            nonzero <- y_draws[y_draws > 0L]
            collected <- c(collected, nonzero)
        }

        if (length(collected) < need) {
            cli_warn(c(
                "Rejection sampling did not produce enough non-zero draws.",
                "i" = "Requested {need}, obtained {length(collected)}.",
                "i" = "Parameters: n={ni}, mu={mui}, kappa={ki}, p0={round(p0, 4)}",
                "i" = "Padding with 1s for {need - length(collected)} missing draw(s)."
            ))
            # Pad with 1 (the smallest valid ZTBB value) to avoid NA
            collected <- c(collected, rep(1L, need - length(collected)))
        }

        result[idx] <- collected[seq_len(need)]
    }

    result
}


# -- 6. compute_p0 -----------------------------------------------------------

#' Compute P(Y = 0) under the Beta-Binomial distribution
#'
#' Computes the probability of observing zero under a
#' \eqn{\text{BetaBin}(n, \mu, \kappa)} distribution using the stable
#' lgamma formula.
#'
#' @details
#' The formula is
#' \deqn{
#'   p_0 = \frac{B(b + n,\; a)}{B(b,\; a)}
#'   = \exp\bigl[\ln\Gamma(b + n) + \ln\Gamma(\kappa)
#'              - \ln\Gamma(b) - \ln\Gamma(\kappa + n)\bigr],
#' }
#' where \eqn{a = \mu\kappa} and \eqn{b = (1 - \mu)\kappa}.
#'
#' **Boundary handling:**
#' \itemize{
#'   \item \eqn{\mu = 0}: Returns 1 (point mass at zero).
#'   \item \eqn{\mu = 1}: Returns 0 (point mass at \eqn{n}).
#'   \item \eqn{n = 0}: Returns 1 (trivially, no trials).
#' }
#'
#' The result is clamped to \eqn{[0, 1]} to guard against floating-point
#' drift.
#'
#' @param n Integer vector of trial sizes (\eqn{n \ge 0}).
#' @param mu Numeric vector of means, each in \eqn{[0, 1]}.
#' @param kappa Numeric vector of concentrations, each \eqn{> 0}.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}.
#'
#' @references
#' Ghosal, S., Ghosh, S., and Moores, M. (2020).
#' \dQuote{Hierarchical beta-binomial models for batch effects in
#' cytometry data.} \emph{Journal of the Royal Statistical Society:
#' Series A}, \strong{183}(4), 1579--1601.
#'
#' @examples
#' # Small n: compare with dbetabinom
#' compute_p0(5, mu = 0.3, kappa = 10)
#' dbetabinom(0, n = 5, mu = 0.3, kappa = 10)
#'
#' # Boundary cases
#' compute_p0(5, mu = 0, kappa = 10)    # 1
#' compute_p0(5, mu = 1, kappa = 10)    # 0
#' compute_p0(0, mu = 0.3, kappa = 10)  # 1
#'
#' @family distributions
#' @export
compute_p0 <- function(n, mu, kappa) {

    assert_integerish(n, lower = 0)
    assert_numeric(mu, lower = 0, upper = 1)
    assert_numeric(kappa, lower = .Machine$double.eps)

    # Recycle
    len <- max(length(n), length(mu), length(kappa))
    n     <- rep_len(as.integer(n), len)
    mu    <- rep_len(mu, len)
    kappa <- rep_len(kappa, len)

    # Allocate
    p0 <- rep(NA_real_, len)

    # Boundary: mu = 0 → point mass at 0
    p0[mu == 0] <- 1

    # Boundary: mu = 1 → point mass at n, so P(Y=0) = 0
    p0[mu == 1] <- 0

    # Boundary: n = 0 → trivially P(Y=0) = 1
    p0[n == 0L] <- 1

    # Main formula for remaining entries
    todo <- is.na(p0)
    if (any(todo)) {
        b_t <- (1 - mu[todo]) * kappa[todo]
        log_p0 <- lgamma(b_t + n[todo]) + lgamma(kappa[todo]) -
            lgamma(b_t) - lgamma(kappa[todo] + n[todo])
        p0[todo] <- exp(log_p0)
    }

    # Clamp to [0, 1]
    pmin(pmax(p0, 0), 1)
}


# -- 7. compute_ztbb_mean ----------------------------------------------------

#' Conditional mean of the zero-truncated Beta-Binomial
#'
#' Returns the conditional mean of the zero-truncated Beta-Binomial,
#' \eqn{E[Y \mid Y > 0] = n \mu / (1 - p_0)}.
#'
#' @details
#' Requires \eqn{n \ge 1} and \eqn{\mu \in (\varepsilon, 1 - \varepsilon)}
#' so that the zero-truncated distribution is well-defined (i.e.,
#' \eqn{p_0 < 1}).
#'
#' @param n Integer vector of trial sizes (\eqn{n \ge 1}).
#' @param mu Numeric vector of means in \eqn{(\varepsilon, 1 - \varepsilon)}.
#' @param kappa Numeric vector of concentrations (\eqn{> 0}).
#'
#' @return A numeric vector of conditional means.
#'
#' @references
#' Ghosal, S., Ghosh, S., and Moores, M. (2020).
#' \dQuote{Hierarchical beta-binomial models for batch effects in
#' cytometry data.} \emph{Journal of the Royal Statistical Society:
#' Series A}, \strong{183}(4), 1579--1601.
#'
#' @examples
#' # Compare analytical mean with empirical mean
#' compute_ztbb_mean(10, mu = 0.3, kappa = 5)
#'
#' set.seed(1)
#' mean(rztbetabinom(10000, n = 10, mu = 0.3, kappa = 5))
#'
#' @family distributions
#' @export
compute_ztbb_mean <- function(n, mu, kappa) {

    assert_integerish(n, lower = 1)
    eps <- .Machine$double.eps
    assert_numeric(mu, lower = eps, upper = 1 - eps)
    assert_numeric(kappa, lower = eps)

    # Recycle
    len <- max(length(n), length(mu), length(kappa))
    n     <- rep_len(as.integer(n), len)
    mu    <- rep_len(mu, len)
    kappa <- rep_len(kappa, len)

    p0 <- compute_p0(n, mu, kappa)
    n * mu / (1 - p0)
}
