# ============================================================================
# hurdle.R — Hurdle Beta-Binomial distribution functions
#
# Two-part model: P(Y = 0) = 1 - q, P(Y = y) = q * f_ZT(y) for y >= 1.
#
# Contents:
#   1. dhurdle_betabinom  — PMF of the hurdle Beta-Binomial
#   2. rhurdle_betabinom  — Random generation from the hurdle Beta-Binomial
#   3. hurdle_mean        — E[Y] under the hurdle model
#   4. hurdle_variance    — Var[Y] under the hurdle model
# ============================================================================


# -- 1. dhurdle_betabinom -----------------------------------------------------

#' Probability mass function of the Hurdle Beta-Binomial
#'
#' Computes the PMF of the two-part hurdle model:
#' \deqn{
#'   P(Y = 0) = 1 - q, \quad
#'   P(Y = y) = q \cdot f_{\text{ZT}}(y \mid n, \mu, \kappa),
#'   \quad y \in \{1, \ldots, n\},
#' }
#' where \eqn{f_{\text{ZT}}} is the zero-truncated Beta-Binomial PMF.
#'
#' @details
#' The structural zero probability \eqn{1 - q} is computed as
#' \code{log1p(-q)} in the log domain for numerical stability when
#' \eqn{q} is near 1.
#'
#' @param y Integer vector of observed counts.
#' @param n Integer vector of trial sizes (\eqn{n \ge 1}).
#' @param q Numeric vector of participation probabilities, each in
#'   \eqn{[0, 1]}.
#' @param mu Numeric vector of intensity means, each in
#'   \eqn{(\varepsilon, 1 - \varepsilon)}.
#' @param kappa Numeric vector of concentrations (\eqn{> 0}).
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
#' # PMF over full support
#' dhurdle_betabinom(0:5, n = 5, q = 0.7, mu = 0.4, kappa = 8)
#'
#' # Should sum to 1
#' sum(dhurdle_betabinom(0:5, n = 5, q = 0.7, mu = 0.4, kappa = 8))
#'
#' # Log scale
#' dhurdle_betabinom(0:5, n = 5, q = 0.7, mu = 0.4, kappa = 8, log = TRUE)
#'
#' @family distributions
#' @export
dhurdle_betabinom <- function(y, n, q, mu, kappa, log = FALSE) {

    assert_flag(log)
    assert_integerish(y)
    assert_integerish(n, lower = 1)
    assert_numeric(q, lower = 0, upper = 1)
    eps <- .Machine$double.eps
    assert_numeric(mu, lower = eps, upper = 1 - eps)
    assert_numeric(kappa, lower = eps)

    # Recycle to common length
    len <- max(length(y), length(n), length(q), length(mu), length(kappa))
    y     <- rep_len(as.numeric(y), len)
    n     <- rep_len(as.numeric(n), len)
    q     <- rep_len(q, len)
    mu    <- rep_len(mu, len)
    kappa <- rep_len(kappa, len)

    # Allocate
    log_f <- rep(NA_real_, len)

    # Mark out-of-support
    invalid <- (y < 0) | (y > n) | (y != round(y))

    # --- y = 0: P(Y = 0) = 1 - q ---
    is_zero <- (!invalid) & (y == 0)
    if (any(is_zero)) {
        log_f[is_zero] <- log1p(-q[is_zero])
    }

    # --- y >= 1: P(Y = y) = q * f_ZT(y) ---
    is_pos <- (!invalid) & (y >= 1)
    if (any(is_pos)) {
        log_zt <- dztbetabinom(y[is_pos], n[is_pos], mu[is_pos],
                               kappa[is_pos], log = TRUE)
        log_f[is_pos] <- log(q[is_pos]) + log_zt
    }

    # Invalid entries
    log_f[invalid] <- -Inf

    if (log) log_f else exp(log_f)
}


# -- 2. rhurdle_betabinom -----------------------------------------------------

#' Random generation from the Hurdle Beta-Binomial
#'
#' Generates random draws from the two-part hurdle model.
#' First draws \eqn{Z \sim \text{Bernoulli}(q)}, then for
#' \eqn{Z = 1} draws \eqn{Y \sim \text{ZT-BetaBin}(n, \mu, \kappa)}.
#'
#' @details
#' All parameter vectors (\code{n}, \code{q}, \code{mu}, \code{kappa})
#' are recycled to length \code{nn}. The function is fully vectorised.
#'
#' @param nn Number of draws to generate (positive integer).
#' @param n Integer vector of trial sizes (\eqn{n \ge 1}).
#' @param q Numeric vector of participation probabilities, each in
#'   \eqn{[0, 1]}.
#' @param mu Numeric vector of intensity means in
#'   \eqn{(\varepsilon, 1 - \varepsilon)}.
#' @param kappa Numeric vector of concentrations (\eqn{> 0}).
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
#' x <- rhurdle_betabinom(1000, n = 10, q = 0.7, mu = 0.3, kappa = 5)
#' mean(x == 0)  # approximately 0.3
#' hist(x, breaks = -0.5:10.5)
#'
#' @family distributions
#' @export
rhurdle_betabinom <- function(nn, n, q, mu, kappa) {

    assert_count(nn)
    if (nn == 0L) return(integer(0L))

    assert_integerish(n, lower = 1)
    assert_numeric(q, lower = 0, upper = 1)
    eps <- .Machine$double.eps
    assert_numeric(mu, lower = eps, upper = 1 - eps)
    assert_numeric(kappa, lower = eps)

    n     <- rep_len(as.integer(n), nn)
    q     <- rep_len(q, nn)
    mu    <- rep_len(mu, nn)
    kappa <- rep_len(kappa, nn)

    # Part 1: z ~ Bernoulli(q)
    z <- rbinom(nn, size = 1L, prob = q)

    # Allocate result: zeros by default
    result <- integer(nn)

    # Part 2: for z = 1, draw from ZT-BetaBin
    active <- which(z == 1L)
    if (length(active) > 0L) {
        result[active] <- rztbetabinom(
            nn    = length(active),
            n     = n[active],
            mu    = mu[active],
            kappa = kappa[active]
        )
    }

    result
}


# -- 3. hurdle_mean -----------------------------------------------------------

#' Expected value under the Hurdle Beta-Binomial
#'
#' Computes the unconditional mean
#' \deqn{
#'   E[Y] = q \cdot \frac{n \mu}{1 - p_0},
#' }
#' where \eqn{p_0 = P(Y = 0)} under the (untruncated) Beta-Binomial.
#'
#' @param n Integer vector of trial sizes (\eqn{n \ge 1}).
#' @param q Numeric vector of participation probabilities, each in
#'   \eqn{[0, 1]}.
#' @param mu Numeric vector of intensity means in
#'   \eqn{(\varepsilon, 1 - \varepsilon)}.
#' @param kappa Numeric vector of concentrations (\eqn{> 0}).
#'
#' @return A numeric vector of expected values.
#'
#' @references
#' Ghosal, S., Ghosh, S., and Moores, M. (2020).
#' \dQuote{Hierarchical beta-binomial models for batch effects in
#' cytometry data.} \emph{Journal of the Royal Statistical Society:
#' Series A}, \strong{183}(4), 1579--1601.
#'
#' @examples
#' hurdle_mean(n = 10, q = 0.7, mu = 0.3, kappa = 5)
#'
#' # q = 0 always gives 0
#' hurdle_mean(n = 10, q = 0, mu = 0.3, kappa = 5)
#'
#' @family distributions
#' @export
hurdle_mean <- function(n, q, mu, kappa) {

    assert_integerish(n, lower = 1)
    assert_numeric(q, lower = 0, upper = 1)
    eps <- .Machine$double.eps
    assert_numeric(mu, lower = eps, upper = 1 - eps)
    assert_numeric(kappa, lower = eps)

    # Recycle
    len <- max(length(n), length(q), length(mu), length(kappa))
    n     <- rep_len(as.integer(n), len)
    q     <- rep_len(q, len)
    mu    <- rep_len(mu, len)
    kappa <- rep_len(kappa, len)

    # Handle q = 0 explicitly
    result <- rep(NA_real_, len)
    result[q == 0] <- 0

    # For q > 0, compute E[Y] = q * n * mu / (1 - p0)
    active <- (q > 0)
    if (any(active)) {
        ztbb_mean <- compute_ztbb_mean(n[active], mu[active], kappa[active])
        result[active] <- q[active] * ztbb_mean
    }

    result
}


# -- 4. hurdle_variance -------------------------------------------------------

#' Variance under the Hurdle Beta-Binomial
#'
#' Computes the unconditional variance via the law of total variance:
#' \deqn{
#'   \text{Var}[Y] = q \cdot V_{\text{ZT}}
#'                  + q (1 - q) \cdot (E_{\text{ZT}})^2,
#' }
#' where \eqn{E_{\text{ZT}} = n\mu / (1 - p_0)} is the zero-truncated
#' mean and \eqn{V_{\text{ZT}}} is the zero-truncated variance.
#'
#' @details
#' The zero-truncated variance is derived from the untruncated moments:
#' \deqn{
#'   V_{\text{ZT}} = \frac{E_{\text{BB}}[Y^2]}{1 - p_0}
#'                 - \left(\frac{E_{\text{BB}}[Y]}{1 - p_0}\right)^2,
#' }
#' where
#' \deqn{
#'   E_{\text{BB}}[Y^2] = \text{Var}_{\text{BB}} + (n\mu)^2
#'   = n\mu(1 - \mu) \frac{n + \kappa}{1 + \kappa} + n^2 \mu^2.
#' }
#'
#' The result is clamped to be non-negative via \code{pmax(., 0)}.
#'
#' @param n Integer vector of trial sizes (\eqn{n \ge 1}).
#' @param q Numeric vector of participation probabilities, each in
#'   \eqn{[0, 1]}.
#' @param mu Numeric vector of intensity means in
#'   \eqn{(\varepsilon, 1 - \varepsilon)}.
#' @param kappa Numeric vector of concentrations (\eqn{> 0}).
#'
#' @return A numeric vector of variances (non-negative).
#'
#' @references
#' Ghosal, S., Ghosh, S., and Moores, M. (2020).
#' \dQuote{Hierarchical beta-binomial models for batch effects in
#' cytometry data.} \emph{Journal of the Royal Statistical Society:
#' Series A}, \strong{183}(4), 1579--1601.
#'
#' @examples
#' hurdle_variance(n = 10, q = 0.7, mu = 0.3, kappa = 5)
#'
#' # q = 0 always gives 0
#' hurdle_variance(n = 10, q = 0, mu = 0.3, kappa = 5)
#'
#' @family distributions
#' @export
hurdle_variance <- function(n, q, mu, kappa) {

    assert_integerish(n, lower = 1)
    assert_numeric(q, lower = 0, upper = 1)
    eps <- .Machine$double.eps
    assert_numeric(mu, lower = eps, upper = 1 - eps)
    assert_numeric(kappa, lower = eps)

    # Recycle
    len <- max(length(n), length(q), length(mu), length(kappa))
    n     <- rep_len(as.integer(n), len)
    q     <- rep_len(q, len)
    mu    <- rep_len(mu, len)
    kappa <- rep_len(kappa, len)

    # Handle q = 0 explicitly
    result <- rep(NA_real_, len)
    result[q == 0] <- 0

    active <- (q > 0)
    if (any(active)) {
        ni <- n[active]
        qi <- q[active]
        mui <- mu[active]
        ki <- kappa[active]

        p0 <- compute_p0(ni, mui, ki)
        one_mp0 <- 1 - p0

        # BB moments
        E_bb  <- ni * mui                              # E_BB[Y]
        V_bb  <- ni * mui * (1 - mui) * (ni + ki) / (1 + ki)  # Var_BB
        E2_bb <- V_bb + E_bb^2                        # E_BB[Y^2]

        # ZT moments
        E_zt  <- E_bb / one_mp0                       # E[Y | Y > 0]
        V_zt  <- E2_bb / one_mp0 - E_zt^2             # Var[Y | Y > 0]

        # Law of total variance
        result[active] <- qi * V_zt + qi * (1 - qi) * E_zt^2
    }

    pmax(result, 0)
}
