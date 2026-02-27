# ============================================================================
# utils.R — Internal utility functions for hurdlebb
#
# Contents:
#   1. inv_logit         — Inverse logit (sigmoid)
#   2. logit             — Logit (log-odds)
#   3. normalize_weights — Rescale survey weights
#   4. log1mexp          — log(1 - exp(-x)), Machler (2012)
#
# All functions are internal (@noRd, not exported).
# ============================================================================


# -- 1. inv_logit -------------------------------------------------------------

#' Inverse logit (sigmoid function)
#'
#' Computes \eqn{1 / (1 + \exp(-x))} via \code{stats::plogis}.
#'
#' @param x Numeric vector.
#' @return Numeric vector in \eqn{(0, 1)}.
#' @noRd
inv_logit <- function(x) {
    plogis(x)
}


# -- 2. logit -----------------------------------------------------------------

#' Logit (log-odds) function
#'
#' Computes \eqn{\log(p / (1 - p))} via \code{stats::qlogis}.
#'
#' @param p Numeric vector in \eqn{(0, 1)}.
#' @return Numeric vector.
#' @noRd
logit <- function(p) {
    qlogis(p)
}


# -- 3. normalize_weights -----------------------------------------------------

#' Normalize survey weights
#'
#' Rescales a positive weight vector so that
#' \eqn{\sum \tilde{w}_i = N} where \eqn{N = \text{length}(w)}.
#'
#' @param w Numeric vector of positive weights.
#' @return Numeric vector of rescaled weights with \code{sum(w) == length(w)}.
#' @noRd
normalize_weights <- function(w) {
    checkmate::assert_numeric(w, lower = .Machine$double.eps,
                              any.missing = FALSE, min.len = 1L)
    w * length(w) / sum(w)
}


# -- 4. log1mexp --------------------------------------------------------------

#' Compute log(1 - exp(-x)) for non-negative x
#'
#' Uses the two-branch algorithm from Machler (2012) for numerical
#' stability:
#' \itemize{
#'   \item \eqn{x \le \log 2}: \code{log(-expm1(-x))}
#'   \item \eqn{x > \log 2}: \code{log1p(-exp(-x))}
#' }
#'
#' @details
#' When \eqn{x = 0}, \eqn{1 - \exp(0) = 0}, so \code{log(0) = -Inf}.
#' \code{NA} values are propagated.
#'
#' @param x Numeric vector of non-negative values.
#' @return Numeric vector of \code{log(1 - exp(-x))}.
#'
#' @references
#' Machler, M. (2012). \dQuote{Accurately computing log(1 - exp(-|a|)):
#' assessed by the Rmpfr package.} Technical report.
#'
#' @noRd
log1mexp <- function(x) {
    out <- rep(NA_real_, length(x))

    # Propagate NA
    not_na <- !is.na(x)

    # x = 0 → log(1 - 1) = log(0) = -Inf
    at_zero <- not_na & (x == 0)
    out[at_zero] <- -Inf

    # Two-branch algorithm for positive x
    small <- not_na & !at_zero & (x <= log(2))
    large <- not_na & !at_zero & (x > log(2))

    if (any(small)) {
        out[small] <- log(-expm1(-x[small]))
    }
    if (any(large)) {
        out[large] <- log1p(-exp(-x[large]))
    }

    out
}
