# ============================================================================
# copula-utils.R --- Gaussian Copula Infrastructure for Synthetic Data
#
# Generates synthetic covariate data that preserves the marginal
# distributions and rank-correlation structure of empirical data via
# the Gaussian copula (Nelsen, 2006; Miratrix, 2025).
#
# Pipeline:
#   1. Convert observed data to z-scores: rank -> uniform -> qnorm
#   2. Estimate the correlation matrix in the z-space
#   3. Generate multivariate normal draws, transform to uniform,
#      apply variable-specific inverse ECDFs
#
# Contents:
#   1. .convert_to_z          --- Rank-based normal-score transformation
#   2. .make_inv_ecdf         --- Inverse ECDF factory (closure)
#   3. .copula_generate       --- MVN + inverse ECDF synthesis
#   4. .extract_copula_params --- One-call extraction from a data frame
#   5. .ensure_pd             --- PD repair helper
#
# All functions are internal (@noRd, not exported).
# ============================================================================


# -- 1. .convert_to_z ---------------------------------------------------------

#' Convert a numeric vector to normal scores via rank transform
#'
#' Implements the first stage of the Gaussian copula: observed values are
#' mapped to their empirical CDF positions (via rank), then transformed to
#' standard normal quantiles.
#'
#' @param y Numeric vector of observed values. May contain \code{NA}s.
#'
#' @return Numeric vector of z-scores (standard normal quantiles), same
#'   length as \code{y}. \code{NA} values are propagated.
#'
#' @details
#' Steps:
#' \enumerate{
#'   \item Remove NAs, work on observed values only.
#'   \item Compute ranks with jittered tie-breaking (see below).
#'   \item Convert to pseudo-uniform: \eqn{u_i = r_i / (N + 1)}, where
#'     \eqn{N} is the number of non-NA observations. The \eqn{N + 1}
#'     denominator (Hazen plotting position) keeps \eqn{u} strictly in
#'     \eqn{(0, 1)}, avoiding \eqn{\pm\infty} from \code{qnorm}.
#'   \item Convert to normal: \eqn{z_i = \Phi^{-1}(u_i)}.
#' }
#'
#' **Tie handling.** Ties are broken by adding tiny uniform noise to the
#' data before ranking. The noise magnitude is \eqn{10^{-10}} times the
#' data range, so it does not alter the empirical distribution shape but
#' does make the z-scores depend on the current RNG state. This is
#' acceptable because the copula itself is a stochastic generator.
#'
#' @section Edge cases:
#' \describe{
#'   \item{Empty vector (length 0)}{Returns \code{numeric(0)}.}
#'   \item{All-NA vector}{Returns all-NA with a cli warning.}
#'   \item{Single non-NA value}{Returns 0 with a cli warning, since
#'     one observation has no ranking context.}
#'   \item{Constant vector (zero variance)}{Returns \code{rep(0, N)}
#'     with a cli warning. A zero-variance variable contributes no
#'     rank information to the copula.}
#' }
#'
#' @noRd
.convert_to_z <- function(y) {

    assert_numeric(y)
    n_total <- length(y)
    if (n_total == 0L) return(numeric(0L))

    # --- Locate non-missing observations ---
    is_obs <- !is.na(y)
    n_obs <- sum(is_obs)

    # --- Edge case: all NA ---
    if (n_obs == 0L) {
        cli_warn(c(
            "All values are {.code NA} in {.fn .convert_to_z}.",
            "i" = "Returning a vector of {n_total} {.code NA}s."
        ))
        return(rep(NA_real_, n_total))
    }

    # --- Edge case: single non-NA observation ---
    if (n_obs == 1L) {
        cli_warn(c(
            "Only 1 non-{.code NA} value in {.fn .convert_to_z}.",
            "i" = "Returning 0 for the single observation (median rank)."
        ))
        z <- rep(NA_real_, n_total)
        z[is_obs] <- 0
        return(z)
    }

    y_obs <- y[is_obs]

    # --- Edge case: constant vector (zero variance) ---
    y_range <- diff(range(y_obs))
    if (y_range < .Machine$double.eps * max(1, abs(max(y_obs)))) {
        cli_warn(c(
            "All non-{.code NA} values are identical in {.fn .convert_to_z}.",
            "i" = "Returning 0 for all {n_obs} observations (no rank info)."
        ))
        z <- rep(NA_real_, n_total)
        z[is_obs] <- 0
        return(z)
    }

    # --- Main path: jittered rank -> uniform -> qnorm ---

    # Break ties with tiny noise. The magnitude is negligible relative to
    # the data range so marginal shapes are preserved. We use uniform noise
    # rather than normal to bound the maximum perturbation.
    jitter_mag <- y_range * 1e-10
    y_jittered <- y_obs + runif(n_obs, min = -jitter_mag, max = jitter_mag)

    # Rank the jittered values (ties virtually eliminated).
    # Use "average" as a safety net in the astronomically unlikely event
    # that two jittered values are still identical.
    ranks <- rank(y_jittered, ties.method = "average")

    # Pseudo-uniform: u = rank / (N + 1) keeps u in (0, 1) strictly.
    # This avoids qnorm(0) = -Inf and qnorm(1) = +Inf.
    u <- ranks / (n_obs + 1L)

    # Standard normal z-scores
    z_obs <- qnorm(u)

    # --- Assemble output, preserving NA positions ---
    z <- rep(NA_real_, n_total)
    z[is_obs] <- z_obs
    z
}


# -- 2. .make_inv_ecdf --------------------------------------------------------

#' Create an inverse empirical CDF function (closure)
#'
#' Builds a self-contained closure that maps uniform \eqn{[0, 1]} quantiles
#' back to the original data scale. Used as the final stage of the Gaussian
#' copula to recover original-scale marginals.
#'
#' The returned function carries its own lookup data (captured in the
#' enclosing environment) and has no external dependencies.
#'
#' @param vals Numeric vector of observed values. \code{NA}s are silently
#'   removed.
#' @param method Character; interpolation method. One of \code{"linear"}
#'   (for continuous variables) or \code{"constant"} (for discrete/integer
#'   variables, including binary). Default is \code{"linear"}.
#' @param n_grid Integer; target number of quantile knots for downsampling.
#'   Used only when the number of unique values exceeds \code{n_grid}.
#'   The min, max, and approximately \code{n_grid - 2} interior quantiles
#'   are kept. Default is 300.
#'
#' @return A function \code{f(u)} that accepts a numeric vector of
#'   probabilities in \eqn{[0, 1]} and returns values on the original
#'   data scale.
#'
#' @section Edge cases:
#' \describe{
#'   \item{No valid values}{Errors via cli.}
#'   \item{Single unique value}{Returns a constant function.}
#'   \item{Binary variable (0/1) with \code{"constant"}}{Returns exactly
#'     0 or 1, with the cutpoint at the empirical proportion of zeros.}
#'   \item{Heavily tied integer data}{Use \code{"constant"} to ensure
#'     output values are always exact observed values.}
#'   \item{Very large input (> n_grid unique values)}{Downsampled to
#'     ~n_grid knots; extreme values always preserved.}
#' }
#'
#' @noRd
.make_inv_ecdf <- function(vals,
                           method = c("linear", "constant"),
                           n_grid = 300L) {

    method <- match.arg(method)
    assert_numeric(vals)
    assert_count(n_grid, positive = TRUE)

    # Remove NAs silently
    vals <- vals[!is.na(vals)]
    n <- length(vals)

    if (n == 0L) {
        cli_abort(c(
            "Cannot build inverse ECDF: no valid (non-{.code NA}) values.",
            "x" = "Input was empty or all {.code NA}."
        ))
    }

    # --- Single unique value: constant function ---
    unique_vals <- sort(unique(vals))
    n_unique <- length(unique_vals)

    if (n_unique == 1L) {
        const_val <- unique_vals[1L]
        return(function(u) {
            assert_numeric(u, lower = 0, upper = 1)
            rep(const_val, length(u))
        })
    }

    # --- Compute ECDF positions for each unique value ---
    sorted_all <- sort(vals)

    ecdf_at <- vapply(unique_vals, function(v) {
        sum(sorted_all <= v) / n
    }, numeric(1L))

    # --- Downsample if too many unique values ---
    if (n_unique > n_grid) {
        # Keep first (min), last (max), and ~(n_grid - 2) interior points
        # chosen at evenly-spaced indices. This preserves the tail extremes
        # exactly while compressing the interior.
        interior_idx <- unique(round(
            seq(2, n_unique - 1L, length.out = max(1L, n_grid - 2L))
        ))
        idx <- unique(c(1L, interior_idx, n_unique))
        ecdf_at     <- ecdf_at[idx]
        unique_vals <- unique_vals[idx]
    }

    # --- Build the inverse CDF function ---
    if (method == "linear") {
        inv_fn <- stats::approxfun(
            x    = ecdf_at,
            y    = unique_vals,
            rule = 2,              # constant extrapolation beyond boundaries
            ties = "ordered"       # ecdf_at is already sorted
        )
    } else {
        # Step function for discrete data. f = 1 gives right-continuous steps
        # (i.e., for u between F(x_i) and F(x_{i+1}), return x_{i+1}),
        # matching the standard quantile function Q(u) = inf{x : F(x) >= u}.
        inv_fn <- stats::approxfun(
            x      = ecdf_at,
            y      = unique_vals,
            method = "constant",
            f      = 1,
            rule   = 2,
            ties   = "ordered"
        )
    }

    # Wrap in a closure that validates input and makes the function
    # self-contained (inv_fn is captured in this environment).
    function(u) {
        assert_numeric(u, lower = 0, upper = 1)
        inv_fn(u)
    }
}


# -- 3. .copula_generate -------------------------------------------------------

#' Generate synthetic data via the Gaussian copula
#'
#' Draws \code{n} observations from a multivariate distribution that
#' preserves the pairwise rank-correlation structure (in z-space) and the
#' marginal distributions (via inverse ECDFs) of the original data.
#'
#' @param n Positive integer; number of observations to generate.
#' @param cor_matrix Numeric \eqn{p \times p} correlation matrix
#'   estimated in z-space. Must be symmetric positive (semi-)definite.
#'   If not PD, repaired automatically via \code{.ensure_pd()}.
#' @param inv_ecdfs Named list of \eqn{p} functions, each mapping
#'   \eqn{[0, 1] \to} original scale.
#' @param seed Optional integer seed for reproducibility. If non-NULL,
#'   the RNG state is saved before, set to \code{seed}, and restored
#'   on exit so the caller's state is unaffected.
#'
#' @return A data frame with \code{n} rows and \code{p} columns, named
#'   according to \code{names(inv_ecdfs)}.
#'
#' @section Edge cases:
#' \describe{
#'   \item{Non-PD cor_matrix}{Repaired via \code{Matrix::nearPD()} with
#'     a cli warning.}
#'   \item{p = 1 (single variable)}{Uses rnorm + pnorm + inv_ecdf
#'     directly, avoiding \code{MASS::mvrnorm} overhead.}
#'   \item{n = 1}{MASS::mvrnorm returns a vector; forced to 1-row
#'     matrix.}
#'   \item{Seed handling}{RNG state saved and restored via
#'     \code{on.exit()}.}
#' }
#'
#' @noRd
.copula_generate <- function(n, cor_matrix, inv_ecdfs, seed = NULL) {

    # --- Input validation ---
    assert_count(n, positive = TRUE)
    assert_matrix(cor_matrix, mode = "numeric")
    p <- ncol(cor_matrix)

    if (nrow(cor_matrix) != p) {
        cli_abort(c(
            "Correlation matrix must be square.",
            "x" = "Got {nrow(cor_matrix)} x {p}."
        ))
    }

    checkmate::assert_list(inv_ecdfs, types = "function", min.len = 1L)

    if (length(inv_ecdfs) != p) {
        cli_abort(c(
            "Dimension mismatch between {.arg cor_matrix} and
             {.arg inv_ecdfs}.",
            "x" = paste0(
                "{.arg cor_matrix} is {p} x {p} but {.arg inv_ecdfs} ",
                "has {length(inv_ecdfs)} element{?s}."
            )
        ))
    }

    if (is.null(names(inv_ecdfs)) || any(names(inv_ecdfs) == "")) {
        cli_abort(c(
            "All elements of {.arg inv_ecdfs} must be named.",
            "i" = "Names define the column names of the output data frame."
        ))
    }

    # --- Seed handling: save, set, restore on exit ---
    if (!is.null(seed)) {
        assert_integerish(seed, len = 1L, any.missing = FALSE)
        if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
            old_seed <- get(".Random.seed", envir = globalenv())
            on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
                    add = TRUE)
        } else {
            # No seed existed before; clean up when done
            on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
        }
        set.seed(as.integer(seed))
    }

    # --- Ensure cor_matrix is PD ---
    cor_matrix <- .ensure_pd(cor_matrix)

    # --- Generate MVN(0, Sigma) draws ---
    if (p == 1L) {
        # Univariate: skip MASS::mvrnorm for efficiency
        z_matrix <- matrix(rnorm(n), ncol = 1L)
    } else {
        z_matrix <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = cor_matrix)
        # When n = 1, MASS::mvrnorm returns a vector; coerce to matrix
        if (!is.matrix(z_matrix)) {
            z_matrix <- matrix(z_matrix, nrow = 1L)
        }
    }

    # --- Convert to uniform via Phi(z) ---
    u_matrix <- pnorm(z_matrix)

    # --- Apply inverse ECDFs column by column ---
    var_names <- names(inv_ecdfs)
    x_list <- vector("list", p)
    names(x_list) <- var_names

    for (j in seq_len(p)) {
        x_list[[j]] <- inv_ecdfs[[j]](u_matrix[, j])
    }

    as.data.frame(x_list, stringsAsFactors = FALSE)
}


# -- 4. .extract_copula_params ------------------------------------------------

#' Extract Gaussian copula parameters from observed data
#'
#' One-call convenience function that takes a data frame and returns
#' all objects needed by \code{.copula_generate}: the z-space correlation
#' matrix and a named list of inverse ECDF functions.
#'
#' @param data Data frame containing the original data.
#' @param var_names Character vector of column names to extract.
#' @param types Character vector specifying variable types. Can be either
#'   (a) a named vector with names matching \code{var_names}, or
#'   (b) an unnamed vector of the same length as \code{var_names} (auto-named).
#'   Values must be one of \code{"continuous"}, \code{"discrete"},
#'   \code{"binary"}, or \code{"count"}. This determines the interpolation
#'   method: \code{"linear"} for continuous, \code{"constant"} for all others.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{cor_matrix}}{Numeric \eqn{p \times p} correlation
#'     matrix estimated from z-scores (complete cases), possibly
#'     PD-repaired.}
#'   \item{\code{inv_ecdfs}}{Named list of \eqn{p} inverse ECDF
#'     closures.}
#'   \item{\code{n_complete}}{Integer; number of complete cases used
#'     for estimation.}
#'   \item{\code{n_dropped}}{Integer; number of rows dropped due to
#'     missing values.}
#'   \item{\code{var_names}}{Character vector of variable names.}
#'   \item{\code{types}}{The types vector (subset to var_names).}
#' }
#'
#' @section Edge cases:
#' \describe{
#'   \item{Missing data}{Complete-case analysis. Warns if > 50\% of
#'     rows are dropped.}
#'   \item{Zero-variance variable}{Warned by \code{.convert_to_z()};
#'     z-scores are all-zero, which may produce singular correlation.}
#'   \item{Resulting non-PD matrix}{Repaired via \code{.ensure_pd()}.}
#'   \item{Types must be named}{Errors clearly if types vector lacks
#'     names or is missing entries.}
#' }
#'
#' @noRd
.extract_copula_params <- function(data, var_names, types) {

    # --- Validation ---
    assert_data_frame(data, min.rows = 2L)
    checkmate::assert_character(var_names, min.len = 1L, any.missing = FALSE,
                                unique = TRUE)
    checkmate::assert_character(types, min.len = 1L, any.missing = FALSE)

    # Accept either named or unnamed types.
    # Unnamed types are auto-named from var_names if same length.
    if (is.null(names(types))) {
        if (length(types) != length(var_names)) {
            cli_abort(c(
                "Unnamed {.arg types} must have the same length as
                 {.arg var_names}.",
                "x" = paste0(
                    "{.arg types} has {length(types)} element{?s} but ",
                    "{.arg var_names} has {length(var_names)}."
                ),
                "i" = "Alternatively, supply a named vector: {.code c(x1 = 'continuous', x2 = 'discrete')}."
            ))
        }
        names(types) <- var_names
    }

    missing_type <- setdiff(var_names, names(types))
    if (length(missing_type) > 0L) {
        cli_abort(c(
            "Missing type specification for variable{?s}:
             {.val {missing_type}}.",
            "i" = paste0(
                "All variables in {.arg var_names} must have a ",
                "corresponding entry in {.arg types}."
            )
        ))
    }

    valid_types <- c("continuous", "binary", "count", "discrete")
    bad_types <- types[var_names][!types[var_names] %in% valid_types]
    if (length(bad_types) > 0L) {
        cli_abort(c(
            "Invalid type{?s}: {.val {unique(bad_types)}}.",
            "i" = "Allowed types: {.val {valid_types}}."
        ))
    }

    # Check columns exist
    missing_cols <- setdiff(var_names, names(data))
    if (length(missing_cols) > 0L) {
        cli_abort(c(
            "Column{?s} not found in {.arg data}: {.val {missing_cols}}."
        ))
    }

    # --- Extract and restrict to complete cases ---
    sub <- data[, var_names, drop = FALSE]
    complete_mask <- stats::complete.cases(sub)
    n_total <- nrow(sub)
    n_complete <- sum(complete_mask)
    n_dropped <- n_total - n_complete

    if (n_complete < 2L) {
        cli_abort(c(
            "Fewer than 2 complete cases after removing {.code NA} rows.",
            "x" = paste0(
                "{n_total} total rows, {n_dropped} dropped, ",
                "{n_complete} remaining."
            ),
            "i" = "Cannot estimate a correlation matrix."
        ))
    }

    if (n_dropped > 0L) {
        pct_dropped <- round(100 * n_dropped / n_total, 1)
        if (pct_dropped > 50) {
            cli_warn(c(
                "{n_dropped} row{?s} ({pct_dropped}%) dropped due to
                 missing values.",
                "!" = "More than half of rows lost. Check data quality."
            ))
        } else {
            cli_inform(c(
                "Dropped {n_dropped} row{?s} with {.code NA} values
                 ({n_complete} complete cases remain)."
            ))
        }
    }

    sub_complete <- sub[complete_mask, , drop = FALSE]

    # --- Check for zero-variance variables (pre-warning) ---
    for (v in var_names) {
        vals <- sub_complete[[v]]
        val_range <- diff(range(vals))
        if (val_range < .Machine$double.eps * max(1, abs(max(vals)))) {
            cli_warn(c(
                "Variable {.val {v}} has zero variance after {.code NA}
                 removal.",
                "i" = paste0(
                    "All {n_complete} values equal {.val {vals[1]}}. ",
                    "z-scores will be all-zero; correlation with this ",
                    "variable is undefined."
                )
            ))
        }
    }

    # --- Convert each column to z-scores ---
    p <- length(var_names)
    z_matrix <- matrix(NA_real_, nrow = n_complete, ncol = p)
    colnames(z_matrix) <- var_names

    for (j in seq_len(p)) {
        z_matrix[, j] <- .convert_to_z(sub_complete[[var_names[j]]])
    }

    # --- Estimate correlation matrix in z-space ---
    cor_matrix <- stats::cor(z_matrix)

    # Enforce exact symmetry (floating-point safety)
    cor_matrix <- (cor_matrix + t(cor_matrix)) / 2

    # Force diagonal to exactly 1
    diag(cor_matrix) <- 1

    rownames(cor_matrix) <- var_names
    colnames(cor_matrix) <- var_names

    # Repair non-PD if necessary
    cor_matrix <- .ensure_pd(cor_matrix)

    # --- Build inverse ECDFs ---
    inv_ecdfs <- vector("list", p)
    names(inv_ecdfs) <- var_names

    for (j in seq_len(p)) {
        vtype <- types[var_names[j]]
        ecdf_method <- if (vtype == "continuous") "linear" else "constant"
        inv_ecdfs[[j]] <- .make_inv_ecdf(
            vals   = sub_complete[[var_names[j]]],
            method = ecdf_method
        )
    }

    # --- Return ---
    list(
        cor_matrix = cor_matrix,
        inv_ecdfs  = inv_ecdfs,
        n_complete = n_complete,
        n_dropped  = n_dropped,
        var_names  = var_names,
        types      = types[var_names]
    )
}


# -- 5. .ensure_pd (internal helper) ------------------------------------------

#' Ensure a matrix is positive definite
#'
#' If all eigenvalues of \code{mat} are positive, returns it unchanged.
#' Otherwise, repairs via \code{Matrix::nearPD()} (Higham's algorithm),
#' forcing unit diagonal and emitting a cli warning.
#'
#' @param mat Symmetric numeric matrix.
#' @return A positive-definite numeric matrix.
#' @noRd
.ensure_pd <- function(mat) {

    eigen_vals <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
    min_eval <- min(eigen_vals)

    # Tolerance: eigenvalues within floating-point noise of zero are okay
    # for nearPD, but truly negative ones need correction.
    if (min_eval > -sqrt(.Machine$double.eps)) {
        return(mat)
    }

    cli_warn(c(
        "Correlation matrix is not positive definite; repairing via
         {.fn Matrix::nearPD}.",
        "i" = paste0(
            "Smallest eigenvalue: {.val {signif(min_eval, 4)}}. ",
            "This can occur with listwise deletion or heavy ties."
        )
    ))

    pd_result <- Matrix::nearPD(mat, corr = TRUE, keepDiag = TRUE,
                                 maxit = 200L)
    mat_fixed <- as.matrix(pd_result$mat)

    # Enforce exact symmetry and unit diagonal
    mat_fixed <- (mat_fixed + t(mat_fixed)) / 2
    diag(mat_fixed) <- 1

    mat_fixed
}
