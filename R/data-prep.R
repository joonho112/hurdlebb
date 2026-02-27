# ============================================================================
# data-prep.R --- Stan Data Preparation for Hurdle Beta-Binomial Models
#
# Transforms an hbb_formula object and data frame into a validated
# Stan-ready data list (S3 class "hbb_data"). Handles design matrix
# construction, standardisation, group indexing, survey weight
# normalisation, and dimensional bookkeeping.
#
# Contents:
#   1. prepare_stan_data    --- Build Stan data list (exported)
#   2. validate_hbb_data    --- Validate a Stan data list (exported)
#   3. print.hbb_data       --- Compact summary print method
#   4. .build_X_matrix      --- Design matrix with optional standardisation
#   5. .build_V_matrix      --- Policy design matrix with optional standardisation
#   6. .build_group_index   --- Sequential integer group index
#   7. .build_survey_design --- Extract and normalise survey variables
#   8. .compute_K           --- Random effects dimension
#   9. .determine_model_type --- Classify model variant
# ============================================================================


# ============================================================================
# 1. prepare_stan_data — Build Stan data list (exported)
# ============================================================================

#' Prepare Stan Data for a Hurdle Beta-Binomial Model
#'
#' Converts a formula specification and data frame into a validated list
#' suitable for passing to Stan. Constructs design matrices with optional
#' standardisation, builds group indices, normalises survey weights, and
#' computes all required dimensions.
#'
#' @param formula An object of class `"hbb_formula"` (from [hbb_formula()]),
#'   or a raw formula that will be parsed via `hbb_formula()`.
#' @param data A data frame containing provider-level variables: the
#'   response, trials, fixed effects, and optionally the grouping variable.
#' @param weights Character string naming the column in `data` containing
#'   survey sampling weights, or `NULL` (default) for unweighted analysis.
#' @param stratum Character string naming the column in `data` containing
#'   sampling stratum identifiers, or `NULL`.
#' @param psu Character string naming the column in `data` containing
#'   primary sampling unit (PSU) identifiers, or `NULL`.
#' @param state_data A data frame of group-level (state-level) variables
#'   for policy moderators. Required if `formula$policy` is not `NULL`.
#'   Must contain a column matching the grouping variable for merging.
#' @param prior A prior specification list, or `NULL` (default) to use
#'   package defaults. (Prior specification API is under development.)
#' @param center Logical. If `TRUE` (default), subtract column means from
#'   non-intercept columns of X and V. The means are stored in the output
#'   for back-transformation.
#' @param scale Logical. If `TRUE` (default), divide non-intercept columns
#'   of X and V by their standard deviations after centering. The SDs are
#'   stored in the output for back-transformation.
#'
#' @return An S3 object of class `"hbb_data"` (a list). See **Value** for
#'   full field descriptions.
#'
#' @section Value fields:
#' \describe{
#'   \item{`N`}{Integer. Number of observations.}
#'   \item{`P`}{Integer. Number of provider covariates including intercept.}
#'   \item{`S`}{Integer. Number of groups (states). 1 if no grouping.}
#'   \item{`Q`}{Integer. Number of policy covariates including intercept.}
#'   \item{`N_pos`}{Integer. Number of positive observations (z == 1).}
#'   \item{`K`}{Integer. Random effects dimension per group.}
#'   \item{`y`}{Integer vector of length N. Response counts.}
#'   \item{`n_trial`}{Integer vector of length N. Trial sizes.}
#'   \item{`z`}{Integer vector of length N. Participation indicators.}
#'   \item{`X`}{Numeric matrix N x P. Design matrix (column 1 = intercept).}
#'   \item{`state`}{Integer vector of length N. Group indices 1..S.}
#'   \item{`V`}{Numeric matrix S x Q, or `NULL` if no groups.}
#'   \item{`idx_pos`}{Integer vector. Row indices where z == 1.}
#'   \item{`w_tilde`}{Numeric vector of normalised weights, or `NULL`.}
#'   \item{`stratum_idx`}{Integer vector of stratum indices, or `NULL`.}
#'   \item{`psu_idx`}{Integer vector of PSU indices, or `NULL`.}
#'   \item{`n_strata`}{Integer. Number of unique strata, or `NULL`.}
#'   \item{`n_psu`}{Integer. Number of unique PSUs, or `NULL`.}
#'   \item{`x_center`}{Named numeric vector of column means, or `NULL`.}
#'   \item{`x_scale`}{Named numeric vector of column SDs, or `NULL`.}
#'   \item{`v_center`}{Named numeric vector of V column means, or `NULL`.}
#'   \item{`v_scale`}{Named numeric vector of V column SDs, or `NULL`.}
#'   \item{`formula`}{The `hbb_formula` object.}
#'   \item{`prior`}{The prior specification (or `NULL`).}
#'   \item{`group_levels`}{Character vector of original group labels.}
#'   \item{`model_type`}{Character. One of `"base"`, `"weighted"`,
#'     `"svc"`, `"svc_weighted"`.}
#' }
#'
#' @examples
#' # Minimal example with the small synthetic dataset
#' data(nsece_synth_small, package = "hurdlebb")
#' f <- hbb_formula(y | trials(n_trial) ~ poverty + urban)
#' d <- prepare_stan_data(f, nsece_synth_small)
#' d
#'
#' # With grouping and weights
#' f2 <- hbb_formula(
#'   y | trials(n_trial) ~ poverty + urban + (1 | state_id)
#' )
#' d2 <- prepare_stan_data(f2, nsece_synth_small, weights = "weight")
#' d2
#'
#' # Full SVC model with policy moderators
#' data(nsece_state_policy, package = "hurdlebb")
#' f3 <- hbb_formula(
#'   y | trials(n_trial) ~ poverty + urban +
#'     (poverty + urban | state_id) +
#'     state_level(mr_pctile)
#' )
#' d3 <- prepare_stan_data(
#'   f3, nsece_synth_small,
#'   weights = "weight", state_data = nsece_state_policy
#' )
#' d3
#'
#' @seealso [hbb_formula()], [validate_hbb_data()]
#' @family data
#' @export
prepare_stan_data <- function(formula,
                              data,
                              weights    = NULL,
                              stratum    = NULL,
                              psu        = NULL,
                              state_data = NULL,
                              prior      = NULL,
                              center     = TRUE,
                              scale      = TRUE) {

    # -- Input validation -------------------------------------------------------

    # Convert raw formula if needed
    if (inherits(formula, "formula") && !inherits(formula, "hbb_formula")) {
        formula <- hbb_formula(formula)
    }
    if (!inherits(formula, "hbb_formula")) {
        cli_abort(c(
            "Expected an {.cls hbb_formula} or a raw {.cls formula}.",
            "x" = "Got {.cls {class(formula)}}.",
            "i" = "Use {.fun hbb_formula} or pass a formula directly."
        ))
    }

    assert_data_frame(data, min.rows = 1L, .var.name = "data")
    assert_flag(center)
    assert_flag(scale)

    if (!is.null(weights)) {
        assert_string(weights, .var.name = "weights")
    }
    if (!is.null(stratum)) {
        assert_string(stratum, .var.name = "stratum")
    }
    if (!is.null(psu)) {
        assert_string(psu, .var.name = "psu")
    }

    # -- Validate formula against data ------------------------------------------
    validate_hbb_formula(formula, data, state_data)

    # -- Outcome vectors --------------------------------------------------------
    y       <- as.integer(data[[formula$response]])
    n_trial <- as.integer(data[[formula$trials]])
    z       <- as.integer(y > 0L)
    N       <- length(y)

    # -- Design matrix X --------------------------------------------------------
    x_result <- .build_X_matrix(formula$fixed, data, center = center,
                                scale = scale)
    X        <- x_result$X
    P        <- ncol(X)
    x_center <- x_result$means
    x_scale  <- x_result$sds

    # -- Group index ------------------------------------------------------------
    grp_result   <- .build_group_index(formula$group, data)
    state_idx    <- grp_result$state_idx
    S            <- grp_result$S
    group_levels <- grp_result$group_levels

    # -- Policy design matrix V -------------------------------------------------
    v_result <- .build_V_matrix(formula$policy, formula$group, state_data,
                                group_levels, S,
                                center = center, scale = scale)
    V        <- v_result$V
    Q        <- v_result$Q
    v_center <- v_result$means
    v_scale  <- v_result$sds

    # -- Positive index ---------------------------------------------------------
    idx_pos <- which(z == 1L)
    N_pos   <- length(idx_pos)

    # -- Survey design ----------------------------------------------------------
    survey <- .build_survey_design(data, weights, stratum, psu)

    # -- Random effects dimension -----------------------------------------------
    K <- .compute_K(formula, P)

    # -- Model type -------------------------------------------------------------
    has_weights <- !is.null(survey$w_tilde)
    model_type  <- .determine_model_type(formula$svc, has_weights)

    # -- Assemble return object -------------------------------------------------
    out <- structure(
        list(
            # Dimensions
            N     = N,
            P     = P,
            S     = S,
            Q     = Q,
            N_pos = N_pos,
            K     = K,
            # Outcome
            y       = y,
            n_trial = n_trial,
            z       = z,
            # Design matrices
            X     = X,
            state = state_idx,
            V     = V,
            # Positive index
            idx_pos = idx_pos,
            # Survey design
            w_tilde     = survey$w_tilde,
            stratum_idx = survey$stratum_idx,
            psu_idx     = survey$psu_idx,
            n_strata    = survey$n_strata,
            n_psu       = survey$n_psu,
            # Metadata for back-transformation
            x_center = x_center,
            x_scale  = x_scale,
            v_center = v_center,
            v_scale  = v_scale,
            # Links
            formula      = formula,
            prior        = prior,
            group_levels = group_levels,
            model_type   = model_type
        ),
        class = "hbb_data"
    )

    # -- Final validation -------------------------------------------------------
    validate_hbb_data(out)

    out
}


# ============================================================================
# 2. validate_hbb_data — Validate a Stan data list (exported)
# ============================================================================

#' Validate an hbb_data Object
#'
#' Performs comprehensive checks on a Stan data list to ensure internal
#' consistency before passing to CmdStan. Verifies field existence,
#' dimensional consistency, value ranges, and cross-field coherence.
#'
#' @param data_list An object of class `"hbb_data"`, as returned by
#'   [prepare_stan_data()].
#'
#' @return Invisibly returns `TRUE` if all checks pass. Throws an
#'   informative error otherwise.
#'
#' @examples
#' data(nsece_synth_small, package = "hurdlebb")
#' f <- hbb_formula(y | trials(n_trial) ~ poverty + urban)
#' d <- prepare_stan_data(f, nsece_synth_small)
#' validate_hbb_data(d)
#'
#' @seealso [prepare_stan_data()]
#' @family data
#' @export
validate_hbb_data <- function(data_list) {

    # -- Check class ------------------------------------------------------------
    if (!inherits(data_list, "hbb_data")) {
        cli_abort(c(
            "Expected an {.cls hbb_data} object.",
            "x" = "Got {.cls {class(data_list)}}.",
            "i" = "Use {.fun prepare_stan_data} to create data for Stan."
        ))
    }

    # -- Required fields --------------------------------------------------------
    required <- c("N", "P", "S", "Q", "N_pos", "K",
                  "y", "n_trial", "z",
                  "X", "state",
                  "idx_pos",
                  "formula", "model_type")
    missing_flds <- setdiff(required, names(data_list))
    if (length(missing_flds) > 0L) {
        cli_abort(c(
            "{qty(length(missing_flds))}Required field{?s} missing from {.cls hbb_data}.",
            "x" = "Missing: {.val {missing_flds}}."
        ))
    }

    N <- data_list$N
    P <- data_list$P
    S <- data_list$S
    Q <- data_list$Q

    # -- Dimensional consistency ------------------------------------------------

    # y, n_trial, z must be length N
    if (length(data_list$y) != N) {
        cli_abort("Length of {.field y} ({length(data_list$y)}) does not match {.field N} ({N}).")
    }
    if (length(data_list$n_trial) != N) {
        cli_abort("Length of {.field n_trial} ({length(data_list$n_trial)}) does not match {.field N} ({N}).")
    }
    if (length(data_list$z) != N) {
        cli_abort("Length of {.field z} ({length(data_list$z)}) does not match {.field N} ({N}).")
    }

    # X must be N x P
    if (!is.matrix(data_list$X)) {
        cli_abort("{.field X} must be a matrix.")
    }
    if (nrow(data_list$X) != N || ncol(data_list$X) != P) {
        cli_abort(
            "{.field X} dimensions ({nrow(data_list$X)} x {ncol(data_list$X)}) do not match N={N}, P={P}."
        )
    }

    # state must be length N
    if (length(data_list$state) != N) {
        cli_abort("Length of {.field state} ({length(data_list$state)}) does not match {.field N} ({N}).")
    }

    # V must be S x Q (if present)
    if (!is.null(data_list$V)) {
        if (!is.matrix(data_list$V)) {
            cli_abort("{.field V} must be a matrix (or NULL).")
        }
        if (nrow(data_list$V) != S || ncol(data_list$V) != Q) {
            cli_abort(
                "{.field V} dimensions ({nrow(data_list$V)} x {ncol(data_list$V)}) do not match S={S}, Q={Q}."
            )
        }
    }

    # N_pos must match idx_pos length
    if (length(data_list$idx_pos) != data_list$N_pos) {
        cli_abort(
            "Length of {.field idx_pos} ({length(data_list$idx_pos)}) does not match {.field N_pos} ({data_list$N_pos})."
        )
    }

    # -- Value ranges -----------------------------------------------------------

    # y >= 0
    if (any(data_list$y < 0L)) {
        cli_abort("{.field y} contains negative values.")
    }

    # n_trial >= 1
    if (any(data_list$n_trial < 1L)) {
        cli_abort("{.field n_trial} contains values < 1.")
    }

    # y <= n_trial
    if (any(data_list$y > data_list$n_trial)) {
        cli_abort("{.field y} exceeds {.field n_trial} for some observations.")
    }

    # z in {0, 1}
    if (!all(data_list$z %in% c(0L, 1L))) {
        cli_abort("{.field z} must contain only 0 and 1.")
    }

    # state in 1..S
    if (any(data_list$state < 1L) || any(data_list$state > S)) {
        cli_abort("{.field state} indices must be in 1..{S}.")
    }

    # -- Cross-consistency checks -----------------------------------------------

    # z == 1 iff y > 0
    expected_z <- as.integer(data_list$y > 0L)
    if (!identical(data_list$z, expected_z)) {
        cli_abort("{.field z} must equal {.code as.integer(y > 0)}.")
    }

    # idx_pos must point to z == 1 rows
    if (data_list$N_pos > 0L) {
        if (!all(data_list$z[data_list$idx_pos] == 1L)) {
            cli_abort("{.field idx_pos} contains indices where {.field z} != 1.")
        }
        # idx_pos should be exactly the indices where z == 1
        expected_idx <- which(data_list$z == 1L)
        if (!identical(data_list$idx_pos, expected_idx)) {
            cli_abort("{.field idx_pos} does not match {.code which(z == 1)}.")
        }
    }

    # -- Survey design consistency ----------------------------------------------
    if (!is.null(data_list$w_tilde)) {
        if (length(data_list$w_tilde) != N) {
            cli_abort(
                "Length of {.field w_tilde} ({length(data_list$w_tilde)}) does not match {.field N} ({N})."
            )
        }
        if (any(data_list$w_tilde <= 0)) {
            cli_abort("{.field w_tilde} must contain strictly positive values.")
        }
    }

    if (!is.null(data_list$stratum_idx)) {
        if (length(data_list$stratum_idx) != N) {
            cli_abort(
                "Length of {.field stratum_idx} ({length(data_list$stratum_idx)}) does not match {.field N} ({N})."
            )
        }
        if (is.null(data_list$n_strata)) {
            cli_abort("{.field n_strata} is required when {.field stratum_idx} is present.")
        }
        if (any(data_list$stratum_idx < 1L) ||
            any(data_list$stratum_idx > data_list$n_strata)) {
            cli_abort("{.field stratum_idx} values must be in 1..{data_list$n_strata}.")
        }
    }

    if (!is.null(data_list$psu_idx)) {
        if (length(data_list$psu_idx) != N) {
            cli_abort(
                "Length of {.field psu_idx} ({length(data_list$psu_idx)}) does not match {.field N} ({N})."
            )
        }
        if (is.null(data_list$n_psu)) {
            cli_abort("{.field n_psu} is required when {.field psu_idx} is present.")
        }
        if (any(data_list$psu_idx < 1L) ||
            any(data_list$psu_idx > data_list$n_psu)) {
            cli_abort("{.field psu_idx} values must be in 1..{data_list$n_psu}.")
        }
    }

    # -- Model type check -------------------------------------------------------
    valid_types <- c("base", "weighted", "svc", "svc_weighted")
    if (!data_list$model_type %in% valid_types) {
        cli_abort(
            "{.field model_type} must be one of {.val {valid_types}}, got {.val {data_list$model_type}}."
        )
    }

    invisible(TRUE)
}


# ============================================================================
# 3. print.hbb_data — Compact summary print method
# ============================================================================

#' @export
print.hbb_data <- function(x, ...) {

    zero_rate <- round(1 - mean(x$z), 3)

    cat("Hurdle Beta-Binomial Data\n")
    cat("-------------------------\n")
    cat("  Observations (N) :", x$N, "\n")
    cat("  Covariates   (P) :", x$P,
        paste0("(intercept + ", x$P - 1L, " predictors)"), "\n")
    cat("  Groups       (S) :", x$S, "\n")
    cat("  Policy vars  (Q) :", x$Q, "\n")
    cat("  Positive obs     :", x$N_pos,
        paste0("(", round(100 * x$N_pos / x$N, 1), "%)"), "\n")
    cat("  Zero rate        :", zero_rate, "\n")
    cat("  RE dimension (K) :", x$K, "\n")
    cat("  Model type       :", x$model_type, "\n")

    # Survey design info
    if (!is.null(x$w_tilde)) {
        cat("  Survey weights   : yes",
            paste0("(range ", round(min(x$w_tilde), 2), "--",
                   round(max(x$w_tilde), 2), ")"), "\n")
    } else {
        cat("  Survey weights   : no\n")
    }
    if (!is.null(x$stratum_idx)) {
        cat("  Strata           :", x$n_strata, "\n")
    }
    if (!is.null(x$psu_idx)) {
        cat("  PSUs             :", x$n_psu, "\n")
    }

    # Standardisation info
    if (!is.null(x$x_center) || !is.null(x$x_scale)) {
        cat("  X standardised   : yes\n")
    } else {
        cat("  X standardised   : no\n")
    }

    invisible(x)
}


# ============================================================================
# 4. .build_X_matrix — Design matrix with optional standardisation
# ============================================================================

#' Build the provider-level design matrix X
#'
#' Constructs the N x P design matrix from fixed effect variable names.
#' Always prepends an intercept column. Optionally centers and/or scales
#' non-intercept columns, recording the transformation parameters.
#'
#' @param fixed Character vector of fixed effect variable names.
#' @param data Data frame.
#' @param center Logical.
#' @param scale Logical.
#' @return A list with `X` (matrix), `means` (named numeric or NULL),
#'   `sds` (named numeric or NULL).
#' @noRd
.build_X_matrix <- function(fixed, data, center, scale) {

    N <- nrow(data)

    # Start with intercept
    X <- matrix(1, nrow = N, ncol = 1L)
    colnames(X) <- "(Intercept)"

    means <- NULL
    sds   <- NULL

    if (length(fixed) > 0L) {
        # Extract columns as matrix
        X_pred <- do.call(cbind, lapply(fixed, function(v) {
            as.numeric(data[[v]])
        }))
        colnames(X_pred) <- fixed

        # Standardise non-intercept columns
        if (center || scale) {
            col_means <- colMeans(X_pred, na.rm = TRUE)
            col_sds   <- apply(X_pred, 2L, sd, na.rm = TRUE)

            # Guard against zero SD (constant columns)
            col_sds[col_sds < .Machine$double.eps] <- 1

            if (center) {
                means <- setNames(col_means, fixed)
                X_pred <- sweep(X_pred, 2L, col_means, "-")
            }
            if (scale) {
                sds <- setNames(col_sds, fixed)
                X_pred <- sweep(X_pred, 2L, col_sds, "/")
            }
        }

        X <- cbind(X, X_pred)
    }

    list(X = X, means = means, sds = sds)
}


# ============================================================================
# 5. .build_V_matrix — Policy design matrix with optional standardisation
# ============================================================================

#' Build the group-level policy design matrix V
#'
#' Constructs the S x Q design matrix from policy moderator variables.
#' Always prepends an intercept column. Returns NULL for V and Q = 1
#' when there are no groups.
#'
#' @param policy Character vector of policy variable names (or NULL).
#' @param group Character. The grouping variable name (or NULL).
#' @param state_data Data frame of group-level data (or NULL).
#' @param group_levels Character vector of original group labels in order.
#' @param S Integer. Number of groups.
#' @param center Logical.
#' @param scale Logical.
#' @return A list with `V` (matrix or NULL), `Q` (integer),
#'   `means` (named numeric or NULL), `sds` (named numeric or NULL).
#' @noRd
.build_V_matrix <- function(policy, group, state_data, group_levels, S,
                            center, scale) {

    means <- NULL
    sds   <- NULL

    # No grouping at all: no V matrix needed
    if (is.null(group)) {
        return(list(V = NULL, Q = 1L, means = NULL, sds = NULL))
    }

    # Groups exist but no policy moderators: intercept-only V
    if (is.null(policy) || is.null(state_data)) {
        V <- matrix(1, nrow = S, ncol = 1L)
        colnames(V) <- "(Intercept)"
        return(list(V = V, Q = 1L, means = NULL, sds = NULL))
    }

    # Merge state_data to match group_levels ordering
    # group_levels is the character vector of original labels in 1..S order
    sd_row_order <- match(group_levels, as.character(state_data[[group]]))
    if (anyNA(sd_row_order)) {
        cli_abort(c(
            "Not all groups found in {.arg state_data}.",
            "i" = "This should have been caught by {.fun validate_hbb_formula}."
        ))
    }
    state_data_ordered <- state_data[sd_row_order, , drop = FALSE]

    # Build predictor columns
    V_pred <- do.call(cbind, lapply(policy, function(v) {
        as.numeric(state_data_ordered[[v]])
    }))
    colnames(V_pred) <- policy

    # Standardise
    if (center || scale) {
        col_means <- colMeans(V_pred, na.rm = TRUE)
        col_sds   <- apply(V_pred, 2L, sd, na.rm = TRUE)
        col_sds[col_sds < .Machine$double.eps] <- 1

        if (center) {
            means <- setNames(col_means, policy)
            V_pred <- sweep(V_pred, 2L, col_means, "-")
        }
        if (scale) {
            sds <- setNames(col_sds, policy)
            V_pred <- sweep(V_pred, 2L, col_sds, "/")
        }
    }

    # Prepend intercept
    V <- cbind(matrix(1, nrow = S, ncol = 1L), V_pred)
    colnames(V) <- c("(Intercept)", policy)
    Q <- ncol(V)

    list(V = V, Q = Q, means = means, sds = sds)
}


# ============================================================================
# 6. .build_group_index — Sequential integer group index
# ============================================================================

#' Create sequential 1..S group indices from a grouping variable
#'
#' @param group Character. Grouping variable name (or NULL).
#' @param data Data frame.
#' @return A list with `state_idx` (integer vector of length N),
#'   `S` (integer), `group_levels` (character vector of original labels).
#' @noRd
.build_group_index <- function(group, data) {

    N <- nrow(data)

    if (is.null(group)) {
        return(list(
            state_idx    = rep(1L, N),
            S            = 1L,
            group_levels = "1"
        ))
    }

    grp_raw <- data[[group]]

    # Convert to factor to get stable ordering
    if (is.factor(grp_raw)) {
        grp_factor <- grp_raw
    } else {
        # Sort levels for consistent ordering
        grp_factor <- factor(grp_raw, levels = sort(unique(grp_raw)))
    }

    state_idx    <- as.integer(grp_factor)
    group_levels <- as.character(levels(grp_factor))
    S            <- length(group_levels)

    list(
        state_idx    = state_idx,
        S            = S,
        group_levels = group_levels
    )
}


# ============================================================================
# 7. .build_survey_design — Extract and normalise survey variables
# ============================================================================

#' Extract and normalise survey design variables
#'
#' @param data Data frame.
#' @param weights Character string column name (or NULL).
#' @param stratum Character string column name (or NULL).
#' @param psu Character string column name (or NULL).
#' @return A list with `w_tilde`, `stratum_idx`, `psu_idx`, `n_strata`,
#'   `n_psu` (all NULL if not specified).
#' @noRd
.build_survey_design <- function(data, weights, stratum, psu) {

    w_tilde     <- NULL
    stratum_idx <- NULL
    psu_idx     <- NULL
    n_strata    <- NULL
    n_psu       <- NULL

    # -- Weights ----------------------------------------------------------------
    if (!is.null(weights)) {
        if (!weights %in% names(data)) {
            cli_abort(c(
                "Weight variable {.val {weights}} not found in {.arg data}.",
                "i" = "Available columns: {.val {head(names(data), 20)}}."
            ))
        }
        w_raw <- data[[weights]]
        if (!is.numeric(w_raw)) {
            cli_abort(c(
                "Weight variable {.val {weights}} must be numeric.",
                "x" = "Got {.cls {class(w_raw)}}."
            ))
        }
        if (anyNA(w_raw)) {
            cli_abort("Weight variable {.val {weights}} contains missing values.")
        }
        if (any(w_raw <= 0)) {
            cli_abort("Weight variable {.val {weights}} must contain strictly positive values.")
        }
        w_tilde <- normalize_weights(w_raw)
    }

    # -- Stratum ----------------------------------------------------------------
    if (!is.null(stratum)) {
        if (!stratum %in% names(data)) {
            cli_abort(c(
                "Stratum variable {.val {stratum}} not found in {.arg data}.",
                "i" = "Available columns: {.val {head(names(data), 20)}}."
            ))
        }
        strat_raw <- data[[stratum]]
        strat_fac <- factor(strat_raw, levels = sort(unique(strat_raw)))
        stratum_idx <- as.integer(strat_fac)
        n_strata    <- nlevels(strat_fac)
    }

    # -- PSU --------------------------------------------------------------------
    if (!is.null(psu)) {
        if (!psu %in% names(data)) {
            cli_abort(c(
                "PSU variable {.val {psu}} not found in {.arg data}.",
                "i" = "Available columns: {.val {head(names(data), 20)}}."
            ))
        }
        psu_raw <- data[[psu]]
        psu_fac <- factor(psu_raw, levels = sort(unique(psu_raw)))
        psu_idx <- as.integer(psu_fac)
        n_psu   <- nlevels(psu_fac)
    }

    list(
        w_tilde     = w_tilde,
        stratum_idx = stratum_idx,
        psu_idx     = psu_idx,
        n_strata    = n_strata,
        n_psu       = n_psu
    )
}


# ============================================================================
# 8. .compute_K — Random effects dimension
# ============================================================================

#' Compute the random effects dimension K
#'
#' Determines the total dimension of the random effects vector per group.
#'
#' @details
#' - No random effects (group is NULL): K = 0
#' - Random intercept only (random == "1"): K = 2 (one per margin)
#' - SVC (random slopes): K = 2 * P (intercept + slopes for both margins)
#'
#' @param formula An `hbb_formula` object.
#' @param P Integer. Number of columns in X (including intercept).
#' @return Integer.
#' @noRd
.compute_K <- function(formula, P) {

    # No grouping means no random effects
    if (is.null(formula$group)) {
        return(0L)
    }

    # Random intercept only
    if (identical(formula$random, "1")) {
        return(2L)
    }

    # SVC: both margins get intercept + all random slopes
    # P already includes the intercept column
    2L * P
}


# ============================================================================
# 9. .determine_model_type — Classify model variant
# ============================================================================

#' Determine the model type string
#'
#' @param svc Logical. Whether state-varying coefficients are present.
#' @param has_weights Logical. Whether survey weights are present.
#' @return Character string: one of "base", "weighted", "svc",
#'   "svc_weighted".
#' @noRd
.determine_model_type <- function(svc, has_weights) {
    if (svc && has_weights) {
        "svc_weighted"
    } else if (svc && !has_weights) {
        "svc"
    } else if (!svc && has_weights) {
        "weighted"
    } else {
        "base"
    }
}
