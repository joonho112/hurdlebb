# ============================================================================
# predict.R --- S3 Predict Method for hbb_fit Objects
#
# Provides out-of-sample and in-sample prediction for Hurdle Beta-Binomial
# models, with optional posterior credible intervals.
#
# Theory:
#   The predicted enrollment proportion for observation i with covariate
#   vector x_i is
#     E[Y_i / n_i] = q(x_i) * g(mu(x_i))
#   where
#     q(x_i)  = logistic(x_i' alpha)           [extensive margin]
#     mu(x_i) = logistic(x_i' beta)            [intensive margin]
#     g(mu)   = mu / (1 - p_0(n, mu, kappa))   [ZT intensity, Eq. 2.5]
#
#   For SVC models with known state s(i), the linear predictor includes
#   state random effects:
#     eta_ext(x_i) = x_i' alpha + x_i' delta_ext[s(i)]
#     eta_int(x_i) = x_i' beta  + x_i' delta_int[s(i)]
#
#   Credible intervals are computed by propagating the full posterior
#   draw matrix through the inverse-link function:
#     theta^{(m)} -> q^{(m)}(x_i) * mu^{(m)}(x_i)
#   and taking empirical quantiles across draws m = 1, ..., M.
#
#   When newdata is provided, the design matrix is constructed by
#   applying the same centering and scaling transformations that were
#   used during training (stored in hbb_data$x_center and x_scale).
#
#
# Existing helpers reused (NOT redefined here):
#   .validate_hbb_fit_methods(x)      from methods.R
#   .validate_level(level)            from cholesky-transform.R
#   .make_stan_param_names(P)         from cholesky-transform.R
#   .build_param_labels(hbb_data)     from sandwich.R
#   .extract_delta_means(fit, S, K)   from sandwich.R
#
# Contents:
#   Exported:
#     1. predict.hbb_fit           --- S3 predict method
#   Internal:
#     2. .predict_build_X          --- Design matrix reconstruction
#     3. .predict_resolve_delta    --- State random-effect resolution
#     4. .predict_compute_point    --- Vectorised point prediction
#     5. .predict_compute_interval --- Credible intervals via draws
# ============================================================================


# ============================================================================
# 1. predict.hbb_fit --- S3 predict method (exported)
# ============================================================================

#' Predictions from a Hurdle Beta-Binomial Model
#'
#' Generates point predictions and optional posterior credible intervals
#' from a fitted hurdle Beta-Binomial model.  Supports both in-sample
#' prediction (using the training data) and out-of-sample prediction
#' (using \code{newdata}).
#'
#' @description
#' Point predictions are computed at the posterior mean of the fixed
#' effects (and state random effects for SVC models):
#' \deqn{
#'   \hat{y}_i = \hat{q}(x_i) \cdot \hat\mu(x_i),
#' }
#' where \eqn{\hat{q}(x) = \mathrm{logistic}(x'\hat\alpha)} and
#' \eqn{\hat\mu(x) = \mathrm{logistic}(x'\hat\beta)}.
#'
#' When \code{interval = "credible"}, the full posterior distribution is
#' propagated through the link function to produce empirical credible
#' intervals.  For each MCMC draw \eqn{m}:
#' \deqn{
#'   \hat{y}_i^{(m)} = q^{(m)}(x_i) \cdot \mu^{(m)}(x_i),
#' }
#' and the \eqn{(\alpha/2,\; 1-\alpha/2)} quantiles across draws give
#' the interval bounds.
#'
#' @section Prediction types:
#' Three types of predictions are available:
#' \describe{
#'   \item{\code{"response"}}{(Default.) The composite enrollment
#'     proportion \eqn{\hat{q}_i \cdot \hat\mu_i \in [0,1]}.  This
#'     represents the unconditional expected proportion
#'     \eqn{E[Y_i/n_i]}.}
#'   \item{\code{"extensive"}}{The participation probability
#'     \eqn{\hat{q}_i = \mathrm{logistic}(x_i'\hat\alpha) \in (0,1)}.}
#'   \item{\code{"intensive"}}{The conditional enrollment share
#'     \eqn{\hat\mu_i = \mathrm{logistic}(x_i'\hat\beta) \in (0,1)}.}
#' }
#'
#' @section Newdata design matrix:
#' When \code{newdata} is provided, the design matrix is constructed by:
#' \enumerate{
#'   \item Extracting the covariate columns named in
#'     \code{object$formula$fixed} from \code{newdata}.
#'   \item Applying the same centering and scaling as the training data,
#'     using the stored \code{hbb_data$x_center} (column means) and
#'     \code{hbb_data$x_scale} (column SDs).  This ensures that
#'     regression coefficients remain on the standardised scale.
#'   \item Prepending an intercept column of ones.
#' }
#' All covariate columns must be present in \code{newdata}.  Missing
#' columns cause an informative error.
#'
#' @section State-varying coefficients (SVC):
#' For SVC models, out-of-sample prediction requires specifying
#' \code{state} (either a column name in \code{newdata} or a vector
#' of state labels).  The state labels must match
#' \code{hbb_data$group_levels} from the training data.  Unknown
#' state labels cause an error.
#'
#' If \code{state = NULL} for an SVC model, predictions use fixed
#' effects only (i.e., the population-average prediction without state
#' random effects), and a warning is issued.
#'
#' @section Credible interval theory:
#' The posterior credible interval for prediction \eqn{i} at level
#' \eqn{1-\alpha} is
#' \deqn{
#'   \bigl[Q_{\alpha/2}(\hat{y}_i^{(1:M)}),\;
#'         Q_{1-\alpha/2}(\hat{y}_i^{(1:M)})\bigr],
#' }
#' where \eqn{Q_p} denotes the empirical \eqn{p}-quantile over \eqn{M}
#' posterior draws.  This is a \emph{conditional} credible interval for
#' the expected proportion \eqn{E[Y_i/n_i | \theta]}, not a predictive
#' interval for \eqn{Y_i} itself (which would additionally account for
#' sampling variability from the Beta-Binomial).
#'
#' For computational efficiency, \code{ndraws} can be used to thin the
#' draws before propagation.  With \eqn{M = 4000} total draws, setting
#' \code{ndraws = 500} reduces computation 8-fold while typically
#' preserving interval accuracy.
#'
#' @param object An object of class \code{"hbb_fit"} returned by
#'   \code{\link{hbb}}.
#' @param newdata An optional data frame containing the covariates for
#'   prediction.  If \code{NULL} (default), in-sample predictions are
#'   returned using the training design matrix.
#' @param type Character string: \code{"response"} (default),
#'   \code{"extensive"}, or \code{"intensive"}.  See \strong{Prediction
#'   types} above.
#' @param interval Character string: \code{"none"} (default) for point
#'   predictions only, or \code{"credible"} for posterior credible
#'   intervals.
#' @param level Numeric scalar in \eqn{(0, 1)}.  Credible level when
#'   \code{interval = "credible"}.  Default is \code{0.95}.
#' @param ndraws Optional positive integer.  If specified, the posterior
#'   draws are thinned to \code{ndraws} before computing intervals.
#'   Ignored when \code{interval = "none"}.  Default is \code{NULL}
#'   (use all draws).
#' @param state Character string or character vector specifying states
#'   for SVC prediction.
#'   \describe{
#'     \item{Character of length 1:}{If \code{newdata} is provided,
#'       interpreted as a column name in \code{newdata}.  Otherwise,
#'       interpreted as a single state label to apply to all
#'       observations.}
#'     \item{Character vector of length \code{nrow(newdata)}:}{Direct
#'       state labels for each observation.}
#'     \item{\code{NULL}:}{No state information.  For SVC models,
#'       predictions use fixed effects only (with a warning).}
#'   }
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with one row per observation.  Columns:
#' \describe{
#'   \item{\code{fit}}{Numeric: point prediction (posterior mean).}
#'   \item{\code{lwr}}{Numeric: lower credible bound.  Present only
#'     when \code{interval = "credible"}.}
#'   \item{\code{upr}}{Numeric: upper credible bound.  Present only
#'     when \code{interval = "credible"}.}
#' }
#'
#' @seealso
#' \code{\link{fitted.hbb_fit}} for in-sample fitted values (simpler
#' interface without newdata support),
#' \code{\link{residuals.hbb_fit}} for residuals,
#' \code{\link{summary.hbb_fit}} for model summary.
#'
#' @references
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B.,
#' Vehtari, A., and Rubin, D. B. (2013).
#' \emph{Bayesian Data Analysis} (3rd ed.).  Chapman and Hall/CRC.
#'
#' Ghosal, R., Ghosh, S. K., and Maiti, T. (2020).
#' Two-part regression models for longitudinal zero-inflated count data.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' \strong{183}(4), 1603--1626.
#'
#' @examples
#' \dontrun{
#' fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
#'
#' # In-sample point predictions
#' pred <- predict(fit)
#'
#' # Out-of-sample with credible intervals
#' new_df <- data.frame(poverty = c(0.1, 0.3, 0.5),
#'                      urban   = c(1, 0, 1))
#' pred_ci <- predict(fit, newdata = new_df,
#'                    interval = "credible", level = 0.95)
#'
#' # Extensive-margin predictions only
#' pred_ext <- predict(fit, type = "extensive")
#'
#' # SVC model with state assignment
#' pred_svc <- predict(fit_svc, newdata = new_df,
#'                     state = c("AL", "CA", "NY"),
#'                     interval = "credible")
#' }
#'
#' @importFrom stats complete.cases
#' @method predict hbb_fit
#' @export
predict.hbb_fit <- function(object,
                              newdata = NULL,
                              type = c("response", "extensive", "intensive"),
                              interval = c("none", "credible"),
                              level = 0.95,
                              ndraws = NULL,
                              state = NULL,
                              ...) {

    # ========================================================================
    # 0. Comprehensive input validation
    # ========================================================================

    .validate_hbb_fit_methods(object)

    hd <- object$hbb_data
    P  <- hd$P
    D  <- 2L * P + 1L

    # -- Type validation ---------------------------------------------------
    type <- tryCatch(
        match.arg(type),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "{.arg type} must be one of {.val response}, \\
                       {.val extensive}, or {.val intensive}.",
                "i" = "Received: {.val {type[1]}}."
            ))
        }
    )

    # -- Interval validation -------------------------------------------------
    interval <- tryCatch(
        match.arg(interval),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "{.arg interval} must be one of {.val none} or \\
                       {.val credible}.",
                "i" = "Received: {.val {interval[1]}}."
            ))
        }
    )

    # -- Level validation (only when interval != "none") ---------------------
    if (interval != "none") {
        .validate_level(level)
    }

    # -- ndraws validation -------------------------------------------------
    if (!is.null(ndraws)) {
        if (!is.numeric(ndraws) || length(ndraws) != 1L ||
            is.na(ndraws) || ndraws < 1L) {
            cli::cli_abort(c(
                "x" = "{.arg ndraws} must be a positive integer.",
                "i" = "Received: {.val {ndraws}}."
            ))
        }
        ndraws <- as.integer(ndraws)
    }

    # -- newdata validation ------------------------------------------------
    if (!is.null(newdata)) {
        if (!is.data.frame(newdata)) {
            cli::cli_abort(c(
                "x" = "{.arg newdata} must be a data.frame.",
                "i" = "Received class: {.cls {class(newdata)}}."
            ))
        }
        if (nrow(newdata) == 0L) {
            cli::cli_abort(c(
                "x" = "{.arg newdata} has zero rows.",
                "i" = "Provide at least one observation."
            ))
        }
        # Check for NAs and warn
        n_na <- sum(!complete.cases(newdata))
        if (n_na > 0L) {
            cli::cli_warn(c(
                "!" = "{n_na} row{?s} in {.arg newdata} contain NA values.",
                "i" = "Predictions will be NA for these rows."
            ))
        }
    }

    # -- State validation --------------------------------------------------
    is_svc <- isTRUE(object$model_type %in% c("svc", "svc_weighted"))

    if (!is.null(state) && !is_svc) {
        cli::cli_abort(c(
            "x" = "{.arg state} can only be used with SVC models.",
            "i" = "Current model type: {.val {object$model_type}}.",
            "i" = "Remove the {.arg state} argument to proceed."
        ))
    }

    # ========================================================================
    # 1. Build design matrix
    # ========================================================================

    X_new <- .predict_build_X(object, newdata)
    N_new <- nrow(X_new)

    # ========================================================================
    # 2. Resolve SVC delta
    # ========================================================================

    delta_info <- .predict_resolve_delta(object, state, N_new, newdata)

    # ========================================================================
    # 3. Point predictions
    # ========================================================================

    # -- Resolve n_trial for ZT correction ------------------------------------
    if (is.null(newdata)) {
        n_trial_pred <- as.integer(object$hbb_data$n_trial)
    } else {
        # Try to extract trials variable from newdata
        trials_var <- object$formula$trials
        if (!is.null(trials_var) && trials_var %in% names(newdata)) {
            n_trial_pred <- as.integer(newdata[[trials_var]])
        } else {
            # Use median from training data as fallback
            n_trial_pred <- as.integer(median(object$hbb_data$n_trial))
            if (type != "extensive") {
                cli::cli_warn(c(
                    "!" = "Trials variable not found in {.arg newdata}.",
                    "i" = "Using median training n_trial = {n_trial_pred} \\
                           for ZT correction.",
                    "i" = "Include the trials column in newdata for exact \\
                           predictions."
                ))
            }
        }
    }

    point_pred <- .predict_compute_point(
        object  = object,
        X_new   = X_new,
        P       = P,
        delta   = delta_info,
        type    = type,
        n_trial = n_trial_pred
    )

    # ========================================================================
    # 4. Credible intervals (if requested)
    # ========================================================================

    if (interval == "credible") {
        ci <- .predict_compute_interval(
            object  = object,
            X_new   = X_new,
            P       = P,
            delta   = delta_info,
            type    = type,
            level   = level,
            ndraws  = ndraws,
            n_trial = n_trial_pred
        )

        result <- data.frame(
            fit = point_pred,
            lwr = ci$lwr,
            upr = ci$upr,
            stringsAsFactors = FALSE,
            row.names = NULL
        )
    } else {
        result <- data.frame(
            fit = point_pred,
            stringsAsFactors = FALSE,
            row.names = NULL
        )
    }

    result
}


# ============================================================================
# 2. .predict_build_X --- Design matrix reconstruction (internal)
# ============================================================================

#' Build the design matrix for prediction
#'
#' If \code{newdata} is \code{NULL}, returns the training design matrix
#' \code{object$hbb_data$X}.  Otherwise, constructs a new design matrix
#' by extracting the covariates named in \code{object$formula$fixed},
#' applying the stored centering and scaling transformations, and
#' prepending an intercept column.
#'
#' @section Centering and scaling:
#' The stored transformations \code{hbb_data$x_center} (named numeric
#' vector of training column means) and \code{hbb_data$x_scale} (named
#' numeric vector of training column SDs) are applied to newdata
#' covariates to ensure the design matrix is on the same scale as the
#' training data.  The centering/scaling formula for covariate \eqn{j}
#' is:
#' \deqn{
#'   x_{ij}^* = \frac{x_{ij} - \bar{x}_j^{\mathrm{train}}}
#'                   {s_j^{\mathrm{train}}},
#' }
#' where \eqn{\bar{x}_j^{\mathrm{train}}} and
#' \eqn{s_j^{\mathrm{train}}} are the mean and SD from the training
#' data.
#'
#' @section Transformation steps:
#' \enumerate{
#'   \item Extract columns named in \code{object$formula$fixed}.
#'   \item If training centred (\code{x_center} is not \code{NULL}),
#'     subtract the \emph{training} column means via \code{sweep()}.
#'   \item If training scaled (\code{x_scale} is not \code{NULL}),
#'     divide by the \emph{training} column SDs via \code{sweep()}.
#'   \item Prepend an intercept column of ones.
#'   \item Validate that \code{ncol(X_new) == P}.
#' }
#'
#' @param object An \code{hbb_fit} object.
#' @param newdata A data frame or \code{NULL}.
#' @return Numeric matrix of dimension \eqn{N_{\mathrm{new}} \times P}
#'   (with intercept column).
#' @keywords internal
.predict_build_X <- function(object, newdata) {

    if (is.null(newdata)) {
        return(object$hbb_data$X)
    }

    hd          <- object$hbb_data
    P           <- hd$P
    fixed_names <- object$formula$fixed

    # ---- 1. Validate that all required columns exist ----
    missing_cols <- setdiff(fixed_names, names(newdata))
    if (length(missing_cols) > 0L) {
        cli::cli_abort(c(
            "x" = "{.arg newdata} is missing required covariate{?s}: \\
                   {.field {missing_cols}}.",
            "i" = "The model formula requires: {.val {fixed_names}}.",
            "i" = "Available columns: {.val {head(names(newdata), 15)}}."
        ))
    }

    N_new <- nrow(newdata)

    # ---- 2. Intercept-only model (short circuit) ----
    if (length(fixed_names) == 0L) {
        X_new <- matrix(1, nrow = N_new, ncol = 1L)
        colnames(X_new) <- "(Intercept)"
        if (ncol(X_new) != P) {
            cli::cli_abort(c(
                "x" = "Design matrix column mismatch: constructed \\
                       {ncol(X_new)} but model has P = {P}."
            ))
        }
        return(X_new)
    }

    # ---- 3. Build raw predictor matrix ----
    X_pred <- do.call(cbind, lapply(fixed_names, function(v) {
        col_data <- newdata[[v]]
        if (!is.numeric(col_data)) {
            cli::cli_abort(c(
                "x" = "Column {.field {v}} in {.arg newdata} must be numeric.",
                "i" = "Got {.cls {class(col_data)}}."
            ))
        }
        as.numeric(col_data)
    }))
    colnames(X_pred) <- fixed_names

    # Check for NA/NaN/Inf
    n_bad <- sum(!is.finite(X_pred))
    if (n_bad > 0L) {
        cli::cli_warn(c(
            "!" = "{n_bad} non-finite value{?s} in newdata predictors.",
            "i" = "Predictions may contain NA."
        ))
    }

    # ---- 4. Apply training centering via sweep ----
    x_center <- hd$x_center
    if (!is.null(x_center)) {
        center_vals <- x_center[fixed_names]
        if (anyNA(center_vals)) {
            # Some variables not covered
            uncovered <- fixed_names[is.na(center_vals)]
            cli::cli_warn(c(
                "!" = "Centering info missing for covariate{?s}: \\
                       {.val {uncovered}}.",
                "i" = "These columns will NOT be centred."
            ))
            center_vals[is.na(center_vals)] <- 0
        }
        X_pred <- sweep(X_pred, 2L, center_vals, "-")
    }

    # ---- 5. Apply training scaling via sweep ----
    x_scale <- hd$x_scale
    if (!is.null(x_scale)) {
        scale_vals <- x_scale[fixed_names]
        if (anyNA(scale_vals)) {
            uncovered <- fixed_names[is.na(scale_vals)]
            cli::cli_warn(c(
                "!" = "Scaling info missing for covariate{?s}: \\
                       {.val {uncovered}}.",
                "i" = "These columns will NOT be scaled."
            ))
            scale_vals[is.na(scale_vals)] <- 1
        }
        # Zero-SD guard
        zero_sd <- which(scale_vals == 0)
        if (length(zero_sd) > 0L) {
            cli::cli_warn(c(
                "!" = "Training-set SD is zero for: \\
                       {.val {fixed_names[zero_sd]}}.",
                "i" = "These columns will be centred but not scaled."
            ))
            scale_vals[zero_sd] <- 1
        }
        X_pred <- sweep(X_pred, 2L, scale_vals, "/")
    }

    # ---- 6. Prepend intercept ----
    X_new <- cbind(matrix(1, nrow = N_new, ncol = 1L), X_pred)
    colnames(X_new) <- c("(Intercept)", fixed_names)

    # ---- 7. Validate dimensions ----
    if (ncol(X_new) != P) {
        cli::cli_abort(c(
            "x" = "Design matrix column mismatch: constructed \\
                   {ncol(X_new)} columns but model has P = {P}.",
            "i" = "Model fixed effects: {.val {fixed_names}}.",
            "i" = "Check that {.arg newdata} matches the training formula."
        ))
    }

    X_new
}


# ============================================================================
# 3. .predict_resolve_delta --- State random-effect resolution (internal)
# ============================================================================

#' Resolve state random effects for prediction
#'
#' For non-SVC models or when \code{state} is \code{NULL}, returns
#' \code{NULL}.  For SVC models with state information, maps state
#' labels to integer indices via \code{hbb_data$group_levels} and
#' returns the posterior mean delta matrix along with the state index
#' vector.
#'
#' @param object An \code{hbb_fit} object.
#' @param state State specification (character vector, column name, or
#'   NULL).
#' @param N_new Number of prediction observations.
#' @param newdata Data frame or NULL.
#' @return \code{NULL} if no SVC adjustment, or a list with:
#'   \describe{
#'     \item{\code{delta_hat}}{S x K matrix of posterior mean deltas.}
#'     \item{\code{state_idx}}{Integer vector of length N_new mapping
#'       observations to state indices.}
#'   }
#' @keywords internal
.predict_resolve_delta <- function(object, state, N_new, newdata = NULL) {

    is_svc <- isTRUE(object$model_type %in% c("svc", "svc_weighted"))

    if (!is_svc) {
        if (!is.null(state)) {
            cli::cli_warn(c(
                "!" = "{.arg state} is ignored for non-SVC models.",
                "i" = "Model type is {.val {object$model_type}}."
            ))
        }
        return(NULL)
    }

    # SVC model but no state provided
    if (is.null(state)) {
        cli::cli_warn(c(
            "!" = "SVC model detected but {.arg state} is NULL.",
            "i" = "Predictions will use fixed effects only \\
                   (population-average).",
            "i" = "Provide {.arg state} for state-specific predictions."
        ))
        return(NULL)
    }

    hd <- object$hbb_data
    P  <- hd$P
    S  <- hd$S
    K  <- 2L * P

    # ---- Resolve state vector ----
    if (is.character(state) && length(state) == 1L && !is.null(newdata)) {
        # Interpret as column name in newdata
        if (!state %in% names(newdata)) {
            cli::cli_abort(c(
                "x" = "Column {.field {state}} not found in {.arg newdata}.",
                "i" = "Available columns: {.val {names(newdata)}}."
            ))
        }
        state_labels <- as.character(newdata[[state]])
    } else if (is.character(state) && length(state) == N_new) {
        state_labels <- state
    } else if (is.character(state) && length(state) == 1L && is.null(newdata)) {
        # Single state label for all in-sample observations
        state_labels <- rep(state, N_new)
    } else {
        cli::cli_abort(c(
            "x" = "{.arg state} must be a column name (length 1) or a \\
                   character vector of length {N_new}.",
            "i" = "Received: length {length(state)}, class {.cls {class(state)}}."
        ))
    }

    # ---- Map labels to integer indices via group_levels ----
    group_levels <- hd$group_levels
    if (is.null(group_levels)) {
        cli::cli_abort(c(
            "x" = "{.field group_levels} is NULL in {.field hbb_data}.",
            "i" = "Cannot map state labels to indices."
        ))
    }

    state_idx <- match(state_labels, group_levels)
    unknown   <- which(is.na(state_idx))
    if (length(unknown) > 0L) {
        bad_labels <- unique(state_labels[unknown])
        n_show     <- min(length(bad_labels), 10L)
        cli::cli_abort(c(
            "x" = "{length(unknown)} observation{?s} have unknown \\
                   state label{?s}.",
            "i" = "Unknown labels: {.val {bad_labels[seq_len(n_show)]}}.",
            "i" = "Known states: {.val {head(group_levels, 10)}} \\
                   ({length(group_levels)} total)."
        ))
    }

    # ---- Extract posterior mean delta matrix (S x K) ----
    delta_hat <- tryCatch(
        .extract_delta_means(object$fit, S, K),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "Failed to extract delta random effects for prediction.",
                "i" = conditionMessage(e)
            ))
        }
    )

    list(
        delta_hat = delta_hat,
        state_idx = state_idx
    )
}


# ============================================================================
# 4. .predict_compute_point --- Point predictions (internal)
# ============================================================================

#' Compute point predictions at the posterior mean
#'
#' Computes predictions using the posterior mean of the fixed-effect
#' parameter vector \eqn{\hat\theta}, optionally incorporating
#' posterior mean state random effects.
#'
#' The computation is vectorised via BLAS matrix multiplication:
#' \deqn{
#'   \hat\eta_{\mathrm{ext}} = X_{\mathrm{new}} \hat\alpha
#'     + X_{\mathrm{new}} \hat\delta_{\mathrm{ext}}[s(\cdot)],
#' }
#' and similarly for the intensive margin, followed by the logistic
#' inverse link.
#'
#' @param object An \code{hbb_fit} object.
#' @param X_new Numeric matrix N_new x P.
#' @param P Integer: number of covariates per margin.
#' @param delta List from \code{.predict_resolve_delta} or NULL.
#' @param type One of "response", "extensive", "intensive".
#' @return Numeric vector of length N_new.
#' @keywords internal
.predict_compute_point <- function(object, X_new, P, delta, type,
                                    n_trial = NULL) {

    D <- 2L * P + 1L
    param_names <- .make_stan_param_names(P)

    # Extract posterior means of fixed effects
    draws <- tryCatch(
        object$fit$draws(variables = param_names, format = "matrix"),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "Failed to extract posterior draws for prediction.",
                "i" = conditionMessage(e)
            ))
        }
    )
    if (!is.matrix(draws)) draws <- as.matrix(draws)
    theta_hat <- colMeans(draws)

    # NaN/Inf guard
    if (any(!is.finite(theta_hat))) {
        n_bad <- sum(!is.finite(theta_hat))
        cli::cli_warn(c(
            "!" = "{n_bad} non-finite posterior mean{?s} detected.",
            "i" = "Predictions may be unreliable."
        ))
    }

    alpha_hat <- theta_hat[seq_len(P)]
    beta_hat  <- theta_hat[(P + 1L):(2L * P)]

    # Linear predictors via BLAS
    eta_ext <- as.numeric(X_new %*% alpha_hat)
    eta_int <- as.numeric(X_new %*% beta_hat)

    # Add SVC random effects if available
    if (!is.null(delta)) {
        delta_hat <- delta$delta_hat
        state_idx <- delta$state_idx

        # delta_ext: columns 1:P; delta_int: columns (P+1):K
        delta_ext <- delta_hat[state_idx, seq_len(P), drop = FALSE]
        delta_int <- delta_hat[state_idx, (P + 1L):(2L * P), drop = FALSE]

        eta_ext <- eta_ext + rowSums(X_new * delta_ext)
        eta_int <- eta_int + rowSums(X_new * delta_int)
    }

    q_hat  <- plogis(eta_ext)
    mu_hat <- plogis(eta_int)

    # --- ZT correction: g(mu) = mu / (1 - p0) --------------------------------
    log_kappa_hat <- theta_hat[D]
    kappa_hat <- pmin(exp(log_kappa_hat), 1e15)

    if (!is.null(n_trial) && type != "extensive") {
        p0_hat   <- compute_p0(as.integer(n_trial), mu_hat, kappa_hat)
        mu_trunc <- mu_hat / pmax(1 - p0_hat, .Machine$double.eps)
    } else {
        mu_trunc <- mu_hat
    }

    pred <- switch(type,
        response  = q_hat * mu_trunc,
        extensive = q_hat,
        intensive = mu_trunc
    )

    # NaN/Inf guard and [0,1] clamping
    n_nonfinite <- sum(!is.finite(pred))
    if (n_nonfinite > 0L) {
        cli::cli_warn(c(
            "!" = "{n_nonfinite} non-finite prediction{?s} detected.",
            "i" = "Setting non-finite predictions to NA."
        ))
        pred[!is.finite(pred)] <- NA_real_
    }
    pred <- pmin(pmax(pred, 0), 1)

    pred
}


# ============================================================================
# 5. .predict_compute_interval --- Posterior credible intervals (internal)
# ============================================================================

#' Compute posterior credible intervals via draw propagation
#'
#' Propagates each MCMC draw through the inverse-link function to
#' obtain draw-level predictions, then computes empirical quantiles
#' for the credible interval.
#'
#' @section Algorithm:
#' \enumerate{
#'   \item Extract the M x D matrix of fixed-effect draws.
#'   \item If \code{ndraws} is specified and smaller than M, subsample
#'     M draws via deterministic evenly-spaced thinning.
#'   \item For each draw \eqn{m}, compute
#'     \eqn{q_i^{(m)} = \mathrm{logistic}(x_i' \alpha^{(m)})} and
#'     \eqn{\mu_i^{(m)} = \mathrm{logistic}(x_i' \beta^{(m)})},
#'     then form the prediction according to \code{type}.
#'   \item Compute \eqn{(\alpha/2,\; 1-\alpha/2)} quantiles across
#'     draws for each observation.
#' }
#'
#' @section Vectorisation:
#' The draw-level computation is vectorised over observations using
#' matrix multiplication: for each draw \eqn{m}, the N_new-vector
#' of predictions is computed via \eqn{X_{\mathrm{new}} \alpha^{(m)}}
#' etc., avoiding explicit loops over observations.
#'
#' @param object An \code{hbb_fit} object.
#' @param X_new Numeric matrix N_new x P.
#' @param P Integer: covariates per margin.
#' @param delta List from \code{.predict_resolve_delta} or NULL.
#' @param type One of "response", "extensive", "intensive".
#' @param level Numeric in (0,1).
#' @param ndraws Integer or NULL.
#' @return List with \code{lwr} (N_new-vector) and \code{upr}
#'   (N_new-vector).
#' @keywords internal
.predict_compute_interval <- function(object, X_new, P, delta, type,
                                       level, ndraws, n_trial = NULL) {

    D <- 2L * P + 1L
    N_new <- nrow(X_new)
    param_names <- .make_stan_param_names(P)

    # Extract full draw matrix (M x D)
    draws <- tryCatch(
        object$fit$draws(variables = param_names, format = "matrix"),
        error = function(e) {
            cli::cli_abort(c(
                "x" = "Failed to extract posterior draws for credible intervals.",
                "i" = conditionMessage(e)
            ))
        }
    )
    if (!is.matrix(draws)) draws <- as.matrix(draws)

    M <- nrow(draws)

    # M < 2 guard
    if (M < 2L) {
        cli::cli_abort(c(
            "x" = "Only {M} MCMC draw{?s} available; need at least 2 \\
                   for credible intervals.",
            "i" = "Check CmdStan fit for sampling failures.",
            "i" = "Use {.code interval = 'none'} for point predictions."
        ))
    }

    # Thin draws if ndraws specified
    if (!is.null(ndraws)) {
        if (ndraws > M) {
            cli::cli_warn(c(
                "!" = "{.arg ndraws} = {ndraws} exceeds available draws \\
                       M = {M}.",
                "i" = "Using all {M} draws."
            ))
            ndraws <- M
        }
        if (ndraws < M) {
            # Deterministic evenly-spaced thinning
            thin_idx <- unique(round(seq(1, M, length.out = ndraws)))
            draws <- draws[thin_idx, , drop = FALSE]
            M <- nrow(draws)
        }
    }

    # Pre-allocate prediction matrix: M x N_new
    pred_mat <- matrix(NA_real_, nrow = M, ncol = N_new)

    # Resolve delta observation-level matrices once
    delta_ext_obs <- NULL
    delta_int_obs <- NULL
    if (!is.null(delta)) {
        delta_hat <- delta$delta_hat
        state_idx <- delta$state_idx
        delta_ext_obs <- delta_hat[state_idx, seq_len(P), drop = FALSE]
        delta_int_obs <- delta_hat[state_idx, (P + 1L):(2L * P), drop = FALSE]
    }

    # Loop over draws, vectorised within each draw
    for (m in seq_len(M)) {
        alpha_m <- draws[m, seq_len(P)]
        beta_m  <- draws[m, (P + 1L):(2L * P)]
        log_kappa_m <- draws[m, D]
        kappa_m <- pmin(exp(log_kappa_m), 1e15)

        eta_ext_m <- as.numeric(X_new %*% alpha_m)
        eta_int_m <- as.numeric(X_new %*% beta_m)

        # Add SVC random effects (posterior mean, conditional interval)
        if (!is.null(delta_ext_obs)) {
            eta_ext_m <- eta_ext_m + rowSums(X_new * delta_ext_obs)
            eta_int_m <- eta_int_m + rowSums(X_new * delta_int_obs)
        }

        q_m  <- plogis(eta_ext_m)
        mu_m <- plogis(eta_int_m)

        # ZT correction
        if (!is.null(n_trial) && type != "extensive") {
            p0_m     <- compute_p0(as.integer(n_trial), mu_m, kappa_m)
            mu_trunc <- mu_m / pmax(1 - p0_m, .Machine$double.eps)
        } else {
            mu_trunc <- mu_m
        }

        pred_mat[m, ] <- switch(type,
            response  = q_m * mu_trunc,
            extensive = q_m,
            intensive = mu_trunc
        )
    }

    # NaN/Inf guard
    n_nonfinite <- sum(!is.finite(pred_mat))
    if (n_nonfinite > 0L) {
        pct <- round(100 * n_nonfinite / length(pred_mat), 1)
        cli::cli_warn(c(
            "!" = "{n_nonfinite} non-finite value{?s} ({pct}%) in \\
                   draw-level predictions.",
            "i" = "These will be treated as NA for quantile computation."
        ))
        pred_mat[!is.finite(pred_mat)] <- NA_real_
    }

    # Clamp to [0, 1]
    pred_mat <- pmin(pmax(pred_mat, 0), 1)

    # Compute quantiles across draws for each observation
    alpha_half <- (1 - level) / 2
    probs <- c(alpha_half, 1 - alpha_half)

    # apply over columns (observations)
    ci_mat <- apply(pred_mat, 2L, quantile, probs = probs, na.rm = TRUE)

    # Guard for N_new == 1 (ensure matrix)
    if (!is.matrix(ci_mat)) {
        ci_mat <- matrix(ci_mat, nrow = 2L, ncol = 1L)
    }

    list(
        lwr = ci_mat[1L, ],
        upr = ci_mat[2L, ]
    )
}
