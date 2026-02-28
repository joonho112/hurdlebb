# ============================================================================
# sandwich.R --- Cluster-Robust Sandwich Variance Estimator for hbb_fit
#
# Implements the survey-design-corrected sandwich variance estimator:
#
#   V_sand = H_obs^{-1} * J_cluster * H_obs^{-1}
#
# where H_obs is the block-diagonal observed Fisher information and
# J_cluster is the cluster-robust "meat" matrix computed from weighted
# score totals aggregated to the PSU level.
#
# Theory:
#   In hierarchical models the MCMC posterior covariance Sigma_MCMC
#   does NOT equal H_obs^{-1} because priors and random effects inflate
#   the marginal posterior variance of fixed effects (prior inflation
#   ratio 632--6552 in the NSECE application). The classic
#   Williams & Savitsky (2021) formula V = Sigma * J * Sigma therefore
#   produces astronomical DER. We instead use explicit H_obs as bread.
#
# Score vectors (D = 2P + 1 fixed effects):
#   s_i = [ score_ext[i, 1:P]  ;   d ell_i / d alpha
#           score_int[i, 1:P]  ;   d ell_i / d beta
#           score_kappa[i]     ]   d ell_i / d log_kappa
#
# Contents:
#   1. sandwich_variance       --- Top-level function (exported)
#   2. compute_score_matrix    --- Extract N x D score matrix (exported)
#   3. compute_H_obs           --- Block-diagonal observed information (exported)
#   4. compute_J_cluster       --- Cluster-robust meat matrix (exported)
#   5. compute_der             --- Design effect ratio (exported)
#   6. print.hbb_sandwich      --- S3 print method (exported)
#   7. .validate_sandwich_fit  --- Internal validation helper
#   8. .extract_posterior_means --- Safe posterior mean extraction
#   9. .extract_delta_means    --- Delta random effects extraction
#  10. .build_param_labels     --- Parameter label construction
#  11. .block_diag             --- Block-diagonal matrix assembly
#  12. .check_pd               --- PD check with diagnostics
#  13. .ridge_regularize       --- Ridge regularisation for near-singular H
# ============================================================================


# ============================================================================
# 1. sandwich_variance --- Top-level orchestrator (exported)
# ============================================================================

#' Cluster-Robust Sandwich Variance Estimator
#'
#' Computes the design-corrected sandwich variance-covariance matrix for
#' the fixed effects of a survey-weighted hurdle Beta-Binomial model.
#' The estimator uses the stratified cluster-robust formula
#'
#' \deqn{V_{\mathrm{sand}} = H_{\mathrm{obs}}^{-1}\, J_{\mathrm{cluster}}\,
#'       H_{\mathrm{obs}}^{-1}}
#'
#' where \eqn{H_{\mathrm{obs}}} is the block-diagonal observed Fisher
#' information (analytic for the extensive margin, empirical information
#' identity for the intensive margin plus dispersion), and
#' \eqn{J_{\mathrm{cluster}}} is the cluster-robust "meat" matrix that
#' accounts for within-PSU correlation and unequal survey weights.
#'
#' @section Why not Sigma_MCMC as bread:
#' In hierarchical models with state random effects, the MCMC posterior
#' covariance of fixed effects absorbs prior and random-effect variance,
#' leading to \eqn{\Sigma_{\mathrm{MCMC}} \gg H_{\mathrm{obs}}^{-1}}.
#' Using \eqn{\Sigma_{\mathrm{MCMC}}} as bread yields astronomical Design
#' Effect Ratios (3000--26000 in the NSECE application). The explicit
#' \eqn{H_{\mathrm{obs}}} resolves this.
#'
#' @section Fixed-effect parameter vector:
#' The \eqn{D = 2P + 1} dimensional parameter vector is ordered as
#' \deqn{\theta = (\alpha_1, \ldots, \alpha_P,\; \beta_1, \ldots,
#'       \beta_P,\; \log\kappa)}
#' where \eqn{\alpha} governs the extensive margin (Bernoulli), \eqn{\beta}
#' governs the intensive margin (zero-truncated Beta-Binomial), and
#' \eqn{\kappa} is the dispersion parameter.
#'
#' @section Design Effect Ratio:
#' The DER for each parameter is
#' \eqn{\mathrm{DER}_p = V_{\mathrm{sand}}[p,p] / H_{\mathrm{obs}}^{-1}[p,p]}.
#' Expected values are 1--5 for typical survey designs.
#'
#' @param fit An object of class `"hbb_fit"`, as returned by [hbb()].
#'   Must be a weighted model (`model_type %in% c("weighted",
#'   "svc_weighted")`), since score generated quantities are only
#'   produced by survey-weighted Stan models.
#'
#' @return An S3 object of class `"hbb_sandwich"` containing:
#' \describe{
#'   \item{`V_sand`}{Numeric \eqn{D \times D} sandwich variance matrix.}
#'   \item{`H_obs`}{Numeric \eqn{D \times D} observed Fisher information
#'     (block-diagonal: extensive \eqn{P \times P} and intensive+kappa
#'     \eqn{(P+1) \times (P+1)}).}
#'   \item{`H_obs_inv`}{Numeric \eqn{D \times D} inverse of `H_obs`,
#'     representing the data-only variance (no design correction).}
#'   \item{`J_cluster`}{Numeric \eqn{D \times D} cluster-robust meat
#'     matrix.}
#'   \item{`Sigma_MCMC`}{Numeric \eqn{D \times D} posterior covariance
#'     of fixed-effect draws (for comparison, not used as bread).}
#'   \item{`scores`}{Numeric \eqn{N \times D} posterior mean score
#'     matrix.}
#'   \item{`DER`}{Named numeric vector of length \eqn{D}. Design Effect
#'     Ratios: \eqn{\mathrm{DER}_p = V_{\mathrm{sand}}[p,p] /
#'     H_{\mathrm{obs}}^{-1}[p,p]}.}
#'   \item{`param_labels`}{Character vector of length \eqn{D} with
#'     human-readable parameter labels.}
#'   \item{`D`}{Integer. Total fixed-effect dimension.}
#'   \item{`N`}{Integer. Number of observations.}
#'   \item{`P`}{Integer. Number of covariates (including intercept).}
#'   \item{`model_type`}{Character. The model type from `fit`.}
#'   \item{`survey_info`}{List with survey design summaries: `n_strata`,
#'     `n_psu`, `df`, `n_singleton`.}
#'   \item{`matrix_diagnostics`}{List of eigenvalue/condition diagnostics
#'     for `H_obs`, `J_cluster`, `V_sand`, and `Sigma_MCMC`.}
#'   \item{`nearPD_applied`}{Logical. Whether nearPD correction was
#'     applied to V_sand.}
#'   \item{`call`}{The matched call.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- hbb(
#'   y | trials(n_trial) ~ poverty + urban,
#'   data = my_data, weights = "weight",
#'   stratum = "vstratum", psu = "vpsu"
#' )
#' sand <- sandwich_variance(fit)
#' print(sand)
#' compute_der(sand)
#' }
#'
#' @references
#' Williams, M. R. and Savitsky, T. D. (2021). Uncertainty estimation
#' for pseudo-Bayesian inference under complex sampling. *International
#' Statistical Review*, **89**(1), 72--107.
#'
#' @seealso [compute_score_matrix()], [compute_H_obs()],
#'   [compute_J_cluster()], [compute_der()]
#' @family sandwich
#' @export
sandwich_variance <- function(fit) {

    cl <- match.call()

    # -- Input validation --------------------------------------------
    .validate_sandwich_fit(fit)

    hbb_data   <- fit$hbb_data
    model_type <- fit$model_type
    N          <- hbb_data$N
    P          <- hbb_data$P
    D          <- 2L * P + 1L

    cli_alert_info(
        "Computing sandwich variance for {.val {model_type}} model \\
         (N={N}, D={D})."
    )

    # -- Step 1: Extract score matrix ------------------------------------------
    cli_alert_info("Step 1/5: Extracting posterior mean scores...")
    scores <- compute_score_matrix(fit)

    # -- Step 2: Compute H_obs ------------------------------------------------
    cli_alert_info("Step 2/5: Computing block-diagonal observed information H_obs...")
    H_result <- compute_H_obs(fit, scores = scores)
    H_obs     <- H_result$H_obs
    H_obs_inv <- H_result$H_obs_inv

    # -- Step 3: Compute J_cluster --------------------------------------------
    cli_alert_info("Step 3/5: Computing cluster-robust meat matrix J_cluster...")
    J_cluster <- compute_J_cluster(scores, hbb_data)

    # -- Step 4: Assemble V_sand = H_obs_inv %*% J_cluster %*% H_obs_inv -----
    cli_alert_info("Step 4/5: Assembling sandwich variance V_sand...")
    V_sand <- H_obs_inv %*% J_cluster %*% H_obs_inv

    # Force exact symmetry (numerical precision)
    V_sand <- (V_sand + t(V_sand)) / 2

    # -- PD check on V_sand with nearPD fallback --------------------
    nearPD_applied <- FALSE
    pd_check_V <- .check_pd(V_sand, label = "V_sand")

    if (!pd_check_V$is_pd) {
        cli_warn(c(
            "!" = "Sandwich variance {.field V_sand} is not positive definite.",
            "i" = "Minimum eigenvalue: {format(pd_check_V$min_eig, digits = 4)}.",
            "i" = "Applying {.fun Matrix::nearPD} correction."
        ))

        rlang::check_installed("Matrix",
            reason = "to apply nearPD correction to non-PD sandwich variance."
        )
        V_sand_corrected <- tryCatch(
            as.matrix(Matrix::nearPD(V_sand, corr = FALSE)$mat),
            error = function(e) {
                cli_abort(c(
                    "x" = "nearPD correction failed for {.field V_sand}.",
                    "i" = conditionMessage(e),
                    "i" = "The sandwich variance matrix may be severely ill-conditioned."
                ))
            }
        )
        V_sand <- V_sand_corrected
        nearPD_applied <- TRUE
    }

    # -- Step 5: Compute Sigma_MCMC for comparison ----------------------------
    cli_alert_info("Step 5/5: Computing Sigma_MCMC for diagnostic comparison...")
    Sigma_MCMC <- .compute_sigma_mcmc(fit$fit, P)

    # -- Build parameter labels -------------------------------------
    param_labels <- .build_param_labels(hbb_data)
    colnames(V_sand) <- rownames(V_sand) <- param_labels
    colnames(H_obs)  <- rownames(H_obs)  <- param_labels
    colnames(H_obs_inv) <- rownames(H_obs_inv) <- param_labels
    colnames(J_cluster) <- rownames(J_cluster) <- param_labels
    colnames(Sigma_MCMC) <- rownames(Sigma_MCMC) <- param_labels
    colnames(scores) <- param_labels

    # -- DER ------------------------------------------------------------------
    DER <- diag(V_sand) / diag(H_obs_inv)
    names(DER) <- param_labels

    # -- Matrix diagnostics -----------------------------------------
    eig_V <- eigen(V_sand, symmetric = TRUE, only.values = TRUE)$values
    eig_H <- eigen(H_obs, symmetric = TRUE, only.values = TRUE)$values
    eig_J <- eigen(J_cluster, symmetric = TRUE, only.values = TRUE)$values
    eig_S <- eigen(Sigma_MCMC, symmetric = TRUE, only.values = TRUE)$values

    matrix_diagnostics <- list(
        V_sand = list(
            min_eig = min(eig_V),
            max_eig = max(eig_V),
            cond    = max(eig_V) / max(min(eig_V), .Machine$double.eps),
            is_pd   = min(eig_V) > 0
        ),
        H_obs = list(
            min_eig = min(eig_H),
            max_eig = max(eig_H),
            cond    = max(eig_H) / max(min(eig_H), .Machine$double.eps),
            is_pd   = min(eig_H) > 0
        ),
        J_cluster = list(
            min_eig = min(eig_J),
            max_eig = max(eig_J),
            is_pd   = min(eig_J) > 0
        ),
        Sigma_MCMC = list(
            min_eig = min(eig_S),
            max_eig = max(eig_S),
            cond    = max(eig_S) / max(min(eig_S), .Machine$double.eps),
            is_pd   = min(eig_S) > 0
        )
    )

    # -- Survey design info ---------------------------------------------------
    strat_idx  <- hbb_data$stratum_idx
    psu_idx    <- hbb_data$psu_idx
    unique_sp  <- unique(paste(strat_idx, psu_idx, sep = ":"))
    n_strata   <- length(unique(strat_idx))
    n_psu      <- length(unique_sp)

    # Count PSUs per stratum for df
    psu_counts <- tapply(
        paste(strat_idx, psu_idx, sep = ":"),
        strat_idx,
        function(x) length(unique(x))
    )
    n_singleton <- sum(psu_counts == 1L)
    df_total    <- sum(pmax(psu_counts - 1L, 0L))

    survey_info <- list(
        n_strata    = n_strata,
        n_psu       = n_psu,
        df          = df_total,
        n_singleton = n_singleton
    )

    # -- Construct return object -----------------------------------------------
    cli_alert_success(
        "Sandwich variance computed. DER range: \\
         [{round(min(DER), 2)}, {round(max(DER), 2)}]."
    )

    structure(
        list(
            V_sand             = V_sand,
            H_obs              = H_obs,
            H_obs_inv          = H_obs_inv,
            J_cluster          = J_cluster,
            Sigma_MCMC         = Sigma_MCMC,
            scores             = scores,
            DER                = DER,
            param_labels       = param_labels,
            D                  = D,
            N                  = N,
            P                  = P,
            model_type         = model_type,
            survey_info        = survey_info,
            matrix_diagnostics = matrix_diagnostics,
            nearPD_applied     = nearPD_applied,
            call               = cl
        ),
        class = "hbb_sandwich"
    )
}


# ============================================================================
# 2. compute_score_matrix --- Extract N x D score matrix (exported)
# ============================================================================

#' Extract Posterior Mean Score Matrix from an hbb_fit
#'
#' Extracts the posterior means of the Stan-generated score variables
#' (`score_ext`, `score_int`, `score_kappa`) and assembles them into
#' an \eqn{N \times D} matrix, where \eqn{D = 2P + 1}:
#' \eqn{(\alpha_1, \ldots, \alpha_P, \beta_1, \ldots, \beta_P,
#' \log\kappa)}.
#'
#' @param fit An object of class `"hbb_fit"`. Must be a weighted
#'   model variant (`"weighted"` or `"svc_weighted"`).
#'
#' @return Numeric matrix of dimension \eqn{N \times D}. Columns are
#'   ordered as (score_ext columns 1 to P, score_int columns 1 to P,
#'   score_kappa).
#'
#' @details
#' Score variables are generated quantities in the weighted/svc_weighted
#' Stan models. CmdStanR names matrix variables as `var[1,1], var[1,2],
#' ..., var[1,P], var[2,1], ...` (first index varies slowest, row-major
#' ordering). Posterior means are computed across all MCMC draws.
#'
#' @seealso [sandwich_variance()]
#' @family sandwich
#' @export
compute_score_matrix <- function(fit) {

    if (!inherits(fit, "hbb_fit")) {
        cli_abort(c(
            "Expected an {.cls hbb_fit} object.",
            "x" = "Got {.cls {class(fit)}}."
        ))
    }

    model_type <- fit$model_type
    if (!model_type %in% c("weighted", "svc_weighted")) {
        cli_abort(c(
            "Score extraction requires a weighted model.",
            "x" = "Model type is {.val {model_type}}.",
            "i" = "Scores are only generated in {.val weighted} or \\
                   {.val svc_weighted} Stan models."
        ))
    }

    N <- fit$hbb_data$N
    P <- fit$hbb_data$P
    D <- 2L * P + 1L

    cmdstan_fit <- fit$fit

    # -- Extract score_ext: N x P (posterior means) ----------------------------
    score_ext_mat <- .extract_posterior_means(
        cmdstan_fit, "score_ext", N, P, label = "score_ext"
    )

    # -- Extract score_int: N x P (posterior means) ----------------------------
    score_int_mat <- .extract_posterior_means(
        cmdstan_fit, "score_int", N, P, label = "score_int"
    )

    # -- Extract score_kappa: length N (posterior means) -----------------------
    score_kappa_mat <- .extract_posterior_means(
        cmdstan_fit, "score_kappa", N, 1L, label = "score_kappa"
    )
    score_kappa_vec <- as.numeric(score_kappa_mat)

    # -- Assemble N x D score matrix -------------------------------------------
    S_mat <- cbind(score_ext_mat, score_int_mat, score_kappa_vec)

    if (nrow(S_mat) != N || ncol(S_mat) != D) {
        cli_abort(c(
            "Score matrix has unexpected dimensions.",
            "x" = "Expected {N} x {D}, got {nrow(S_mat)} x {ncol(S_mat)}.",
            "i" = "Check that the Stan model generated quantities match \\
                   the data dimensions."
        ))
    }

    S_mat
}


# ============================================================================
# 3. compute_H_obs --- Block-diagonal observed information (exported)
# ============================================================================

#' Compute the Block-Diagonal Observed Information Matrix
#'
#' Computes the observed information matrix \eqn{H_{\mathrm{obs}}} for
#' the fixed effects of a hurdle Beta-Binomial model. Uses a
#' block-diagonal structure:
#'
#' * **Extensive margin** (alpha block, \eqn{P \times P}): analytic
#'   Fisher information from logistic regression.
#' * **Intensive + kappa block** (\eqn{(P+1) \times (P+1)}): empirical
#'   information identity using posterior mean scores.
#'
#' @section Extensive margin (H_ext):
#' For the base weighted model:
#' \deqn{H_{\mathrm{ext}} = \sum_{i=1}^N \tilde{w}_i\, q_i (1-q_i)\,
#'   X_i X_i^\top}
#' where \eqn{q_i = \mathrm{logit}^{-1}(X_i^\top \hat{\alpha})}.
#'
#' For SVC models, the linear predictor includes state random effects:
#' \eqn{q_i = \mathrm{logit}^{-1}(X_i^\top \hat{\alpha} + X_i^\top
#' \hat{\delta}^{\mathrm{ext}}_{s[i]})}.
#'
#' @section Intensive margin (H_int):
#' The zero-truncated Beta-Binomial Hessian is analytically complex,
#' so we use the information identity:
#' \deqn{H_{\mathrm{int}} = \sum_{i=1}^N \tilde{w}_i\,
#'   s_i^{\mathrm{int}} (s_i^{\mathrm{int}})^\top}
#' where \eqn{s_i^{\mathrm{int}} = (s_{\beta,i}, s_{\kappa,i})^\top}
#' is the \eqn{(P+1)}-dimensional intensive score at the posterior mean.
#'
#' @param fit An object of class `"hbb_fit"`.
#' @param scores Optional \eqn{N \times D} score matrix (as from
#'   [compute_score_matrix()]). If `NULL`, computed internally.
#'
#' @return A list with components:
#' \describe{
#'   \item{`H_obs`}{Numeric matrix (\eqn{D \times D}). Block-diagonal
#'     observed information.}
#'   \item{`H_obs_inv`}{Numeric matrix (\eqn{D \times D}). Inverse of
#'     H_obs (the bread matrix).}
#'   \item{`H_ext`}{Numeric matrix (\eqn{P \times P}). Extensive-margin
#'     Fisher information block.}
#'   \item{`H_int`}{Numeric matrix (\eqn{(P+1) \times (P+1)}).
#'     Intensive + kappa block from information identity.}
#'   \item{`ridge_applied`}{Logical. Whether ridge regularisation was
#'     needed.}
#' }
#'
#' @seealso [sandwich_variance()], [compute_score_matrix()]
#' @family sandwich
#' @export
compute_H_obs <- function(fit, scores = NULL) {

    if (!inherits(fit, "hbb_fit")) {
        cli_abort(c(
            "Expected an {.cls hbb_fit} object.",
            "x" = "Got {.cls {class(fit)}}."
        ))
    }

    hbb_data   <- fit$hbb_data
    model_type <- fit$model_type
    N          <- hbb_data$N
    P          <- hbb_data$P
    D          <- 2L * P + 1L
    P_int      <- P + 1L   # beta[1:P] + log_kappa

    # Get scores if not provided
    if (is.null(scores)) {
        scores <- compute_score_matrix(fit)
    }

    checkmate::assert_matrix(scores, nrows = N, ncols = D,
                             .var.name = "scores")

    # -- Extract design matrix and weights ------------------------------------
    X_mat   <- hbb_data$X           # N x P
    w_tilde <- hbb_data$w_tilde     # length N

    if (is.null(w_tilde)) {
        cli_abort(c(
            "Survey weights {.field w_tilde} are required for H_obs.",
            "i" = "This model was fit without weights."
        ))
    }

    # =====================================================================
    # H_ext: Analytic logistic Fisher information (P x P)
    # H_ext = X' diag(w * q * (1-q)) X = crossprod(sqrt(W_ext) * X)
    # =====================================================================

    # Extract posterior means of alpha (P-vector)
    alpha_draws <- fit$fit$draws(
        variables = paste0("alpha[", seq_len(P), "]"),
        format = "draws_matrix"
    )
    alpha_hat <- colMeans(alpha_draws)

    # Compute linear predictor for extensive margin
    eta_ext <- as.numeric(X_mat %*% alpha_hat)  # N-vector

    # For SVC models: add state random effects delta_ext (vectorised)
    is_svc <- model_type %in% c("svc", "svc_weighted")
    if (is_svc) {
        S <- hbb_data$S
        K <- 2L * P
        delta_hat <- .extract_delta_means(fit$fit, S, K)

        # delta_ext: first P columns of delta
        delta_ext_per_obs <- delta_hat[hbb_data$state, seq_len(P), drop = FALSE]
        eta_ext <- eta_ext + rowSums(X_mat * delta_ext_per_obs)
    }

    q_vec <- plogis(eta_ext)                    # N-vector

    # Diagnostic: track extreme q_i values
    n_extreme_q <- sum(q_vec < 1e-10 | q_vec > (1 - 1e-10))
    if (n_extreme_q > 0L) {
        cli_warn(c(
            "!" = "{n_extreme_q} observation{?s} with extreme q_i \\
                   (near 0 or 1).",
            "i" = "H_ext may be near-singular due to perfect separation.",
            "i" = "Ridge regularisation will be applied if needed."
        ))
    }

    # Vectorised computation
    W_ext <- w_tilde * q_vec * (1 - q_vec)      # N-vector
    X_scaled <- sqrt(W_ext) * X_mat              # N x P
    H_ext <- crossprod(X_scaled)                  # P x P

    # =====================================================================
    # H_int: Empirical information identity ((P+1) x (P+1))
    # H_int = sum_i w_i * s_int_full_i s_int_full_i'
    #       = crossprod(sqrt(w) * S_int_full)
    # =====================================================================

    # scores columns: [alpha_1..P, beta_1..P, log_kappa]
    # Intensive block = columns (P+1):(2P+1)
    S_int_full <- scores[, (P + 1L):D, drop = FALSE]   # N x (P+1)

    S_int_scaled <- sqrt(w_tilde) * S_int_full           # N x (P+1)
    H_int <- crossprod(S_int_scaled)                      # (P+1) x (P+1)

    # =====================================================================
    # Assemble block-diagonal H_obs (D x D)
    # =====================================================================

    H_obs <- .block_diag(H_ext, H_int)

    # Force symmetry
    H_obs <- (H_obs + t(H_obs)) / 2

    # =====================================================================
    # PD check and ridge regularisation
    # =====================================================================

    ridge_applied <- FALSE
    pd_check <- .check_pd(H_obs, label = "H_obs")

    if (!pd_check$is_pd) {
        cli_warn(c(
            "!" = "Observed information {.field H_obs} is not positive definite.",
            "i" = "Minimum eigenvalue: {format(pd_check$min_eig, digits = 4)}.",
            "i" = "Applying ridge regularisation."
        ))
        H_obs <- .ridge_regularize(H_obs)
        ridge_applied <- TRUE
    }

    # Condition number warning
    if (pd_check$cond > 1e10) {
        cli_warn(c(
            "!" = "H_obs has high condition number: \\
                   {format(pd_check$cond, digits = 4, scientific = TRUE)}.",
            "i" = "Sandwich variance estimates may be numerically unstable."
        ))
    }

    # =====================================================================
    # Invert H_obs
    # =====================================================================

    H_obs_inv <- tryCatch(
        solve(H_obs),
        error = function(e) {
            cli_abort(c(
                "x" = "Failed to invert {.field H_obs}.",
                "i" = conditionMessage(e),
                "i" = "H_obs may be singular. Check for collinear predictors \\
                       or extreme q_i values."
            ))
        }
    )

    # Force symmetry on inverse
    H_obs_inv <- (H_obs_inv + t(H_obs_inv)) / 2

    list(
        H_obs         = H_obs,
        H_obs_inv     = H_obs_inv,
        H_ext         = H_ext,
        H_int         = H_int,
        ridge_applied = ridge_applied
    )
}


# ============================================================================
# 4. compute_J_cluster --- Cluster-robust meat matrix (exported)
# ============================================================================

#' Compute the Cluster-Robust Meat Matrix
#'
#' Aggregates weighted scores to the stratum--PSU level and computes
#' the Taylor linearisation variance estimator with degrees-of-freedom
#' (Bessel) correction at the stratum level.
#'
#' @param scores Numeric matrix of dimension \eqn{N \times D}. The
#'   posterior mean score matrix from [compute_score_matrix()].
#' @param hbb_data An S3 object of class `"hbb_data"` (from
#'   [prepare_stan_data()]) or a list containing the fields
#'   `stratum_idx`, `psu_idx`, `w_tilde`, and `N`.
#'
#' @return Numeric matrix of dimension \eqn{D \times D}. The
#'   cluster-robust meat matrix (positive semi-definite).
#'
#' @details
#' ## Algorithm
#'
#' For each stratum \eqn{h} with \eqn{C_h} PSUs:
#' 1. Compute weighted score totals per PSU:
#'    \eqn{s_{hc} = \sum_{i \in \mathrm{PSU}(h,c)} \tilde{w}_i s_i}
#' 2. Compute stratum mean: \eqn{\bar{s}_h = C_h^{-1} \sum_c s_{hc}}
#' 3. Center: \eqn{\delta_{hc} = s_{hc} - \bar{s}_h}
#' 4. Accumulate with FPC: \eqn{J_h = \frac{C_h}{C_h - 1}
#'    \sum_c \delta_{hc} \delta_{hc}^\top}
#'
#' Uses `rowsum()` for efficient PSU-level aggregation.
#' Singleton strata (\eqn{C_h = 1}) are skipped with a warning.
#'
#' @seealso [sandwich_variance()]
#' @family sandwich
#' @export
compute_J_cluster <- function(scores, hbb_data) {

    # -- Input validation ------------------------------------------------------
    N <- nrow(scores)
    D <- ncol(scores)

    stratum_idx <- hbb_data$stratum_idx
    psu_idx     <- hbb_data$psu_idx
    w_tilde     <- hbb_data$w_tilde

    if (is.null(stratum_idx)) {
        cli_abort(c(
            "{.field stratum_idx} is required for cluster-robust variance.",
            "i" = "Provide {.arg stratum} when calling {.fun hbb}."
        ))
    }
    if (is.null(psu_idx)) {
        cli_abort(c(
            "{.field psu_idx} is required for cluster-robust variance.",
            "i" = "Provide {.arg psu} when calling {.fun hbb}."
        ))
    }

    if (length(stratum_idx) != N || length(psu_idx) != N) {
        cli_abort(c(
            "Dimension mismatch between scores ({N} rows) and survey indices.",
            "x" = "stratum_idx has {length(stratum_idx)}, psu_idx has \\
                   {length(psu_idx)}."
        ))
    }

    # Default to uniform weights if NULL
    if (is.null(w_tilde)) {
        w_tilde <- rep(1, N)
    }

    if (length(w_tilde) != N) {
        cli_abort(
            "Length of {.field w_tilde} ({length(w_tilde)}) does not \\
             match score rows ({N})."
        )
    }

    # -- Weighted scores: W_scores[i, d] = w_tilde[i] * scores[i, d] ---------
    W_scores <- w_tilde * scores  # N x D (recycling w_tilde column-wise)

    # -- Aggregate weighted scores to PSU level via rowsum ----------
    psu_group <- interaction(stratum_idx, psu_idx, drop = TRUE)
    psu_totals <- rowsum(W_scores, group = psu_group, reorder = FALSE)
    # psu_totals: n_psu_unique x D matrix

    # Map each PSU to its stratum
    psu_group_levels <- levels(psu_group)
    first_idx <- match(psu_group_levels, as.character(psu_group))
    psu_strata <- stratum_idx[first_idx]

    # -- Count PSUs per stratum -----------------------------------------------
    unique_strata <- sort(unique(psu_strata))
    C_by_stratum <- tabulate(match(psu_strata, unique_strata))

    # -- Warn about singletons --------------------------------------
    n_singleton <- sum(C_by_stratum == 1L)
    if (n_singleton > 0L) {
        singleton_ids <- unique_strata[C_by_stratum == 1L]
        id_display <- paste(head(singleton_ids, 10), collapse = ", ")
        if (n_singleton > 10L) {
            id_display <- paste0(id_display, " ... and ",
                                 n_singleton - 10L, " more")
        }
        cli_warn(c(
            "!" = paste0(n_singleton,
                         " singleton stratum/strata (C_h = 1) detected."),
            "i" = "Singleton strata are skipped because within-stratum \\
                   variance is inestimable.",
            "i" = paste0("Stratum IDs: ", id_display, ".")
        ))
    }

    # -- Accumulate J_cluster stratum by stratum ------------------------------
    J_cluster <- matrix(0, D, D)

    for (h_idx in seq_along(unique_strata)) {
        C_h <- C_by_stratum[h_idx]

        # Skip singleton strata (division by zero)
        if (C_h < 2L) next

        # Indices of PSU-level totals belonging to this stratum
        in_stratum <- which(psu_strata == unique_strata[h_idx])

        # C_h x D sub-matrix of PSU score totals
        s_hc <- psu_totals[in_stratum, , drop = FALSE]

        # Center within stratum
        s_bar_h <- colMeans(s_hc)
        delta_hc <- sweep(s_hc, 2L, s_bar_h, FUN = "-")

        # Accumulate: J += (C_h / (C_h - 1)) * crossprod(delta_hc)
        # Note: C_h/(C_h-1) is the Bessel (df) correction, NOT an FPC
        bessel <- C_h / (C_h - 1L)
        J_cluster <- J_cluster + bessel * crossprod(delta_hc)
    }

    # Force symmetry
    J_cluster <- (J_cluster + t(J_cluster)) / 2

    J_cluster
}


# ============================================================================
# 5. compute_der --- Design Effect Ratio (exported)
# ============================================================================

#' Compute Design Effect Ratios from a Sandwich Variance Object
#'
#' Returns the Design Effect Ratio (DER) for each fixed-effect
#' parameter:
#' \deqn{\mathrm{DER}_p = V_{\mathrm{sand}}[p,p] /
#'       H_{\mathrm{obs}}^{-1}[p,p]}
#'
#' DER > 1 means the survey design inflates variance beyond what a
#' naive (non-survey-corrected) analysis would produce. Expected
#' values are typically 1--5 for stratified cluster designs (Kish
#' DEFF ~ 3--4).
#'
#' @param sandwich An object of class `"hbb_sandwich"` from
#'   [sandwich_variance()].
#'
#' @return Named numeric vector of length \eqn{D}.
#'
#' @details
#' The following warnings are issued:
#' * DER < 0.5: survey design appears to reduce variance substantially,
#'   which is unusual and may signal a problem.
#' * DER > 20: very large design effect, possibly indicating model
#'   misspecification or extreme weight variability.
#'
#' @seealso [sandwich_variance()]
#' @family sandwich
#' @export
compute_der <- function(sandwich) {

    if (!inherits(sandwich, "hbb_sandwich")) {
        cli_abort(c(
            "Expected an {.cls hbb_sandwich} object.",
            "x" = "Got {.cls {class(sandwich)}}."
        ))
    }

    V_diag     <- diag(sandwich$V_sand)
    H_inv_diag <- diag(sandwich$H_obs_inv)

    DER <- V_diag / H_inv_diag
    names(DER) <- sandwich$param_labels

    # -- Warn on extreme DER values ---------------------------------
    low_idx  <- which(DER < 0.5)
    high_idx <- which(DER > 20)

    if (length(low_idx) > 0L) {
        low_labels <- sandwich$param_labels[low_idx]
        low_vals   <- round(DER[low_idx], 3)
        cli_warn(c(
            "!" = "Unusually low DER (< 0.5) for {length(low_idx)} \\
                   parameter{?s}.",
            "i" = "Parameter{?s}: {.val {low_labels}}.",
            "i" = "DER value{?s}: {low_vals}."
        ))
    }

    if (length(high_idx) > 0L) {
        high_labels <- sandwich$param_labels[high_idx]
        high_vals   <- round(DER[high_idx], 3)
        cli_warn(c(
            "!" = "Unusually high DER (> 20) for {length(high_idx)} \\
                   parameter{?s}.",
            "i" = "Parameter{?s}: {.val {high_labels}}.",
            "i" = "DER value{?s}: {high_vals}."
        ))
    }

    DER
}


# ============================================================================
# 6. print.hbb_sandwich --- S3 print method (exported)
# ============================================================================

#' Print Method for hbb_sandwich Objects
#'
#' Displays a compact summary of the sandwich variance estimator,
#' including the design effect ratios (DER), survey design information,
#' and matrix condition diagnostics. Wrapped in `tryCatch` so that
#' printing never fails even on malformed objects.
#'
#' @param x An object of class `"hbb_sandwich"`.
#' @param ... Currently unused.
#'
#' @return Invisibly returns `x`.
#'
#' @seealso [sandwich_variance()], [compute_der()]
#' @family sandwich
#' @method print hbb_sandwich
#' @export
print.hbb_sandwich <- function(x, ...) {

    tryCatch({

        cat("Hurdle Beta-Binomial Sandwich Variance Estimator\n")
        cat("=================================================\n\n")

        # -- Model info -------------------------------------------------------
        cat("  Model type       :", x$model_type %||% "(unknown)", "\n")
        cat("  Observations (N) :", x$N %||% NA, "\n")
        cat("  Parameters   (D) :", x$D %||% NA, "\n")

        # -- Survey design info -----------------------------------------------
        si <- x$survey_info
        if (!is.null(si)) {
            cat("\nSurvey design:\n")
            cat("  Strata           :", si$n_strata %||% NA, "\n")
            cat("  PSUs             :", si$n_psu %||% NA, "\n")
            cat("  Degrees of freedom:", si$df %||% NA, "\n")
            if (!is.null(si$n_singleton) && si$n_singleton > 0L) {
                cat("  WARNING:", si$n_singleton,
                    "singleton strata skipped\n")
            }
        }

        # -- DER table --------------------------------------------------------
        der <- x$DER
        labels <- x$param_labels

        if (!is.null(der) && !is.null(labels)) {
            D <- length(der)
            sand_sd <- tryCatch(
                sqrt(diag(x$V_sand)),
                error = function(e) rep(NA_real_, D)
            )
            h_inv_sd <- tryCatch(
                sqrt(diag(x$H_obs_inv)),
                error = function(e) rep(NA_real_, D)
            )
            sigma_sd <- tryCatch(
                sqrt(diag(x$Sigma_MCMC)),
                error = function(e) rep(NA_real_, D)
            )

            cat("\nDesign Effect Ratios (DER = V_sand / H_obs_inv):\n")
            cat(sprintf("  %-28s %10s %10s %10s %8s %8s\n",
                        "Parameter", "H_inv SD",
                        "Sand. SD", "Sigma SD", "DER", "Flag"))
            cat(sprintf("  %s\n", paste(rep("-", 80), collapse = "")))

            for (d in seq_len(D)) {
                flag <- ""
                if (!is.na(der[d])) {
                    if (der[d] < 0.5) flag <- "[LOW]"
                    else if (der[d] > 20) flag <- "[HIGH]"
                    else if (der[d] > 5) flag <- "[*]"
                }
                cat(sprintf("  %-28s %10.6f %10.6f %10.6f %8.3f %8s\n",
                            labels[d],
                            h_inv_sd[d],
                            sand_sd[d],
                            sigma_sd[d],
                            der[d],
                            flag))
            }

            cat(sprintf(
                "\n  DER summary: min=%.3f, median=%.3f, mean=%.3f, max=%.3f\n",
                min(der), median(der), mean(der), max(der)
            ))
        }

        # -- Matrix diagnostics -----------------------------------------------
        md <- x$matrix_diagnostics
        if (!is.null(md)) {
            cat("\nMatrix diagnostics:\n")

            .print_diag <- function(name, entry) {
                if (!is.null(entry)) {
                    pd_label <- if (isTRUE(entry$is_pd)) "YES"
                                else if (!is.null(entry$min_eig) &&
                                         entry$min_eig >= -1e-10) "PSD"
                                else "NO"
                    cond_str <- if (!is.null(entry$cond)) {
                        sprintf(", cond=%.2e", entry$cond)
                    } else ""
                    cat(sprintf("  %-12s: PD=%s%s\n", name, pd_label,
                                cond_str))
                }
            }

            .print_diag("V_sand", md$V_sand)
            .print_diag("H_obs", md$H_obs)
            .print_diag("J_cluster", md$J_cluster)
            .print_diag("Sigma_MCMC", md$Sigma_MCMC)
        }

        if (isTRUE(x$nearPD_applied)) {
            cat("\n  NOTE: Matrix::nearPD correction was applied to V_sand.\n")
        }

        cat("\n")

    }, error = function(e) {
        cat("Hurdle Beta-Binomial Sandwich Variance Estimator\n")
        cat("  (print failed: ", conditionMessage(e), ")\n")
    })

    invisible(x)
}


# ============================================================================
# 7. .validate_sandwich_fit --- Internal validation
# ============================================================================

#' Validate that an hbb_fit is suitable for sandwich variance estimation
#'
#' Checks: (1) class, (2) model_type is weighted, (3) fit$fit is not
#' NULL, (4) hbb_data has stratum_idx and psu_idx.
#'
#' @param fit An R object to validate.
#' @return Invisibly returns `TRUE`.
#' @noRd
.validate_sandwich_fit <- function(fit) {

    # Class check
    if (!inherits(fit, "hbb_fit")) {
        cli_abort(c(
            "Expected an {.cls hbb_fit} object.",
            "x" = "Got {.cls {class(fit)}}.",
            "i" = "Use {.fun hbb} to fit a hurdle Beta-Binomial model first."
        ))
    }

    # Model type check
    model_type <- fit$model_type
    if (!model_type %in% c("weighted", "svc_weighted")) {
        cli_abort(c(
            "Sandwich variance requires a survey-weighted model.",
            "x" = "Model type is {.val {model_type}}.",
            "i" = "Only {.val weighted} and {.val svc_weighted} Stan models \\
                   include generated-quantity score variables.",
            "i" = "Refit with {.code weights = \"your_weight_column\"} in \\
                   {.fun hbb}."
        ))
    }

    # CmdStanMCMC object check
    if (is.null(fit$fit)) {
        cli_abort(c(
            "The {.field fit} slot is NULL.",
            "i" = "The CmdStanMCMC object is required for score extraction."
        ))
    }

    # hbb_data check
    hbb_data <- fit$hbb_data
    if (is.null(hbb_data)) {
        cli_abort(c(
            "The {.field hbb_data} slot is missing.",
            "i" = "This should not happen with a properly constructed \\
                   {.cls hbb_fit}."
        ))
    }

    # Survey design check
    if (is.null(hbb_data$stratum_idx)) {
        cli_abort(c(
            "{.field stratum_idx} not found in {.field hbb_data}.",
            "i" = "Provide the {.arg stratum} argument when calling {.fun hbb}."
        ))
    }

    if (is.null(hbb_data$psu_idx)) {
        cli_abort(c(
            "{.field psu_idx} not found in {.field hbb_data}.",
            "i" = "Provide the {.arg psu} argument when calling {.fun hbb}."
        ))
    }

    invisible(TRUE)
}


# ============================================================================
# 8. .extract_posterior_means --- Safe posterior mean extraction
# ============================================================================

#' Safely extract posterior means for a named block of Stan GQ variables
#'
#' Extracts draws for a variable block (e.g., `"score_ext"`) from a
#' CmdStanMCMC object and computes column means. Returns an `n_rows`
#' x `n_cols` matrix.
#'
#' @param cmdstan_fit A CmdStanMCMC object.
#' @param var_name Character. Base variable name (e.g., "score_ext").
#' @param n_rows Integer. Expected number of rows (N).
#' @param n_cols Integer. Expected number of columns (P or 1).
#' @param label Character. Label for error messages.
#' @return Numeric matrix of dimension `n_rows` x `n_cols`.
#' @noRd
.extract_posterior_means <- function(cmdstan_fit, var_name,
                                     n_rows, n_cols, label = var_name) {

    draws <- tryCatch(
        cmdstan_fit$draws(variables = var_name, format = "draws_matrix"),
        error = function(e) {
            cli_abort(c(
                "Failed to extract {.val {label}} draws from the Stan fit.",
                "x" = conditionMessage(e),
                "i" = "Ensure the Stan model includes {.val {var_name}} \\
                       as a generated quantity."
            ))
        }
    )

    if (is.null(draws) || ncol(draws) == 0L) {
        cli_abort(c(
            "No draws found for {.val {label}}.",
            "i" = "Expected {n_rows * n_cols} elements ({n_rows} x {n_cols})."
        ))
    }

    # Column means across MCMC draws
    col_means <- colMeans(draws)

    expected_len <- n_rows * n_cols
    actual_len   <- length(col_means)

    if (actual_len != expected_len) {
        cli_abort(c(
            "Dimension mismatch for {.val {label}} posterior means.",
            "x" = "Expected {expected_len} elements ({n_rows} x {n_cols}), \\
                   got {actual_len}."
        ))
    }

    # Reshape into matrix
    # CmdStanR names matrix variables as: var[1,1], var[1,2], ..., var[1,P],
    # var[2,1], ... (first index varies slowest = row-major ordering).
    # Therefore byrow = TRUE is correct.
    matrix(col_means, nrow = n_rows, ncol = n_cols, byrow = TRUE)
}


# ============================================================================
# 9. .extract_delta_means --- Delta random effects extraction
# ============================================================================

#' Extract posterior means of delta (state random effects)
#'
#' @param cmdstan_fit A CmdStanMCMC object.
#' @param S Integer. Number of states.
#' @param K Integer. Total RE dimension (2*P).
#' @return An S x K numeric matrix of posterior means.
#' @noRd
.extract_delta_means <- function(cmdstan_fit, S, K) {
    # Build parameter names: delta[1,1], delta[1,2], ..., delta[S,K]
    delta_names <- character(S * K)
    idx <- 1L
    for (s in seq_len(S)) {
        for (k in seq_len(K)) {
            delta_names[idx] <- sprintf("delta[%d,%d]", s, k)
            idx <- idx + 1L
        }
    }
    delta_summary <- cmdstan_fit$summary(variables = delta_names)
    delta_vec <- delta_summary$mean
    matrix(delta_vec, nrow = S, ncol = K, byrow = TRUE)
}


# ============================================================================
# 10. .build_param_labels --- Parameter label construction
# ============================================================================

#' Build human-readable parameter labels from hbb_data
#'
#' Constructs labels of the form `alpha_name`, `beta_name`, `log_kappa`
#' using the column names of the design matrix X.
#'
#' @param hbb_data An hbb_data object.
#' @return Character vector of length D = 2P + 1.
#' @noRd
.build_param_labels <- function(hbb_data) {

    P <- hbb_data$P
    X <- hbb_data$X

    # Use column names of X if available, else generic names
    if (!is.null(colnames(X))) {
        covariate_names <- colnames(X)
    } else {
        covariate_names <- paste0("x", seq_len(P))
    }

    c(
        paste0("alpha_", covariate_names),
        paste0("beta_", covariate_names),
        "log_kappa"
    )
}


# ============================================================================
# 11. .block_diag --- Block-diagonal matrix assembly
# ============================================================================

#' Assemble two square matrices into a block-diagonal matrix
#'
#' @param A Numeric matrix (first block).
#' @param B Numeric matrix (second block).
#' @return Numeric matrix of dimension (nrow(A)+nrow(B)) x (ncol(A)+ncol(B)).
#' @noRd
.block_diag <- function(A, B) {
    nA <- nrow(A)
    nB <- nrow(B)
    D  <- nA + nB

    out <- matrix(0, nrow = D, ncol = D)
    out[seq_len(nA), seq_len(nA)] <- A
    out[(nA + 1L):D, (nA + 1L):D] <- B
    out
}


# ============================================================================
# 12. .check_pd --- Positive-definiteness check
# ============================================================================

#' Check if a symmetric matrix is positive definite
#'
#' @param mat Numeric matrix.
#' @param label Character label for messages.
#' @return A list with `is_pd` (logical), `min_eig`, `max_eig`, `cond`.
#' @noRd
.check_pd <- function(mat, label = "matrix") {

    eig_vals <- tryCatch(
        eigen(mat, symmetric = TRUE, only.values = TRUE)$values,
        error = function(e) {
            cli_warn(c(
                "!" = "Eigenvalue computation failed for {.field {label}}.",
                "i" = conditionMessage(e)
            ))
            return(NA_real_)
        }
    )

    if (anyNA(eig_vals)) {
        return(list(
            is_pd   = FALSE,
            min_eig = NA_real_,
            max_eig = NA_real_,
            cond    = Inf
        ))
    }

    min_eig <- min(eig_vals)
    max_eig <- max(eig_vals)
    cond_num <- max_eig / max(min_eig, .Machine$double.eps)

    list(
        is_pd   = min_eig > 0,
        min_eig = min_eig,
        max_eig = max_eig,
        cond    = cond_num
    )
}


# ============================================================================
# 13. .ridge_regularize --- Ridge regularisation
# ============================================================================

#' Apply ridge regularisation to make a matrix positive definite
#'
#' Adds `max(|min_eig| * 2, eps * ||mat||_F) * I` to the matrix.
#'
#' @param mat Symmetric numeric matrix.
#' @param eps Numeric. Minimum ridge relative to Frobenius norm.
#'   Default `1e-6`.
#' @return The regularised matrix (guaranteed positive definite).
#' @noRd
.ridge_regularize <- function(mat, eps = 1e-6) {

    D <- nrow(mat)
    eig_vals <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
    min_eig <- min(eig_vals)

    # Ridge magnitude: max of twice the negative eigenvalue and a
    # fraction of the Frobenius norm
    frob_norm <- sqrt(sum(mat^2))
    ridge <- max(abs(min_eig) * 2, eps * frob_norm)

    mat_reg <- mat + ridge * diag(D)

    # Verify PD
    new_min <- min(eigen(mat_reg, symmetric = TRUE, only.values = TRUE)$values)
    if (new_min <= 0) {
        cli_warn(c(
            "!" = "Ridge regularisation did not fully restore positive \\
                   definiteness.",
            "i" = "Minimum eigenvalue after ridge: \\
                   {format(new_min, digits = 4)}."
        ))
    }

    mat_reg
}


# ============================================================================
# 14. .compute_sigma_mcmc --- Posterior covariance of fixed effects
# ============================================================================

#' Compute MCMC posterior covariance of fixed-effect draws
#'
#' Extracts the posterior draws for alpha(1..P), beta(1..P), and
#' log_kappa, and computes their sample covariance matrix. This is
#' included for diagnostic comparison only --- it should NOT be used
#' as bread in the sandwich.
#'
#' @param cmdstan_fit A CmdStanMCMC object.
#' @param P Integer. Number of covariates.
#' @return A D x D covariance matrix.
#' @noRd
.compute_sigma_mcmc <- function(cmdstan_fit, P) {

    D <- 2L * P + 1L
    param_names <- c(
        paste0("alpha[", seq_len(P), "]"),
        paste0("beta[", seq_len(P), "]"),
        "log_kappa"
    )

    draws_mat <- tryCatch(
        cmdstan_fit$draws(variables = param_names, format = "draws_matrix"),
        error = function(e) {
            cli_warn(c(
                "!" = "Failed to extract fixed-effect draws for Sigma_MCMC.",
                "i" = conditionMessage(e)
            ))
            return(NULL)
        }
    )

    if (is.null(draws_mat)) {
        return(diag(D) * NA_real_)
    }

    cov(draws_mat)
}
