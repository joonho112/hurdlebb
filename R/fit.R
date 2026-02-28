# ============================================================================
# fit.R --- Main Fitting Function for Hurdle Beta-Binomial Models
#
# The primary user-facing function for fitting two-part hurdle
# Beta-Binomial models via Stan/CmdStan. Orchestrates formula parsing,
# data preparation, Stan data assembly, model compilation, MCMC
# sampling, post-sampling diagnostic checks, and construction of the
# hbb_fit S3 return object.
#
# Contents:
#   1. hbb                   --- Main fitting function (exported)
#   2. .build_stan_input     --- Assemble Stan-ready data list (internal)
#   3. .validate_cmdstan_fit --- Post-sampling MCMC diagnostics (internal)
#   4. .model_name_from_type --- Map model_type to Stan model name (internal)
# ============================================================================


# ============================================================================
# 1. hbb --- Main fitting function (exported)
# ============================================================================

#' Fit a Hurdle Beta-Binomial Model
#'
#' Fits a two-part hurdle Beta-Binomial model using Hamiltonian Monte
#' Carlo (HMC) via CmdStan. The model comprises an extensive margin
#' (Bernoulli for whether a center serves infant/toddler children) and
#' an intensive margin (zero-truncated Beta-Binomial for the enrollment
#' share among servers). Both margins share the same design matrix with
#' separate coefficient vectors.
#'
#' @param formula A two-sided formula of the form
#'   `y | trials(n) ~ predictors`, or an object of class
#'   [hbb_formula]. See [hbb_formula()] for the full syntax including
#'   random effects and policy moderators.
#' @param data A data frame containing provider-level variables: the
#'   response, trials, fixed effects, and optionally the grouping
#'   variable.
#' @param weights Character string naming the column in `data`
#'   containing survey sampling weights, or `NULL` (default) for
#'   unweighted analysis.
#' @param stratum Character string naming the column in `data`
#'   containing sampling stratum identifiers, or `NULL`.
#' @param psu Character string naming the column in `data` containing
#'   primary sampling unit (PSU) identifiers, or `NULL`.
#' @param state_data A data frame of group-level (state-level) variables
#'   for policy moderators. Required if the formula includes
#'   `state_level()` terms. Must contain a column matching the grouping
#'   variable for merging.
#' @param prior A prior specification object of class [hbb_prior], or
#'   `NULL` to use [default_prior()]. See [hbb_prior()] for
#'   customisation.
#' @param chains Integer. Number of Markov chains. Default `4`.
#' @param iter_warmup Integer. Number of warmup iterations per chain.
#'   Default `1000`.
#' @param iter_sampling Integer. Number of post-warmup iterations per
#'   chain. Default `1000`.
#' @param adapt_delta Numeric in `(0, 1)`. Target average proposal
#'   acceptance probability during warmup adaptation. Higher values
#'   reduce divergent transitions at the cost of slower sampling.
#'   Default `0.95`.
#' @param max_treedepth Integer. Maximum tree depth for the NUTS
#'   sampler. Default `12`.
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#' @param refresh Integer or `NULL`. How often to print progress (every
#'   `refresh` iterations). Set to `0` to suppress progress output.
#'   `NULL` uses the CmdStan default.
#' @param cpp_options Named list of C++ compilation options passed to
#'   `cmdstanr::cmdstan_model()`. For example,
#'   `list(stan_threads = TRUE)` enables within-chain threading.
#'   Default `list()`.
#' @param ... Additional arguments passed to the CmdStanModel
#'   `$sample()` method.
#'
#' @return An S3 object of class `"hbb_fit"`. See [is.hbb_fit()] and
#'   [print.hbb_fit()] for methods. The object contains:
#' \describe{
#'   \item{`fit`}{The `CmdStanMCMC` object from cmdstanr.}
#'   \item{`stan_data`}{The list passed to Stan's `$sample()`.}
#'   \item{`hbb_data`}{The full `hbb_data` object for downstream
#'     methods.}
#'   \item{`formula`}{The `hbb_formula` object.}
#'   \item{`prior`}{The `hbb_prior` object used.}
#'   \item{`model_type`}{Character: one of `"base"`, `"weighted"`,
#'     `"svc"`, `"svc_weighted"`.}
#'   \item{`model_name`}{Character: the Stan model file name (e.g.,
#'     `"hbb_base"`).}
#'   \item{`call`}{The matched call.}
#'   \item{`elapsed`}{Numeric: wall-clock time in seconds.}
#' }
#'
#' @details
#' ## Model variants
#'
#' The Stan model selected depends on the formula and whether survey
#' weights are provided:
#' \tabular{lll}{
#'   **Model** \tab **Formula pattern** \tab **Weights** \cr
#'   `hbb_base`         \tab `y | trials(n) ~ x1 + x2` \tab No \cr
#'   `hbb_weighted`     \tab same \tab Yes \cr
#'   `hbb_svc`          \tab `... + (x1 | group)` \tab No \cr
#'   `hbb_svc_weighted` \tab `... + (x1 | group)` \tab Yes \cr
#' }
#'
#' Note that random-intercept-only models (`(1 | state_id)`) currently
#' use the `base` or `weighted` variants, as the random intercept is
#' absorbed into the SVC framework only when random slopes are present.
#'
#' ## Workflow
#'
#' `hbb()` performs the following steps:
#' 1. Parses the formula (if a raw formula rather than an `hbb_formula`).
#' 2. Validates the prior specification and MCMC arguments.
#' 3. Prepares Stan data via [prepare_stan_data()].
#' 4. Maps the prepared data and priors to the Stan data block.
#' 5. Compiles the appropriate Stan model (cached after first use).
#' 6. Runs MCMC sampling via CmdStan.
#' 7. Checks for MCMC pathologies (divergences, treedepth, E-BFMI).
#' 8. Returns an `hbb_fit` object for downstream analysis.
#'
#' ## Compilation
#'
#' Stan models are compiled on first use and cached persistently. See
#' [hbb_compile()] for details on the caching strategy and explicit
#' pre-compilation.
#'
#' @examples
#' \dontrun{
#' data(nsece_synth_small, package = "hurdlebb")
#'
#' # Base model (no weights, no random effects)
#' fit1 <- hbb(
#'   y | trials(n_trial) ~ poverty + urban,
#'   data = nsece_synth_small,
#'   chains = 2, iter_warmup = 500, iter_sampling = 500
#' )
#' print(fit1)
#'
#' # Weighted model
#' fit2 <- hbb(
#'   y | trials(n_trial) ~ poverty + urban,
#'   data = nsece_synth_small,
#'   weights = "weight",
#'   chains = 2, iter_warmup = 500, iter_sampling = 500
#' )
#'
#' # SVC model with policy moderators
#' data(nsece_state_policy, package = "hurdlebb")
#' fit3 <- hbb(
#'   y | trials(n_trial) ~ poverty + urban +
#'     (poverty + urban | state_id) +
#'     state_level(mr_pctile),
#'   data = nsece_synth_small,
#'   state_data = nsece_state_policy,
#'   weights = "weight",
#'   chains = 2, iter_warmup = 500, iter_sampling = 500
#' )
#' }
#'
#' @seealso [hbb_formula()], [prepare_stan_data()], [hbb_prior()],
#'   [hbb_compile()], [is.hbb_fit()]
#' @family fitting
#' @export
hbb <- function(formula,
                data,
                weights       = NULL,
                stratum       = NULL,
                psu           = NULL,
                state_data    = NULL,
                prior         = default_prior(),
                chains        = 4L,
                iter_warmup   = 1000L,
                iter_sampling = 1000L,
                adapt_delta   = 0.95,
                max_treedepth = 12L,
                seed          = NULL,
                refresh       = NULL,
                cpp_options   = list(),
                ...) {

    # -- Capture the call for reproducibility ----------------------------------
    cl <- match.call()

    # =========================================================================
    # INPUT VALIDATION
    # =========================================================================

    # -- Formula: accept hbb_formula or raw formula ----------------------------
    if (inherits(formula, "formula") && !inherits(formula, "hbb_formula")) {
        formula_obj <- hbb_formula(formula)
    } else if (inherits(formula, "hbb_formula")) {
        formula_obj <- formula
    } else {
        cli_abort(c(
            "Expected an {.cls hbb_formula} or a raw {.cls formula}.",
            "x" = "Got {.cls {class(formula)}}.",
            "i" = "Use brms-style syntax: {.code y | trials(n) ~ x1 + x2}."
        ))
    }

    # -- Data: must be a data frame -------------------------------------------
    assert_data_frame(data, min.rows = 1L, .var.name = "data")

    # -- Prior: NULL -> default, must be hbb_prior ----------------------------
    if (is.null(prior)) {
        prior <- default_prior()
    }
    if (!inherits(prior, "hbb_prior")) {
        cli_abort(c(
            "Expected a {.cls hbb_prior} object for {.arg prior}.",
            "x" = "Got {.cls {class(prior)}}.",
            "i" = "Use {.fun default_prior} or {.fun hbb_prior} to create one."
        ))
    }

    # -- MCMC control parameters ----------------------------------------------
    assert_count(chains, positive = TRUE, .var.name = "chains")
    assert_count(iter_warmup, positive = TRUE, .var.name = "iter_warmup")
    assert_count(iter_sampling, positive = TRUE, .var.name = "iter_sampling")
    assert_number(adapt_delta, lower = 0, upper = 1,
                  .var.name = "adapt_delta")
    assert_count(max_treedepth, positive = TRUE,
                 .var.name = "max_treedepth")

    if (!is.null(seed)) {
        assert_integerish(seed, len = 1L, lower = 1L,
                          .var.name = "seed")
        seed <- as.integer(seed)
    }
    if (!is.null(refresh)) {
        assert_count(refresh, .var.name = "refresh")
        refresh <- as.integer(refresh)
    }

    # -- Optional string arguments --------------------------------------------
    if (!is.null(weights)) {
        assert_string(weights, .var.name = "weights")
    }
    if (!is.null(stratum)) {
        assert_string(stratum, .var.name = "stratum")
    }
    if (!is.null(psu)) {
        assert_string(psu, .var.name = "psu")
    }

    # -- C++ options: must be a list ------------------------------------------
    checkmate::assert_list(cpp_options, .var.name = "cpp_options")

    # -- CmdStan availability -------------------------------------------------
    check_installed(
        "cmdstanr",
        reason = paste(
            "to fit Stan models.",
            "Install from {.url https://mc-stan.org/cmdstanr/}."
        )
    )

    # =========================================================================
    # DATA PREPARATION
    # =========================================================================

    cli_alert_info("Preparing data...")

    hbb_data <- prepare_stan_data(
        formula    = formula_obj,
        data       = data,
        weights    = weights,
        stratum    = stratum,
        psu        = psu,
        state_data = state_data,
        prior      = prior,
        center     = TRUE,
        scale      = TRUE
    )

    # =========================================================================
    # STAN DATA ASSEMBLY
    # =========================================================================

    stan_input <- .build_stan_input(hbb_data, prior)

    # =========================================================================
    # MODEL COMPILATION
    # =========================================================================

    model_type <- hbb_data$model_type
    model_name <- .model_name_from_type(model_type)

    cli_alert_info("Using Stan model {.val {model_name}}.")
    mod <- .hbb_model(model_name, cpp_options = cpp_options)

    # =========================================================================
    # MCMC SAMPLING
    # =========================================================================

    cli_alert_info(
        "Sampling: {chains} chain{?s}, {iter_warmup} warmup + \\
         {iter_sampling} sampling iterations each."
    )

    # Assemble sampler arguments
    n_cores <- parallel::detectCores(logical = FALSE)
    if (is.na(n_cores) || is.null(n_cores)) n_cores <- 1L

    sample_args <- list(
        data            = stan_input,
        chains          = as.integer(chains),
        parallel_chains = min(as.integer(chains), n_cores),
        iter_warmup     = as.integer(iter_warmup),
        iter_sampling   = as.integer(iter_sampling),
        adapt_delta     = adapt_delta,
        max_treedepth   = as.integer(max_treedepth)
    )

    # Add optional arguments only if non-NULL
    if (!is.null(seed))    sample_args$seed    <- seed
    if (!is.null(refresh)) sample_args$refresh <- refresh

    # Merge user-supplied ... arguments
    dots <- list(...)
    if (length(dots) > 0L) {
        sample_args <- c(sample_args, dots)
    }

    # Time the sampling
    t0 <- proc.time()

    cmdstan_fit <- tryCatch(
        do.call(mod$sample, sample_args),
        error = function(e) {
            cli_abort(
                c(
                    "Stan sampling failed for model {.val {model_name}}.",
                    "x" = conditionMessage(e),
                    "i" = "Check your data and model specification.",
                    "i" = "Try increasing {.arg adapt_delta} or \\
                           {.arg max_treedepth} if divergences occurred."
                ),
                parent = e
            )
        }
    )

    elapsed <- (proc.time() - t0)[["elapsed"]]

    # =========================================================================
    # POST-SAMPLING DIAGNOSTICS
    # =========================================================================

    .validate_cmdstan_fit(cmdstan_fit, adapt_delta, max_treedepth)

    # =========================================================================
    # CONSTRUCT RETURN OBJECT
    # =========================================================================

    cli_alert_success(
        "Model {.val {model_name}} fitted in {round(elapsed, 1)}s."
    )

    structure(
        list(
            fit        = cmdstan_fit,
            stan_data  = stan_input,
            hbb_data   = hbb_data,
            formula    = formula_obj,
            prior      = prior,
            model_type = model_type,
            model_name = model_name,
            call       = cl,
            elapsed    = elapsed
        ),
        class = "hbb_fit"
    )
}


# ============================================================================
# 2. .build_stan_input --- Assemble Stan-ready data list (internal)
# ============================================================================

#' Build the Stan input data list from an hbb_data object and prior
#'
#' Maps hbb_data fields and prior hyperparameters to the exact variable
#' names expected by the Stan model code. The set of variables included
#' depends on the model type. Includes critical pre-flight checks for
#' data integrity before passing to CmdStan.
#'
#' @param hbb_data An `hbb_data` object from [prepare_stan_data()].
#' @param prior An `hbb_prior` object.
#' @return A named list suitable for passing to `CmdStanModel$sample()`.
#' @noRd
.build_stan_input <- function(hbb_data, prior) {

    model_type  <- hbb_data$model_type
    is_weighted <- model_type %in% c("weighted", "svc_weighted")
    is_svc      <- model_type %in% c("svc", "svc_weighted")
    N           <- hbb_data$N
    P           <- hbb_data$P

    # -- Base fields (ALL models) ---------------------------------------------
    stan_input <- list(
        N                = N,
        P                = P,
        y                = hbb_data$y,
        n_trial          = hbb_data$n_trial,
        z                = hbb_data$z,
        X                = hbb_data$X,
        prior_alpha_mean = prior$alpha$mean,
        prior_alpha_sd   = prior$alpha$sd,
        prior_beta_mean  = prior$beta$mean,
        prior_beta_sd    = prior$beta$sd,
        prior_kappa_mean = prior$log_kappa$mean,
        prior_kappa_sd   = prior$log_kappa$sd
    )

    # -- Weighted models: add normalised weights ------------------------------
    if (is_weighted) {
        stan_input$w_tilde <- hbb_data$w_tilde
    }

    # -- SVC models: add group structure and hierarchical priors --------------
    if (is_svc) {
        stan_input$S              <- hbb_data$S
        stan_input$Q              <- hbb_data$Q
        stan_input$state          <- hbb_data$state
        stan_input$v_state        <- hbb_data$V   # R field V -> Stan v_state
        stan_input$prior_gamma_sd <- prior$gamma$sd
        stan_input$prior_tau_sd   <- prior$tau$sd
        stan_input$prior_lkj_eta  <- prior$lkj_eta
    }

    # =========================================================================
    # PRE-FLIGHT CHECKS (defense-in-depth)
    # =========================================================================

    # Check 1: Column 1 of X must be all 1s (intercept)
    intercept_err <- sum(abs(stan_input$X[, 1] - 1))
    if (intercept_err > .Machine$double.eps * N) {
        cli_abort(c(
            "Design matrix {.field X} column 1 must be all 1s (intercept).",
            "x" = "Sum of |X[,1] - 1| = {round(intercept_err, 6)}.",
            "i" = "This is an internal error in data preparation."
        ))
    }

    # Check 2: All n_trial >= 1
    if (any(stan_input$n_trial < 1L)) {
        cli_abort(c(
            "All {.field n_trial} values must be >= 1.",
            "x" = "Found {sum(stan_input$n_trial < 1L)} value{?s} below 1.",
            "i" = "This is an internal error in data preparation."
        ))
    }

    # Check 3: Weighted models -- w_tilde > 0 and sums to N
    if (is_weighted) {
        if (!all(stan_input$w_tilde > 0)) {
            n_bad <- sum(stan_input$w_tilde <= 0)
            cli_abort(c(
                "All normalized weights {.field w_tilde} must be strictly positive.",
                "x" = "Found {n_bad} non-positive weight{?s}."
            ))
        }
        w_sum_diff <- abs(sum(stan_input$w_tilde) - N)
        if (w_sum_diff > 0.01) {
            cli_abort(c(
                "Normalized weights {.field w_tilde} must sum to N = {N}.",
                "x" = "sum(w_tilde) = {round(sum(stan_input$w_tilde), 4)}, \\
                       |diff| = {round(w_sum_diff, 6)}.",
                "i" = "This is an internal error in weight normalization."
            ))
        }
    }

    # Check 4: SVC models -- state indices valid, v_state dimensions correct
    if (is_svc) {
        if (!all(stan_input$state >= 1L & stan_input$state <= stan_input$S)) {
            cli_abort(c(
                "State indices must be in 1..{stan_input$S}.",
                "x" = "Range: [{min(stan_input$state)}, {max(stan_input$state)}].",
                "i" = "This is an internal error in group indexing."
            ))
        }
        if (!is.matrix(stan_input$v_state)) {
            cli_abort(c(
                "{.field v_state} must be a matrix.",
                "x" = "Got {.cls {class(stan_input$v_state)}}."
            ))
        }
        if (nrow(stan_input$v_state) != stan_input$S ||
            ncol(stan_input$v_state) != stan_input$Q) {
            cli_abort(c(
                "{.field v_state} dimensions must be S x Q = \\
                 {stan_input$S} x {stan_input$Q}.",
                "x" = "Got {nrow(stan_input$v_state)} x \\
                       {ncol(stan_input$v_state)}."
            ))
        }
    }

    stan_input
}


# ============================================================================
# 3. .validate_cmdstan_fit --- Post-sampling MCMC diagnostics (internal)
# ============================================================================

#' Post-sampling diagnostic checks
#'
#' Extracts MCMC diagnostics from a CmdStanMCMC object and issues
#' warnings for divergences, max treedepth hits, and low E-BFMI.
#' All extractions are wrapped in tryCatch to be robust against
#' unusual CmdStan output formats.
#'
#' @param cmdstan_fit A `CmdStanMCMC` object from `$sample()`.
#' @param adapt_delta Numeric. Current adapt_delta value (for message).
#' @param max_treedepth Integer. Current max_treedepth (for message).
#' @return Invisibly returns `NULL`.
#' @noRd
.validate_cmdstan_fit <- function(cmdstan_fit,
                                  adapt_delta = NULL,
                                  max_treedepth = NULL) {

    diag <- tryCatch(
        cmdstan_fit$diagnostic_summary(quiet = TRUE),
        error = function(e) NULL
    )

    if (is.null(diag)) return(invisible(NULL))

    # -- Divergent transitions ------------------------------------------------
    n_div <- tryCatch(sum(diag$num_divergent), error = function(e) 0L)
    if (n_div > 0L) {
        msg <- "{n_div} divergent transition{?s} after warmup."
        hints <- c(
            "i" = "See {.url https://mc-stan.org/misc/warnings.html#divergent-transitions}."
        )
        if (!is.null(adapt_delta)) {
            hints <- c(
                "i" = "Consider increasing {.arg adapt_delta} \\
                       (currently {adapt_delta}).",
                hints
            )
        }
        cli_warn(c("!" = msg, hints))
    }

    # -- Max treedepth hits ---------------------------------------------------
    n_max_td <- tryCatch(sum(diag$num_max_treedepth), error = function(e) 0L)
    if (n_max_td > 0L) {
        msg <- "{n_max_td} transition{?s} hit max treedepth."
        hints <- c(
            "i" = "See {.url https://mc-stan.org/misc/warnings.html#maximum-treedepth}."
        )
        if (!is.null(max_treedepth)) {
            hints <- c(
                "i" = "Consider increasing {.arg max_treedepth} \\
                       (currently {max_treedepth}).",
                hints
            )
        }
        cli_warn(c("!" = msg, hints))
    }

    # -- Low E-BFMI -----------------------------------------------------------
    ebfmi <- tryCatch(diag$ebfmi, error = function(e) NULL)
    if (!is.null(ebfmi) && is.numeric(ebfmi)) {
        low_chains <- which(ebfmi < 0.2)
        if (length(low_chains) > 0L) {
            cli_warn(c(
                "!" = "Low E-BFMI for chain{?s} {low_chains} \\
                       (E-BFMI = {round(ebfmi[low_chains], 3)}).",
                "i" = "This may indicate a funnel or other pathology.",
                "i" = "Consider reparameterization or stronger priors."
            ))
        }
    }

    invisible(NULL)
}


# ============================================================================
# 4. .model_name_from_type --- Map model_type to Stan model name (internal)
# ============================================================================

#' Map model_type string to Stan model file name
#'
#' @param model_type Character scalar. One of `"base"`, `"weighted"`,
#'   `"svc"`, `"svc_weighted"`.
#' @return Character scalar. The Stan model file stem (e.g.,
#'   `"hbb_base"`).
#' @noRd
.model_name_from_type <- function(model_type) {
    valid_types <- c("base", "weighted", "svc", "svc_weighted")
    if (!model_type %in% valid_types) {
        cli_abort(c(
            "Unknown {.field model_type}: {.val {model_type}}.",
            "i" = "Expected one of: {.val {valid_types}}."
        ))
    }
    paste0("hbb_", model_type)
}
