# ============================================================================
# data-raw/precompute_vignettes.R
#
# Purpose:  Precompute all vignette results so they render without CmdStan.
#           Each vignette loads saved .rds files from inst/extdata/ instead
#           of fitting models at build time.
#
# Estimated runtime: 60-90 minutes depending on hardware (4 model fits,
#                    4 chains each, plus post-estimation).
#
# How to run:
#   cd hurdlebb/
#   Rscript data-raw/precompute_vignettes.R
#
# Prerequisites:
#   - CmdStan installed and configured (cmdstanr)
#   - hurdlebb package installed or loaded via devtools::load_all()
#   - Working directory should be the hurdlebb package root
#
# Outputs (saved to inst/extdata/):
#   vig01_results.rds     — V1 Getting Started: summary, coefs, fitted, draws
#   vig01_ppc.rds         — V1 Getting Started: posterior predictive checks
#   vig01_loo.rds         — V1 Getting Started: LOO-CV
#   vig02_results.rds     — V2 Survey Design: weighted model summary + Wald
#   vig02_sandwich.rds    — V2 Survey Design: sandwich variance object
#   vig02_wald_ci.rds     — V2 Survey Design: Wald CIs
#   vig02_loo_compare.rds — V2 Survey Design: LOO comparison (base vs weighted)
#   vig03_results.rds     — V3 SVC: summary, tau, correlation, random effects
#   vig03_loo.rds         — V3 SVC: LOO-CV
#   vig04_results.rds     — V4 Policy: summary, gamma table, tau comparison
#   vig04_loo.rds         — V4 Policy: LOO-CV
#   vig05_ame.rds         — V5 AME: hbb_ame object
#   vig05_ame_decomp.rds  — V5 AME: decomposition data frame
#   vig05_results.rds     — V5 AME: print text + poverty draws for density plots
# ============================================================================

library(hurdlebb)
library(posterior)
library(loo)

# -- Output directory ---------------------------------------------------------
out_dir <- file.path("inst", "extdata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -- Load datasets ------------------------------------------------------------
data(nsece_synth_small, package = "hurdlebb")
data(nsece_synth, package = "hurdlebb")
data(nsece_state_policy, package = "hurdlebb")

# -- Timing -------------------------------------------------------------------
t_start <- proc.time()
cat("=== Precompute vignettes started at", format(Sys.time()), "===\n\n")


# ============================================================================
# Section 1: Vignette 01 — Getting Started (base model on nsece_synth_small)
# ============================================================================
cat("--- Section 1: Vignette 01 (Getting Started) ---\n")

cat("  Fitting base model on nsece_synth_small ...\n")
fit_base <- tryCatch(
    hbb(
        y | trials(n_trial) ~ poverty + urban,
        data = nsece_synth_small,
        chains = 4,
        seed = 42
    ),
    error = function(e) {
        stop("Failed to fit base model: ", conditionMessage(e))
    }
)
cat("  Base model fit complete.\n")

cat("  Extracting vig01 results ...\n")
vig01_results <- tryCatch({
    # S3 summary
    summary_obj  <- summary(fit_base)
    summary_text <- capture.output(summary(fit_base))

    # Coefficients
    coef_both <- coef(fit_base, margin = "both")
    coef_ext  <- coef(fit_base, margin = "extensive")
    coef_int  <- coef(fit_base, margin = "intensive")

    # Observation count
    n_obs <- nobs(fit_base)

    # Fitted values
    fitted_response <- fitted(fit_base, type = "response")
    fitted_ext      <- fitted(fit_base, type = "extensive")
    fitted_int      <- fitted(fit_base, type = "intensive")

    # Residuals
    residuals_response <- residuals(fit_base, type = "response")
    residuals_pearson  <- residuals(fit_base, type = "pearson")

    # MCMC draws for trace/density plots (P = 3: intercept + poverty + urban)
    draw_vars <- c(
        "alpha[1]", "alpha[2]", "alpha[3]",
        "beta[1]",  "beta[2]",  "beta[3]",
        "log_kappa"
    )
    draws_df <- posterior::as_draws_df(fit_base$fit$draws(draw_vars))

    # Print output
    print_text <- capture.output(print(fit_base))

    list(
        summary_obj        = summary_obj,
        summary_text       = summary_text,
        coef_both          = coef_both,
        coef_ext           = coef_ext,
        coef_int           = coef_int,
        nobs               = n_obs,
        fitted_response    = fitted_response,
        fitted_ext         = fitted_ext,
        fitted_int         = fitted_int,
        residuals_response = residuals_response,
        residuals_pearson  = residuals_pearson,
        draws_df           = draws_df,
        print_text         = print_text
    )
}, error = function(e) {
    stop("Failed to extract vig01 results: ", conditionMessage(e))
})

saveRDS(vig01_results, file.path(out_dir, "vig01_results.rds"))
cat("  Saved vig01_results.rds\n")

cat("  Computing PPC ...\n")
vig01_ppc <- tryCatch(
    ppc(fit_base),
    error = function(e) {
        stop("Failed to compute PPC: ", conditionMessage(e))
    }
)
saveRDS(vig01_ppc, file.path(out_dir, "vig01_ppc.rds"))
cat("  Saved vig01_ppc.rds\n")

cat("  Computing LOO ...\n")
vig01_loo <- tryCatch(
    loo(fit_base),
    error = function(e) {
        stop("Failed to compute LOO: ", conditionMessage(e))
    }
)
saveRDS(vig01_loo, file.path(out_dir, "vig01_loo.rds"))
cat("  Saved vig01_loo.rds\n")

cat("  Section 1 complete.\n\n")


# ============================================================================
# Section 2: Vignette 02 — Survey Design (weighted model on nsece_synth_small)
# ============================================================================
cat("--- Section 2: Vignette 02 (Survey Design) ---\n")

cat("  Fitting weighted model on nsece_synth_small ...\n")
fit_weighted <- tryCatch(
    hbb(
        y | trials(n_trial) ~ poverty + urban,
        data = nsece_synth_small,
        weights = "weight",
        stratum = "stratum",
        psu = "psu",
        chains = 4,
        seed = 42
    ),
    error = function(e) {
        stop("Failed to fit weighted model: ", conditionMessage(e))
    }
)
cat("  Weighted model fit complete.\n")

cat("  Computing sandwich variance ...\n")
sandwich <- tryCatch(
    sandwich_variance(fit_weighted),
    error = function(e) {
        stop("Failed to compute sandwich variance: ", conditionMessage(e))
    }
)

cat("  Computing Wald CIs ...\n")
coef_weighted <- coef(fit_weighted)
wald_ci <- tryCatch(
    compute_wald_ci(coef_weighted, sandwich$V_sand, level = 0.95),
    error = function(e) {
        stop("Failed to compute Wald CIs: ", conditionMessage(e))
    }
)

cat("  Extracting vig02 results ...\n")
vig02_results <- tryCatch({
    summary_obj  <- summary(fit_weighted)
    summary_text <- capture.output(summary(fit_weighted))
    print_text   <- capture.output(print(fit_weighted))
    coef_both    <- coef(fit_weighted, margin = "both")

    # MCMC draws for the weighted model
    draw_vars <- c(
        "alpha[1]", "alpha[2]", "alpha[3]",
        "beta[1]",  "beta[2]",  "beta[3]",
        "log_kappa"
    )
    draws_df <- posterior::as_draws_df(fit_weighted$fit$draws(draw_vars))

    # Summary with Wald CIs (sandwich-adjusted)
    summary_wald_obj  <- summary(fit_weighted, sandwich = sandwich)
    summary_wald_text <- capture.output(summary(fit_weighted, sandwich = sandwich))

    list(
        summary_obj       = summary_obj,
        summary_text      = summary_text,
        print_text        = print_text,
        coef_both         = coef_both,
        draws_df          = draws_df,
        summary_wald_obj  = summary_wald_obj,
        summary_wald_text = summary_wald_text
    )
}, error = function(e) {
    stop("Failed to extract vig02 results: ", conditionMessage(e))
})

saveRDS(vig02_results, file.path(out_dir, "vig02_results.rds"))
cat("  Saved vig02_results.rds\n")

saveRDS(sandwich, file.path(out_dir, "vig02_sandwich.rds"))
cat("  Saved vig02_sandwich.rds\n")

saveRDS(wald_ci, file.path(out_dir, "vig02_wald_ci.rds"))
cat("  Saved vig02_wald_ci.rds\n")

cat("  Computing LOO comparison (base vs weighted) ...\n")
vig02_loo_compare <- tryCatch(
    hbb_loo_compare(base = fit_base, weighted = fit_weighted),
    error = function(e) {
        stop("Failed to compute LOO comparison: ", conditionMessage(e))
    }
)
saveRDS(vig02_loo_compare, file.path(out_dir, "vig02_loo_compare.rds"))
cat("  Saved vig02_loo_compare.rds\n")

cat("  Section 2 complete.\n\n")


# ============================================================================
# Section 3: Vignette 03 — State-Varying Coefficients (SVC on nsece_synth)
# ============================================================================
cat("--- Section 3: Vignette 03 (State-Varying Coefficients) ---\n")
cat("  NOTE: This fit may take 30-45 minutes.\n")

tryCatch({
    cat("  Fitting SVC model on nsece_synth (N =", nrow(nsece_synth), ") ...\n")
    cat("  Start time:", format(Sys.time(), "%H:%M:%S"), "\n")
    t3_start <- Sys.time()

    fit_svc <- hbb(
        y | trials(n_trial) ~ poverty + urban + (poverty + urban | state_id),
        data = nsece_synth,
        chains = 4,
        seed = 42,
        adapt_delta = 0.97
    )

    t3_elapsed <- difftime(Sys.time(), t3_start, units = "mins")
    cat("  SVC model fit complete in", round(as.numeric(t3_elapsed), 1), "min.\n")

    # -- Summary --
    cat("  Computing summary ...\n")
    summary_svc      <- summary(fit_svc)
    summary_text_svc <- capture.output(summary(fit_svc))
    print_text_svc   <- capture.output(print(fit_svc))
    coef_both_svc    <- coef(fit_svc, margin = "both")

    # -- Tau estimates (random effect SDs) --
    # P = 3, so 2P = 6 tau parameters
    cat("  Extracting tau ...\n")
    P  <- fit_svc$hbb_data$P
    K  <- 2L * P  # 6
    tau_vars    <- paste0("tau[", seq_len(K), "]")
    tau_summary <- fit_svc$fit$summary(variables = tau_vars)

    cov_names <- colnames(fit_svc$hbb_data$X)
    if (is.null(cov_names)) cov_names <- paste0("x", seq_len(P))
    tau_names <- c(
        paste0("alpha_", cov_names),
        paste0("beta_",  cov_names)
    )

    tau_table <- data.frame(
        parameter = tau_names,
        margin    = rep(c("extensive", "intensive"), each = P),
        covariate = rep(cov_names, 2),
        estimate  = tau_summary$mean,
        median    = tau_summary$median,
        sd        = tau_summary$sd,
        ci_lower  = tau_summary$q5,
        ci_upper  = tau_summary$q95,
        rhat      = tau_summary$rhat,
        ess_bulk  = tau_summary$ess_bulk,
        ess_tail  = tau_summary$ess_tail,
        stringsAsFactors = FALSE
    )

    # -- Correlation matrix from L_Omega --
    cat("  Computing correlation matrix from L_Omega ...\n")
    L_mat <- matrix(0, nrow = K, ncol = K)
    for (i in seq_len(K)) {
        for (j in seq_len(i)) {
            var_name <- paste0("L_Omega[", i, ",", j, "]")
            draws_ij <- fit_svc$fit$draws(variables = var_name, format = "matrix")
            L_mat[i, j] <- mean(draws_ij)
        }
    }
    corr_matrix <- L_mat %*% t(L_mat)
    rownames(corr_matrix) <- tau_names
    colnames(corr_matrix) <- tau_names

    # -- State random effects delta[s, k] --
    cat("  Extracting state random effects (delta) ...\n")
    S <- fit_svc$stan_data$S

    re_list <- vector("list", S * K)
    idx <- 0
    for (s in seq_len(S)) {
        for (k in seq_len(K)) {
            idx <- idx + 1
            var_name <- paste0("delta[", s, ",", k, "]")
            draws_sk <- fit_svc$fit$draws(variables = var_name, format = "matrix")
            re_list[[idx]] <- data.frame(
                state_id  = s,
                param_idx = k,
                param     = tau_names[k],
                margin    = if (k <= P) "extensive" else "intensive",
                covariate = cov_names[((k - 1L) %% P) + 1L],
                estimate  = mean(draws_sk),
                sd        = sd(draws_sk),
                ci_lower  = quantile(draws_sk, 0.025),
                ci_upper  = quantile(draws_sk, 0.975),
                stringsAsFactors = FALSE
            )
        }
        if (s %% 10 == 0) cat("    ... processed state", s, "of", S, "\n")
    }
    random_effects_df <- do.call(rbind, re_list)
    rownames(random_effects_df) <- NULL

    # -- MCMC draws for trace plots --
    cat("  Extracting MCMC draws for trace plots ...\n")
    trace_vars <- c(
        paste0("alpha[", seq_len(P), "]"),
        paste0("beta[",  seq_len(P), "]"),
        "log_kappa",
        paste0("tau[", seq_len(K), "]")
    )
    draws_df_svc <- posterior::as_draws_df(fit_svc$fit$draws(trace_vars))

    # -- Save results --
    vig03_results <- list(
        summary_obj       = summary_svc,
        summary_text      = summary_text_svc,
        print_text        = print_text_svc,
        coef_both         = coef_both_svc,
        tau_table         = tau_table,
        corr_matrix       = corr_matrix,
        random_effects_df = random_effects_df,
        draws_df          = draws_df_svc
    )

    saveRDS(vig03_results, file.path(out_dir, "vig03_results.rds"))
    cat("  Saved vig03_results.rds\n")

    cat("  Computing LOO ...\n")
    vig03_loo <- loo(fit_svc)
    saveRDS(vig03_loo, file.path(out_dir, "vig03_loo.rds"))
    cat("  Saved vig03_loo.rds\n")

    cat("  Section 3 complete.\n\n")

}, error = function(e) {
    cat("  ERROR in Section 3:", conditionMessage(e), "\n")
    cat("  Skipping remainder of Section 3.\n\n")
})


# ============================================================================
# Section 4: Vignette 04 — Policy Moderators (SVC + policy on nsece_synth)
# ============================================================================
cat("--- Section 4: Vignette 04 (Policy Moderators) ---\n")
cat("  NOTE: This fit may take 45-60 minutes.\n")

tryCatch({
    cat("  Fitting policy moderator model on nsece_synth ...\n")
    cat("  Start time:", format(Sys.time(), "%H:%M:%S"), "\n")
    t4_start <- Sys.time()

    fit_policy <- hbb(
        y | trials(n_trial) ~ poverty + urban +
            (poverty + urban | state_id) +
            state_level(mr_pctile + tiered_reim + it_addon),
        data = nsece_synth,
        state_data = nsece_state_policy,
        chains = 4,
        seed = 42,
        adapt_delta = 0.97
    )

    t4_elapsed <- difftime(Sys.time(), t4_start, units = "mins")
    cat("  Policy model fit complete in", round(as.numeric(t4_elapsed), 1), "min.\n")

    # -- Summary --
    cat("  Computing summary ...\n")
    summary_policy      <- summary(fit_policy)
    summary_text_policy <- capture.output(summary(fit_policy))
    print_text_policy   <- capture.output(print(fit_policy))
    coef_both_policy    <- coef(fit_policy, margin = "both")

    # -- Gamma table (policy moderator coefficients) --
    # Q = 4 (intercept + 3 policies), 2P = 6
    cat("  Extracting Gamma ...\n")
    P  <- fit_policy$hbb_data$P
    K  <- 2L * P  # 6
    Q  <- fit_policy$stan_data$Q  # 4

    # Gamma is matrix[K, Q] in Stan: rows = random effect dims, cols = policies
    gamma_vars    <- paste0("Gamma[", rep(seq_len(K), Q), ",",
                            rep(seq_len(Q), each = K), "]")
    gamma_summary <- fit_policy$fit$summary(variables = gamma_vars)

    policy_names <- c("intercept", "mr_pctile", "tiered_reim", "it_addon")
    cov_names <- colnames(fit_policy$hbb_data$X)
    if (is.null(cov_names)) cov_names <- paste0("x", seq_len(P))
    param_names <- c(
        paste0("alpha_", cov_names),
        paste0("beta_",  cov_names)
    )

    # Gamma[k,q]: iterate k (1..K) within each q (1..Q)
    gamma_table <- data.frame(
        variable  = gamma_summary$variable,
        policy    = rep(policy_names, each = K),
        parameter = rep(param_names, Q),
        margin    = rep(rep(c("extensive", "intensive"), each = P), Q),
        estimate  = gamma_summary$mean,
        median    = gamma_summary$median,
        sd        = gamma_summary$sd,
        ci_lower  = gamma_summary$q5,
        ci_upper  = gamma_summary$q95,
        rhat      = gamma_summary$rhat,
        ess_bulk  = gamma_summary$ess_bulk,
        stringsAsFactors = FALSE
    )

    # -- Tau estimates (with policy moderators) --
    cat("  Extracting tau ...\n")
    tau_vars_pol    <- paste0("tau[", seq_len(K), "]")
    tau_summary_pol <- fit_policy$fit$summary(variables = tau_vars_pol)

    tau_table_policy <- data.frame(
        parameter = param_names,
        margin    = rep(c("extensive", "intensive"), each = P),
        covariate = rep(cov_names, 2),
        estimate  = tau_summary_pol$mean,
        median    = tau_summary_pol$median,
        sd        = tau_summary_pol$sd,
        ci_lower  = tau_summary_pol$q5,
        ci_upper  = tau_summary_pol$q95,
        rhat      = tau_summary_pol$rhat,
        ess_bulk  = tau_summary_pol$ess_bulk,
        ess_tail  = tau_summary_pol$ess_tail,
        stringsAsFactors = FALSE
    )

    # -- Tau comparison: SVC (V3) vs Policy (V4) --
    cat("  Building tau comparison ...\n")
    tau_comparison <- NULL
    vig03_path <- file.path(out_dir, "vig03_results.rds")
    if (exists("vig03_results") || file.exists(vig03_path)) {
        if (!exists("vig03_results")) {
            vig03_results <- readRDS(vig03_path)
        }
        tau_v3 <- vig03_results$tau_table

        tau_comparison <- data.frame(
            parameter           = param_names,
            margin              = rep(c("extensive", "intensive"), each = P),
            covariate           = rep(cov_names, 2),
            tau_svc             = tau_v3$estimate,
            tau_svc_ci_lower    = tau_v3$ci_lower,
            tau_svc_ci_upper    = tau_v3$ci_upper,
            tau_policy          = tau_table_policy$estimate,
            tau_policy_ci_lower = tau_table_policy$ci_lower,
            tau_policy_ci_upper = tau_table_policy$ci_upper,
            reduction_pct       = round(
                100 * (tau_v3$estimate - tau_table_policy$estimate) /
                    tau_v3$estimate, 1
            ),
            stringsAsFactors = FALSE
        )
    } else {
        cat("  WARNING: vig03_results not available; tau comparison skipped.\n")
    }

    # -- Save results --
    vig04_results <- list(
        summary_obj    = summary_policy,
        summary_text   = summary_text_policy,
        print_text     = print_text_policy,
        coef_both      = coef_both_policy,
        gamma_table    = gamma_table,
        tau_table      = tau_table_policy,
        tau_comparison = tau_comparison
    )

    saveRDS(vig04_results, file.path(out_dir, "vig04_results.rds"))
    cat("  Saved vig04_results.rds\n")

    cat("  Computing LOO ...\n")
    vig04_loo <- loo(fit_policy)
    saveRDS(vig04_loo, file.path(out_dir, "vig04_loo.rds"))
    cat("  Saved vig04_loo.rds\n")

    cat("  Section 4 complete.\n\n")

}, error = function(e) {
    cat("  ERROR in Section 4:", conditionMessage(e), "\n")
    cat("  Skipping remainder of Section 4.\n\n")
})


# ============================================================================
# Section 5: Vignette 05 — Marginal Effects (reuses fit_svc from Section 3)
# ============================================================================
cat("--- Section 5: Vignette 05 (Marginal Effects) ---\n")

tryCatch({
    if (!exists("fit_svc")) {
        stop("fit_svc not found. Section 3 must complete successfully first.")
    }

    cat("  Computing average marginal effects from SVC model ...\n")
    t5_start <- Sys.time()

    # -- AME computation --
    ame_result <- ame(fit_svc, level = 0.95)

    t5_elapsed <- difftime(Sys.time(), t5_start, units = "mins")
    cat("  AME computation complete in", round(as.numeric(t5_elapsed), 1), "min.\n")

    # -- Decomposition --
    cat("  Computing AME decomposition ...\n")
    decomp <- ame_decomposition(ame_result)

    # -- Capture print output --
    ame_print_text <- capture.output(print(ame_result))

    # -- Extract poverty AME draws for density plots --
    # ext_ame_draws and int_ame_draws are M x P matrices
    # Column 1 = intercept, Column 2 = poverty, Column 3 = urban (for P=3)
    cat("  Extracting poverty AME draws for density plots ...\n")
    poverty_idx <- 2  # poverty is the 2nd covariate
    poverty_ext_draws <- ame_result$ext_ame_draws[, poverty_idx]
    poverty_int_draws <- ame_result$int_ame_draws[, poverty_idx]
    poverty_total_draws <- ame_result$total_ame_draws[, poverty_idx]

    # -- Save --
    saveRDS(ame_result, file.path(out_dir, "vig05_ame.rds"))
    cat("  Saved vig05_ame.rds\n")

    saveRDS(decomp, file.path(out_dir, "vig05_ame_decomp.rds"))
    cat("  Saved vig05_ame_decomp.rds\n")

    vig05_results <- list(
        ame_print_text      = ame_print_text,
        decomp_df           = decomp,
        poverty_ext_draws   = poverty_ext_draws,
        poverty_int_draws   = poverty_int_draws,
        poverty_total_draws = poverty_total_draws
    )
    saveRDS(vig05_results, file.path(out_dir, "vig05_results.rds"))
    cat("  Saved vig05_results.rds\n")

    cat("  Section 5 complete.\n\n")

}, error = function(e) {
    cat("  ERROR in Section 5:", conditionMessage(e), "\n")
    cat("  Skipping remainder of Section 5.\n\n")
})


# ============================================================================
# Closing: Summary
# ============================================================================
cat("================================================================\n")
cat("  Precomputation Complete\n")
cat("================================================================\n\n")

total_elapsed <- (proc.time() - t_start)[["elapsed"]]
cat("Total elapsed time:", round(total_elapsed / 60, 1), "minutes\n\n")

cat("Generated files:\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

output_files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(output_files) > 0) {
    file_info <- file.info(output_files)
    for (i in seq_along(output_files)) {
        fname <- basename(output_files[i])
        fsize <- file_info$size[i]
        if (fsize >= 1024 * 1024) {
            size_str <- sprintf("%.1f MB", fsize / (1024 * 1024))
        } else {
            size_str <- sprintf("%.1f KB", fsize / 1024)
        }
        cat(sprintf("  %-30s  %s\n", fname, size_str))
    }
    cat(paste(rep("-", 60), collapse = ""), "\n")
    cat(sprintf("  Total: %d files\n", length(output_files)))
} else {
    cat("  No .rds files found. Check error messages above.\n")
}

cat("\nDone!\n")
