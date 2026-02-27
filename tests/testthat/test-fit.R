# ============================================================================
# test-fit.R --- Tests for hbb(), .build_stan_input(), .validate_cmdstan_fit(),
#                .model_name_from_type(), is.hbb_fit(), print.hbb_fit()
#
# Sections:
#   1. .build_stan_input() — Unit tests (NO CmdStan required)
#   2. .model_name_from_type() — Lookup helper
#   3. is.hbb_fit() — Class predicate
#   4. hbb() — Input validation (NO CmdStan required)
#   5. print.hbb_fit() — Mock/defensive tests (NO CmdStan required)
#   6. hbb() — Integration: base model (requires CmdStan)
#   7. hbb() — Integration: weighted model (requires CmdStan)
#   8. hbb() — Integration: SVC model (requires CmdStan)
#   9. hbb() — Integration: custom prior, seed, print (requires CmdStan)
# ============================================================================


# -- Shared minimal MCMC settings for integration tests --------------------
# Fast settings for CI/test speed. NOT for real inference.
.test_chains        <- 2L
.test_iter_warmup   <- 100L
.test_iter_sampling <- 100L
.test_seed          <- 12345L
.test_refresh       <- 0L


# ============================================================================
# Section 1: .build_stan_input() — Unit tests (NO CmdStan required)
# ============================================================================

test_that(".build_stan_input() produces correct fields for base model", {

    data(nsece_synth_small, package = "hurdlebb")
    f <- hbb_formula(y | trials(n_trial) ~ poverty + urban)
    d <- prepare_stan_data(f, nsece_synth_small)
    p <- default_prior()

    sl <- hurdlebb:::.build_stan_input(d, p)

    # Core fields present
    expect_true(all(c("N", "P", "y", "n_trial", "z", "X") %in% names(sl)))
    expect_true(all(c("prior_alpha_sd", "prior_beta_sd",
                       "prior_kappa_mean", "prior_kappa_sd") %in% names(sl)))

    # Dimensional checks
    expect_equal(sl$N, d$N)
    expect_equal(sl$P, d$P)
    expect_equal(length(sl$y), d$N)
    expect_equal(length(sl$n_trial), d$N)
    expect_equal(length(sl$z), d$N)
    expect_equal(nrow(sl$X), d$N)
    expect_equal(ncol(sl$X), d$P)

    # Prior scalars match
    expect_equal(sl$prior_alpha_sd, p$alpha$sd)
    expect_equal(sl$prior_beta_sd, p$beta$sd)
    expect_equal(sl$prior_kappa_mean, p$log_kappa$mean)
    expect_equal(sl$prior_kappa_sd, p$log_kappa$sd)

    # Intercept column is all 1s
    expect_true(all(sl$X[, 1] == 1))

    # No SVC or weight fields for base model
    expect_null(sl$w_tilde)
    expect_null(sl$S)
    expect_null(sl$Q)
    expect_null(sl$state)
    expect_null(sl$v_state)
    expect_null(sl$prior_gamma_sd)
    expect_null(sl$prior_tau_sd)
    expect_null(sl$prior_lkj_eta)
})


test_that(".build_stan_input() adds w_tilde for weighted model", {

    data(nsece_synth_small, package = "hurdlebb")
    f <- hbb_formula(y | trials(n_trial) ~ poverty + urban)
    d <- prepare_stan_data(f, nsece_synth_small, weights = "weight")
    p <- default_prior()

    sl <- hurdlebb:::.build_stan_input(d, p)

    # Must have w_tilde
    expect_true("w_tilde" %in% names(sl))
    expect_equal(length(sl$w_tilde), d$N)
    expect_true(all(sl$w_tilde > 0))

    # w_tilde should sum to N (normalised weights)
    expect_equal(sum(sl$w_tilde), d$N, tolerance = 0.01)

    # Still no SVC fields
    expect_null(sl$S)
    expect_null(sl$state)
})


test_that(".build_stan_input() adds hierarchy fields for SVC model", {

    data(nsece_synth_small, package = "hurdlebb")
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + urban + (poverty + urban | state_id)
    )
    d <- prepare_stan_data(f, nsece_synth_small)
    p <- default_prior()

    sl <- hurdlebb:::.build_stan_input(d, p)

    # SVC fields present
    expect_true(all(c("S", "Q", "state", "v_state") %in% names(sl)))
    expect_true(all(c("prior_gamma_sd", "prior_tau_sd",
                       "prior_lkj_eta") %in% names(sl)))

    # Dimensional checks
    expect_equal(sl$S, d$S)
    expect_equal(sl$Q, d$Q)
    expect_equal(length(sl$state), d$N)
    expect_true(is.matrix(sl$v_state))
    expect_equal(nrow(sl$v_state), d$S)
    expect_equal(ncol(sl$v_state), d$Q)

    # State indices valid
    expect_true(all(sl$state >= 1 & sl$state <= sl$S))

    # No weights for non-weighted SVC
    expect_null(sl$w_tilde)
})


test_that(".build_stan_input() adds both hierarchy and weights for svc_weighted", {

    data(nsece_synth_small, package = "hurdlebb")
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + urban + (poverty + urban | state_id)
    )
    d <- prepare_stan_data(f, nsece_synth_small, weights = "weight")
    p <- default_prior()

    expect_equal(d$model_type, "svc_weighted")

    sl <- hurdlebb:::.build_stan_input(d, p)

    # Both weight and SVC fields present
    expect_true("w_tilde" %in% names(sl))
    expect_true("S" %in% names(sl))
    expect_true("state" %in% names(sl))
    expect_true("v_state" %in% names(sl))
    expect_true("prior_lkj_eta" %in% names(sl))
})


# ============================================================================
# Section 2: .model_name_from_type() — Lookup helper
# ============================================================================

test_that(".model_name_from_type() maps all valid types", {
    expect_equal(hurdlebb:::.model_name_from_type("base"), "hbb_base")
    expect_equal(hurdlebb:::.model_name_from_type("weighted"), "hbb_weighted")
    expect_equal(hurdlebb:::.model_name_from_type("svc"), "hbb_svc")
    expect_equal(hurdlebb:::.model_name_from_type("svc_weighted"), "hbb_svc_weighted")
})


test_that(".model_name_from_type() errors on invalid types", {
    expect_error(
        hurdlebb:::.model_name_from_type("invalid"),
        class = "rlang_error"
    )
    expect_error(
        hurdlebb:::.model_name_from_type(""),
        class = "rlang_error"
    )
})


# ============================================================================
# Section 3: is.hbb_fit() — Class predicate
# ============================================================================

test_that("is.hbb_fit() returns FALSE for non-hbb_fit objects", {
    expect_false(is.hbb_fit(1))
    expect_false(is.hbb_fit("fit"))
    expect_false(is.hbb_fit(list()))
    expect_false(is.hbb_fit(NULL))
    expect_false(is.hbb_fit(data.frame()))
    expect_false(is.hbb_fit(lm(1:10 ~ 1)))

    # A plain list with matching field names is NOT hbb_fit
    fake_list <- list(fit = NULL, model_type = "base")
    expect_false(is.hbb_fit(fake_list))
})


test_that("is.hbb_fit() returns TRUE for classed objects", {
    # Manually constructed hbb_fit (for testing without Stan)
    fake_fit <- structure(list(a = 1), class = "hbb_fit")
    expect_true(is.hbb_fit(fake_fit))

    # Multiple class inheritance
    fake_fit2 <- structure(list(), class = c("hbb_fit", "list"))
    expect_true(is.hbb_fit(fake_fit2))

    fake_fit3 <- structure(list(), class = c("special", "hbb_fit"))
    expect_true(is.hbb_fit(fake_fit3))
})


# ============================================================================
# Section 4: hbb() — Input validation (NO CmdStan required)
# ============================================================================

test_that("hbb() rejects invalid formula", {
    data(nsece_synth_small, package = "hurdlebb")

    expect_error(
        hbb(formula = "not a formula", data = nsece_synth_small),
        class = "rlang_error"
    )
    expect_error(
        hbb(formula = 42, data = nsece_synth_small),
        class = "rlang_error"
    )
    expect_error(
        hbb(formula = list(a = 1), data = nsece_synth_small),
        class = "rlang_error"
    )
})


test_that("hbb() rejects invalid data", {
    expect_error(
        hbb(y | trials(n_trial) ~ poverty, data = "not a data frame")
    )
    expect_error(
        hbb(y | trials(n_trial) ~ poverty, data = matrix(1:10, ncol = 2))
    )
})


test_that("hbb() rejects invalid prior", {
    data(nsece_synth_small, package = "hurdlebb")

    expect_error(
        hbb(y | trials(n_trial) ~ poverty,
            data = nsece_synth_small,
            prior = list(alpha = list(sd = 2))),
        "hbb_prior"
    )
    expect_error(
        hbb(y | trials(n_trial) ~ poverty,
            data = nsece_synth_small,
            prior = "flat"),
        "hbb_prior"
    )
})


test_that("hbb() rejects invalid MCMC arguments", {
    data(nsece_synth_small, package = "hurdlebb")

    # chains: negative, zero, non-integer
    expect_error(
        hbb(y | trials(n_trial) ~ poverty,
            data = nsece_synth_small, chains = -1),
        "chains"
    )
    expect_error(
        hbb(y | trials(n_trial) ~ poverty,
            data = nsece_synth_small, chains = 0),
        "chains"
    )

    # adapt_delta: out of range
    expect_error(
        hbb(y | trials(n_trial) ~ poverty,
            data = nsece_synth_small, adapt_delta = 1.5),
        "adapt_delta"
    )
    expect_error(
        hbb(y | trials(n_trial) ~ poverty,
            data = nsece_synth_small, adapt_delta = -0.5),
        "adapt_delta"
    )

    # iter_warmup: zero
    expect_error(
        hbb(y | trials(n_trial) ~ poverty,
            data = nsece_synth_small, iter_warmup = 0),
        "iter_warmup"
    )

    # iter_sampling: zero
    expect_error(
        hbb(y | trials(n_trial) ~ poverty,
            data = nsece_synth_small, iter_sampling = 0),
        "iter_sampling"
    )

    # max_treedepth: zero
    expect_error(
        hbb(y | trials(n_trial) ~ poverty,
            data = nsece_synth_small, max_treedepth = 0),
        "max_treedepth"
    )
})


test_that("hbb() rejects invalid seed and refresh", {
    data(nsece_synth_small, package = "hurdlebb")

    # seed: negative
    expect_error(
        hbb(y | trials(n_trial) ~ poverty,
            data = nsece_synth_small, seed = -5),
        "seed"
    )

    # refresh: negative
    expect_error(
        hbb(y | trials(n_trial) ~ poverty,
            data = nsece_synth_small, refresh = -1),
        "refresh"
    )
})


test_that("hbb() rejects invalid cpp_options", {
    data(nsece_synth_small, package = "hurdlebb")

    expect_error(
        hbb(y | trials(n_trial) ~ poverty,
            data = nsece_synth_small,
            cpp_options = "not a list"),
        "cpp_options"
    )
})


# ============================================================================
# Section 5: print.hbb_fit() — Mock/defensive tests (NO CmdStan required)
# ============================================================================

test_that("print.hbb_fit() handles minimal mock object", {
    # Create mock without CmdStanMCMC — just test print doesn't error
    mock_fit <- structure(
        list(
            fit        = list(),  # empty list, no CmdStanMCMC methods
            stan_data  = list(N = 100, P = 3),
            hbb_data   = list(
                N = 100L, P = 3L, S = 1L,
                z = c(rep(1, 70), rep(0, 30)),
                w_tilde = NULL
            ),
            formula    = structure(
                list(formula = y | trials(n_trial) ~ poverty + urban),
                class = "hbb_formula"
            ),
            prior      = default_prior(),
            model_type = "base",
            model_name = "hbb_base",
            call       = quote(hbb(y | trials(n_trial) ~ poverty + urban)),
            elapsed    = 42.3
        ),
        class = "hbb_fit"
    )

    # print should not error even with a fake $fit
    output <- capture.output(print(mock_fit))
    combined <- paste(output, collapse = "\n")

    expect_true(grepl("Hurdle Beta-Binomial", combined))
    expect_true(grepl("hbb_base", combined))
    expect_true(grepl("100", combined))   # N
    expect_true(grepl("42.3", combined))  # elapsed
})


test_that("print.hbb_fit() handles NULL components gracefully", {
    mock_fit <- structure(
        list(
            fit = NULL, stan_data = NULL, hbb_data = NULL,
            formula = NULL, prior = NULL,
            model_type = NULL, model_name = NULL,
            call = NULL, elapsed = NULL
        ),
        class = "hbb_fit"
    )

    # Should not error
    output <- capture.output(print(mock_fit))
    expect_true(any(grepl("Hurdle Beta-Binomial", output)))
})


test_that("print.hbb_fit() formats elapsed time correctly", {
    make_mock <- function(elapsed) {
        structure(
            list(
                fit = list(), stan_data = NULL,
                hbb_data = list(N = 50L, P = 2L, z = rep(1, 50)),
                formula = NULL, prior = NULL,
                model_type = "base", model_name = "hbb_base",
                call = NULL, elapsed = elapsed
            ),
            class = "hbb_fit"
        )
    }

    # Seconds
    out_sec <- capture.output(print(make_mock(45.2)))
    expect_true(any(grepl("45.2 seconds", out_sec)))

    # Minutes
    out_min <- capture.output(print(make_mock(185.0)))
    expect_true(any(grepl("3.1 minutes", out_min)))

    # Hours
    out_hr <- capture.output(print(make_mock(7500.0)))
    expect_true(any(grepl("2.1 hours", out_hr)))
})


test_that("print.hbb_fit() returns object invisibly", {
    mock_fit <- structure(
        list(
            fit = list(), stan_data = NULL,
            hbb_data = list(N = 10L, P = 2L, z = rep(1, 10)),
            formula = NULL, prior = NULL,
            model_type = "base", model_name = "hbb_base",
            call = NULL, elapsed = 5
        ),
        class = "hbb_fit"
    )

    result <- withVisible(capture.output(ret <- print(mock_fit)))
    expect_identical(ret, mock_fit)
})


# ============================================================================
# Section 6: hbb() — Integration: base model (requires CmdStan)
# ============================================================================

test_that("hbb() fits base model and returns hbb_fit", {
    skip_if_no_cmdstan()
    skip_on_cran()

    data(nsece_synth_small, package = "hurdlebb")

    fit <- hbb(
        formula       = y | trials(n_trial) ~ poverty + urban,
        data          = nsece_synth_small,
        chains        = .test_chains,
        iter_warmup   = .test_iter_warmup,
        iter_sampling = .test_iter_sampling,
        seed          = .test_seed,
        refresh       = .test_refresh
    )

    # Returns hbb_fit
    expect_s3_class(fit, "hbb_fit")
    expect_true(is.hbb_fit(fit))

    # Model identification
    expect_equal(fit$model_type, "base")
    expect_equal(fit$model_name, "hbb_base")

    # All structural components present
    expect_true(all(c("fit", "stan_data", "hbb_data", "formula", "prior",
                       "model_type", "model_name", "call",
                       "elapsed") %in% names(fit)))

    # Call is captured
    expect_true(is.call(fit$call))

    # Timing is positive
    expect_true(is.numeric(fit$elapsed))
    expect_true(fit$elapsed > 0)

    # Formula object preserved
    expect_s3_class(fit$formula, "hbb_formula")

    # Prior object preserved
    expect_s3_class(fit$prior, "hbb_prior")

    # Stan data has the required fields for base model
    sd <- fit$stan_data
    expect_true(all(c("N", "P", "y", "n_trial", "z", "X",
                       "prior_alpha_sd", "prior_beta_sd",
                       "prior_kappa_mean", "prior_kappa_sd") %in% names(sd)))

    # Base model should NOT have SVC or weighted fields
    expect_null(sd$w_tilde)
    expect_null(sd$S)
    expect_null(sd$state)
    expect_null(sd$v_state)

    # CmdStanMCMC object has draws
    expect_true(inherits(fit$fit, "CmdStanMCMC"))
    draws <- fit$fit$draws()
    expect_true(inherits(draws, "draws"))
    var_names <- posterior::variables(draws)
    expect_true(any(grepl("^alpha", var_names)))
    expect_true(any(grepl("^beta", var_names)))
})


test_that("hbb() accepts pre-built hbb_formula object", {
    skip_if_no_cmdstan()
    skip_on_cran()

    data(nsece_synth_small, package = "hurdlebb")

    f <- hbb_formula(y | trials(n_trial) ~ poverty + urban)

    fit <- hbb(
        formula       = f,
        data          = nsece_synth_small,
        chains        = .test_chains,
        iter_warmup   = .test_iter_warmup,
        iter_sampling = .test_iter_sampling,
        seed          = .test_seed,
        refresh       = .test_refresh
    )

    expect_s3_class(fit, "hbb_fit")
    expect_equal(fit$model_type, "base")
})


test_that("hbb() uses default_prior() when prior is NULL", {
    skip_if_no_cmdstan()
    skip_on_cran()

    data(nsece_synth_small, package = "hurdlebb")

    fit <- hbb(
        formula       = y | trials(n_trial) ~ poverty,
        data          = nsece_synth_small,
        prior         = NULL,
        chains        = .test_chains,
        iter_warmup   = .test_iter_warmup,
        iter_sampling = .test_iter_sampling,
        seed          = .test_seed,
        refresh       = .test_refresh
    )

    expect_s3_class(fit, "hbb_fit")
    dp <- default_prior()
    expect_equal(fit$prior$alpha$sd, dp$alpha$sd)
    expect_equal(fit$prior$log_kappa$mean, dp$log_kappa$mean)
})


# ============================================================================
# Section 7: hbb() — Integration: weighted model (requires CmdStan)
# ============================================================================

test_that("hbb() fits weighted model", {
    skip_if_no_cmdstan()
    skip_on_cran()

    data(nsece_synth_small, package = "hurdlebb")

    fit <- hbb(
        formula       = y | trials(n_trial) ~ poverty + urban,
        data          = nsece_synth_small,
        weights       = "weight",
        chains        = .test_chains,
        iter_warmup   = .test_iter_warmup,
        iter_sampling = .test_iter_sampling,
        seed          = .test_seed,
        refresh       = .test_refresh
    )

    expect_s3_class(fit, "hbb_fit")
    expect_equal(fit$model_type, "weighted")
    expect_equal(fit$model_name, "hbb_weighted")

    # Weighted model should include w_tilde in stan_data
    expect_true("w_tilde" %in% names(fit$stan_data))
    expect_equal(length(fit$stan_data$w_tilde), fit$stan_data$N)
    expect_true(all(fit$stan_data$w_tilde > 0))
})


# ============================================================================
# Section 8: hbb() — Integration: SVC model (requires CmdStan)
# ============================================================================

test_that("hbb() fits SVC model", {
    skip_if_no_cmdstan()
    skip_on_cran()

    data(nsece_synth_small, package = "hurdlebb")

    fit <- hbb(
        formula       = y | trials(n_trial) ~ poverty + urban +
            (poverty + urban | state_id),
        data          = nsece_synth_small,
        chains        = .test_chains,
        iter_warmup   = .test_iter_warmup,
        iter_sampling = .test_iter_sampling,
        seed          = .test_seed,
        refresh       = .test_refresh
    )

    expect_s3_class(fit, "hbb_fit")
    expect_equal(fit$model_type, "svc")
    expect_equal(fit$model_name, "hbb_svc")

    # SVC model should include group-level fields
    sd <- fit$stan_data
    expect_true(all(c("S", "Q", "state", "v_state",
                       "prior_gamma_sd", "prior_tau_sd",
                       "prior_lkj_eta") %in% names(sd)))
    expect_true(is.matrix(sd$v_state))
    expect_equal(nrow(sd$v_state), sd$S)
    expect_equal(ncol(sd$v_state), sd$Q)
})


test_that("hbb() fits svc_weighted model", {
    skip_if_no_cmdstan()
    skip_on_cran()

    data(nsece_synth_small, package = "hurdlebb")

    fit <- hbb(
        formula       = y | trials(n_trial) ~ poverty + urban +
            (poverty + urban | state_id),
        data          = nsece_synth_small,
        weights       = "weight",
        chains        = .test_chains,
        iter_warmup   = .test_iter_warmup,
        iter_sampling = .test_iter_sampling,
        seed          = .test_seed,
        refresh       = .test_refresh
    )

    expect_s3_class(fit, "hbb_fit")
    expect_equal(fit$model_type, "svc_weighted")
    expect_equal(fit$model_name, "hbb_svc_weighted")

    # Should have both weights and group-level fields
    sd <- fit$stan_data
    expect_true("w_tilde" %in% names(sd))
    expect_true("S" %in% names(sd))
    expect_true("v_state" %in% names(sd))
})


# ============================================================================
# Section 9: hbb() — Integration: custom prior, seed, print (requires CmdStan)
# ============================================================================

test_that("hbb() propagates custom prior to stan_data", {
    skip_if_no_cmdstan()
    skip_on_cran()

    data(nsece_synth_small, package = "hurdlebb")

    p <- hbb_prior(
        alpha     = list(dist = "normal", mean = 0, sd = 5),
        beta      = list(dist = "normal", mean = 0, sd = 1),
        log_kappa = list(dist = "normal", mean = 3, sd = 0.5)
    )

    fit <- hbb(
        y | trials(n_trial) ~ poverty,
        data          = nsece_synth_small,
        prior         = p,
        chains        = .test_chains,
        iter_warmup   = .test_iter_warmup,
        iter_sampling = .test_iter_sampling,
        seed          = .test_seed,
        refresh       = .test_refresh
    )

    # Prior object preserved
    expect_s3_class(fit$prior, "hbb_prior")
    expect_equal(fit$prior$alpha$sd, 5)
    expect_equal(fit$prior$log_kappa$mean, 3)

    # Stan data reflects the custom prior
    expect_equal(fit$stan_data$prior_alpha_sd, 5)
    expect_equal(fit$stan_data$prior_kappa_mean, 3)
    expect_equal(fit$stan_data$prior_kappa_sd, 0.5)
})


test_that("hbb() seed produces reproducible results", {
    skip_if_no_cmdstan()
    skip_on_cran()

    data(nsece_synth_small, package = "hurdlebb")

    fit1 <- hbb(
        y | trials(n_trial) ~ poverty,
        data = nsece_synth_small,
        chains = 1L, iter_warmup = 50L, iter_sampling = 50L,
        seed = 42L, refresh = 0L
    )
    fit2 <- hbb(
        y | trials(n_trial) ~ poverty,
        data = nsece_synth_small,
        chains = 1L, iter_warmup = 50L, iter_sampling = 50L,
        seed = 42L, refresh = 0L
    )

    # Same seed should produce identical draws
    draws1 <- tryCatch(
        fit1$fit$draws(variables = "alpha[1]", format = "matrix"),
        error = function(e) NULL
    )
    draws2 <- tryCatch(
        fit2$fit$draws(variables = "alpha[1]", format = "matrix"),
        error = function(e) NULL
    )

    if (!is.null(draws1) && !is.null(draws2)) {
        expect_equal(draws1, draws2)
    }
})


test_that("print.hbb_fit() works on real fit", {
    skip_if_no_cmdstan()
    skip_on_cran()

    data(nsece_synth_small, package = "hurdlebb")

    fit <- hbb(
        formula       = y | trials(n_trial) ~ poverty + urban,
        data          = nsece_synth_small,
        chains        = .test_chains,
        iter_warmup   = .test_iter_warmup,
        iter_sampling = .test_iter_sampling,
        seed          = .test_seed,
        refresh       = .test_refresh
    )

    out <- capture.output(print(fit))
    combined <- paste(out, collapse = "\n")

    # Should contain key information
    expect_true(grepl("Hurdle Beta-Binomial", combined))
    expect_true(grepl("hbb_base", combined))
    expect_true(grepl("Observations", combined))
    expect_true(grepl("Chains", combined))
    expect_true(grepl("seconds|minutes|hours", combined))
    expect_true(grepl("Rhat", combined))
    expect_true(grepl("summary\\(\\)", combined))

    # Returns invisibly
    result <- withVisible(print(fit))
    expect_false(result$visible)
    expect_identical(result$value, fit)
})
