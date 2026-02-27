# ============================================================================
# test-hbb-fit-class.R --- testthat tests for hbb-fit-class.R
#
# Covers:
#   1. is.hbb_fit (lines 110-112 in the uncovered report)
#   2. print.hbb_fit (lines 146,173,176-179,183-186,190-195,206-223)
#      - Total draws computation
#      - Diagnostics section (divergences, max treedepth, E-BFMI)
#      - Rhat/ESS section
#      - SVC/policy section
# ============================================================================


# ============================================================================
# SECTION A: is.hbb_fit
# ============================================================================

test_that("A1. is.hbb_fit returns TRUE for hbb_fit objects", {
    obj <- structure(list(), class = "hbb_fit")
    expect_true(is.hbb_fit(obj))
})

test_that("A2. is.hbb_fit returns FALSE for non-hbb_fit objects", {
    expect_false(is.hbb_fit(1))
    expect_false(is.hbb_fit(list()))
    expect_false(is.hbb_fit(NULL))
    expect_false(is.hbb_fit("hbb_fit"))
    expect_false(is.hbb_fit(data.frame()))
})

test_that("A3. is.hbb_fit returns TRUE for objects with multiple classes", {
    obj <- structure(list(), class = c("hbb_fit", "list"))
    expect_true(is.hbb_fit(obj))
})


# ============================================================================
# SECTION B: print.hbb_fit - Basic and Model Type Labels
# ============================================================================

# Helper: minimal mock hbb_fit for print testing
make_print_mock <- function(
    model_type = "base",
    model_name = "hbb_test",
    N = 50L, P = 3L,
    elapsed = 30,
    n_chains = 4L,
    iter_warmup = 500L,
    iter_sampling = 1000L,
    has_diagnostics = FALSE,
    has_summary = FALSE,
    svc = FALSE,
    S = NULL, Q = NULL, K = NULL,
    formula_obj = y ~ x1 + x2
) {
    hbb_data <- list(
        N = N,
        P = P,
        z = rbinom(N, 1, 0.65),
        y = sample(0:10, N, replace = TRUE),
        n_trial = rep(20L, N)
    )
    if (svc && !is.null(S)) {
        hbb_data$S <- S
        if (!is.null(Q)) hbb_data$Q <- Q
        if (!is.null(K)) hbb_data$K <- K
    }

    mock_fit <- list(
        num_chains = function() n_chains,
        metadata = function() {
            list(
                iter_warmup   = iter_warmup,
                iter_sampling = iter_sampling
            )
        }
    )

    if (has_diagnostics) {
        mock_fit$diagnostic_summary <- function(quiet = TRUE) {
            list(
                num_divergent       = c(0L, 0L, 0L, 0L),
                num_max_treedepth   = c(0L, 0L, 0L, 0L),
                ebfmi               = c(0.85, 0.92, 0.88, 0.90)
            )
        }
    }

    if (has_summary) {
        mock_fit$summary <- function() {
            data.frame(
                rhat     = c(1.001, 1.002, 1.003),
                ess_bulk = c(800, 1200, 600),
                ess_tail = c(500, 700, 400)
            )
        }
    }

    structure(
        list(
            fit        = mock_fit,
            hbb_data   = hbb_data,
            model_type = model_type,
            model_name = model_name,
            elapsed    = elapsed,
            formula    = list(formula = formula_obj)
        ),
        class = "hbb_fit"
    )
}


test_that("B1. print.hbb_fit produces basic output for base model", {
    fit <- make_print_mock(model_type = "base")
    out <- capture.output(print(fit))
    expect_true(any(grepl("Hurdle Beta-Binomial", out)))
    expect_true(any(grepl("Base", out)))
    expect_true(any(grepl("hbb_test", out)))
})

test_that("B2. print.hbb_fit shows weighted model type", {
    fit <- make_print_mock(model_type = "weighted")
    out <- capture.output(print(fit))
    expect_true(any(grepl("Weighted", out)))
})

test_that("B3. print.hbb_fit shows SVC model type and group dimensions", {
    fit <- make_print_mock(
        model_type = "svc",
        svc = TRUE,
        S = 51L,
        Q = 4L,
        K = 6L
    )
    out <- capture.output(print(fit))
    expect_true(any(grepl("SVC", out)))
    expect_true(any(grepl("51", out)))
    expect_true(any(grepl("4", out)))
    expect_true(any(grepl("6", out)))
})

test_that("B4. print.hbb_fit shows SVC weighted model type", {
    fit <- make_print_mock(
        model_type = "svc_weighted",
        svc = TRUE,
        S = 51L
    )
    out <- capture.output(print(fit))
    expect_true(any(grepl("SVC.*survey-weighted", out)))
})

test_that("B5. print.hbb_fit handles unknown model type", {
    fit <- make_print_mock(model_type = "custom_model")
    out <- capture.output(print(fit))
    expect_true(any(grepl("custom_model", out)))
})


# ============================================================================
# SECTION C: print.hbb_fit - MCMC and Total Draws
# ============================================================================

test_that("C1. print.hbb_fit shows total draws = chains * iter_sampling", {
    fit <- make_print_mock(n_chains = 4L, iter_sampling = 1000L)
    out <- capture.output(print(fit))
    expect_true(any(grepl("4000", out)))
    expect_true(any(grepl("Chains.*4", out)))
})

test_that("C2. print.hbb_fit shows (unknown) when metadata fails", {
    fit <- make_print_mock()
    fit$fit$metadata <- function() stop("metadata unavailable")
    out <- capture.output(print(fit))
    expect_true(any(grepl("unknown", out)))
})

test_that("C3. print.hbb_fit handles elapsed time in seconds", {
    fit <- make_print_mock(elapsed = 45)
    out <- capture.output(print(fit))
    expect_true(any(grepl("45.0 seconds", out)))
})

test_that("C4. print.hbb_fit handles elapsed time in minutes", {
    fit <- make_print_mock(elapsed = 300)
    out <- capture.output(print(fit))
    expect_true(any(grepl("5.0 minutes", out)))
})

test_that("C5. print.hbb_fit handles elapsed time in hours", {
    fit <- make_print_mock(elapsed = 7200)
    out <- capture.output(print(fit))
    expect_true(any(grepl("2.0 hours", out)))
})


# ============================================================================
# SECTION D: print.hbb_fit - Diagnostics Section
# ============================================================================

test_that("D1. print.hbb_fit shows diagnostics with OK status", {
    fit <- make_print_mock(has_diagnostics = TRUE)
    out <- capture.output(print(fit))
    expect_true(any(grepl("Diagnostics", out)))
    expect_true(any(grepl("Divergent.*0.*OK", out)))
    expect_true(any(grepl("Max treedepth.*0.*OK", out)))
    expect_true(any(grepl("E-BFMI", out)))
})

test_that("D2. print.hbb_fit shows WARNING for divergences > 0", {
    fit <- make_print_mock()
    fit$fit$diagnostic_summary <- function(quiet = TRUE) {
        list(
            num_divergent       = c(5L, 0L, 2L, 0L),
            num_max_treedepth   = c(0L, 0L, 0L, 0L),
            ebfmi               = c(0.85, 0.92, 0.88, 0.90)
        )
    }
    out <- capture.output(print(fit))
    expect_true(any(grepl("7.*WARNING", out)))
})

test_that("D3. print.hbb_fit shows WARNING for max treedepth hits > 0", {
    fit <- make_print_mock()
    fit$fit$diagnostic_summary <- function(quiet = TRUE) {
        list(
            num_divergent       = c(0L, 0L, 0L, 0L),
            num_max_treedepth   = c(10L, 5L, 0L, 0L),
            ebfmi               = c(0.85, 0.92, 0.88, 0.90)
        )
    }
    out <- capture.output(print(fit))
    expect_true(any(grepl("15.*WARNING", out)))
})

test_that("D4. print.hbb_fit shows WARNING for low E-BFMI (< 0.2)", {
    fit <- make_print_mock()
    fit$fit$diagnostic_summary <- function(quiet = TRUE) {
        list(
            num_divergent       = c(0L, 0L, 0L, 0L),
            num_max_treedepth   = c(0L, 0L, 0L, 0L),
            ebfmi               = c(0.05, 0.10, 0.88, 0.90)
        )
    }
    out <- capture.output(print(fit))
    expect_true(any(grepl("WARNING.*Low E-BFMI", out)))
})

test_that("D5. print.hbb_fit handles diagnostics failure gracefully", {
    fit <- make_print_mock()
    fit$fit$diagnostic_summary <- function(quiet = TRUE) {
        stop("diagnostics unavailable")
    }
    # Should not error, just skip the section
    out <- capture.output(print(fit))
    expect_true(any(grepl("Hurdle Beta-Binomial", out)))
    expect_false(any(grepl("Diagnostics", out)))
})


# ============================================================================
# SECTION E: print.hbb_fit - Rhat and ESS Section
# ============================================================================

test_that("E1. print.hbb_fit shows Rhat and ESS when summary available", {
    fit <- make_print_mock(has_summary = TRUE)
    out <- capture.output(print(fit))
    expect_true(any(grepl("Max Rhat", out)))
    expect_true(any(grepl("Min bulk ESS", out)))
    expect_true(any(grepl("Min tail ESS", out)))
})

test_that("E2. print.hbb_fit shows WARNING for high Rhat (> 1.01)", {
    fit <- make_print_mock()
    fit$fit$summary <- function() {
        data.frame(
            rhat     = c(1.05, 1.02, 1.001),
            ess_bulk = c(800, 1200, 600),
            ess_tail = c(500, 700, 400)
        )
    }
    out <- capture.output(print(fit))
    expect_true(any(grepl("1.05.*WARNING", out)))
})

test_that("E3. print.hbb_fit shows LOW for small ESS (< 400)", {
    fit <- make_print_mock()
    fit$fit$summary <- function() {
        data.frame(
            rhat     = c(1.001, 1.002, 1.001),
            ess_bulk = c(100, 200, 150),
            ess_tail = c(50, 60, 70)
        )
    }
    out <- capture.output(print(fit))
    expect_true(any(grepl("Min bulk ESS.*LOW", out)))
    expect_true(any(grepl("Min tail ESS.*LOW", out)))
})

test_that("E4. print.hbb_fit handles missing ess_tail column", {
    fit <- make_print_mock()
    fit$fit$summary <- function() {
        data.frame(
            rhat     = c(1.001, 1.002),
            ess_bulk = c(800, 1200)
        )
    }
    out <- capture.output(print(fit))
    expect_true(any(grepl("Max Rhat", out)))
    expect_true(any(grepl("Min bulk ESS", out)))
    # No tail ESS line
    expect_false(any(grepl("Min tail ESS", out)))
})

test_that("E5. print.hbb_fit skips summary section when summary() fails", {
    fit <- make_print_mock()
    fit$fit$summary <- function() {
        stop("summary unavailable")
    }
    out <- capture.output(print(fit))
    expect_false(any(grepl("Max Rhat", out)))
})


# ============================================================================
# SECTION F: print.hbb_fit - Data Unavailable Path
# ============================================================================

test_that("F1. print.hbb_fit shows (data summary unavailable) when hbb_data NULL", {
    fit <- make_print_mock()
    fit$hbb_data <- NULL
    out <- capture.output(print(fit))
    expect_true(any(grepl("data summary unavailable", out)))
})


# ============================================================================
# SECTION G: print.hbb_fit - Num Chains Failure
# ============================================================================

test_that("G1. print.hbb_fit shows (unknown) when num_chains fails", {
    fit <- make_print_mock()
    fit$fit$num_chains <- function() stop("no chains")
    out <- capture.output(print(fit))
    expect_true(any(grepl("Chains.*unknown", out)))
})
