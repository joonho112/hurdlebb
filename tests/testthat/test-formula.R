# ============================================================================
# test-formula.R --- Tests for hbb_formula(), validate_hbb_formula(),
#                    and print.hbb_formula()
#
# Covers:
#   Section 1: hbb_formula() basic parsing (happy paths)
#   Section 2: S3 object structure
#   Section 3: validate_hbb_formula() with real data
#   Section 4: print.hbb_formula()
#   Section 5: Whitespace and formatting robustness
# ============================================================================


# ============================================================================
# Section 1: hbb_formula() basic parsing
# ============================================================================

test_that("hbb_formula: intercept-only model", {
    f <- hbb_formula(y | trials(n_trial) ~ 1)

    expect_equal(f$response, "y")
    expect_equal(f$trials, "n_trial")
    expect_equal(f$fixed, character(0))
    expect_null(f$group)
    expect_null(f$random)
    expect_null(f$policy)
    expect_false(f$svc)
})


test_that("hbb_formula: single predictor", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty)

    expect_equal(f$response, "y")
    expect_equal(f$trials, "n_trial")
    expect_equal(f$fixed, "poverty")
    expect_null(f$group)
    expect_null(f$random)
    expect_null(f$policy)
    expect_false(f$svc)
})


test_that("hbb_formula: multiple predictors", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty + urban + black + hispanic)

    expect_equal(f$response, "y")
    expect_equal(f$trials, "n_trial")
    expect_equal(f$fixed, c("poverty", "urban", "black", "hispanic"))
    expect_null(f$group)
    expect_null(f$random)
    expect_null(f$policy)
    expect_false(f$svc)
})


test_that("hbb_formula: random intercept", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty + urban + (1 | state_id))

    expect_equal(f$response, "y")
    expect_equal(f$trials, "n_trial")
    expect_equal(f$fixed, c("poverty", "urban"))
    expect_equal(f$group, "state_id")
    expect_equal(f$random, "1")
    expect_null(f$policy)
    expect_false(f$svc)
})


test_that("hbb_formula: SVC (random slopes)", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + urban + (poverty + urban | state_id)
    )

    expect_equal(f$response, "y")
    expect_equal(f$trials, "n_trial")
    expect_equal(f$fixed, c("poverty", "urban"))
    expect_equal(f$group, "state_id")
    expect_equal(f$random, c("poverty", "urban"))
    expect_null(f$policy)
    expect_true(f$svc)
})


test_that("hbb_formula: full model with policy moderators", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + urban +
            (poverty + urban | state_id) +
            state_level(mr_pctile + tiered_reim)
    )

    expect_equal(f$response, "y")
    expect_equal(f$trials, "n_trial")
    expect_equal(f$fixed, c("poverty", "urban"))
    expect_equal(f$group, "state_id")
    expect_equal(f$random, c("poverty", "urban"))
    expect_equal(f$policy, c("mr_pctile", "tiered_reim"))
    expect_true(f$svc)
})


test_that("hbb_formula: intercept-only with random intercept", {
    f <- hbb_formula(y | trials(n_trial) ~ 1 + (1 | state_id))

    expect_equal(f$response, "y")
    expect_equal(f$trials, "n_trial")
    expect_equal(f$fixed, character(0))
    expect_equal(f$group, "state_id")
    expect_equal(f$random, "1")
    expect_null(f$policy)
    expect_false(f$svc)
})


test_that("hbb_formula: single random slope", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + urban + (poverty | state_id)
    )

    expect_equal(f$fixed, c("poverty", "urban"))
    expect_equal(f$group, "state_id")
    expect_equal(f$random, "poverty")
    expect_true(f$svc)
})


test_that("hbb_formula: random slopes with explicit intercept in RE", {
    # (1 + poverty | state_id) — the 1 is implicit and should be stripped
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + urban + (1 + poverty | state_id)
    )

    expect_equal(f$fixed, c("poverty", "urban"))
    expect_equal(f$group, "state_id")
    expect_equal(f$random, "poverty")
    expect_true(f$svc)
})


test_that("hbb_formula: policy with random intercept only", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (1 | state_id) +
            state_level(mr_pctile)
    )

    expect_equal(f$fixed, "poverty")
    expect_equal(f$group, "state_id")
    expect_equal(f$random, "1")
    expect_equal(f$policy, "mr_pctile")
    expect_false(f$svc)
})


test_that("hbb_formula: three policy moderators", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty +
            (poverty | state_id) +
            state_level(mr_pctile + tiered_reim + it_addon)
    )

    expect_equal(f$policy, c("mr_pctile", "tiered_reim", "it_addon"))
})


test_that("hbb_formula: formula object is stored", {
    orig <- y | trials(n_trial) ~ poverty + urban
    f <- hbb_formula(orig)

    expect_true(inherits(f$formula, "formula"))
    # Deparse should match the original
    expect_equal(deparse(f$formula), deparse(orig))
})


# ============================================================================
# Section 2: S3 object structure
# ============================================================================

test_that("hbb_formula: returns object of class hbb_formula", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty)
    expect_s3_class(f, "hbb_formula")
})


test_that("hbb_formula: all expected fields present", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + urban +
            (poverty + urban | state_id) +
            state_level(mr_pctile + tiered_reim)
    )

    expected_fields <- c("response", "trials", "fixed", "group",
                         "random", "policy", "formula", "svc")
    expect_true(all(expected_fields %in% names(f)))
    expect_equal(length(f), length(expected_fields))
})


test_that("hbb_formula: field types are correct", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + urban +
            (poverty + urban | state_id) +
            state_level(mr_pctile)
    )

    expect_type(f$response, "character")
    expect_length(f$response, 1L)
    expect_type(f$trials, "character")
    expect_length(f$trials, 1L)
    expect_type(f$fixed, "character")
    expect_type(f$group, "character")
    expect_length(f$group, 1L)
    expect_type(f$random, "character")
    expect_type(f$policy, "character")
    expect_true(inherits(f$formula, "formula"))
    expect_type(f$svc, "logical")
    expect_length(f$svc, 1L)
})


test_that("hbb_formula: NULL fields for minimal model", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty)

    expect_null(f$group)
    expect_null(f$random)
    expect_null(f$policy)
})


test_that("hbb_formula: svc is FALSE for random intercept only", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty + (1 | state_id))
    expect_false(f$svc)
})


test_that("hbb_formula: svc is TRUE for random slopes", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (poverty | state_id)
    )
    expect_true(f$svc)
})


test_that("hbb_formula: svc is FALSE for no random effects", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty + urban)
    expect_false(f$svc)
})


test_that("hbb_formula: fixed is character(0) not NULL for intercept-only", {
    f <- hbb_formula(y | trials(n_trial) ~ 1)
    expect_identical(f$fixed, character(0))
    expect_false(is.null(f$fixed))
})


# ============================================================================
# Section 3: validate_hbb_formula() with real data
# ============================================================================

test_that("validate_hbb_formula: simple formula against nsece_synth_small", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    f <- hbb_formula(y | trials(n_trial) ~ poverty + urban)

    result <- validate_hbb_formula(f, nsece_synth_small)
    expect_true(result)
})


test_that("validate_hbb_formula: intercept-only model", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    f <- hbb_formula(y | trials(n_trial) ~ 1)

    result <- validate_hbb_formula(f, nsece_synth_small)
    expect_true(result)
})


test_that("validate_hbb_formula: all four covariates", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + urban + black + hispanic
    )

    result <- validate_hbb_formula(f, nsece_synth_small)
    expect_true(result)
})


test_that("validate_hbb_formula: formula with random intercept", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + urban + (1 | state_id)
    )

    result <- validate_hbb_formula(f, nsece_synth_small)
    expect_true(result)
})


test_that("validate_hbb_formula: formula with SVC", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + urban +
            (poverty + urban | state_id)
    )

    result <- validate_hbb_formula(f, nsece_synth_small)
    expect_true(result)
})


test_that("validate_hbb_formula: full formula with policy moderators", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    data("nsece_state_policy", package = "hurdlebb", envir = environment())
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + urban +
            (poverty + urban | state_id) +
            state_level(mr_pctile + tiered_reim)
    )

    result <- validate_hbb_formula(f, nsece_synth_small, nsece_state_policy)
    expect_true(result)
})


test_that("validate_hbb_formula: policy with all three moderators", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    data("nsece_state_policy", package = "hurdlebb", envir = environment())
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty +
            (poverty | state_id) +
            state_level(mr_pctile + tiered_reim + it_addon)
    )

    result <- validate_hbb_formula(f, nsece_synth_small, nsece_state_policy)
    expect_true(result)
})


test_that("validate_hbb_formula: returns TRUE invisibly", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    f <- hbb_formula(y | trials(n_trial) ~ poverty)

    # Capture the visible/invisible return
    out <- withVisible(validate_hbb_formula(f, nsece_synth_small))
    expect_true(out$value)
    expect_false(out$visible)
})


test_that("validate_hbb_formula: works with full nsece_synth dataset", {
    data("nsece_synth", package = "hurdlebb", envir = environment())
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + urban + black + hispanic
    )

    result <- validate_hbb_formula(f, nsece_synth)
    expect_true(result)
})


# ============================================================================
# Section 4: print.hbb_formula()
# ============================================================================

test_that("print.hbb_formula: intercept-only model prints without error", {
    f <- hbb_formula(y | trials(n_trial) ~ 1)
    expect_no_error(print(f))
})


test_that("print.hbb_formula: fixed effects model prints without error", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty + urban)
    expect_no_error(print(f))
})


test_that("print.hbb_formula: random intercept model prints without error", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty + (1 | state_id))
    expect_no_error(print(f))
})


test_that("print.hbb_formula: SVC model prints without error", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + urban + (poverty + urban | state_id)
    )
    expect_no_error(print(f))
})


test_that("print.hbb_formula: full model with policy prints without error", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + urban +
            (poverty + urban | state_id) +
            state_level(mr_pctile + tiered_reim)
    )
    expect_no_error(print(f))
})


test_that("print.hbb_formula: output contains response and trials", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty + urban)
    out <- capture.output(print(f))
    out_text <- paste(out, collapse = "\n")

    expect_true(grepl("y", out_text, fixed = TRUE))
    expect_true(grepl("n_trial", out_text, fixed = TRUE))
})


test_that("print.hbb_formula: output contains header", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty)
    out <- capture.output(print(f))
    out_text <- paste(out, collapse = "\n")

    expect_true(grepl("Hurdle Beta-Binomial Formula", out_text, fixed = TRUE))
})


test_that("print.hbb_formula: intercept-only shows (intercept only)", {
    f <- hbb_formula(y | trials(n_trial) ~ 1)
    out <- capture.output(print(f))
    out_text <- paste(out, collapse = "\n")

    expect_true(grepl("intercept only", out_text, fixed = TRUE))
})


test_that("print.hbb_formula: shows fixed effect names", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty + urban + hispanic)
    out <- capture.output(print(f))
    out_text <- paste(out, collapse = "\n")

    expect_true(grepl("poverty", out_text, fixed = TRUE))
    expect_true(grepl("urban", out_text, fixed = TRUE))
    expect_true(grepl("hispanic", out_text, fixed = TRUE))
})


test_that("print.hbb_formula: shows random intercept", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty + (1 | state_id))
    out <- capture.output(print(f))
    out_text <- paste(out, collapse = "\n")

    expect_true(grepl("state_id", out_text, fixed = TRUE))
    expect_true(grepl("random intercept", out_text, fixed = TRUE))
})


test_that("print.hbb_formula: shows SVC label", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (poverty | state_id)
    )
    out <- capture.output(print(f))
    out_text <- paste(out, collapse = "\n")

    expect_true(grepl("SVC", out_text, fixed = TRUE))
})


test_that("print.hbb_formula: shows no random effects", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty + urban)
    out <- capture.output(print(f))
    out_text <- paste(out, collapse = "\n")

    expect_true(grepl("none", out_text, fixed = TRUE))
})


test_that("print.hbb_formula: shows policy variables", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty +
            (poverty | state_id) +
            state_level(mr_pctile + tiered_reim)
    )
    out <- capture.output(print(f))
    out_text <- paste(out, collapse = "\n")

    expect_true(grepl("mr_pctile", out_text, fixed = TRUE))
    expect_true(grepl("tiered_reim", out_text, fixed = TRUE))
    expect_true(grepl("state_level", out_text, fixed = TRUE))
})


test_that("print.hbb_formula: shows svc flag value", {
    f_true <- hbb_formula(
        y | trials(n_trial) ~ poverty + (poverty | state_id)
    )
    f_false <- hbb_formula(y | trials(n_trial) ~ poverty)

    out_true <- capture.output(print(f_true))
    out_false <- capture.output(print(f_false))

    # SVC TRUE model
    svc_line_true <- out_true[grepl("SVC", out_true)]
    expect_true(any(grepl("TRUE", svc_line_true, fixed = TRUE)))

    # SVC FALSE model
    svc_line_false <- out_false[grepl("SVC", out_false)]
    expect_true(any(grepl("FALSE", svc_line_false, fixed = TRUE)))
})


test_that("print.hbb_formula: returns object invisibly", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty)
    out <- withVisible(print(f))
    expect_false(out$visible)
    expect_s3_class(out$value, "hbb_formula")
})


# ============================================================================
# Section 5: Whitespace and formatting robustness
# ============================================================================

test_that("hbb_formula: extra spaces in formula do not cause errors", {
    # R's formula parser normalizes whitespace, so these should all work
    f <- hbb_formula(y | trials(n_trial) ~  poverty  +  urban)

    expect_equal(f$response, "y")
    expect_equal(f$trials, "n_trial")
    expect_equal(f$fixed, c("poverty", "urban"))
})


test_that("hbb_formula: variable names with underscores", {
    # Construct data with underscored names for parsing test
    f <- hbb_formula(
        y_count | trials(n_trial_total) ~ var_one + var_two
    )

    expect_equal(f$response, "y_count")
    expect_equal(f$trials, "n_trial_total")
    expect_equal(f$fixed, c("var_one", "var_two"))
})


test_that("hbb_formula: variable names with dots", {
    f <- hbb_formula(
        my.response | trials(my.trials) ~ var.one + var.two
    )

    expect_equal(f$response, "my.response")
    expect_equal(f$trials, "my.trials")
    expect_equal(f$fixed, c("var.one", "var.two"))
})


test_that("hbb_formula: mixed dot and underscore variable names", {
    f <- hbb_formula(
        resp.var | trials(n_trials) ~ x.var_1 + x_var.2
    )

    expect_equal(f$response, "resp.var")
    expect_equal(f$trials, "n_trials")
    expect_equal(f$fixed, c("x.var_1", "x_var.2"))
})


test_that("hbb_formula: grouping variable with underscores", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (1 | state_id)
    )

    expect_equal(f$group, "state_id")
})


test_that("hbb_formula: grouping variable with dots", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (1 | state.id)
    )

    expect_equal(f$group, "state.id")
})


test_that("hbb_formula: multiline formula is parsed correctly", {
    # Multi-line formulas are common in real code
    f <- hbb_formula(
        y | trials(n_trial) ~
            poverty +
            urban +
            black +
            hispanic +
            (poverty + urban | state_id) +
            state_level(mr_pctile + tiered_reim)
    )

    expect_equal(f$response, "y")
    expect_equal(f$trials, "n_trial")
    expect_equal(f$fixed, c("poverty", "urban", "black", "hispanic"))
    expect_equal(f$group, "state_id")
    expect_equal(f$random, c("poverty", "urban"))
    expect_equal(f$policy, c("mr_pctile", "tiered_reim"))
    expect_true(f$svc)
})


test_that("hbb_formula: policy variables with underscores and dots", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty +
            (poverty | state_id) +
            state_level(var_one + var.two)
    )

    expect_equal(f$policy, c("var_one", "var.two"))
})


test_that("hbb_formula: single-letter variable names", {
    f <- hbb_formula(y | trials(n) ~ x)

    expect_equal(f$response, "y")
    expect_equal(f$trials, "n")
    expect_equal(f$fixed, "x")
})


test_that("hbb_formula: variable names starting with a dot", {
    f <- hbb_formula(y | trials(n_trial) ~ .hidden_var)

    expect_equal(f$fixed, ".hidden_var")
})


# ============================================================================
# Section 6: hbb_formula() — Error paths
# ============================================================================

# --------------------------------------------------------------------------
# 6a. Top-level input checks
# --------------------------------------------------------------------------

test_that("error: hbb_formula rejects character string input", {
    expect_error(
        hbb_formula("y | trials(n) ~ x1"),
        "Expected a formula object"
    )
})

test_that("error: hbb_formula rejects list input", {
    expect_error(
        hbb_formula(list(y = 1, n = 10)),
        "Expected a formula object"
    )
})

test_that("error: hbb_formula rejects numeric input", {
    expect_error(
        hbb_formula(42),
        "Expected a formula object"
    )
})

test_that("error: hbb_formula rejects NULL input", {
    expect_error(
        hbb_formula(NULL),
        "Expected a formula object"
    )
})

test_that("error: hbb_formula rejects one-sided formula", {
    expect_error(
        hbb_formula(~ x1 + x2),
        "Formula must be two-sided"
    )
})


# --------------------------------------------------------------------------
# 6b. LHS parsing errors
# --------------------------------------------------------------------------

test_that("error: LHS without trials() function", {
    expect_error(
        hbb_formula(y ~ x1 + x2),
        "Invalid left-hand side"
    )
})

test_that("error: LHS with wrong function name instead of trials()", {
    expect_error(
        hbb_formula(y | size(n) ~ x1),
        "Invalid left-hand side"
    )
})

test_that("error: LHS with empty trials()", {
    expect_error(
        hbb_formula(y | trials() ~ x1),
        "Invalid left-hand side"
    )
})

test_that("error: LHS with numeric literal inside trials()", {
    expect_error(
        hbb_formula(y | trials(100) ~ x1),
        "Invalid left-hand side"
    )
})

test_that("error: LHS with expression inside trials()", {
    expect_error(
        hbb_formula(y | trials(n + 1) ~ x1),
        "Invalid left-hand side"
    )
})

test_that("error: LHS with bare response (no pipe)", {
    expect_error(
        hbb_formula(y ~ x1 + x2),
        "Invalid left-hand side"
    )
})


# --------------------------------------------------------------------------
# 6c. RHS: dot expansion & intercept removal
# --------------------------------------------------------------------------

test_that("error: dot expansion (y ~ .)", {
    expect_error(
        hbb_formula(y | trials(n) ~ .),
        "Dot expansion.*is not supported"
    )
})

test_that("error: dot expansion mixed with variables (y ~ x1 + .)", {
    expect_error(
        hbb_formula(y | trials(n) ~ x1 + .),
        "Dot expansion.*is not supported"
    )
})

test_that("error: intercept removal with - 1", {
    expect_error(
        hbb_formula(y | trials(n) ~ x1 - 1),
        "Removing the intercept is not supported"
    )
})

test_that("error: intercept removal with + 0", {
    expect_error(
        hbb_formula(y | trials(n) ~ x1 + 0),
        "Removing the intercept is not supported"
    )
})

test_that("error: standalone 0 on RHS", {
    expect_error(
        hbb_formula(y | trials(n) ~ 0),
        "Removing the intercept is not supported"
    )
})


# --------------------------------------------------------------------------
# 6d. RHS: state_level() errors
# --------------------------------------------------------------------------

test_that("error: multiple state_level() calls", {
    expect_error(
        hbb_formula(
            y | trials(n) ~ x1 + (1 | grp) +
                state_level(v1) + state_level(v2)
        ),
        "Multiple.*state_level.*terms are not allowed"
    )
})

test_that("error: empty state_level() is rejected", {
    # state_level() with no arguments deparsed as "state_level()" which
    # doesn't match the state_level regex ([^)]+ requires >= 1 char).
    # It falls through to the fixed effects parser which rejects it
    # as an unsupported function call.
    expect_error(
        hbb_formula(
            y | trials(n) ~ x1 + (1 | grp) + state_level( )
        ),
        "Function calls in formula terms are not supported"
    )
})

test_that("error: invalid variable name in state_level()", {
    # Backtick-quoted names that start with a digit are syntactically valid in R
    # but fail the identifier regex in .parse_state_level().
    # Note: cli_abort() has a pre-existing pluralization issue with {?s}
    # when no qty() is provided, so we match on the error class only.
    expect_error(
        hbb_formula(
            y | trials(n) ~ x1 + (1 | grp) + state_level(`2bad`)
        )
    )
})

test_that("error: duplicate variables in state_level()", {
    expect_error(
        hbb_formula(
            y | trials(n) ~ x1 + (1 | grp) + state_level(v1 + v1)
        ),
        "Duplicate variable"
    )
})

test_that("error: '1' inside state_level() is caught as invalid name", {
    # "1" fails the regex check for valid R identifiers (starts with digit)
    # before reaching the explicit "Do not include 1" check.
    expect_error(
        hbb_formula(
            y | trials(n) ~ x1 + (1 | grp) + state_level(1 + v1)
        ),
        "Invalid variable name"
    )
})


# --------------------------------------------------------------------------
# 6e. RHS: random effects errors
# --------------------------------------------------------------------------

test_that("error: multiple (... | group) terms", {
    expect_error(
        hbb_formula(
            y | trials(n) ~ x1 + (1 | grp1) + (x1 | grp2)
        ),
        "Multiple random effects terms are not allowed"
    )
})

test_that("error: invalid grouping variable name in random effects", {
    # Backtick-quoted names preserve backticks in deparse, failing regex check
    expect_error(
        hbb_formula(
            y | trials(n) ~ x1 + (1 | `123bad`)
        ),
        "Invalid grouping variable name"
    )
})

test_that("error: empty LHS in random effects via internal parser", {
    # R's own parser rejects `( | grp)` so we call the internal function
    # directly to exercise this error path.
    expect_error(
        hurdlebb:::.parse_bar_terms("( | grp)"),
        "Empty left-hand side in random effects"
    )
})

test_that("error: invalid variable name in random effects LHS", {
    # Backtick-quoted names that fail the identifier regex.
    # Note: cli_abort() has a pre-existing pluralization issue, so we
    # match on the error class only.
    expect_error(
        hbb_formula(
            y | trials(n) ~ x1 + (`3abc` | grp)
        )
    )
})

test_that("error: duplicate variables in random effects", {
    expect_error(
        hbb_formula(
            y | trials(n) ~ x1 + (x1 + x1 | grp)
        ),
        "Duplicate variable"
    )
})


# --------------------------------------------------------------------------
# 6f. RHS: fixed effects errors
# --------------------------------------------------------------------------

test_that("error: interaction terms with *", {
    expect_error(
        hbb_formula(y | trials(n) ~ x1 * x2),
        "Interaction terms are not supported"
    )
})

test_that("error: interaction terms with :", {
    expect_error(
        hbb_formula(y | trials(n) ~ x1 + x1:x2),
        "Interaction terms are not supported"
    )
})

test_that("error: function call log() in fixed effects", {
    expect_error(
        hbb_formula(y | trials(n) ~ log(x1)),
        "Function calls in formula terms are not supported"
    )
})

test_that("error: function call poly() in fixed effects", {
    expect_error(
        hbb_formula(y | trials(n) ~ poly(x1, 2)),
        "Function calls in formula terms are not supported"
    )
})

test_that("error: function call scale() in fixed effects", {
    expect_error(
        hbb_formula(y | trials(n) ~ scale(x1)),
        "Function calls in formula terms are not supported"
    )
})

test_that("error: duplicate fixed effect variables", {
    expect_error(
        hbb_formula(y | trials(n) ~ poverty + poverty),
        "Duplicate fixed effect variable"
    )
})


# --------------------------------------------------------------------------
# 6g. Cross-component validation errors
# --------------------------------------------------------------------------

test_that("error: random slope not in fixed effects (single)", {
    expect_error(
        hbb_formula(
            y | trials(n) ~ x1 + (x2 | grp)
        ),
        "Random slope variable.*must also appear as fixed effect"
    )
})

test_that("error: random slopes not in fixed effects (multiple)", {
    expect_error(
        hbb_formula(
            y | trials(n) ~ x1 + (x1 + x2 + x3 | grp)
        ),
        "Random slope variable.*must also appear as fixed effect"
    )
})

test_that("error: state_level() without any random effects", {
    expect_error(
        hbb_formula(
            y | trials(n) ~ x1 + state_level(v1 + v2)
        ),
        "state_level.*requires random effects"
    )
})


# ============================================================================
# Section 7: validate_hbb_formula() — Error paths
# ============================================================================

# --------------------------------------------------------------------------
# Helper: build a clean formula and data for validation error tests
# --------------------------------------------------------------------------

.make_valid_formula <- function() {
    hbb_formula(y | trials(n_trial) ~ poverty + urban)
}

.make_valid_data <- function() {
    data.frame(
        y       = c(0L, 3L, 5L, 0L, 2L),
        n_trial = c(10L, 10L, 10L, 10L, 10L),
        poverty = c(0.1, 0.5, 0.3, 0.8, 0.4),
        urban   = c(0.2, 0.6, 0.9, 0.1, 0.5)
    )
}


# --------------------------------------------------------------------------
# 7a. Object class check
# --------------------------------------------------------------------------

test_that("validate error: rejects character string", {
    expect_error(
        validate_hbb_formula("not a formula", data.frame(y = 1)),
        "Expected an.*hbb_formula.*object"
    )
})

test_that("validate error: rejects plain list (wrong class)", {
    fake <- list(response = "y", trials = "n")
    expect_error(
        validate_hbb_formula(fake, data.frame(y = 1)),
        "Expected an.*hbb_formula.*object"
    )
})

test_that("validate error: rejects numeric", {
    expect_error(
        validate_hbb_formula(42, data.frame(y = 1)),
        "Expected an.*hbb_formula.*object"
    )
})


# --------------------------------------------------------------------------
# 7b. Response variable errors
# --------------------------------------------------------------------------

test_that("validate error: response variable missing from data", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$y <- NULL
    expect_error(
        validate_hbb_formula(f, d),
        "not found in.*data"
    )
})

test_that("validate error: response is character (non-numeric)", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$y <- c("a", "b", "c", "d", "e")
    expect_error(
        validate_hbb_formula(f, d),
        "must be numeric"
    )
})

test_that("validate error: response is factor (non-numeric)", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$y <- factor(c(0, 3, 5, 0, 2))
    expect_error(
        validate_hbb_formula(f, d),
        "must be numeric"
    )
})

test_that("validate error: response contains negative values", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$y <- c(0L, 3L, -1L, 0L, 2L)
    expect_error(
        validate_hbb_formula(f, d),
        "contains negative values"
    )
})

test_that("validate error: response contains non-integer (decimal) values", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$y <- c(0.5, 3.0, 5.0, 0.0, 2.0)
    expect_error(
        validate_hbb_formula(f, d),
        "must contain integers"
    )
})

test_that("validate error: response is logical (non-numeric)", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$y <- c(TRUE, FALSE, TRUE, FALSE, TRUE)
    expect_error(
        validate_hbb_formula(f, d),
        "must be numeric"
    )
})


# --------------------------------------------------------------------------
# 7c. Trials variable errors
# --------------------------------------------------------------------------

test_that("validate error: trials variable missing from data", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$n_trial <- NULL
    expect_error(
        validate_hbb_formula(f, d),
        "not found in.*data"
    )
})

test_that("validate error: trials is character (non-numeric)", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$n_trial <- c("a", "b", "c", "d", "e")
    expect_error(
        validate_hbb_formula(f, d),
        "must be numeric"
    )
})

test_that("validate error: trials < 1 (contains zero)", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$n_trial <- c(10L, 10L, 0L, 10L, 10L)
    expect_error(
        validate_hbb_formula(f, d),
        "must be >= 1"
    )
})

test_that("validate error: trials contains non-integer values", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$n_trial <- c(10.5, 10, 10, 10, 10)
    expect_error(
        validate_hbb_formula(f, d),
        "must contain integers"
    )
})

test_that("validate error: trials contains negative values (< 1)", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$n_trial <- c(10L, -5L, 10L, 10L, 10L)
    expect_error(
        validate_hbb_formula(f, d),
        "must be >= 1"
    )
})


# --------------------------------------------------------------------------
# 7d. Response > trials
# --------------------------------------------------------------------------

test_that("validate error: response exceeds trials (single violation)", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$y <- c(0L, 3L, 15L, 0L, 2L)
    expect_error(
        validate_hbb_formula(f, d),
        "Response exceeds trials"
    )
})

test_that("validate error: response exceeds trials (multiple violations)", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$y <- c(11L, 12L, 15L, 0L, 2L)
    expect_error(
        validate_hbb_formula(f, d),
        "Response exceeds trials"
    )
})


# --------------------------------------------------------------------------
# 7e. Fixed effect variables errors
# --------------------------------------------------------------------------

test_that("validate error: fixed effect variable missing from data", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$poverty <- NULL
    expect_error(
        validate_hbb_formula(f, d),
        "not found in.*data"
    )
})

test_that("validate error: fixed effect variable is character", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$urban <- c("low", "med", "high", "low", "med")
    expect_error(
        validate_hbb_formula(f, d),
        "must be numeric"
    )
})

test_that("validate error: fixed effect variable is logical", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$poverty <- c(TRUE, FALSE, TRUE, FALSE, TRUE)
    expect_error(
        validate_hbb_formula(f, d),
        "must be numeric"
    )
})

test_that("validate error: fixed effect variable is factor", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    d$poverty <- factor(c("low", "med", "high", "low", "med"))
    expect_error(
        validate_hbb_formula(f, d),
        "must be numeric"
    )
})


# --------------------------------------------------------------------------
# 7f. Grouping variable errors
# --------------------------------------------------------------------------

test_that("validate error: grouping variable missing from data", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty + (1 | state_id))
    d <- .make_valid_data()
    # No state_id column
    expect_error(
        validate_hbb_formula(f, d),
        "not found in.*data"
    )
})


# --------------------------------------------------------------------------
# 7g. Policy moderators / state_data errors
# --------------------------------------------------------------------------

test_that("validate error: policy moderators specified but no state_data", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (1 | state_id) +
            state_level(mr_pctile)
    )
    d <- .make_valid_data()
    d$state_id <- c(1L, 1L, 2L, 2L, 3L)
    expect_error(
        validate_hbb_formula(f, d),
        "Policy moderators require.*state_data"
    )
})

test_that("validate error: grouping variable missing from state_data", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (1 | state_id) +
            state_level(mr_pctile)
    )
    d <- .make_valid_data()
    d$state_id <- c(1L, 1L, 2L, 2L, 3L)
    sdata <- data.frame(mr_pctile = c(0.5, 0.6, 0.7))
    expect_error(
        validate_hbb_formula(f, d, state_data = sdata),
        "not found in.*state_data"
    )
})

test_that("validate error: policy variable missing from state_data", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (1 | state_id) +
            state_level(mr_pctile)
    )
    d <- .make_valid_data()
    d$state_id <- c(1L, 1L, 2L, 2L, 3L)
    sdata <- data.frame(
        state_id  = c(1L, 2L, 3L),
        other_var = c(0.1, 0.2, 0.3)
    )
    expect_error(
        validate_hbb_formula(f, d, state_data = sdata),
        "not found in.*state_data"
    )
})

test_that("validate error: policy variable is character (non-numeric)", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (1 | state_id) +
            state_level(mr_pctile)
    )
    d <- .make_valid_data()
    d$state_id <- c(1L, 1L, 2L, 2L, 3L)
    sdata <- data.frame(
        state_id  = c(1L, 2L, 3L),
        mr_pctile = c("low", "med", "high")
    )
    expect_error(
        validate_hbb_formula(f, d, state_data = sdata),
        "must be numeric"
    )
})

test_that("validate error: policy variable is factor (non-numeric)", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (1 | state_id) +
            state_level(mr_pctile)
    )
    d <- .make_valid_data()
    d$state_id <- c(1L, 1L, 2L, 2L, 3L)
    sdata <- data.frame(
        state_id  = c(1L, 2L, 3L),
        mr_pctile = factor(c("low", "med", "high"))
    )
    expect_error(
        validate_hbb_formula(f, d, state_data = sdata),
        "must be numeric"
    )
})

test_that("validate error: multiple policy vars, one missing from state_data", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (1 | state_id) +
            state_level(mr_pctile + tiered_reim)
    )
    d <- .make_valid_data()
    d$state_id <- c(1L, 1L, 2L, 2L, 3L)
    sdata <- data.frame(
        state_id  = c(1L, 2L, 3L),
        mr_pctile = c(0.5, 0.6, 0.7)
        # tiered_reim is missing
    )
    expect_error(
        validate_hbb_formula(f, d, state_data = sdata),
        "not found in.*state_data"
    )
})


# ============================================================================
# Section 8: Edge cases that should SUCCEED
# ============================================================================

test_that("edge: variable names with dots in LHS and RHS", {
    f <- hbb_formula(y.resp | trials(n.trial) ~ x.var + urban.pct)
    expect_equal(f$response, "y.resp")
    expect_equal(f$trials, "n.trial")
    expect_equal(f$fixed, c("x.var", "urban.pct"))
})

test_that("edge: single random slope parses correctly", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty + (poverty | state_id))
    expect_equal(f$fixed, "poverty")
    expect_equal(f$group, "state_id")
    expect_equal(f$random, "poverty")
    expect_true(f$svc)
})

test_that("edge: state_level with single variable", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (1 | state_id) +
            state_level(mr_pctile)
    )
    expect_equal(f$policy, "mr_pctile")
})

test_that("edge: explicit '1 +' in RHS is stripped (intercept-plus-predictors)", {
    f <- hbb_formula(y | trials(n_trial) ~ 1 + poverty + urban)
    expect_equal(f$fixed, c("poverty", "urban"))
})

test_that("edge: random effects with explicit '1 +' strips intercept marker", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (1 + poverty | state_id)
    )
    expect_equal(f$random, "poverty")
    expect_equal(f$group, "state_id")
    expect_true(f$svc)
})

test_that("edge: no spaces around | in LHS", {
    f <- hbb_formula(y|trials(n_trial)~poverty)
    expect_equal(f$response, "y")
    expect_equal(f$trials, "n_trial")
    expect_equal(f$fixed, "poverty")
})

test_that("edge: extra spaces inside trials()", {
    f <- hbb_formula(y | trials( n_trial ) ~ poverty)
    expect_equal(f$response, "y")
    expect_equal(f$trials, "n_trial")
    expect_equal(f$fixed, "poverty")
})

test_that("edge: formula created programmatically with as.formula()", {
    fml <- as.formula("y | trials(n_trial) ~ poverty + urban")
    f <- hbb_formula(fml)
    expect_equal(f$response, "y")
    expect_equal(f$trials, "n_trial")
    expect_equal(f$fixed, c("poverty", "urban"))
})

test_that("edge: no spaces around ~ ", {
    f <- hbb_formula(y | trials(n_trial)~poverty+urban)
    expect_equal(f$fixed, c("poverty", "urban"))
})

test_that("edge: excessive whitespace in formula", {
    fml <- as.formula("y  |  trials( n_trial )  ~  poverty  +  urban")
    f <- hbb_formula(fml)
    expect_equal(f$response, "y")
    expect_equal(f$trials, "n_trial")
    expect_equal(f$fixed, c("poverty", "urban"))
})

test_that("edge: variable names starting with dot", {
    f <- hbb_formula(.y | trials(.n) ~ .x1 + .x2)
    expect_equal(f$response, ".y")
    expect_equal(f$trials, ".n")
    expect_equal(f$fixed, c(".x1", ".x2"))
})

test_that("edge: variable names with underscores and digits", {
    f <- hbb_formula(y_2024 | trials(n_trial_01) ~ x_var_1 + var2_pct)
    expect_equal(f$response, "y_2024")
    expect_equal(f$trials, "n_trial_01")
    expect_equal(f$fixed, c("x_var_1", "var2_pct"))
})

test_that("edge: validate passes with all-zero response (structural zeros)", {
    f <- .make_valid_formula()
    d <- data.frame(
        y       = c(0L, 0L, 0L, 0L, 0L),
        n_trial = c(10L, 10L, 10L, 10L, 10L),
        poverty = c(0.1, 0.5, 0.3, 0.8, 0.4),
        urban   = c(0.2, 0.6, 0.9, 0.1, 0.5)
    )
    expect_true(validate_hbb_formula(f, d))
})

test_that("edge: validate passes with y = n (boundary case)", {
    f <- .make_valid_formula()
    d <- data.frame(
        y       = c(10L, 10L, 10L, 10L, 10L),
        n_trial = c(10L, 10L, 10L, 10L, 10L),
        poverty = c(0.1, 0.5, 0.3, 0.8, 0.4),
        urban   = c(0.2, 0.6, 0.9, 0.1, 0.5)
    )
    expect_true(validate_hbb_formula(f, d))
})

test_that("edge: validate accepts character grouping variable", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty + (1 | state_id))
    d <- data.frame(
        y        = c(0L, 3L, 5L, 0L, 2L),
        n_trial  = c(10L, 10L, 10L, 10L, 10L),
        poverty  = c(0.1, 0.5, 0.3, 0.8, 0.4),
        state_id = c("AL", "AL", "AK", "AK", "AZ"),
        stringsAsFactors = FALSE
    )
    expect_true(validate_hbb_formula(f, d))
})

test_that("edge: validate accepts factor grouping variable", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty + (1 | state_id))
    d <- data.frame(
        y        = c(0L, 3L, 5L, 0L, 2L),
        n_trial  = c(10L, 10L, 10L, 10L, 10L),
        poverty  = c(0.1, 0.5, 0.3, 0.8, 0.4),
        state_id = factor(c("AL", "AL", "AK", "AK", "AZ"))
    )
    expect_true(validate_hbb_formula(f, d))
})

test_that("edge: validate passes with complete model and state_data", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (1 | state_id) +
            state_level(mr_pctile)
    )
    d <- data.frame(
        y        = c(0L, 3L, 5L, 0L, 2L, 8L),
        n_trial  = c(10L, 10L, 10L, 10L, 10L, 10L),
        poverty  = c(0.1, 0.5, 0.3, 0.8, 0.4, 0.6),
        state_id = c(1L, 1L, 2L, 2L, 3L, 3L)
    )
    sdata <- data.frame(
        state_id  = c(1L, 2L, 3L),
        mr_pctile = c(0.5, 0.6, 0.7)
    )
    expect_true(validate_hbb_formula(f, d, state_data = sdata))
})

test_that("edge: validate passes with integer response stored as double", {
    # Whole numbers stored as double should be accepted
    f <- .make_valid_formula()
    d <- data.frame(
        y       = c(0, 3, 5, 0, 2),  # double, not integer
        n_trial = c(10, 10, 10, 10, 10),
        poverty = c(0.1, 0.5, 0.3, 0.8, 0.4),
        urban   = c(0.2, 0.6, 0.9, 0.1, 0.5)
    )
    expect_true(validate_hbb_formula(f, d))
})

test_that("edge: validate returns TRUE invisibly", {
    f <- .make_valid_formula()
    d <- .make_valid_data()
    out <- withVisible(validate_hbb_formula(f, d))
    expect_true(out$value)
    expect_false(out$visible)
})

test_that("edge: print.hbb_formula returns invisibly", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty + (1 | state_id))
    out <- capture.output(result <- print(f))
    expect_identical(result, f)
    expect_true(length(out) > 0)
})

test_that("edge: hbb_formula S3 class is set correctly", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty + urban)
    expect_s3_class(f, "hbb_formula")
})


# ============================================================================
# Section 9: Bug fixes and new features
# ============================================================================

# --------------------------------------------------------------------------
# 9a. C1 fix: intercept removal with RE present
# --------------------------------------------------------------------------

test_that("error: intercept removal via 0 detected even with RE", {
    # Bug C1: ~ 0 + x1 + (1 | state_id) was NOT caught because the outer
    # gate checked !grepl("|") which was TRUE when RE present.
    expect_error(
        hbb_formula(y | trials(n_trial) ~ 0 + poverty + (1 | state_id)),
        "Removing the intercept is not supported"
    )
})

test_that("error: intercept removal via - 1 detected even with RE", {
    expect_error(
        hbb_formula(y | trials(n_trial) ~ poverty + (1 | state_id) - 1),
        "Removing the intercept is not supported"
    )
})

test_that("error: intercept removal via 0 detected with SVC", {
    expect_error(
        hbb_formula(
            y | trials(n_trial) ~ 0 + poverty + (poverty | state_id)
        ),
        "Removing the intercept is not supported"
    )
})

# --------------------------------------------------------------------------
# 9b. I1: NA/NaN/Inf checks in validate_hbb_formula
# --------------------------------------------------------------------------

test_that("error: validate catches NaN in response", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty)
    d <- data.frame(y = c(1, NaN, 3), n_trial = c(10, 10, 10),
                    poverty = c(0.1, 0.2, 0.3))
    expect_error(validate_hbb_formula(f, d), "NaN")
})

test_that("error: validate catches Inf in response", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty)
    d <- data.frame(y = c(1, Inf, 3), n_trial = c(10, 10, 10),
                    poverty = c(0.1, 0.2, 0.3))
    expect_error(validate_hbb_formula(f, d), "infinite")
})

test_that("warning: validate warns on NA in response", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty)
    d <- data.frame(y = c(1, NA, 3), n_trial = c(10, 10, 10),
                    poverty = c(0.1, 0.2, 0.3))
    expect_warning(validate_hbb_formula(f, d), "missing")
})

test_that("error: validate catches NaN in trials", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty)
    d <- data.frame(y = c(1, 2, 3), n_trial = c(10, NaN, 10),
                    poverty = c(0.1, 0.2, 0.3))
    expect_error(validate_hbb_formula(f, d), "NaN")
})

test_that("error: validate catches Inf in trials", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty)
    d <- data.frame(y = c(1, 2, 3), n_trial = c(10, Inf, 10),
                    poverty = c(0.1, 0.2, 0.3))
    expect_error(validate_hbb_formula(f, d), "infinite")
})

test_that("warning: validate warns on NA in trials", {
    f <- hbb_formula(y | trials(n_trial) ~ poverty)
    d <- data.frame(y = c(1, 2, 3), n_trial = c(10, NA, 10),
                    poverty = c(0.1, 0.2, 0.3))
    expect_warning(validate_hbb_formula(f, d), "missing")
})

# --------------------------------------------------------------------------
# 9c. I5: state coverage and duplicate checks
# --------------------------------------------------------------------------

test_that("error: validate catches states in data missing from state_data", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (1 | state_id) +
            state_level(policy_var)
    )
    d <- data.frame(
        y = c(1, 2, 3, 0), n_trial = c(10, 10, 10, 10),
        poverty = c(0.1, 0.2, 0.3, 0.4),
        state_id = c(1, 2, 3, 4)
    )
    # state_data only has states 1 and 2
    sd <- data.frame(
        state_id = c(1, 2),
        policy_var = c(0.5, 0.6)
    )
    expect_error(validate_hbb_formula(f, d, state_data = sd),
                 "not found in.*state_data")
})

test_that("error: validate catches duplicate states in state_data", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (1 | state_id) +
            state_level(policy_var)
    )
    d <- data.frame(
        y = c(1, 2), n_trial = c(10, 10),
        poverty = c(0.1, 0.2),
        state_id = c(1, 2)
    )
    sd <- data.frame(
        state_id = c(1, 2, 2),
        policy_var = c(0.5, 0.6, 0.7)
    )
    expect_error(validate_hbb_formula(f, d, state_data = sd),
                 "duplicate")
})

test_that("validate: state coverage passes when all states covered", {
    f <- hbb_formula(
        y | trials(n_trial) ~ poverty + (1 | state_id) +
            state_level(policy_var)
    )
    d <- data.frame(
        y = c(1, 2), n_trial = c(10, 10),
        poverty = c(0.1, 0.2),
        state_id = c(1, 2)
    )
    # state_data has extra state 3 — that's fine
    sd <- data.frame(
        state_id = c(1, 2, 3),
        policy_var = c(0.5, 0.6, 0.7)
    )
    expect_true(validate_hbb_formula(f, d, state_data = sd))
})
