# ============================================================================
# test-data-prep.R --- Tests for prepare_stan_data(), validate_hbb_data(),
#                      and print.hbb_data()
#
# Covers:
#   Section 1: prepare_stan_data() basic functionality
#   Section 2: S3 structure and model_type
#   Section 3: Standardization
#   Section 4: Dimensional consistency
#   Section 5: Value constraints
#   Section 6: validate_hbb_data()
#   Section 7: Error paths
#   Section 8: print.hbb_data()
# ============================================================================


# -- Load package data once for all tests ------------------------------------
# Each test_that block loads its own copy via envir = environment() to ensure
# isolation, but we document the shared data dependency here.


# ============================================================================
# Section 1: prepare_stan_data() basic functionality
# ============================================================================

test_that("prepare_stan_data: intercept-only model", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ 1,
        data    = nsece_synth_small
    )

    expect_s3_class(d, "hbb_data")
    expect_equal(d$P, 1L)
    # X should be N x 1 matrix of all 1s (intercept only)
    expect_equal(ncol(d$X), 1L)
    expect_equal(nrow(d$X), d$N)
    expect_true(all(d$X[, 1] == 1))
    expect_equal(d$model_type, "base")
})


test_that("prepare_stan_data: fixed effects model", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban,
        data    = nsece_synth_small
    )

    expect_s3_class(d, "hbb_data")
    # Intercept + poverty + urban = 3
    expect_equal(d$P, 3L)
    expect_equal(ncol(d$X), 3L)
    expect_equal(nrow(d$X), d$N)
    # First column should always be intercept (all 1s)
    expect_true(all(d$X[, 1] == 1))
    # Remaining columns should be standardized by default
    for (j in 2:d$P) {
        expect_equal(mean(d$X[, j]), 0, tolerance = 1e-10,
                     label = paste("X column", j, "mean"))
        expect_equal(sd(d$X[, j]), 1, tolerance = 1e-10,
                     label = paste("X column", j, "sd"))
    }
})


test_that("prepare_stan_data: random intercept model", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + (1 | state_id),
        data    = nsece_synth_small
    )

    expect_s3_class(d, "hbb_data")
    expect_equal(d$S, 51L)
    # K = 2 * 1 = 2 for random intercept (one per margin: alpha, beta)
    expect_equal(d$K, 2L)
    # State indices should be integers in 1..S
    expect_true(all(d$state >= 1L))
    expect_true(all(d$state <= d$S))
    expect_equal(length(unique(d$state)), d$S)
})


test_that("prepare_stan_data: SVC model", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban +
            (poverty + urban | state_id),
        data    = nsece_synth_small
    )

    expect_s3_class(d, "hbb_data")
    # P = intercept + poverty + urban = 3
    expect_equal(d$P, 3L)
    expect_equal(d$S, 51L)
    # K = 2 * P = 6 (SVC: each of P coefficients has random effect in both margins)
    expect_equal(d$K, 2L * d$P)
})


test_that("prepare_stan_data: full model with policy moderators", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    data("nsece_state_policy", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban +
            (poverty + urban | state_id) +
            state_level(mr_pctile),
        data       = nsece_synth_small,
        state_data = nsece_state_policy
    )

    expect_s3_class(d, "hbb_data")
    # Q = intercept + mr_pctile = 2
    expect_equal(d$Q, 2L)
    # V should be S x Q
    expect_equal(nrow(d$V), d$S)
    expect_equal(ncol(d$V), d$Q)
    # V first column should be intercept (all 1s)
    expect_true(all(d$V[, 1] == 1))
})


test_that("prepare_stan_data: with survey weights", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban,
        data    = nsece_synth_small,
        weights = "weight"
    )

    expect_s3_class(d, "hbb_data")
    expect_equal(d$model_type, "weighted")
    # w_tilde should be present and normalized: sum = N
    expect_equal(length(d$w_tilde), d$N)
    expect_true(all(d$w_tilde > 0))
    expect_equal(sum(d$w_tilde), d$N, tolerance = 1e-10)
})


test_that("prepare_stan_data: with stratum and psu", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small,
        weights = "weight",
        stratum = "stratum",
        psu     = "psu"
    )

    expect_s3_class(d, "hbb_data")
    # stratum_idx should be sequential integers
    expect_true(all(d$stratum_idx >= 1L))
    expect_true(all(d$stratum_idx <= d$n_strata))
    expect_equal(length(d$stratum_idx), d$N)
    # psu_idx should be sequential integers
    expect_true(all(d$psu_idx >= 1L))
    expect_true(all(d$psu_idx <= d$n_psu))
    expect_equal(length(d$psu_idx), d$N)
    # n_strata and n_psu should be positive integers
    expect_true(d$n_strata >= 1L)
    expect_true(d$n_psu >= 1L)
    expect_true(d$n_psu >= d$n_strata)
})


test_that("prepare_stan_data: raw formula auto-conversion", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    # Passing a raw formula (not an hbb_formula object) should work via
    # auto-conversion inside prepare_stan_data
    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban,
        data    = nsece_synth_small
    )

    expect_s3_class(d, "hbb_data")
    expect_equal(d$P, 3L)
    expect_equal(d$N, nrow(nsece_synth_small))
})


# ============================================================================
# Section 2: S3 structure and model_type
# ============================================================================

test_that("prepare_stan_data: class is 'hbb_data'", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small
    )

    expect_s3_class(d, "hbb_data")
    expect_true(inherits(d, "hbb_data"))
})


test_that("prepare_stan_data: all expected fields present", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    data("nsece_state_policy", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban +
            (poverty + urban | state_id) +
            state_level(mr_pctile),
        data       = nsece_synth_small,
        weights    = "weight",
        stratum    = "stratum",
        psu        = "psu",
        state_data = nsece_state_policy
    )

    expected_fields <- c(
        "N", "P", "S", "Q", "N_pos", "K",
        "y", "n_trial", "z",
        "X", "state", "V", "idx_pos",
        "w_tilde", "stratum_idx", "psu_idx", "n_strata", "n_psu",
        "x_center", "x_scale", "v_center", "v_scale",
        "formula", "prior", "group_levels", "model_type"
    )
    for (field in expected_fields) {
        expect_true(field %in% names(d),
                    info = paste("Missing field:", field))
    }
})


test_that("prepare_stan_data: model_type = 'base'", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + (1 | state_id),
        data    = nsece_synth_small
    )

    expect_equal(d$model_type, "base")
})


test_that("prepare_stan_data: model_type = 'weighted'", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + (1 | state_id),
        data    = nsece_synth_small,
        weights = "weight"
    )

    expect_equal(d$model_type, "weighted")
})


test_that("prepare_stan_data: model_type = 'svc'", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban +
            (poverty + urban | state_id),
        data    = nsece_synth_small
    )

    expect_equal(d$model_type, "svc")
})


test_that("prepare_stan_data: model_type = 'svc_weighted'", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban +
            (poverty + urban | state_id),
        data    = nsece_synth_small,
        weights = "weight"
    )

    expect_equal(d$model_type, "svc_weighted")
})


# ============================================================================
# Section 3: Standardization
# ============================================================================

test_that("prepare_stan_data: center=TRUE, scale=TRUE standardizes covariates", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban + black + hispanic,
        data    = nsece_synth_small,
        center  = TRUE,
        scale   = TRUE
    )

    # Non-intercept columns should have mean ~0 and sd ~1
    for (j in 2:d$P) {
        expect_equal(mean(d$X[, j]), 0, tolerance = 1e-10,
                     label = paste("Column", j, "mean"))
        expect_equal(sd(d$X[, j]), 1, tolerance = 1e-10,
                     label = paste("Column", j, "sd"))
    }
})


test_that("prepare_stan_data: center=TRUE, scale=FALSE centers only", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban,
        data    = nsece_synth_small,
        center  = TRUE,
        scale   = FALSE
    )

    # Non-intercept columns should have mean ~0 but NOT sd ~1
    for (j in 2:d$P) {
        expect_equal(mean(d$X[, j]), 0, tolerance = 1e-10,
                     label = paste("Column", j, "mean"))
    }
    # sd should reflect the original data scale (generally not 1)
    # poverty has sd >> 1 in percentage units
    expect_true(sd(d$X[, 2]) > 1,
                info = "Poverty column should retain original scale")
})


test_that("prepare_stan_data: center=FALSE, scale=FALSE preserves raw values", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban,
        data    = nsece_synth_small,
        center  = FALSE,
        scale   = FALSE
    )

    # Non-intercept columns should match the raw data
    expect_equal(d$X[, 2], nsece_synth_small$poverty, tolerance = 1e-14)
    expect_equal(d$X[, 3], nsece_synth_small$urban, tolerance = 1e-14)
})


test_that("prepare_stan_data: x_center and x_scale stored correctly", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban,
        data    = nsece_synth_small,
        center  = TRUE,
        scale   = TRUE
    )

    # x_center should store the column means of the original data
    expect_equal(unname(d$x_center[1]), mean(nsece_synth_small$poverty),
                 tolerance = 1e-10)
    expect_equal(unname(d$x_center[2]), mean(nsece_synth_small$urban),
                 tolerance = 1e-10)

    # x_scale should store the column SDs of the original data
    expect_equal(unname(d$x_scale[1]), sd(nsece_synth_small$poverty),
                 tolerance = 1e-10)
    expect_equal(unname(d$x_scale[2]), sd(nsece_synth_small$urban),
                 tolerance = 1e-10)
})


test_that("prepare_stan_data: intercept column is NOT standardized", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban,
        data    = nsece_synth_small,
        center  = TRUE,
        scale   = TRUE
    )

    # Column 1 (intercept) must always be 1
    expect_true(all(d$X[, 1] == 1))
})


test_that("prepare_stan_data: V matrix policy vars standardized when continuous", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    data("nsece_state_policy", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban +
            (poverty + urban | state_id) +
            state_level(mr_pctile),
        data       = nsece_synth_small,
        state_data = nsece_state_policy,
        center     = TRUE,
        scale      = TRUE
    )

    # V column 1 is intercept (all 1s)
    expect_true(all(d$V[, 1] == 1))
    # V column 2 (mr_pctile) should be centered and scaled
    # mr_pctile is already approximately standardized, so after
    # centering/scaling it should be close to mean=0, sd=1
    expect_equal(mean(d$V[, 2]), 0, tolerance = 1e-10)
    expect_equal(sd(d$V[, 2]), 1, tolerance = 1e-10)
})


# ============================================================================
# Section 4: Dimensional consistency
# ============================================================================

test_that("prepare_stan_data: X dimensions are N x P", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban + black,
        data    = nsece_synth_small
    )

    expect_equal(nrow(d$X), d$N)
    expect_equal(ncol(d$X), d$P)
})


test_that("prepare_stan_data: V dimensions are S x Q", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    data("nsece_state_policy", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty +
            (poverty | state_id) +
            state_level(mr_pctile + tiered_reim),
        data       = nsece_synth_small,
        state_data = nsece_state_policy
    )

    expect_equal(nrow(d$V), d$S)
    expect_equal(ncol(d$V), d$Q)
    # Q = intercept + mr_pctile + tiered_reim = 3
    expect_equal(d$Q, 3L)
})


test_that("prepare_stan_data: y, n_trial, z lengths equal N", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small
    )

    expect_equal(length(d$y), d$N)
    expect_equal(length(d$n_trial), d$N)
    expect_equal(length(d$z), d$N)
})


test_that("prepare_stan_data: state length is N, range 1..S", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + (1 | state_id),
        data    = nsece_synth_small
    )

    expect_equal(length(d$state), d$N)
    expect_true(all(d$state >= 1L))
    expect_true(all(d$state <= d$S))
})


test_that("prepare_stan_data: idx_pos length equals N_pos", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small
    )

    expect_equal(length(d$idx_pos), d$N_pos)
    expect_equal(d$N_pos, sum(d$z == 1))
})


test_that("prepare_stan_data: w_tilde length is N when weights present", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small,
        weights = "weight"
    )

    expect_equal(length(d$w_tilde), d$N)
})


# ============================================================================
# Section 5: Value constraints
# ============================================================================

test_that("prepare_stan_data: y >= 0", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small
    )

    expect_true(all(d$y >= 0))
})


test_that("prepare_stan_data: n_trial >= 1", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small
    )

    expect_true(all(d$n_trial >= 1))
})


test_that("prepare_stan_data: y <= n_trial", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small
    )

    expect_true(all(d$y <= d$n_trial))
})


test_that("prepare_stan_data: z = I(y > 0)", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small
    )

    expect_equal(d$z, as.integer(d$y > 0))
})


test_that("prepare_stan_data: state indices are in 1..S", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + (1 | state_id),
        data    = nsece_synth_small
    )

    expect_true(all(d$state %in% 1:d$S))
})


test_that("prepare_stan_data: idx_pos points to z=1 rows", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small
    )

    # idx_pos should be integer indices where z == 1
    expect_equal(sort(d$idx_pos), sort(which(d$z == 1)))
    # All pointed-to rows should have z=1
    expect_true(all(d$z[d$idx_pos] == 1))
})


test_that("prepare_stan_data: w_tilde > 0 and sum approx N", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small,
        weights = "weight"
    )

    expect_true(all(d$w_tilde > 0))
    expect_equal(sum(d$w_tilde), d$N, tolerance = 1e-10)
})


# ============================================================================
# Section 6: validate_hbb_data()
# ============================================================================

test_that("validate_hbb_data: passes on valid hbb_data", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban +
            (1 | state_id),
        data    = nsece_synth_small,
        weights = "weight"
    )

    expect_true(validate_hbb_data(d))
})


test_that("validate_hbb_data: errors on missing fields", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small
    )

    # Remove a required field
    d_broken <- d
    d_broken$N <- NULL

    expect_error(validate_hbb_data(d_broken))
})


test_that("validate_hbb_data: errors on dimensional mismatch", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small
    )

    # Corrupt X dimensions
    d_broken <- d
    d_broken$X <- d_broken$X[1:10, , drop = FALSE]

    expect_error(validate_hbb_data(d_broken))
})


test_that("validate_hbb_data: errors on value violations", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small
    )

    # Corrupt y to have negative values
    d_broken <- d
    d_broken$y[1] <- -1L

    expect_error(validate_hbb_data(d_broken))
})


# ============================================================================
# Section 7: Error paths
# ============================================================================

test_that("prepare_stan_data: errors on non-data-frame input", {
    expect_error(
        prepare_stan_data(
            formula = y | trials(n_trial) ~ poverty,
            data    = list(y = 1:10, n_trial = rep(20, 10), poverty = rnorm(10))
        ),
        regexp = "data\\.frame|data frame"
    )
})


test_that("prepare_stan_data: errors on missing columns", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    # Specify a predictor that does not exist
    expect_error(
        prepare_stan_data(
            formula = y | trials(n_trial) ~ nonexistent_var,
            data    = nsece_synth_small
        ),
        regexp = "nonexistent_var|not found"
    )
})


test_that("prepare_stan_data: errors on invalid weights column name", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    expect_error(
        prepare_stan_data(
            formula = y | trials(n_trial) ~ poverty,
            data    = nsece_synth_small,
            weights = "nonexistent_weight"
        ),
        regexp = "nonexistent_weight|not found|weight"
    )
})


test_that("prepare_stan_data: errors on negative weights", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    bad_data <- nsece_synth_small
    bad_data$weight[1] <- -1

    expect_error(
        prepare_stan_data(
            formula = y | trials(n_trial) ~ poverty,
            data    = bad_data,
            weights = "weight"
        ),
        regexp = "negative|positive|weight"
    )
})


test_that("prepare_stan_data: errors when state_data lacks group variable", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    data("nsece_state_policy", package = "hurdlebb", envir = environment())

    # Remove state_id from state_data
    bad_state_data <- nsece_state_policy
    bad_state_data$state_id <- NULL

    expect_error(
        prepare_stan_data(
            formula = y | trials(n_trial) ~ poverty +
                (1 | state_id) +
                state_level(mr_pctile),
            data       = nsece_synth_small,
            state_data = bad_state_data
        ),
        regexp = "state_id|not found|group"
    )
})


# ============================================================================
# Section 8: print.hbb_data()
# ============================================================================

test_that("print.hbb_data: prints without error", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban,
        data    = nsece_synth_small,
        weights = "weight"
    )

    expect_output(print(d))
})


test_that("print.hbb_data: output contains key info", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())

    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty + urban,
        data    = nsece_synth_small,
        weights = "weight"
    )

    out <- capture.output(print(d))
    combined <- paste(out, collapse = "\n")

    # Should mention N, P, S or equivalent dimension info
    expect_true(grepl(as.character(d$N), combined),
                info = "N not found in print output")
    expect_true(grepl(as.character(d$P), combined),
                info = "P not found in print output")
    # Should mention model_type
    expect_true(grepl(d$model_type, combined, ignore.case = TRUE),
                info = "model_type not found in print output")
})


# ============================================================================
# Section 9: Additional Coverage Tests --- validate_hbb_data error paths
# ============================================================================

# Helper: build a minimal valid hbb_data object for corruption testing
.make_valid_hbb_data <- function(N = 20, P = 2, S = 3) {
    set.seed(42)
    y       <- c(rep(0L, 8), sample(1:10, N - 8, replace = TRUE))
    n_trial <- rep(20L, N)
    z       <- as.integer(y > 0L)
    X       <- cbind(1, matrix(rnorm(N * (P - 1)), N, P - 1))
    colnames(X) <- c("(Intercept)", paste0("x", seq_len(P - 1)))
    state   <- rep(seq_len(S), length.out = N)

    f <- hbb_formula(y | trials(n_trial) ~ x1)

    structure(
        list(
            N = N, P = P, S = S, Q = 1L, N_pos = sum(z), K = 0L,
            y = y, n_trial = n_trial, z = z,
            X = X, state = state, V = NULL,
            idx_pos = which(z == 1L),
            w_tilde = NULL, stratum_idx = NULL, psu_idx = NULL,
            n_strata = NULL, n_psu = NULL,
            x_center = NULL, x_scale = NULL,
            v_center = NULL, v_scale = NULL,
            formula = f, prior = NULL,
            group_levels = as.character(seq_len(S)),
            model_type = "base"
        ),
        class = "hbb_data"
    )
}


test_that("validate_hbb_data: errors on non-hbb_data class (line 279-283)", {
    bad <- list(N = 10, P = 2)
    expect_error(validate_hbb_data(bad), "hbb_data")
})


test_that("validate_hbb_data: errors on y length mismatch (line 309)", {
    d <- .make_valid_hbb_data()
    d$y <- d$y[1:5]   # corrupt y length
    expect_error(validate_hbb_data(d), "y")
})


test_that("validate_hbb_data: errors on n_trial length mismatch (line 312)", {
    d <- .make_valid_hbb_data()
    d$n_trial <- d$n_trial[1:5]
    expect_error(validate_hbb_data(d), "n_trial")
})


test_that("validate_hbb_data: errors on z length mismatch (line 315)", {
    d <- .make_valid_hbb_data()
    d$z <- d$z[1:5]
    expect_error(validate_hbb_data(d), "z")
})


test_that("validate_hbb_data: errors when X is not a matrix (line 320)", {
    d <- .make_valid_hbb_data()
    d$X <- as.data.frame(d$X)
    expect_error(validate_hbb_data(d), "matrix")
})


test_that("validate_hbb_data: errors on state length mismatch (line 330)", {
    d <- .make_valid_hbb_data()
    d$state <- d$state[1:5]
    expect_error(validate_hbb_data(d), "state")
})


test_that("validate_hbb_data: errors when V is not a matrix (line 336)", {
    d <- .make_valid_hbb_data(S = 3)
    d$V <- data.frame(x = 1:3)   # not NULL, but not a matrix
    expect_error(validate_hbb_data(d), "matrix")
})


test_that("validate_hbb_data: errors on V dimension mismatch (line 339-341)", {
    d <- .make_valid_hbb_data(S = 3)
    d$V <- matrix(1, nrow = 5, ncol = 1)   # wrong nrow
    d$Q <- 1L
    expect_error(validate_hbb_data(d), "V")
})


test_that("validate_hbb_data: errors on idx_pos / N_pos mismatch (line 347-349)", {
    d <- .make_valid_hbb_data()
    d$N_pos <- d$N_pos + 5L   # corrupt
    expect_error(validate_hbb_data(d), "idx_pos")
})


test_that("validate_hbb_data: errors on n_trial < 1 (line 361)", {
    d <- .make_valid_hbb_data()
    d$n_trial[1] <- 0L
    expect_error(validate_hbb_data(d), "n_trial")
})


test_that("validate_hbb_data: errors on y > n_trial (line 366)", {
    d <- .make_valid_hbb_data()
    d$y[1]       <- 25L
    d$n_trial[1] <- 20L
    d$z[1]       <- 1L
    # Recalculate idx_pos and N_pos to keep them consistent
    d$idx_pos <- which(d$z == 1L)
    d$N_pos   <- length(d$idx_pos)
    expect_error(validate_hbb_data(d), "exceeds")
})


test_that("validate_hbb_data: errors on z not in {0,1} (line 371)", {
    d <- .make_valid_hbb_data()
    d$z[1] <- 2L
    expect_error(validate_hbb_data(d), "z")
})


test_that("validate_hbb_data: errors on state out of range (line 376)", {
    d <- .make_valid_hbb_data(S = 3)
    d$state[1] <- 99L
    expect_error(validate_hbb_data(d), "state")
})


test_that("validate_hbb_data: errors on z != (y > 0) (line 384)", {
    d <- .make_valid_hbb_data()
    # Set y > 0 but z = 0 (inconsistent)
    d$y[1] <- 5L
    d$z[1] <- 0L
    expect_error(validate_hbb_data(d), "z")
})


test_that("validate_hbb_data: errors on idx_pos pointing to z=0 (line 390)", {
    d <- .make_valid_hbb_data()
    # Corrupt idx_pos to include a zero-row
    zero_idx <- which(d$z == 0L)[1]
    if (!is.na(zero_idx)) {
        d$idx_pos <- c(d$idx_pos, zero_idx)
        d$N_pos   <- length(d$idx_pos)
    }
    expect_error(validate_hbb_data(d), "idx_pos|z")
})


test_that("validate_hbb_data: errors on idx_pos not matching which(z==1) (line 395)", {
    d <- .make_valid_hbb_data()
    # Reverse idx_pos (correct values, wrong order)
    d$idx_pos <- rev(d$idx_pos)
    expect_error(validate_hbb_data(d), "idx_pos")
})


test_that("validate_hbb_data: errors on w_tilde length mismatch (line 402-404)", {
    d <- .make_valid_hbb_data()
    d$w_tilde <- rep(1, 5)   # wrong length
    expect_error(validate_hbb_data(d), "w_tilde")
})


test_that("validate_hbb_data: errors on w_tilde <= 0 (line 407)", {
    d <- .make_valid_hbb_data()
    d$w_tilde <- rep(1, d$N)
    d$w_tilde[1] <- -0.5
    expect_error(validate_hbb_data(d), "w_tilde")
})


test_that("validate_hbb_data: errors on stratum_idx length mismatch (line 413-415)", {
    d <- .make_valid_hbb_data()
    d$stratum_idx <- c(1L, 1L, 2L)  # wrong length
    d$n_strata    <- 2L
    expect_error(validate_hbb_data(d), "stratum_idx")
})


test_that("validate_hbb_data: errors when n_strata is NULL but stratum_idx present (line 418)", {
    d <- .make_valid_hbb_data()
    d$stratum_idx <- rep(1L, d$N)
    d$n_strata    <- NULL
    expect_error(validate_hbb_data(d), "n_strata")
})


test_that("validate_hbb_data: errors on stratum_idx out of range (line 422)", {
    d <- .make_valid_hbb_data()
    d$stratum_idx    <- rep(1L, d$N)
    d$stratum_idx[1] <- 99L
    d$n_strata       <- 2L
    expect_error(validate_hbb_data(d), "stratum_idx")
})


test_that("validate_hbb_data: errors on psu_idx length mismatch (line 428-430)", {
    d <- .make_valid_hbb_data()
    d$psu_idx <- c(1L, 2L, 3L)  # wrong length
    d$n_psu   <- 3L
    expect_error(validate_hbb_data(d), "psu_idx")
})


test_that("validate_hbb_data: errors when n_psu is NULL but psu_idx present (line 433)", {
    d <- .make_valid_hbb_data()
    d$psu_idx <- rep(1L, d$N)
    d$n_psu   <- NULL
    expect_error(validate_hbb_data(d), "n_psu")
})


test_that("validate_hbb_data: errors on psu_idx out of range (line 437)", {
    d <- .make_valid_hbb_data()
    d$psu_idx    <- rep(1L, d$N)
    d$psu_idx[1] <- 99L
    d$n_psu      <- 2L
    expect_error(validate_hbb_data(d), "psu_idx")
})


test_that("validate_hbb_data: errors on invalid model_type (line 444-446)", {
    d <- .make_valid_hbb_data()
    d$model_type <- "invalid_type"
    expect_error(validate_hbb_data(d), "model_type")
})


# ============================================================================
# Section 10: prepare_stan_data formula class error path (line 136-140)
# ============================================================================

test_that("prepare_stan_data: errors on non-formula non-hbb_formula input (line 136-140)", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    expect_error(
        prepare_stan_data(formula = 42, data = nsece_synth_small),
        "hbb_formula|formula"
    )
    expect_error(
        prepare_stan_data(formula = "not_a_formula", data = nsece_synth_small),
        "hbb_formula|formula"
    )
})


# ============================================================================
# Section 11: print.hbb_data additional branches (lines 481, 484, 487, 494)
# ============================================================================

test_that("print.hbb_data: no-weights branch prints 'no' (line 481)", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small
    )
    out <- capture.output(print(d))
    combined <- paste(out, collapse = "\n")
    expect_true(grepl("Survey weights.*no", combined))
})


test_that("print.hbb_data: strata and PSU lines printed (lines 484, 487)", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small,
        weights = "weight",
        stratum = "stratum",
        psu     = "psu"
    )
    out <- capture.output(print(d))
    combined <- paste(out, collapse = "\n")
    expect_true(grepl("Strata", combined))
    expect_true(grepl("PSU", combined))
})


test_that("print.hbb_data: X standardised = no when center/scale FALSE (line 494)", {
    data("nsece_synth_small", package = "hurdlebb", envir = environment())
    d <- prepare_stan_data(
        formula = y | trials(n_trial) ~ poverty,
        data    = nsece_synth_small,
        center  = FALSE,
        scale   = FALSE
    )
    out <- capture.output(print(d))
    combined <- paste(out, collapse = "\n")
    expect_true(grepl("X standardised.*no", combined))
})


# ============================================================================
# Section 12: .build_V_matrix error path (lines 603-606)
# ============================================================================

test_that(".build_V_matrix: errors when groups not found in state_data (lines 603-606)", {
    # Directly test the internal helper
    # group_levels has a level not in state_data
    state_data <- data.frame(grp = c("A", "B"), mr = c(1, 2))
    expect_error(
        hurdlebb:::.build_V_matrix(
            policy = "mr", group = "grp",
            state_data = state_data,
            group_levels = c("A", "B", "C"),  # "C" is not in state_data
            S = 3L, center = TRUE, scale = TRUE
        ),
        "Not all groups"
    )
})


# ============================================================================
# Section 13: .build_group_index factor input (line 668)
# ============================================================================

test_that(".build_group_index: handles pre-existing factor grouping variable (line 668)", {
    df <- data.frame(
        x  = 1:6,
        grp = factor(c("B", "A", "C", "B", "A", "C"), levels = c("A", "B", "C"))
    )
    result <- hurdlebb:::.build_group_index("grp", df)
    # Factor levels should be preserved
    expect_equal(result$group_levels, c("A", "B", "C"))
    expect_equal(result$S, 3L)
    # row 1 is "B" = level 2, row 2 is "A" = level 1, etc.
    expect_equal(result$state_idx, c(2L, 1L, 3L, 2L, 1L, 3L))
})


# ============================================================================
# Section 14: .build_survey_design error paths (lines 717-720, 723, 734-736)
# ============================================================================

test_that(".build_survey_design: errors on missing weight column (line 710-713)", {
    df <- data.frame(x = 1:5)
    expect_error(
        hurdlebb:::.build_survey_design(df, weights = "nonexistent",
                                         stratum = NULL, psu = NULL),
        "not found"
    )
})


test_that(".build_survey_design: errors on non-numeric weights (lines 717-720)", {
    df <- data.frame(x = 1:5, wt = c("a", "b", "c", "d", "e"))
    expect_error(
        hurdlebb:::.build_survey_design(df, weights = "wt",
                                         stratum = NULL, psu = NULL),
        "numeric"
    )
})


test_that(".build_survey_design: errors on NA weights (line 723)", {
    df <- data.frame(x = 1:5, wt = c(1, 2, NA, 4, 5))
    expect_error(
        hurdlebb:::.build_survey_design(df, weights = "wt",
                                         stratum = NULL, psu = NULL),
        "missing"
    )
})


test_that(".build_survey_design: errors on missing stratum column (lines 734-736)", {
    df <- data.frame(x = 1:5)
    expect_error(
        hurdlebb:::.build_survey_design(df, weights = NULL,
                                         stratum = "nonexistent", psu = NULL),
        "not found"
    )
})


test_that(".build_survey_design: errors on missing PSU column", {
    df <- data.frame(x = 1:5)
    expect_error(
        hurdlebb:::.build_survey_design(df, weights = NULL,
                                         stratum = NULL, psu = "nonexistent"),
        "not found"
    )
})
