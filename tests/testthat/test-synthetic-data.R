# ============================================================================
# test-synthetic-data.R --- Tests for Synthetic NSECE Datasets
#
# Covers:
#   A. nsece_synth         — dimensions, types, constraints, statistics
#   B. nsece_synth_small   — dimensions, types, constraints, state coverage
#   C. nsece_state_policy  — dimensions, types, values
#   D. Cross-dataset consistency
# ============================================================================


# ============================================================================
# A. nsece_synth — Main synthetic dataset
# ============================================================================

test_that("nsece_synth loads without error", {
  expect_no_error(data("nsece_synth", package = "hurdlebb", envir = environment()))
  expect_true(exists("nsece_synth", envir = environment()))
})


test_that("nsece_synth has correct dimensions", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  # Approximately 6,785 rows (allow some flexibility for generation)
  expect_true(nrow(nsece_synth) >= 6000)
  expect_true(nrow(nsece_synth) <= 7500)
  # Exactly 13 columns

  expect_equal(ncol(nsece_synth), 13L)
})


test_that("nsece_synth has correct column names", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  expected_names <- c("provider_id", "state_id", "y", "n_trial", "z",
                      "it_share", "poverty", "urban", "black", "hispanic",
                      "weight", "stratum", "psu")
  expect_equal(sort(names(nsece_synth)), sort(expected_names))
})


test_that("nsece_synth has correct column types", {
  data("nsece_synth", package = "hurdlebb", envir = environment())

  # Integer columns
  expect_true(is.integer(nsece_synth$provider_id) ||
              is.numeric(nsece_synth$provider_id))
  expect_true(is.integer(nsece_synth$state_id) ||
              is.numeric(nsece_synth$state_id))
  expect_true(is.integer(nsece_synth$y) ||
              is.numeric(nsece_synth$y))
  expect_true(is.integer(nsece_synth$n_trial) ||
              is.numeric(nsece_synth$n_trial))
  expect_true(is.integer(nsece_synth$z))
  expect_true(is.integer(nsece_synth$stratum) ||
              is.numeric(nsece_synth$stratum))
  expect_true(is.integer(nsece_synth$psu) ||
              is.numeric(nsece_synth$psu))

  # Numeric columns
  expect_true(is.numeric(nsece_synth$it_share))
  expect_true(is.numeric(nsece_synth$poverty))
  expect_true(is.numeric(nsece_synth$urban))
  expect_true(is.numeric(nsece_synth$black))
  expect_true(is.numeric(nsece_synth$hispanic))
  expect_true(is.numeric(nsece_synth$weight))
})


test_that("nsece_synth has no missing values", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  for (col in names(nsece_synth)) {
    expect_false(anyNA(nsece_synth[[col]]),
                 info = paste("Column", col, "has NAs"))
  }
})


test_that("nsece_synth: y is non-negative and bounded by n_trial", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  expect_true(all(nsece_synth$y >= 0))
  expect_true(all(nsece_synth$y <= nsece_synth$n_trial))
})


test_that("nsece_synth: n_trial >= 1", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  expect_true(all(nsece_synth$n_trial >= 1))
})


test_that("nsece_synth: z == I(y > 0)", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  expect_equal(nsece_synth$z, as.integer(nsece_synth$y > 0))
})


test_that("nsece_synth: it_share is consistent with y / n_trial", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  expected_share <- nsece_synth$y / nsece_synth$n_trial
  expect_equal(nsece_synth$it_share, expected_share,
               tolerance = 1e-10)
})


test_that("nsece_synth: it_share is zero when z is zero", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  zero_mask <- nsece_synth$z == 0L
  if (any(zero_mask)) {
    expect_true(all(nsece_synth$it_share[zero_mask] == 0))
  }
})


test_that("nsece_synth: weights are strictly positive", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  expect_true(all(nsece_synth$weight > 0))
})


test_that("nsece_synth: state_id in 1:51", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  expect_true(all(nsece_synth$state_id %in% 1:51))
})


test_that("nsece_synth: exactly 51 unique state_ids", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  expect_equal(length(unique(nsece_synth$state_id)), 51L)
})


test_that("nsece_synth: z is binary (0 or 1 only)", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  expect_true(all(nsece_synth$z %in% c(0L, 1L)))
})


test_that("nsece_synth: provider_id values are unique", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  expect_equal(length(unique(nsece_synth$provider_id)), nrow(nsece_synth))
})


test_that("nsece_synth: y values are whole numbers", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  expect_equal(nsece_synth$y, round(nsece_synth$y))
})


test_that("nsece_synth: n_trial values are whole numbers", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  expect_equal(nsece_synth$n_trial, round(nsece_synth$n_trial))
})


# --- Statistical properties ---

test_that("nsece_synth: zero rate is between 0.25 and 0.45", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  zero_rate <- 1 - mean(nsece_synth$z)
  expect_true(zero_rate >= 0.25,
              info = paste("Zero rate too low:", round(zero_rate, 3)))
  expect_true(zero_rate <= 0.45,
              info = paste("Zero rate too high:", round(zero_rate, 3)))
})


test_that("nsece_synth: mean IT share among servers between 0.35 and 0.60", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  servers <- nsece_synth$z == 1L
  expect_true(any(servers), info = "No servers found")
  mean_share <- mean(nsece_synth$it_share[servers])
  expect_true(mean_share >= 0.35,
              info = paste("Mean IT share too low:", round(mean_share, 3)))
  expect_true(mean_share <= 0.60,
              info = paste("Mean IT share too high:", round(mean_share, 3)))
})


test_that("nsece_synth: multiple strata (> 10)", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  n_strata <- length(unique(nsece_synth$stratum))
  expect_true(n_strata > 10,
              info = paste("Only", n_strata, "strata"))
})


test_that("nsece_synth: multiple PSUs per stratum", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  psu_per_stratum <- tapply(nsece_synth$psu, nsece_synth$stratum,
                            function(x) length(unique(x)))
  # Most strata should have > 1 PSU (allow a few singletons)
  expect_true(median(psu_per_stratum) > 1,
              info = paste("Median PSUs per stratum:",
                           median(psu_per_stratum)))
})


test_that("nsece_synth: covariates have reasonable ranges", {
  data("nsece_synth", package = "hurdlebb", envir = environment())

  # Poverty: percentage, should be roughly 0-60
  expect_true(min(nsece_synth$poverty) >= 0)
  expect_true(max(nsece_synth$poverty) <= 100)

  # Urban: percentage, 0-100
  expect_true(min(nsece_synth$urban) >= 0)
  expect_true(max(nsece_synth$urban) <= 100)

  # Black: percentage, 0-100
  expect_true(min(nsece_synth$black) >= 0)
  expect_true(max(nsece_synth$black) <= 100)

  # Hispanic: percentage, 0-100
  expect_true(min(nsece_synth$hispanic) >= 0)
  expect_true(max(nsece_synth$hispanic) <= 100)
})


test_that("nsece_synth: covariates have nonzero variance", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  for (col in c("poverty", "urban", "black", "hispanic")) {
    expect_true(var(nsece_synth[[col]]) > 0,
                info = paste("Column", col, "has zero variance"))
  }
})


# ============================================================================
# B. nsece_synth_small — Small stratified subsample
# ============================================================================

test_that("nsece_synth_small loads without error", {
  expect_no_error(data("nsece_synth_small", package = "hurdlebb",
                       envir = environment()))
  expect_true(exists("nsece_synth_small", envir = environment()))
})


test_that("nsece_synth_small has correct dimensions", {
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  # Between 400 and 600 rows
  expect_true(nrow(nsece_synth_small) >= 400)
  expect_true(nrow(nsece_synth_small) <= 600)
  # Same 13 columns as nsece_synth
  expect_equal(ncol(nsece_synth_small), 13L)
})


test_that("nsece_synth_small has same column names as nsece_synth", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  expect_equal(sort(names(nsece_synth_small)), sort(names(nsece_synth)))
})


test_that("nsece_synth_small has no missing values", {
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  for (col in names(nsece_synth_small)) {
    expect_false(anyNA(nsece_synth_small[[col]]),
                 info = paste("Column", col, "has NAs"))
  }
})


test_that("nsece_synth_small: y is non-negative and bounded by n_trial", {
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  expect_true(all(nsece_synth_small$y >= 0))
  expect_true(all(nsece_synth_small$y <= nsece_synth_small$n_trial))
})


test_that("nsece_synth_small: n_trial >= 1", {
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  expect_true(all(nsece_synth_small$n_trial >= 1))
})


test_that("nsece_synth_small: z == I(y > 0)", {
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  expect_equal(nsece_synth_small$z, as.integer(nsece_synth_small$y > 0))
})


test_that("nsece_synth_small: it_share is consistent with y / n_trial", {
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  expected_share <- nsece_synth_small$y / nsece_synth_small$n_trial
  expect_equal(nsece_synth_small$it_share, expected_share,
               tolerance = 1e-10)
})


test_that("nsece_synth_small: weights are strictly positive", {
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  expect_true(all(nsece_synth_small$weight > 0))
})


test_that("nsece_synth_small: state_id in 1:51", {
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  expect_true(all(nsece_synth_small$state_id %in% 1:51))
})


test_that("nsece_synth_small: all 51 state_ids present", {
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  expect_equal(sort(unique(nsece_synth_small$state_id)), 1:51)
})


test_that("nsece_synth_small: z is binary (0 or 1 only)", {
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  expect_true(all(nsece_synth_small$z %in% c(0L, 1L)))
})


test_that("nsece_synth_small: provider_id values are unique", {
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  expect_equal(length(unique(nsece_synth_small$provider_id)),
               nrow(nsece_synth_small))
})


test_that("nsece_synth_small: covariates have reasonable ranges", {
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  for (col in c("poverty", "urban", "black", "hispanic")) {
    expect_true(min(nsece_synth_small[[col]]) >= 0,
                info = paste(col, "has negative values"))
    expect_true(max(nsece_synth_small[[col]]) <= 100,
                info = paste(col, "exceeds 100"))
  }
})


# ============================================================================
# C. nsece_state_policy — State-level policy data
# ============================================================================

test_that("nsece_state_policy loads without error", {
  expect_no_error(data("nsece_state_policy", package = "hurdlebb",
                       envir = environment()))
  expect_true(exists("nsece_state_policy", envir = environment()))
})


test_that("nsece_state_policy has correct dimensions", {
  data("nsece_state_policy", package = "hurdlebb", envir = environment())
  expect_equal(nrow(nsece_state_policy), 51L)
  expect_equal(ncol(nsece_state_policy), 5L)
})


test_that("nsece_state_policy has correct column names", {
  data("nsece_state_policy", package = "hurdlebb", envir = environment())
  expected_names <- c("state_id", "state_name", "mr_pctile",
                      "tiered_reim", "it_addon")
  expect_equal(sort(names(nsece_state_policy)), sort(expected_names))
})


test_that("nsece_state_policy: state_id is 1:51", {
  data("nsece_state_policy", package = "hurdlebb", envir = environment())
  expect_equal(sort(nsece_state_policy$state_id), 1:51)
})


test_that("nsece_state_policy: state_name is character and unique", {
  data("nsece_state_policy", package = "hurdlebb", envir = environment())
  expect_true(is.character(nsece_state_policy$state_name))
  expect_equal(length(unique(nsece_state_policy$state_name)), 51L)
})


test_that("nsece_state_policy: mr_pctile is numeric", {
  data("nsece_state_policy", package = "hurdlebb", envir = environment())
  expect_true(is.numeric(nsece_state_policy$mr_pctile))
  expect_false(anyNA(nsece_state_policy$mr_pctile))
})


test_that("nsece_state_policy: mr_pctile is approximately standardized", {
  data("nsece_state_policy", package = "hurdlebb", envir = environment())
  # Mean near 0 (within 0.5 of zero)
  expect_true(abs(mean(nsece_state_policy$mr_pctile)) < 0.5,
              info = paste("Mean:", round(mean(nsece_state_policy$mr_pctile), 3)))
  # SD near 1 (within 0.5 of one)
  expect_true(abs(sd(nsece_state_policy$mr_pctile) - 1) < 0.5,
              info = paste("SD:", round(sd(nsece_state_policy$mr_pctile), 3)))
})


test_that("nsece_state_policy: tiered_reim is binary (0 or 1)", {
  data("nsece_state_policy", package = "hurdlebb", envir = environment())
  expect_true(all(nsece_state_policy$tiered_reim %in% c(0L, 1L)))
})


test_that("nsece_state_policy: it_addon is binary (0 or 1)", {
  data("nsece_state_policy", package = "hurdlebb", envir = environment())
  expect_true(all(nsece_state_policy$it_addon %in% c(0L, 1L)))
})


test_that("nsece_state_policy: tiered_reim prevalence is plausible", {
  data("nsece_state_policy", package = "hurdlebb", envir = environment())
  # Approximately 84% have tiered reimbursement (allow 60-100%)
  pct <- mean(nsece_state_policy$tiered_reim)
  expect_true(pct >= 0.60,
              info = paste("Tiered reim pct too low:", round(pct, 3)))
  expect_true(pct <= 1.00,
              info = paste("Tiered reim pct too high:", round(pct, 3)))
})


test_that("nsece_state_policy: it_addon prevalence is plausible", {
  data("nsece_state_policy", package = "hurdlebb", envir = environment())
  # Approximately 24% have IT add-on (allow 10-50%)
  pct <- mean(nsece_state_policy$it_addon)
  expect_true(pct >= 0.10,
              info = paste("IT add-on pct too low:", round(pct, 3)))
  expect_true(pct <= 0.50,
              info = paste("IT add-on pct too high:", round(pct, 3)))
})


test_that("nsece_state_policy: no missing values", {
  data("nsece_state_policy", package = "hurdlebb", envir = environment())
  for (col in names(nsece_state_policy)) {
    expect_false(anyNA(nsece_state_policy[[col]]),
                 info = paste("Column", col, "has NAs"))
  }
})


# ============================================================================
# D. Cross-dataset consistency
# ============================================================================

test_that("all state_ids in nsece_synth appear in nsece_state_policy", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  data("nsece_state_policy", package = "hurdlebb", envir = environment())
  synth_states <- unique(nsece_synth$state_id)
  policy_states <- nsece_state_policy$state_id
  missing <- setdiff(synth_states, policy_states)
  expect_equal(length(missing), 0L,
               info = paste("Missing state_ids in policy:",
                            paste(missing, collapse = ", ")))
})


test_that("all state_ids in nsece_synth_small are in nsece_synth", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  small_states <- unique(nsece_synth_small$state_id)
  full_states <- unique(nsece_synth$state_id)
  expect_true(all(small_states %in% full_states))
})


test_that("all state_ids in nsece_synth_small appear in nsece_state_policy", {
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  data("nsece_state_policy", package = "hurdlebb", envir = environment())
  small_states <- unique(nsece_synth_small$state_id)
  policy_states <- nsece_state_policy$state_id
  expect_true(all(small_states %in% policy_states))
})


test_that("nsece_synth and nsece_synth_small have identical column types", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  for (col in names(nsece_synth)) {
    expect_equal(class(nsece_synth[[col]]), class(nsece_synth_small[[col]]),
                 info = paste("Column", col, "type mismatch between",
                              "nsece_synth and nsece_synth_small"))
  }
})


test_that("nsece_synth_small is strictly smaller than nsece_synth", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  data("nsece_synth_small", package = "hurdlebb", envir = environment())
  expect_true(nrow(nsece_synth_small) < nrow(nsece_synth))
})


test_that("nsece_state_policy merges cleanly with nsece_synth", {
  data("nsece_synth", package = "hurdlebb", envir = environment())
  data("nsece_state_policy", package = "hurdlebb", envir = environment())
  merged <- merge(nsece_synth, nsece_state_policy, by = "state_id",
                  all.x = TRUE)
  # All rows should survive the merge (no NAs introduced)
  expect_equal(nrow(merged), nrow(nsece_synth))
  expect_false(anyNA(merged$mr_pctile),
               info = "NAs in mr_pctile after merge")
  expect_false(anyNA(merged$tiered_reim),
               info = "NAs in tiered_reim after merge")
  expect_false(anyNA(merged$it_addon),
               info = "NAs in it_addon after merge")
})
