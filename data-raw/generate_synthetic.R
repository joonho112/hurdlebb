# ============================================================================
# generate_synthetic.R --- Generate Synthetic Datasets for hurdlebb Package
#
# This script creates three synthetic datasets from real NSECE 2019 data
# using the Gaussian copula infrastructure in R/copula-utils.R.
#
# Datasets produced:
#   1. nsece_synth         (6785 rows) --- Full synthetic provider dataset
#   2. nsece_synth_small   ( ~500 rows) --- Stratified subsample
#   3. nsece_state_policy  (  51 rows) --- State-level policy indicators
#
# The NSECE restricted-use microdata cannot be shipped. Instead, we:
#   (a) Preserve marginal distributions and rank correlations via copula,
#   (b) Generate hurdle beta-binomial outcomes with known true parameters,
#   (c) Construct a synthetic survey design (strata + PSUs).
#
# Usage (from package root):
#   devtools::load_all()
#   source("data-raw/generate_synthetic.R")
#
# Author: Package build script
# Date:   2026-02-25
# ============================================================================


# ============================================================================
# 0. Setup
# ============================================================================

cat("=== generate_synthetic.R: Starting ===\n\n")

# Load the package in development mode so internal functions are accessible
devtools::load_all()

set.seed(20260225)

# Project root (hurdlebb/) relative to the main it-enrollment repo
project_root <- normalizePath(file.path("..", "dev", "stan", "output"),
                              mustWork = TRUE)

cat("  Project data root:", project_root, "\n\n")


# ============================================================================
# 1. Load Real NSECE Data
# ============================================================================

cat("--- Step 1: Loading real NSECE data ---\n")

analysis_data <- readRDS(file.path(project_root, "analysis_data.rds"))
cat("  analysis_data:", nrow(analysis_data), "rows x",
    ncol(analysis_data), "cols\n")

empirical_cal <- readRDS(file.path(project_root, "simulation",
                                   "empirical_calibration.rds"))
cat("  empirical_calibration: zero_rate =", empirical_cal$zero_rate,
    ", mean_it_share =", empirical_cal$mean_it_share, "\n")

sim_config <- readRDS(file.path(project_root, "simulation", "sim_config.rds"))
true_params <- sim_config$true_params
cat("  sim_config: alpha =", true_params$alpha, "\n")
cat("              beta  =", true_params$beta, "\n")
cat("              kappa =", true_params$kappa, "\n")
cat("              tau   =", true_params$tau, "\n")
cat("              rho   =", true_params$rho, "\n")

stan_data <- readRDS(file.path(project_root, "stan_data.rds"))
V_matrix <- stan_data$V
cat("  V matrix:", nrow(V_matrix), "x", ncol(V_matrix),
    "  (state-level policy)\n\n")


# ============================================================================
# 2. Extract Copula Parameters for Covariates
# ============================================================================

cat("--- Step 2: Extracting copula parameters ---\n")

# Rename columns to short names for the copula extraction
copula_df <- data.frame(
    poverty  = analysis_data$comm_pct_poverty_num,
    urban    = analysis_data$pct_urban,
    black    = analysis_data$comm_pct_black_num,
    hispanic = analysis_data$comm_pct_hisp_num,
    n_trial  = analysis_data$n_trial,
    weight   = analysis_data$weight
)

copula_var_names <- c("poverty", "urban", "black", "hispanic",
                      "n_trial", "weight")
copula_types <- c(
    poverty  = "continuous",
    urban    = "continuous",
    black    = "continuous",
    hispanic = "continuous",
    n_trial  = "count",       # discrete integer
    weight   = "continuous"
)

copula_params <- .extract_copula_params(
    data      = copula_df,
    var_names = copula_var_names,
    types     = copula_types
)

cat("  n_complete:", copula_params$n_complete,
    " n_dropped:", copula_params$n_dropped, "\n")
cat("  Copula correlation matrix (z-space):\n")
print(round(copula_params$cor_matrix, 3))
cat("\n")


# ============================================================================
# 3. Generate Synthetic Covariates (N = 6785)
# ============================================================================

cat("--- Step 3: Generating synthetic covariates via Gaussian copula ---\n")

N <- nrow(analysis_data)  # 6785
S <- 51L                  # number of states

synth_covariates <- .copula_generate(
    n          = N,
    cor_matrix = copula_params$cor_matrix,
    inv_ecdfs  = copula_params$inv_ecdfs,
    seed       = 20260225
)

cat("  Generated:", nrow(synth_covariates), "rows x",
    ncol(synth_covariates), "cols\n")

# Force n_trial to positive integer (at least 1)
synth_covariates$n_trial <- pmax(1L, as.integer(round(synth_covariates$n_trial)))

# Force weight to be strictly positive
synth_covariates$weight <- pmax(synth_covariates$weight, 0.1)

cat("  n_trial range:", range(synth_covariates$n_trial),
    "  median:", median(synth_covariates$n_trial), "\n")
cat("  weight  range:", round(range(synth_covariates$weight), 2),
    "  mean:", round(mean(synth_covariates$weight), 2), "\n\n")


# ============================================================================
# 4. State Assignment
# ============================================================================

cat("--- Step 4: Assigning states ---\n")

state_sizes <- empirical_cal$state_sizes   # integer vector, length 51
stopifnot(length(state_sizes) == S)
stopifnot(sum(state_sizes) == N)

# Create state_id vector: rep(1, n1), rep(2, n2), ..., rep(51, n51)
# then shuffle randomly
state_id_ordered <- rep(seq_len(S), times = state_sizes)
state_id <- sample(state_id_ordered)  # random permutation

cat("  State sizes (first 10):", head(state_sizes, 10), "\n")
cat("  Unique states after shuffle:", length(unique(state_id)), "\n\n")


# ============================================================================
# 5. Survey Design Generation (Strata + PSUs)
# ============================================================================

cat("--- Step 5: Generating synthetic survey design ---\n")

# Target: ~30 strata, ~415 PSUs total
# Strategy: group states into strata based on sorted population sizes,
# then create PSUs within each stratum.

n_strata <- 30L

# Sort states by descending size; assign to strata in round-robin
state_order <- order(state_sizes, decreasing = TRUE)
stratum_of_state <- integer(S)
for (i in seq_along(state_order)) {
    stratum_of_state[state_order[i]] <- ((i - 1L) %% n_strata) + 1L
}

# Within each stratum, create PSUs. Target ~14 PSUs per stratum.
# We assign PSUs proportional to the number of providers per stratum.
n_psu_target <- 415L

# Count providers per stratum
stratum_of_provider <- stratum_of_state[state_id]
providers_per_stratum <- tabulate(stratum_of_provider, nbins = n_strata)

# Allocate PSUs proportional to stratum size, minimum 2 per stratum
psu_alloc_raw <- providers_per_stratum / sum(providers_per_stratum) * n_psu_target
psu_alloc <- pmax(2L, as.integer(round(psu_alloc_raw)))

# Create PSU IDs within each stratum
psu_id <- integer(N)
psu_counter <- 0L

for (s in seq_len(n_strata)) {
    idx <- which(stratum_of_provider == s)
    n_psu_s <- psu_alloc[s]
    # Assign providers to PSUs sequentially (after random ordering within stratum)
    idx_shuffled <- sample(idx)
    psu_labels <- rep_len(seq_len(n_psu_s), length(idx_shuffled))
    psu_id[idx_shuffled] <- psu_labels + psu_counter
    psu_counter <- psu_counter + n_psu_s
}

cat("  Number of strata:", n_strata, "\n")
cat("  Total PSUs:", length(unique(psu_id)), "\n")
cat("  PSU ID range:", range(psu_id), "\n\n")


# ============================================================================
# 6. Standardize Covariates and Build Design Matrix
# ============================================================================

cat("--- Step 6: Standardizing covariates ---\n")

poverty_raw  <- synth_covariates$poverty
urban_raw    <- synth_covariates$urban
black_raw    <- synth_covariates$black
hispanic_raw <- synth_covariates$hispanic

poverty_std  <- scale(poverty_raw)[, 1]
urban_std    <- scale(urban_raw)[, 1]
black_std    <- scale(black_raw)[, 1]
hispanic_std <- scale(hispanic_raw)[, 1]

# Design matrix: intercept + 4 standardized covariates
X <- cbind(1, poverty_std, urban_std, black_std, hispanic_std)
colnames(X) <- c("intercept", "poverty", "urban", "black", "hispanic")

cat("  X dimensions:", nrow(X), "x", ncol(X), "\n")
cat("  Covariate means (should be ~0):",
    round(colMeans(X[, -1]), 4), "\n")
cat("  Covariate SDs   (should be ~1):",
    round(apply(X[, -1], 2, sd), 4), "\n\n")


# ============================================================================
# 7. State-Level Policy Data (V matrix)
# ============================================================================

cat("--- Step 7: Preparing state-level policy data ---\n")

# The V matrix from stan_data contains public policy indicators:
# intercept, MR_pctile (standardized), TieredReim (0/1), ITaddon (0/1).
# These are aggregate state-level policy data, not restricted provider data.

nsece_state_policy <- data.frame(
    state_id    = seq_len(S),
    state_name  = sprintf("State_%02d", seq_len(S)),
    mr_pctile   = V_matrix[, "MR_pctile"],
    tiered_reim = as.integer(V_matrix[, "TieredReim"]),
    it_addon    = as.integer(V_matrix[, "ITaddon"])
)

cat("  States with TieredReim = 1:",
    sum(nsece_state_policy$tiered_reim), "of", S, "\n")
cat("  States with ITaddon = 1:",
    sum(nsece_state_policy$it_addon), "of", S, "\n")
cat("  MR_pctile range:", round(range(nsece_state_policy$mr_pctile), 3), "\n\n")


# ============================================================================
# 8. Generate Outcomes via Hurdle Beta-Binomial
# ============================================================================

cat("--- Step 8: Generating outcomes via Hurdle BB ---\n")

# Extract true parameters
alpha <- true_params$alpha  # length 5 (intercept + 4 covariates)
beta  <- true_params$beta   # length 5
kappa <- true_params$kappa  # scalar
tau   <- true_params$tau    # length 2 (ext, int)
rho   <- true_params$rho    # scalar

# --- 8a. State random effects ---
Omega <- matrix(c(1, rho, rho, 1), 2, 2)
Sigma_delta <- diag(tau) %*% Omega %*% diag(tau)

# Draw 51 state-level random effects (bivariate)
delta <- MASS::mvrnorm(n = S, mu = c(0, 0), Sigma = Sigma_delta)
colnames(delta) <- c("delta_ext", "delta_int")

cat("  State RE delta_ext: mean =", round(mean(delta[, 1]), 3),
    ", sd =", round(sd(delta[, 1]), 3), "\n")
cat("  State RE delta_int: mean =", round(mean(delta[, 2]), 3),
    ", sd =", round(sd(delta[, 2]), 3), "\n")

# --- 8b. Linear predictors ---
eta_ext <- as.numeric(X %*% alpha) + delta[state_id, 1]
eta_int <- as.numeric(X %*% beta)  + delta[state_id, 2]

# --- 8c. Transform to probability scale ---
q  <- plogis(eta_ext)   # participation probability
mu <- plogis(eta_int)   # intensity mean (BB)

cat("  q  (participation): mean =", round(mean(q), 3),
    ", range =", round(range(q), 3), "\n")
cat("  mu (intensity):     mean =", round(mean(mu), 3),
    ", range =", round(range(mu), 3), "\n")

# --- 8d. Sample outcomes ---
n_trial_vec <- synth_covariates$n_trial

# Part 1: participation indicator
z <- rbinom(N, size = 1L, prob = q)

# Part 2: for z == 1, draw from ZT-BetaBin
y <- integer(N)
active <- which(z == 1L)
cat("  Active (z=1):", length(active), "of", N, "\n")

if (length(active) > 0L) {
    y[active] <- rztbetabinom(
        nn    = length(active),
        n     = n_trial_vec[active],
        mu    = mu[active],
        kappa = kappa
    )
}

# Compute IT share
it_share <- ifelse(z == 1L, y / n_trial_vec, 0)

cat("  Zero rate:          ", round(mean(z == 0), 3),
    "  (target ~0.353)\n")
cat("  Mean IT share (z=1):", round(mean(it_share[z == 1L]), 3),
    "  (target ~0.478)\n")
cat("  y range:", range(y), "\n\n")


# ============================================================================
# 9. Assemble nsece_synth (Full Dataset)
# ============================================================================

cat("--- Step 9: Assembling nsece_synth ---\n")

nsece_synth <- data.frame(
    provider_id = seq_len(N),
    state_id    = state_id,
    y           = y,
    n_trial     = n_trial_vec,
    z           = z,
    it_share    = it_share,
    poverty     = poverty_raw,
    urban       = urban_raw,
    black       = black_raw,
    hispanic    = hispanic_raw,
    weight      = synth_covariates$weight,
    stratum     = stratum_of_provider,
    psu         = psu_id
)

cat("  Dimensions:", nrow(nsece_synth), "x", ncol(nsece_synth), "\n\n")


# ============================================================================
# 10. Create nsece_synth_small (Stratified Subsample, ~500 rows)
# ============================================================================

cat("--- Step 10: Creating nsece_synth_small (~500 rows) ---\n")

# Stratified sampling by state: draw proportional to state size,
# with minimum 2 per state, targeting ~500 total.
target_small <- 500L

# Proportional allocation
small_alloc_raw <- (state_sizes / sum(state_sizes)) * target_small
small_alloc <- pmax(2L, as.integer(round(small_alloc_raw)))

# Adjust to be close to target
while (sum(small_alloc) > target_small + 20L) {
    # Reduce from the largest allocations
    max_idx <- which.max(small_alloc)
    small_alloc[max_idx] <- small_alloc[max_idx] - 1L
}

# Sample from each state
small_idx <- integer(0)
for (s in seq_len(S)) {
    state_rows <- which(nsece_synth$state_id == s)
    n_draw <- min(small_alloc[s], length(state_rows))
    if (n_draw > 0L) {
        drawn <- sample(state_rows, size = n_draw, replace = FALSE)
        small_idx <- c(small_idx, drawn)
    }
}

nsece_synth_small <- nsece_synth[small_idx, ]
# Re-index provider_id for the small dataset
nsece_synth_small$provider_id <- seq_len(nrow(nsece_synth_small))
rownames(nsece_synth_small) <- NULL

cat("  Dimensions:", nrow(nsece_synth_small), "x",
    ncol(nsece_synth_small), "\n")
cat("  States represented:", length(unique(nsece_synth_small$state_id)), "\n")
cat("  Zero rate:", round(mean(nsece_synth_small$z == 0), 3), "\n\n")


# ============================================================================
# 11. Validation Checks
# ============================================================================

cat("--- Step 11: Validation checks ---\n\n")

# --- 11a. Outcome consistency ---
cat("  [Check 1] No negative y .............. ")
stopifnot(all(nsece_synth$y >= 0))
cat("PASS\n")

cat("  [Check 2] y <= n_trial ............... ")
stopifnot(all(nsece_synth$y <= nsece_synth$n_trial))
cat("PASS\n")

cat("  [Check 3] z == I(y > 0) .............. ")
stopifnot(all(nsece_synth$z == as.integer(nsece_synth$y > 0)))
cat("PASS\n")

cat("  [Check 4] it_share in [0, 1] ......... ")
stopifnot(all(nsece_synth$it_share >= 0 & nsece_synth$it_share <= 1))
cat("PASS\n")

# --- 11b. Marginal properties ---
obs_zero_rate <- mean(nsece_synth$z == 0)
obs_it_share  <- mean(nsece_synth$it_share[nsece_synth$z == 1])
obs_ntrial_med <- median(nsece_synth$n_trial)

cat("\n  [Check 5] Zero rate:          ", round(obs_zero_rate, 3),
    "  (target ~0.353)\n")
cat("  [Check 6] Mean IT share (z=1):", round(obs_it_share, 3),
    "  (target ~0.478)\n")
cat("  [Check 7] n_trial median:     ", obs_ntrial_med,
    "  (target ~48)\n")

# --- 11c. Spearman correlation preservation ---
cat("\n  [Check 8] Spearman correlation preservation:\n")

orig_covs <- data.frame(
    poverty  = analysis_data$comm_pct_poverty_num,
    urban    = analysis_data$pct_urban,
    black    = analysis_data$comm_pct_black_num,
    hispanic = analysis_data$comm_pct_hisp_num
)
synth_covs <- data.frame(
    poverty  = nsece_synth$poverty,
    urban    = nsece_synth$urban,
    black    = nsece_synth$black,
    hispanic = nsece_synth$hispanic
)

cor_orig  <- cor(orig_covs,  method = "spearman")
cor_synth <- cor(synth_covs, method = "spearman")
cor_diff  <- abs(cor_orig - cor_synth)

cat("    Max |Spearman difference|:", round(max(cor_diff), 4))
if (max(cor_diff) < 0.10) {
    cat("  PASS (< 0.10)\n")
} else {
    cat("  WARNING (>= 0.10)\n")
}

# Print the correlation comparison
cat("\n    Original Spearman correlations:\n")
print(round(cor_orig, 3))
cat("\n    Synthetic Spearman correlations:\n")
print(round(cor_synth, 3))

# --- 11d. Dataset integrity ---
cat("\n  [Check 9]  nsece_synth: all cols present .. ")
expected_cols <- c("provider_id", "state_id", "y", "n_trial", "z",
                   "it_share", "poverty", "urban", "black", "hispanic",
                   "weight", "stratum", "psu")
stopifnot(all(expected_cols %in% names(nsece_synth)))
cat("PASS\n")

cat("  [Check 10] nsece_state_policy: 51 rows .... ")
stopifnot(nrow(nsece_state_policy) == 51)
cat("PASS\n")

cat("  [Check 11] nsece_synth_small: 51 states ... ")
stopifnot(length(unique(nsece_synth_small$state_id)) == 51)
cat("PASS\n")

cat("  [Check 12] No NA in nsece_synth ........... ")
stopifnot(!any(is.na(nsece_synth)))
cat("PASS\n")

cat("  [Check 13] No NA in nsece_state_policy .... ")
stopifnot(!any(is.na(nsece_state_policy)))
cat("PASS\n")

cat("\n")


# ============================================================================
# 12. Save Datasets
# ============================================================================

cat("--- Step 12: Saving datasets via usethis::use_data() ---\n")

usethis::use_data(nsece_synth, overwrite = TRUE)
usethis::use_data(nsece_synth_small, overwrite = TRUE)
usethis::use_data(nsece_state_policy, overwrite = TRUE)


# ============================================================================
# 13. Summary
# ============================================================================

cat("\n=== generate_synthetic.R: Complete ===\n\n")

cat("Datasets saved to data/:\n")
cat("  1. nsece_synth        :", nrow(nsece_synth), "rows x",
    ncol(nsece_synth), "cols\n")
cat("  2. nsece_synth_small  :", nrow(nsece_synth_small), "rows x",
    ncol(nsece_synth_small), "cols\n")
cat("  3. nsece_state_policy :", nrow(nsece_state_policy), "rows x",
    ncol(nsece_state_policy), "cols\n")

cat("\nTrue parameters (for testing / verification):\n")
cat("  alpha:", true_params$alpha, "\n")
cat("  beta :", true_params$beta, "\n")
cat("  kappa:", true_params$kappa, "\n")
cat("  tau  :", true_params$tau, "\n")
cat("  rho  :", true_params$rho, "\n")

cat("\nKey statistics:\n")
cat("  Zero rate:           ", round(obs_zero_rate, 3), "\n")
cat("  Mean IT share (z=1): ", round(obs_it_share, 3), "\n")
cat("  n_trial median:      ", obs_ntrial_med, "\n")
cat("  N states:            ", S, "\n")
cat("  N strata:            ", n_strata, "\n")
cat("  N PSUs:              ", length(unique(psu_id)), "\n")
