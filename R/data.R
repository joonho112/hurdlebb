# ============================================================================
# data.R --- Roxygen2 Documentation for Synthetic Datasets
#
# Three datasets:
#   1. nsece_synth         — Main synthetic NSECE 2019 dataset (~6,785 rows)
#   2. nsece_synth_small   — Stratified subsample (~500 rows)
#   3. nsece_state_policy  — State-level policy variables (51 rows)
#
# All datasets are lazy-loaded from data/*.rda files.
# ============================================================================


# -- 1. nsece_synth -----------------------------------------------------------

#' Synthetic NSECE 2019 Center-Based Provider Data
#'
#' A synthetic dataset mimicking the structure and statistical properties
#' of the 2019 National Survey of Early Care and Education (NSECE)
#' center-based provider data. Generated via Gaussian copula to preserve
#' marginal distributions and rank-correlation structure of the original
#' restricted-use data.
#'
#' @format A data frame with approximately 6,785 rows and 13 columns:
#' \describe{
#'   \item{provider_id}{Integer. Unique provider identifier (1 to N).}
#'   \item{state_id}{Integer. State identifier (1 to 51, including DC).}
#'   \item{y}{Integer. Infant/toddler (IT) enrollment count. Zero for
#'     centers not serving IT (\code{z = 0}).}
#'   \item{n_trial}{Integer. Total enrollment for children ages 0--5.
#'     Always \eqn{\ge 1}.}
#'   \item{z}{Integer. Participation indicator: 1 if the center serves
#'     any IT children, 0 otherwise. Equal to \code{I(y > 0)}.}
#'   \item{it_share}{Numeric. IT enrollment share \code{y / n_trial}.
#'     Zero when \code{z = 0}.}
#'   \item{poverty}{Numeric. Community poverty rate (percent below
#'     poverty line). Range approximately 2--54.}
#'   \item{urban}{Numeric. Community urbanization rate (percent urban
#'     population). Range 0--100, heavily right-skewed.}
#'   \item{black}{Numeric. Community percent Black population.
#'     Range 0--97.}
#'   \item{hispanic}{Numeric. Community percent Hispanic population.
#'     Range 0--98.}
#'   \item{weight}{Numeric. Survey sampling weight. Positive values;
#'     larger weights indicate the provider represents more units in the
#'     population.}
#'   \item{stratum}{Integer. Sampling stratum identifier.}
#'   \item{psu}{Integer. Primary sampling unit (PSU) identifier, nested
#'     within strata.}
#' }
#'
#' @details
#' The data-generating process uses a two-part hurdle Beta-Binomial model:
#' \describe{
#'   \item{Part 1 (Extensive margin)}{Whether a center serves IT children
#'     at all, modeled as \eqn{z_i \sim \textrm{Bernoulli}(q_i)} where
#'     \eqn{\textrm{logit}(q_i) = X_i \alpha + \delta_{1,s_i}}.}
#'   \item{Part 2 (Intensive margin)}{Among servers, the IT enrollment
#'     count follows a zero-truncated Beta-Binomial:
#'     \eqn{y_i \mid z_i = 1 \sim \textrm{ZT-BetaBin}(n_i, \mu_i, \kappa)}
#'     where \eqn{\textrm{logit}(\mu_i) = X_i \beta + \delta_{2,s_i}}.}
#' }
#'
#' The structural zero rate is approximately 35 percent, and the mean IT share
#' among servers is approximately 48 percent. These match the empirical moments
#' of the original NSECE 2019 data.
#'
#' @source
#' Generated from the NSECE 2019 restricted-use data via Gaussian copula
#' (Miratrix, 2025) with hurdle Beta-Binomial outcomes calibrated to the
#' empirical moments of the original data. See
#' \code{data-raw/generate_synthetic.R} for the full generation script.
#'
#' @references
#' National Survey of Early Care and Education (NSECE), 2019.
#' U.
#' S. Department of Health and Human Services.
#'
#' @seealso
#' \code{\link{nsece_synth_small}} for a smaller version suitable for
#' quick examples and tests, \code{\link{nsece_state_policy}} for
#' state-level policy variables.
#'
#' @examples
#' data(nsece_synth)
#' head(nsece_synth)
#'
#' # Basic summary
#' cat("N =", nrow(nsece_synth), "\n")
#' cat("Zero rate:", round(1 - mean(nsece_synth$z), 3), "\n")
#' cat("IT share (servers):",
#'     round(mean(nsece_synth$it_share[nsece_synth$z == 1]), 3), "\n")
#'
#' # State distribution
#' cat("States:", length(unique(nsece_synth$state_id)), "\n")
#'
"nsece_synth"


# -- 2. nsece_synth_small -----------------------------------------------------

#' Small Synthetic NSECE 2019 Dataset
#'
#' A stratified subsample of \code{\link{nsece_synth}} with approximately
#' 500 rows. Designed for quick demonstrations, unit tests, and package
#' examples where computational speed matters more than statistical power.
#'
#' @format A data frame with approximately 500 rows and the same 13
#'   columns as \code{\link{nsece_synth}}:
#' \describe{
#'   \item{provider_id}{Integer. Unique provider identifier (re-indexed
#'     1 to N for this subset).}
#'   \item{state_id}{Integer. State identifier (1 to 51, including DC).}
#'   \item{y}{Integer. Infant/toddler (IT) enrollment count.}
#'   \item{n_trial}{Integer. Total enrollment for children ages 0--5.}
#'   \item{z}{Integer. Participation indicator: 1 if \code{y > 0}, 0
#'     otherwise.}
#'   \item{it_share}{Numeric. IT enrollment share \code{y / n_trial}.}
#'   \item{poverty}{Numeric. Community poverty rate.}
#'   \item{urban}{Numeric. Community urbanization rate.}
#'   \item{black}{Numeric. Community percent Black population.}
#'   \item{hispanic}{Numeric. Community percent Hispanic population.}
#'   \item{weight}{Numeric. Survey sampling weight.}
#'   \item{stratum}{Integer. Sampling stratum identifier.}
#'   \item{psu}{Integer. Primary sampling unit (PSU) identifier.}
#' }
#'
#' @details
#' The subsample is drawn via stratified sampling on \code{state_id} to
#' preserve the state size distribution of the full dataset. All 51
#' state identifiers are represented. The structural zero rate and
#' covariate distributions approximately match the full dataset, though
#' with more sampling variability due to the smaller size.
#'
#' This dataset is \strong{not suitable for model fitting}. With only
#' \eqn{\sim 10} observations per state on average, state-varying
#' coefficient models will not converge reliably. Use
#' \code{\link{nsece_synth}} for any serious modeling.
#'
#' @source
#' Stratified subsample of \code{\link{nsece_synth}}. See
#' \code{data-raw/generate_synthetic.R} for the generation script.
#'
#' @seealso
#' \code{\link{nsece_synth}} for the full dataset,
#' \code{\link{nsece_state_policy}} for state-level policy variables.
#'
#' @examples
#' data(nsece_synth_small)
#'
#' # Quick overview
#' cat("N =", nrow(nsece_synth_small), "\n")
#' cat("States:", length(unique(nsece_synth_small$state_id)), "\n")
#' cat("Zero rate:", round(1 - mean(nsece_synth_small$z), 3), "\n")
#'
#' # Useful for quick tests and demonstrations
#' table(nsece_synth_small$z)
#'
"nsece_synth_small"


# -- 3. nsece_state_policy ----------------------------------------------------

#' State-Level Childcare Policy Variables
#'
#' State-level policy indicators for all 51 jurisdictions (50 states + DC).
#' Used as cross-level moderators in the hurdle Beta-Binomial model with
#' state-varying coefficients (Module D).
#'
#' @format A data frame with 51 rows and 5 columns:
#' \describe{
#'   \item{state_id}{Integer. State identifier matching
#'     \code{nsece_synth$state_id} (1 to 51).}
#'   \item{state_name}{Character. Generic state label
#'     (\code{"State_01"} through \code{"State_51"}). Names are
#'     anonymized to prevent identification of actual states.}
#'   \item{mr_pctile}{Numeric. CCDF market rate percentile
#'     (standardized, mean \eqn{\approx 0}, SD \eqn{\approx 1}).
#'     Higher values indicate more generous reimbursement rates.}
#'   \item{tiered_reim}{Integer. 1 if the state uses tiered
#'     reimbursement (quality-differentiated CCDF rates), 0 otherwise.
#'     Approximately 84 percent of states have tiered reimbursement.}
#'   \item{it_addon}{Integer. 1 if the state provides an IT add-on
#'     payment (supplemental rate for infant/toddler care), 0 otherwise.
#'     Approximately 24 percent of states have IT add-ons.}
#' }
#'
#' @details
#' These policy variables serve as cross-level moderators in the
#' state-varying coefficient (SVC) model. They enter the model as
#' predictors of the state random effects:
#' \deqn{\delta_s = \Gamma v_s + \eta_s}
#' where \eqn{v_s = (1, \texttt{mr\_pctile}_s, \texttt{tiered\_reim}_s,
#' \texttt{it\_addon}_s)'} and \eqn{\Gamma} captures the cross-level
#' interaction effects.
#'
#' The \code{mr_pctile} variable is standardized (centered and scaled)
#' to improve MCMC sampling and interpretability. The binary indicators
#' \code{tiered_reim} and \code{it_addon} are left uncentered.
#'
#' @source
#' Synthetic values calibrated to the distribution of actual U.S.
#' childcare subsidy policies. See the CCDF Policies Database
#' maintained by the Urban Institute for original policy data.
#' Generated in \code{data-raw/generate_synthetic.R}.
#'
#' @references
#' National Survey of Early Care and Education (NSECE), 2019.
#' U.S. Department of Health and Human Services.
#'
#' @seealso
#' \code{\link{nsece_synth}} for the provider-level data that can be
#' merged with this table via \code{state_id}.
#'
#' @examples
#' data(nsece_state_policy)
#' head(nsece_state_policy)
#'
#' # Policy summary
#' cat("Tiered reimbursement:",
#'     round(mean(nsece_state_policy$tiered_reim), 2), "\n")
#' cat("IT add-on:",
#'     round(mean(nsece_state_policy$it_addon), 2), "\n")
#'
#' # Merge with provider data
#' data(nsece_synth)
#' merged <- merge(nsece_synth, nsece_state_policy, by = "state_id")
#' cat("Merged rows:", nrow(merged), "\n")
#'
"nsece_state_policy"
