#' @keywords internal
"_PACKAGE"

# ============================================================================
# Namespace imports
#
# Philosophy: @importFrom for functions used frequently throughout the
# package (avoids repeated pkg:: lookups at runtime). Functions used
# in only one file or rarely should use explicit pkg::fun() notation.
#
# Categories:
#   ALWAYS imported: rlang, cli, checkmate, digest, stats (high-frequency)
#   :: notation: posterior, loo, Matrix, MASS (used in specific modules)
#   :: notation: cmdstanr (Suggests, not Imports)
# ============================================================================

## usethis namespace: start

# -- rlang: Core error handling & tidy eval ----------------------------------
# Used pervasively across nearly every file.
#' @importFrom rlang check_installed abort warn inform caller_env .data

# -- cli: User-facing messages -----------------------------------------------
# cli_alert_* used in stan-cache.R, fit.R, sandwich.R.
# cli_abort/warn for structured errors (used throughout).
#' @importFrom cli cli_alert_info cli_alert_success cli_alert_warning
#' @importFrom cli cli_abort cli_warn cli_inform qty
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done

# -- checkmate: Fast input validation ----------------------------------------
# Used in public-facing functions for argument checking.
#' @importFrom checkmate assert_flag assert_string assert_choice assert_numeric
#' @importFrom checkmate assert_number assert_integerish assert_matrix
#' @importFrom checkmate assert_data_frame assert_count

# -- digest: MD5 hashing for Stan cache -------------------------------------
#' @importFrom digest digest

# -- stats: Base R statistical functions -------------------------------------
# Used across distributions.R, data-prep.R, marginal-effects.R, etc.
#' @importFrom stats model.matrix plogis qlogis qnorm pnorm
#' @importFrom stats rbinom rnorm runif dbinom rbeta
#' @importFrom stats setNames sd dnorm var cov median quantile
#' @importFrom stats residuals fitted predict coef vcov nobs

# -- utils: Utility functions -------------------------------------------------
# head() used in sandwich.R for singleton strata display.
#' @importFrom utils head

## usethis namespace: end
NULL
