# ============================================================================
# helper-skip.R — Test skip helpers for hurdlebb
#
# Loaded automatically by testthat before any test file runs.
# See ?testthat::test_dir for the helper-*.R loading convention.
# ============================================================================


#' Skip test if CmdStan is not installed
#'
#' Tests that require Stan compilation or model fitting should call this
#' at the top of each `test_that()` block. Checks three things:
#'
#' 1. The `cmdstanr` R package is installed.
#' 2. CmdStan itself is installed and findable via `cmdstanr::cmdstan_path()`.
#' 3. The CmdStan binary actually exists at the reported path.
#'
#' @noRd
skip_if_no_cmdstan <- function() {
  testthat::skip_if_not_installed("cmdstanr")

  path <- tryCatch(
    cmdstanr::cmdstan_path(),
    error = function(e) NULL
  )
  testthat::skip_if(
    is.null(path),
    message = "CmdStan is not installed. Install via cmdstanr::install_cmdstan()."
  )

  # Belt-and-suspenders: verify the path actually exists on disk.
  # Protects against stale cmdstanr config pointing to a deleted installation.
  testthat::skip_if(
    !dir.exists(path),
    message = paste0(
      "CmdStan path reported as '", path,
      "' but directory does not exist. Reinstall via cmdstanr::install_cmdstan()."
    )
  )
}


#' Skip test if hurdlebb Stan models have not been compiled
#'
#' For integration tests that need a pre-compiled Stan model but should
#' not spend time compiling during a quick test run. Checks that at least
#' one compiled executable exists in the cache directory.
#'
#' @noRd
skip_if_no_hbb_models <- function() {
  skip_if_no_cmdstan()

  cache_dir <- tools::R_user_dir("hurdlebb", "cache")
  testthat::skip_if(
    !dir.exists(cache_dir),
    message = "hurdlebb model cache not found. Run hbb_compile() first."
  )

  # Check for at least one executable
  model_names <- c("hbb_base", "hbb_weighted", "hbb_svc", "hbb_svc_weighted")
  has_exe <- vapply(model_names, function(nm) {
    base <- file.path(cache_dir, nm)
    file.exists(base) || file.exists(paste0(base, ".exe"))
  }, logical(1L))

  testthat::skip_if(
    !any(has_exe),
    message = "No compiled hurdlebb Stan models found. Run hbb_compile() first."
  )
}


#' Skip test if running on CRAN (convenience alias)
#'
#' Equivalent to `testthat::skip_on_cran()`. Provided here so that test
#' files can use a consistent `skip_*` pattern from a single file.
#'
#' @noRd
skip_on_cran_or_ci <- function() {
  testthat::skip_on_cran()
  # Also skip on CI if tests are too slow (> 60s compile times)
  testthat::skip_if(
    nzchar(Sys.getenv("CI")) && !nzchar(Sys.getenv("HBB_RUN_SLOW_TESTS")),
    message = "Skipping slow test on CI. Set HBB_RUN_SLOW_TESTS=true to enable."
  )
}


#' Skip test if Stan file(s) are not present
#'
#' For tests that need to read Stan source files but do not need CmdStan
#' (e.g., hash computation tests, source parsing tests).
#'
#' @param model_name Name of the model to check. If `NULL`, checks that
#'   the `inst/stan/` directory exists at all.
#' @noRd
skip_if_no_stan_files <- function(model_name = NULL) {
  if (is.null(model_name)) {
    stan_dir <- system.file("stan", package = "hurdlebb")
    testthat::skip_if(
      !nzchar(stan_dir) || !dir.exists(stan_dir),
      message = "inst/stan/ directory not found. Is the package loaded?"
    )
  } else {
    f <- system.file("stan", paste0(model_name, ".stan"), package = "hurdlebb")
    testthat::skip_if(
      !nzchar(f),
      message = paste0("Stan file '", model_name, ".stan' not found in inst/stan/.")
    )
  }
}
