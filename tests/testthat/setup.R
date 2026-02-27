# ============================================================================
# setup.R — Test environment configuration for hurdlebb
#
# This file runs BEFORE any test file. Key responsibilities:
# 1. Redirect Stan model cache to a temp directory (CRAN policy compliance:
#    R CMD check must not write to user directories).
# 2. Set up any test-wide fixtures.
#
# See: https://testthat.r-lib.org/articles/test-fixtures.html
# ============================================================================

# -- Redirect cache to temp directory during tests ----------------------------
# This ensures R CMD check never modifies user-level cache directories.
# The teardown_env() scoping ensures cleanup after all tests finish.
if (requireNamespace("withr", quietly = TRUE)) {
  withr::local_envvar(
    R_USER_CACHE_DIR = tempfile("hurdlebb_test_cache_"),
    .local_envir = teardown_env()
  )
}
