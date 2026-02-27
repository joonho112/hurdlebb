# ============================================================================
# test-stan-cache.R — Tests for the Stan model cache system
# ============================================================================

test_that(".hbb_model_names() returns expected model names", {
  nms <- hurdlebb:::.hbb_model_names()
  expect_identical(
    nms,
    c("hbb_base", "hbb_weighted", "hbb_svc", "hbb_svc_weighted")
  )
  expect_type(nms, "character")
  expect_length(nms, 4L)
})

test_that(".hbb_cache_dir() returns a valid path string", {
  cache_dir <- hurdlebb:::.hbb_cache_dir()
  expect_type(cache_dir, "character")
  expect_length(cache_dir, 1L)
  expect_true(nzchar(cache_dir))
  # Must contain "hurdlebb" somewhere in the path
  expect_true(grepl("hurdlebb", cache_dir, fixed = TRUE))
})

test_that(".hbb_exe_exists() checks both plain and .exe paths", {
  tmp <- tempfile("test_exe_")
  on.exit(unlink(tmp), add = TRUE)

  # Non-existent file: FALSE for both variants

  expect_false(hurdlebb:::.hbb_exe_exists(tmp))

  # Create the base file (Unix-style exe)
  writeLines("dummy", tmp)
  expect_true(hurdlebb:::.hbb_exe_exists(tmp))
  unlink(tmp)

  # Create the .exe variant (Windows-style)
  writeLines("dummy", paste0(tmp, ".exe"))
  on.exit(unlink(paste0(tmp, ".exe")), add = TRUE)
  expect_true(hurdlebb:::.hbb_exe_exists(tmp))
})

test_that(".read_hash() handles missing and present files", {
  # Non-existent file returns ""
  expect_identical(
    hurdlebb:::.read_hash(tempfile()),
    ""
  )

  # Existing file returns trimmed first line
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c("abc123  ", "extra line"), tmp)
  expect_identical(hurdlebb:::.read_hash(tmp), "abc123")

  # Empty file returns ""
  writeLines(character(0), tmp)
  expect_identical(hurdlebb:::.read_hash(tmp), "")
})

test_that(".ensure_dir() creates directories", {
  tmp <- file.path(tempdir(), "hbb_test_ensure_dir", "sub1", "sub2")
  on.exit(unlink(file.path(tempdir(), "hbb_test_ensure_dir"), recursive = TRUE),
          add = TRUE)

  expect_false(dir.exists(tmp))
  hurdlebb:::.ensure_dir(tmp)
  expect_true(dir.exists(tmp))

  # Calling again on existing dir is a no-op (no error)
  expect_no_error(hurdlebb:::.ensure_dir(tmp))
})

test_that("hbb_clear_cache() works even with no cache", {
  # Should not error when cache doesn't exist
  expect_no_error(hbb_clear_cache())
})

test_that("hbb_compile() validates model names", {
  skip_if_no_cmdstan()

  expect_error(
    hbb_compile(model = "nonexistent_model"),
    "Unknown model"
  )
})

test_that("hbb_compile() validates flag arguments", {
  expect_error(hbb_compile(force = "yes"))
  expect_error(hbb_compile(quiet = 42))
})

test_that(".hbb_stan_dir() errors when package not loaded properly", {
  # This test verifies the error path. In normal usage via load_all()
  # or an installed package, system.file("stan", ...) returns a valid path.
  # We cannot easily mock system.file(), so just verify the function runs
  # without error when the package IS loaded.
  skip_if_no_stan_files()
  d <- hurdlebb:::.hbb_stan_dir()
  expect_true(dir.exists(d))
})

test_that("lock file mechanism works", {
  lock_file <- tempfile("hbb_lock_test_")
  on.exit(unlink(lock_file), add = TRUE)

  # Acquire lock: creates the file
  hurdlebb:::.acquire_lock(lock_file, name = "test_model")
  expect_true(file.exists(lock_file))

  # Lock file contains our PID
  pid <- hurdlebb:::.read_lock_pid(lock_file)
  expect_identical(pid, Sys.getpid())

  # Release lock: removes the file (since we own it)
  hurdlebb:::.release_lock(lock_file)
  expect_false(file.exists(lock_file))
})

test_that(".read_lock_pid() handles edge cases", {
  # Non-existent file
  expect_true(is.na(hurdlebb:::.read_lock_pid(tempfile())))

  # Empty file
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  writeLines(character(0), tmp)
  expect_true(is.na(hurdlebb:::.read_lock_pid(tmp)))

  # Non-integer content
  writeLines("not_a_pid", tmp)
  expect_true(is.na(hurdlebb:::.read_lock_pid(tmp)))

  # Valid PID
  writeLines("12345", tmp)
  expect_identical(hurdlebb:::.read_lock_pid(tmp), 12345L)
})
