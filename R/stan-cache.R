# ============================================================================
# stan-cache.R — Stan Model Compilation & Persistent Cache Management
#
# Strategy: Runtime compilation with persistent user-level cache.
# Stan files ship as text in inst/stan/. On first use (or when source
# changes), they are compiled via cmdstanr and the executable is cached
# in tools::R_user_dir("hurdlebb", "cache"). Subsequent calls reuse the
# cached executable unless the Stan source hash has changed.
#
# Works with both devtools::load_all() and installed packages because
# system.file() is shimmed by devtools to resolve inst/ paths.
#
# Precedent: brms, rethinking, EpiNow2 all use runtime compilation.
# ============================================================================

# -- Package-level environment for session cache (NOT persisted) -------------
# Maps model names to CmdStanModel objects within the current R session.
# Cleared on package reload. The persistent cache on disk is separate.
.hbb_cache <- new.env(parent = emptyenv())


# ============================================================================
# Public: hbb_compile()
# ============================================================================

#' Compile Stan Models for hurdlebb
#'
#' Compiles (or retrieves from cache) the Stan models shipped with the
#' package. Compiled executables are stored in a persistent user-level
#' cache directory so that compilation happens only once per Stan source
#' version.
#'
#' You can call this function explicitly, or let the main fitting
#' function trigger compilation automatically on first use.
#'
#' @param model Character vector of model names to compile, or `NULL`
#'   (default) to compile all four models. Valid names: `"hbb_base"`,
#'   `"hbb_weighted"`, `"hbb_svc"`, `"hbb_svc_weighted"`.
#' @param force Logical. If `TRUE`, recompile even when the cached
#'   executable is up-to-date. Default `FALSE`.
#' @param quiet Logical. If `TRUE`, suppress informational messages
#'   during compilation. Default `FALSE`.
#' @param cpp_options Named list of C++ compilation options passed to
#'   `cmdstanr::cmdstan_model()`. For example,
#'   `list(stan_threads = TRUE)` enables within-chain threading.
#'   Default `list()`.
#' @param stanc_options Named list of stanc3 transpiler options.
#'   Default `list()`.
#'
#' @return Invisibly returns a named list of compiled `CmdStanModel`
#'   objects.
#'
#' @details
#' ## Cache location
#' Compiled executables are stored under
#' `tools::R_user_dir("hurdlebb", "cache")`, following R >= 4.0
#' conventions. Each model produces two files in the cache directory:
#' the compiled executable and a `.md5` file recording the MD5 digest
#' of the Stan source at compilation time.
#'
#' ## Hash-based invalidation
#' When a Stan source file changes (e.g., after a package update), the
#' MD5 hash no longer matches the cached `.md5` file, triggering
#' automatic recompilation.
#'
#' ## Thread safety
#' A simple file-based lock mechanism prevents two concurrent R sessions
#' from compiling the same model simultaneously. If a lock file is
#' detected, the function polls and waits (with a configurable timeout).
#' Stale locks from dead processes are broken automatically.
#'
#' @seealso \code{\link{hbb_clear_cache}} to remove all cached executables.
#'
#' @export
#' @examples
#' \dontrun{
#' # Compile all models (one-time operation)
#' hbb_compile()
#'
#' # Compile only the base model with threading support
#' hbb_compile("hbb_base", cpp_options = list(stan_threads = TRUE))
#'
#' # Force recompilation after CmdStan upgrade
#' hbb_compile(force = TRUE)
#' }
hbb_compile <- function(model = NULL,
                        force = FALSE,
                        quiet = FALSE,
                        cpp_options = list(),
                        stanc_options = list()) {
  # -- Input validation -------------------------------------------------------
  assert_flag(force)
  assert_flag(quiet)
  checkmate::assert_list(cpp_options, names = "unique")
  checkmate::assert_list(stanc_options, names = "unique")

  all_models <- .hbb_model_names()

  if (is.null(model)) {
    model <- all_models
  } else {
    checkmate::assert_character(model, min.len = 1L, any.missing = FALSE)
    bad <- setdiff(model, all_models)
    if (length(bad) > 0L) {
      valid <- .hbb_model_names()
      cli::cli_abort(c(
        "Unknown model name{?s}: {.val {bad}}.",
        "i" = "Valid names: {.val {valid}}."
      ))
    }
  }

  # -- Check cmdstanr is available -------------------------------------------
  check_installed("cmdstanr", reason = "to compile and run Stan models.")

  # -- Check CmdStan installation (not just the R interface) -----------------
  .assert_cmdstan_available()

  # -- Create cache directory ------------------------------------------------
  cache_dir <- .hbb_cache_dir()
  .ensure_dir(cache_dir)

  # -- Compile each model ----------------------------------------------------
  compiled <- stats::setNames(vector("list", length(model)), model)

  for (nm in model) {
    compiled[[nm]] <- .hbb_compile_one(
      name          = nm,
      cache_dir     = cache_dir,
      force         = force,
      quiet         = quiet,
      cpp_options   = cpp_options,
      stanc_options = stanc_options
    )
  }

  invisible(compiled)
}


# ============================================================================
# Public: hbb_clear_cache()
# ============================================================================

#' Clear Compiled Stan Model Cache
#'
#' Removes all cached compiled Stan executables and hash files from the
#' persistent user-level cache directory. Also clears the in-memory
#' session cache so that the next call to a fitting function will trigger
#' fresh compilation.
#'
#' @return Invisibly returns `TRUE` if the disk cache existed and was
#'   removed, `FALSE` if no cache was found.
#'
#' @seealso \code{\link{hbb_compile}} to compile Stan models.
#'
#' @export
#' @examples
#' \dontrun{
#' hbb_clear_cache()
#' }
hbb_clear_cache <- function() {
  # Clear in-memory session cache first (always succeeds)
  rm(list = ls(envir = .hbb_cache), envir = .hbb_cache)

  cache_dir <- .hbb_cache_dir()
  if (dir.exists(cache_dir)) {
    unlink(cache_dir, recursive = TRUE)
    cli::cli_alert_success(
      "Cleared hurdlebb Stan model cache at {.path {cache_dir}}."
    )
    invisible(TRUE)
  } else {
    cli::cli_alert_info("No cache directory found at {.path {cache_dir}}.")
    invisible(FALSE)
  }
}


# ============================================================================
# Internal: .hbb_model() — Get a compiled model, auto-compiling if needed
# ============================================================================

#' Retrieve a Compiled Stan Model (internal)
#'
#' Returns a compiled `CmdStanModel` for the requested model name. If not
#' in the session cache (or if the cached exe is gone or the source hash
#' has changed), triggers compilation.
#'
#' @param name A single model name (character scalar).
#' @param cpp_options Named list of C++ options (forwarded to compilation).
#' @return A `CmdStanModel` object.
#' @noRd
.hbb_model <- function(name, cpp_options = list()) {
  assert_string(name)
  assert_choice(name, choices = .hbb_model_names())


  # Fast path: already compiled this session
  if (exists(name, envir = .hbb_cache, inherits = FALSE)) {
    mod <- get(name, envir = .hbb_cache, inherits = FALSE)

    # Guard 1: exe file might have been deleted externally
    exe_ok <- tryCatch(
      file.exists(mod$exe_file()),
      error = function(e) FALSE
    )
    if (!exe_ok) {
      if (!identical(cpp_options, list())) {
        # Cached model may have been compiled with different cpp_options
        # but exe is gone, so recompile.
      }
    } else {
      # Guard 2: source might have changed mid-session (package update)
      stan_file <- .hbb_stan_file(name)
      current_hash <- digest(file = stan_file, algo = "md5")
      hash_file <- file.path(.hbb_cache_dir(), paste0(name, ".md5"))
      cached_hash <- .read_hash(hash_file)

      if (identical(current_hash, cached_hash)) {
        return(mod)
      }
      # Hash mismatch => source changed, fall through to recompile
    }
  }

  # Slow path: compile (or retrieve from persistent cache)
  cache_dir <- .hbb_cache_dir()
  .ensure_dir(cache_dir)

  mod <- .hbb_compile_one(
    name          = name,
    cache_dir     = cache_dir,
    force         = FALSE,
    quiet         = FALSE,
    cpp_options   = cpp_options,
    stanc_options = list()
  )

  mod
}


# ============================================================================
# Internal helpers
# ============================================================================

#' Known model names
#' @return Character vector of valid model names.
#' @noRd
.hbb_model_names <- function() {
  c("hbb_base", "hbb_weighted", "hbb_svc", "hbb_svc_weighted")
}

#' Persistent cache directory path
#' @return Character scalar path.
#' @noRd
.hbb_cache_dir <- function() {
  tools::R_user_dir("hurdlebb", "cache")
}


#' Locate the stan/ directory (works under load_all and installed)
#'
#' @return Absolute path to the stan/ directory.
#' @noRd
.hbb_stan_dir <- function() {
  d <- system.file("stan", package = "hurdlebb")
  if (!nzchar(d) || !dir.exists(d)) {
    abort(
      c(
        "Cannot locate Stan model files for {.pkg hurdlebb}.",
        "i" = "Expected directory: {.path inst/stan/} in the package source.",
        "i" = "If developing, try {.code devtools::load_all()} first."
      ),
      call = NULL
    )
  }
  d
}


#' Resolve path to a single Stan source file
#'
#' @param name Model name without `.stan` extension.
#' @return Absolute path to the `.stan` file.
#' @noRd
.hbb_stan_file <- function(name) {
  f <- file.path(.hbb_stan_dir(), paste0(name, ".stan"))
  if (!file.exists(f)) {
    abort(
      c(
        "Stan file not found: {.path {f}}.",
        "i" = "The package ships {.val {paste0(.hbb_model_names(), '.stan')}}.",
        "i" = "Check that {.path inst/stan/} is complete."
      ),
      call = NULL
    )
  }
  f
}


#' Ensure a directory exists; abort on failure
#' @param dir Path to directory.
#' @noRd
.ensure_dir <- function(dir) {
  if (!dir.exists(dir)) {
    ok <- dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    # dir.create can return FALSE but dir might exist from a race
    if (!ok && !dir.exists(dir)) {
      cli::cli_abort(
        "Cannot create directory {.path {dir}}. Check file permissions."
      )
    }
  }
}


#' Verify CmdStan (not just cmdstanr) is installed
#' @noRd
.assert_cmdstan_available <- function() {
  path <- tryCatch(
    cmdstanr::cmdstan_path(),
    error = function(e) NULL
  )
  if (is.null(path)) {
    cli::cli_abort(c(
      "CmdStan is not installed.",
      "i" = "Run {.code cmdstanr::install_cmdstan()} to install it.",
      "i" = "See {.url https://mc-stan.org/cmdstanr/articles/cmdstanr.html} \\
             for setup instructions."
    ))
  }
}


#' Read a hash file safely
#'
#' @param hash_file Path to the `.md5` file.
#' @return Character scalar (the hash), or `""` if the file cannot be read.
#' @noRd
.read_hash <- function(hash_file) {
  if (!file.exists(hash_file)) return("")
  tryCatch(
    {
      lines <- readLines(hash_file, n = 1L, warn = FALSE)
      if (length(lines) == 0L) "" else trimws(lines[[1L]])
    },
    error = function(e) ""
  )
}


#' Check whether a compiled executable exists (cross-platform)
#'
#' On Windows, CmdStan appends `.exe` to compiled models; on Unix the
#' executable has no extension.
#'
#' @param base_path Path without OS-specific extension.
#' @return Logical.
#' @noRd
.hbb_exe_exists <- function(base_path) {
  file.exists(base_path) || file.exists(paste0(base_path, ".exe"))
}


#' Return the actual exe path (cross-platform)
#'
#' @param base_path Path without OS-specific extension.
#' @return The path that exists, or `NULL` if neither variant exists.
#' @noRd
.hbb_exe_path <- function(base_path) {
  if (file.exists(base_path)) return(base_path)
  with_exe <- paste0(base_path, ".exe")
  if (file.exists(with_exe)) return(with_exe)
  NULL
}


#' Clean up partial compilation artifacts
#'
#' Removes hash and any partial exe / intermediate files left behind
#' by a failed compilation.
#'
#' @param name Model name.
#' @param cache_dir Cache directory path.
#' @noRd
.cleanup_partial <- function(name, cache_dir) {
  extensions <- c("", ".exe", ".hpp", ".o", ".md5", ".lock")
  for (ext in extensions) {
    f <- file.path(cache_dir, paste0(name, ext))
    if (file.exists(f)) unlink(f)
  }
}


# ============================================================================
# Internal: File-lock mechanism for concurrent session safety
# ============================================================================

#' Acquire a file lock (advisory, non-blocking with polling)
#'
#' Creates a lock file containing PID and timestamp. If a lock already
#' exists, waits up to `timeout_secs` polling every `poll_secs`.
#' Automatically breaks stale locks (> `stale_secs` old) and locks
#' from dead processes.
#'
#' @param lock_file Path to the lock file.
#' @param name Model name (for messages).
#' @param timeout_secs Maximum seconds to wait.
#' @param stale_secs Locks older than this are considered stale.
#' @param poll_secs Polling interval.
#' @noRd
.acquire_lock <- function(lock_file, name,
                          timeout_secs = 300,
                          stale_secs = 600,
                          poll_secs = 2) {
  elapsed <- 0
  warned <- FALSE

  while (file.exists(lock_file)) {
    # Is the lock stale by age?
    lock_age <- as.numeric(
      difftime(Sys.time(), file.mtime(lock_file), units = "secs")
    )
    if (lock_age > stale_secs) {
      cli::cli_alert_warning(
        "Breaking stale lock for {.val {name}} (age: {round(lock_age)}s)."
      )
      unlink(lock_file)
      break
    }

    # Is the owning process dead? (Unix only)
    if (.Platform$OS.type == "unix") {
      lock_pid <- .read_lock_pid(lock_file)
      if (!is.na(lock_pid) && !.pid_alive(lock_pid)) {
        cli::cli_alert_warning(
          "Breaking orphaned lock for {.val {name}} (PID {lock_pid} is dead)."
        )
        unlink(lock_file)
        break
      }
    }

    # Timeout check
    if (elapsed >= timeout_secs) {
      cli::cli_abort(c(
        "Timed out waiting {timeout_secs}s for lock on {.val {name}}.",
        "i" = "Another R session may be compiling this model.",
        "i" = "If no other session is compiling, delete {.path {lock_file}} \\
               and retry."
      ))
    }

    if (!warned) {
      cli::cli_alert_info(
        "Waiting for another session to finish compiling {.val {name}}..."
      )
      warned <- TRUE
    }

    Sys.sleep(poll_secs)
    elapsed <- elapsed + poll_secs
  }

  # Write lock file: PID on line 1, timestamp on line 2
  tryCatch(
    writeLines(
      c(as.character(Sys.getpid()), format(Sys.time())),
      lock_file
    ),
    error = function(e) {
      # Non-fatal: we can still compile, just without lock protection
      cli::cli_alert_warning(
        "Could not create lock file: {conditionMessage(e)}"
      )
    }
  )
}


#' Release a file lock (only if we own it)
#' @param lock_file Path to the lock file.
#' @noRd
.release_lock <- function(lock_file) {
  if (!file.exists(lock_file)) return(invisible(NULL))

  lock_pid <- .read_lock_pid(lock_file)
  if (is.na(lock_pid) || lock_pid == Sys.getpid()) {
    unlink(lock_file)
  }
  invisible(NULL)
}


#' Read the PID from a lock file
#' @param lock_file Path to lock file.
#' @return Integer PID, or NA if unreadable.
#' @noRd
.read_lock_pid <- function(lock_file) {
  if (!file.exists(lock_file)) return(NA_integer_)
  tryCatch(
    {
      lines <- readLines(lock_file, n = 1L, warn = FALSE)
      if (length(lines) == 0L) return(NA_integer_)
      suppressWarnings(as.integer(lines[[1L]]))
    },
    error = function(e) NA_integer_
  )
}


#' Check if a PID is alive (Unix only)
#' @param pid Integer PID.
#' @return Logical.
#' @noRd
.pid_alive <- function(pid) {
  tryCatch(
    {
      # signal 0 tests existence without actually signalling
      res <- system2("kill", args = c("-0", as.character(pid)),
                     stdout = FALSE, stderr = FALSE)
      res == 0L
    },
    error = function(e) TRUE  # assume alive if we cannot check
  )
}


# ============================================================================
# Internal: .hbb_compile_one() — Compile/load a single Stan model
# ============================================================================

#' Compile or load a single Stan model from cache (internal)
#'
#' Handles hash checking, lock-file coordination, compilation, and session
#' caching for one model.
#'
#' @param name Model name (e.g., `"hbb_base"`).
#' @param cache_dir Path to the persistent cache directory.
#' @param force Logical. Force recompilation.
#' @param quiet Logical. Suppress messages.
#' @param cpp_options List of C++ options.
#' @param stanc_options List of stanc3 options.
#' @return A `CmdStanModel` object (also stored in `.hbb_cache`).
#' @noRd
.hbb_compile_one <- function(name,
                             cache_dir,
                             force,
                             quiet,
                             cpp_options,
                             stanc_options) {
  stan_file <- .hbb_stan_file(name)
  hash_file <- file.path(cache_dir, paste0(name, ".md5"))
  exe_base  <- file.path(cache_dir, name)
  lock_file <- file.path(cache_dir, paste0(name, ".lock"))

  # Compute current source hash
  current_hash <- digest(file = stan_file, algo = "md5")
  cached_hash  <- .read_hash(hash_file)

  # -- Fast path: cache hit --------------------------------------------------
  needs_compile <- force
  if (!needs_compile) {
    if (!identical(current_hash, cached_hash)) {
      needs_compile <- TRUE
    } else if (!.hbb_exe_exists(exe_base)) {
      # Hash matches but exe is gone (deleted externally)
      needs_compile <- TRUE
    }
  }

  if (!needs_compile) {
    # Load from cache
    if (!quiet) {
      cli::cli_alert_info("Model {.val {name}} loaded from cache.")
    }
    mod <- .load_cached_model(stan_file, exe_base)
    if (!is.null(mod)) {
      assign(name, mod, envir = .hbb_cache)
      return(mod)
    }
    # Loading failed (corrupted exe?) — fall through to recompile
    if (!quiet) {
      cli::cli_alert_warning(
        "Cached executable for {.val {name}} could not be loaded. Recompiling."
      )
    }
  }

  # -- Slow path: compile ----------------------------------------------------
  if (!quiet) {
    cli::cli_alert_info("Compiling Stan model {.val {name}}...")
  }

  # Acquire lock (prevents concurrent compilation of the same model)
  .acquire_lock(lock_file, name = name)
  on.exit(.release_lock(lock_file), add = TRUE)

  # Re-check cache after acquiring lock: another session may have compiled
  # while we were waiting.
  if (!force) {
    recheck_hash <- .read_hash(hash_file)
    if (identical(current_hash, recheck_hash) && .hbb_exe_exists(exe_base)) {
      mod <- .load_cached_model(stan_file, exe_base)
      if (!is.null(mod)) {
        if (!quiet) {
          cli::cli_alert_success(
            "Model {.val {name}} was compiled by another session."
          )
        }
        assign(name, mod, envir = .hbb_cache)
        return(mod)
      }
    }
  }

  # Actually compile
  mod <- tryCatch(
    {
      cmdstanr::cmdstan_model(
        stan_file      = stan_file,
        dir            = cache_dir,
        compile        = TRUE,
        quiet          = quiet,
        cpp_options    = cpp_options,
        stanc_options  = stanc_options
      )
    },
    error = function(e) {
      .cleanup_partial(name, cache_dir)
      cli::cli_abort(
        c(
          "Stan compilation failed for {.val {name}}.",
          "x" = conditionMessage(e),
          "i" = "Verify CmdStan is installed: {.code cmdstanr::cmdstan_path()}",
          "i" = "Try {.code cmdstanr::install_cmdstan()} to reinstall."
        ),
        parent = e
      )
    }
  )

  # Verify that the executable was actually produced
  if (!.hbb_exe_exists(exe_base)) {
    .cleanup_partial(name, cache_dir)
    cli::cli_abort(
      "Compilation of {.val {name}} appeared to succeed, \\
       but no executable was found in {.path {cache_dir}}. \\
       Try {.code hbb_clear_cache()} and then {.code hbb_compile(force = TRUE)}."
    )
  }

  # Write hash ONLY after successful compilation + exe verification
  tryCatch(
    writeLines(current_hash, hash_file),
    error = function(e) {
      cli::cli_alert_warning(
        "Could not write hash file: {conditionMessage(e)}"
      )
    }
  )

  if (!quiet) {
    cli::cli_alert_success("Model {.val {name}} compiled successfully.")
  }

  assign(name, mod, envir = .hbb_cache)
  mod
}


#' Safely load a CmdStanModel from a cached executable
#'
#' @param stan_file Path to the Stan source.
#' @param exe_base Base path to the executable (without OS extension).
#' @return A `CmdStanModel`, or `NULL` on failure.
#' @noRd
.load_cached_model <- function(stan_file, exe_base) {
  exe_path <- .hbb_exe_path(exe_base)
  if (is.null(exe_path)) return(NULL)

  tryCatch(
    cmdstanr::cmdstan_model(
      stan_file = stan_file,
      exe_file  = exe_path,
      compile   = FALSE
    ),
    error = function(e) NULL
  )
}
