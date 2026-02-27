# Compile Stan Models for hurdlebb

Compiles (or retrieves from cache) the Stan models shipped with the
package. Compiled executables are stored in a persistent user-level
cache directory so that compilation happens only once per Stan source
version.

## Usage

``` r
hbb_compile(
  model = NULL,
  force = FALSE,
  quiet = FALSE,
  cpp_options = list(),
  stanc_options = list()
)
```

## Arguments

- model:

  Character vector of model names to compile, or `NULL` (default) to
  compile all four models. Valid names: `"hbb_base"`, `"hbb_weighted"`,
  `"hbb_svc"`, `"hbb_svc_weighted"`.

- force:

  Logical. If `TRUE`, recompile even when the cached executable is
  up-to-date. Default `FALSE`.

- quiet:

  Logical. If `TRUE`, suppress informational messages during
  compilation. Default `FALSE`.

- cpp_options:

  Named list of C++ compilation options passed to
  [`cmdstanr::cmdstan_model()`](https://mc-stan.org/cmdstanr/reference/cmdstan_model.html).
  For example, `list(stan_threads = TRUE)` enables within-chain
  threading. Default [`list()`](https://rdrr.io/r/base/list.html).

- stanc_options:

  Named list of stanc3 transpiler options. Default
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

Invisibly returns a named list of compiled `CmdStanModel` objects.

## Details

You can call this function explicitly, or let the main fitting function
trigger compilation automatically on first use.

### Cache location

Compiled executables are stored under
`tools::R_user_dir("hurdlebb", "cache")`, following R \>= 4.0
conventions. Each model produces two files in the cache directory: the
compiled executable and a `.md5` file recording the MD5 digest of the
Stan source at compilation time.

### Hash-based invalidation

When a Stan source file changes (e.g., after a package update), the MD5
hash no longer matches the cached `.md5` file, triggering automatic
recompilation.

### Thread safety

A simple file-based lock mechanism prevents two concurrent R sessions
from compiling the same model simultaneously. If a lock file is
detected, the function polls and waits (with a configurable timeout).
Stale locks from dead processes are broken automatically.

## See also

[`hbb_clear_cache`](https://joonho112.github.io/hurdlebb/reference/hbb_clear_cache.md)
to remove all cached executables.

## Examples

``` r
if (FALSE) { # \dontrun{
# Compile all models (one-time operation)
hbb_compile()

# Compile only the base model with threading support
hbb_compile("hbb_base", cpp_options = list(stan_threads = TRUE))

# Force recompilation after CmdStan upgrade
hbb_compile(force = TRUE)
} # }
```
