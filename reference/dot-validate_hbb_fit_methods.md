# Validate that an object is a well-formed hbb_fit for S3 methods

Performs comprehensive checks: class verification, presence of
`hbb_data` with required fields (`N`, `P`, `X`), dimensional consistency
of the design matrix, type checks on `N` and `P`, and presence of the
CmdStanMCMC `fit` slot. Aborts with structured `cli_abort` messages on
failure.

## Usage

``` r
.validate_hbb_fit_methods(x)
```

## Arguments

- x:

  An object to validate.

## Value

Invisible `TRUE`; aborts on failure.
