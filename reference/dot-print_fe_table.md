# Print a subset of the fixed-effects table with aligned columns

Formats and prints a fixed-effects data.frame with right-aligned numeric
columns using `sprintf`. Used internally by `print.summary.hbb_fit` to
display the extensive and intensive margin tables separately.

## Usage

``` r
.print_fe_table(fe_sub, digits)
```

## Arguments

- fe_sub:

  A data.frame (subset of fixed_effects from `summary.hbb_fit`).

- digits:

  Integer; number of decimal places for numeric output.

## Value

Invisible `NULL` (side effect: prints to console).
