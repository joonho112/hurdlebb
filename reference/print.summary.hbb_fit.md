# Print Method for summary.hbb_fit Objects

Displays a structured summary of the hurdle Beta-Binomial model,
including model information, fixed-effect estimates separated by margin
with confidence intervals, dispersion parameter, MCMC diagnostics, and
design effect ratios (when sandwich SEs are used).

## Usage

``` r
# S3 method for class 'summary.hbb_fit'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `"summary.hbb_fit"` returned by
  [`summary.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/summary.hbb_fit.md).

- digits:

  Integer; number of significant digits for numeric output. Default is
  3.

- ...:

  Additional arguments (currently unused).

## Value

The object `x`, invisibly.

## Details

The entire method is wrapped in a top-level `tryCatch` so that it never
crashes, even if the summary object is corrupted or incomplete.

## See also

[`summary.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/summary.hbb_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
s <- summary(fit)
print(s, digits = 4)
} # }
```
