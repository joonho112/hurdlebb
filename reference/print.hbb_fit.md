# Print Method for hbb_fit Objects

Displays a compact summary of a fitted hurdle Beta-Binomial model,
including the model type, formula, data dimensions, MCMC configuration,
elapsed time, and basic diagnostic indicators. Does **not** print
parameter estimates; use
[`summary()`](https://rdrr.io/r/base/summary.html) for that.

## Usage

``` r
# S3 method for class 'hbb_fit'
print(x, ...)
```

## Arguments

- x:

  An object of class `"hbb_fit"`, as returned by
  [`hbb()`](https://joonho112.github.io/hurdlebb/reference/hbb.md).

- ...:

  Currently unused; included for S3 method consistency.

## Value

Invisibly returns `x`.

## Details

All diagnostic extractions are defensive (wrapped in `tryCatch`) so
printing never fails, even if the underlying CmdStan output files have
been moved or the fit object is incomplete.

## See also

[`hbb()`](https://joonho112.github.io/hurdlebb/reference/hbb.md),
[`is.hbb_fit()`](https://joonho112.github.io/hurdlebb/reference/is.hbb_fit.md)

Other fitting:
[`hbb()`](https://joonho112.github.io/hurdlebb/reference/hbb.md),
[`is.hbb_fit()`](https://joonho112.github.io/hurdlebb/reference/is.hbb_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- hbb(y | trials(n_trial) ~ poverty, data = my_data)
print(fit)
} # }
```
