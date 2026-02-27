# Number of Observations in an hbb_fit

Extracts the number of provider-level observations \\N\\ used to fit the
hurdle Beta-Binomial model.

## Usage

``` r
# S3 method for class 'hbb_fit'
nobs(object, ...)
```

## Arguments

- object:

  An object of class `"hbb_fit"` returned by
  [`hbb`](https://joonho112.github.io/hurdlebb/reference/hbb.md).

- ...:

  Currently unused; included for S3 method consistency.

## Value

An integer scalar giving the number of observations.

## See also

[`hbb`](https://joonho112.github.io/hurdlebb/reference/hbb.md),
[`coef.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/coef.hbb_fit.md),
[`summary.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/summary.hbb_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
nobs(fit)
} # }
```
