# Print Method for hbb_ppc Objects

Displays a compact numerical summary of the posterior predictive check,
showing the observed test statistic, the posterior predictive mean and
median, the credible interval, and whether the observed value falls
within the interval.

## Usage

``` r
# S3 method for class 'hbb_ppc'
print(x, digits = 4L, ...)
```

## Arguments

- x:

  An object of class `"hbb_ppc"` returned by
  [`ppc`](https://joonho112.github.io/hurdlebb/reference/ppc.md).

- digits:

  Integer; number of decimal places for numeric output. Default `4`.

- ...:

  Additional arguments (currently unused).

## Value

The object `x`, invisibly.

## See also

[`ppc`](https://joonho112.github.io/hurdlebb/reference/ppc.md),
[`plot.hbb_ppc`](https://joonho112.github.io/hurdlebb/reference/plot.hbb_ppc.md)

Other model-checking:
[`hbb_loo_compare()`](https://joonho112.github.io/hurdlebb/reference/hbb_loo_compare.md),
[`loo.hbb_fit()`](https://joonho112.github.io/hurdlebb/reference/loo.hbb_fit.md),
[`plot.hbb_ppc()`](https://joonho112.github.io/hurdlebb/reference/plot.hbb_ppc.md),
[`ppc()`](https://joonho112.github.io/hurdlebb/reference/ppc.md),
[`print.hbb_loo_compare()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_loo_compare.md)
