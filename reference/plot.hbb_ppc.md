# Plot Method for hbb_ppc Objects

Each panel shows:

- A histogram of the \\M\\ posterior predictive draws of the test
  statistic \\T(\mathbf{y}^{\mathrm{rep}})\\.

- A solid red vertical line at the observed value \\T(\mathbf{y})\\,
  annotated with its numeric value.

- Dashed green vertical lines at the lower and upper bounds of the
  \\(1-\alpha)\\ posterior predictive interval.

A pass/fail indicator appears in the panel subtitle.

## Usage

``` r
# S3 method for class 'hbb_ppc'
plot(x, type = NULL, ...)
```

## Arguments

- x:

  An object of class `"hbb_ppc"` returned by
  [`ppc`](https://joonho112.github.io/hurdlebb/reference/ppc.md).

- type:

  Character string or `NULL`; which statistic to plot. If `NULL`
  (default), uses `x$type`. Otherwise one of `"both"`, `"zero_rate"`,
  `"it_share"`.

- ...:

  Additional arguments (currently unused).

## Value

A `ggplot` object (single-statistic) or a `patchwork` object
(`type = "both"`), returned invisibly.

## Details

Produces a histogram of the posterior predictive draws for each
requested test statistic, overlaid with a vertical line at the observed
value. When `type = "both"`, the two panels are combined using
patchwork.

## ggplot2 Requirement

This function requires ggplot2 (and patchwork for `type = "both"`). Both
are in `Suggests`; they will be checked via
[`rlang::check_installed()`](https://rlang.r-lib.org/reference/is_installed.html)
and the user will be prompted to install them if absent.

## References

Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A.
(2019). Visualisation in Bayesian workflow. *Journal of the Royal
Statistical Society: Series A*, **182**(2), 389–402.

## See also

[`ppc`](https://joonho112.github.io/hurdlebb/reference/ppc.md),
[`print.hbb_ppc`](https://joonho112.github.io/hurdlebb/reference/print.hbb_ppc.md)

Other model-checking:
[`hbb_loo_compare()`](https://joonho112.github.io/hurdlebb/reference/hbb_loo_compare.md),
[`loo.hbb_fit()`](https://joonho112.github.io/hurdlebb/reference/loo.hbb_fit.md),
[`ppc()`](https://joonho112.github.io/hurdlebb/reference/ppc.md),
[`print.hbb_loo_compare()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_loo_compare.md),
[`print.hbb_ppc()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_ppc.md)
