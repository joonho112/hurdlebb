# Print Method for hbb_loo_compare Objects

The display is structured into three sections:

1.  A model ranking table with \\\mathrm{ELPD}\\,
    \\\Delta\mathrm{ELPD}\\, \\\mathrm{SE}(\Delta)\\, and
    \\\mathrm{LOOIC}\\.

2.  A pairwise z-ratio table for consecutive model pairs.

3.  A Pareto-k diagnostic summary (counts by category).

## Usage

``` r
# S3 method for class 'hbb_loo_compare'
print(x, digits = 1L, ...)
```

## Arguments

- x:

  An object of class `"hbb_loo_compare"` returned by
  [`hbb_loo_compare`](https://joonho112.github.io/hurdlebb/reference/hbb_loo_compare.md).

- digits:

  Integer; decimal places for numeric values. Default is `1`.

- ...:

  Additional arguments (currently unused).

## Value

The object `x`, invisibly.

## Details

Displays a formatted LOO comparison table including ELPD differences,
standard errors, and pairwise z-statistics.

## See also

[`hbb_loo_compare`](https://joonho112.github.io/hurdlebb/reference/hbb_loo_compare.md),
[`loo.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/loo.hbb_fit.md)

Other model-checking:
[`hbb_loo_compare()`](https://joonho112.github.io/hurdlebb/reference/hbb_loo_compare.md),
[`loo.hbb_fit()`](https://joonho112.github.io/hurdlebb/reference/loo.hbb_fit.md),
[`plot.hbb_ppc()`](https://joonho112.github.io/hurdlebb/reference/plot.hbb_ppc.md),
[`ppc()`](https://joonho112.github.io/hurdlebb/reference/ppc.md),
[`print.hbb_ppc()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_ppc.md)
