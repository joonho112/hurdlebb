# Extract the AME Decomposition Table

Extracts the decomposition table from an `hbb_ame` object, showing the
extensive and intensive components of the Average Marginal Effect for
each non-intercept covariate.

## Usage

``` r
ame_decomposition(ame_result)
```

## Arguments

- ame_result:

  An object of class `"hbb_ame"` returned by
  [`ame`](https://joonho112.github.io/hurdlebb/reference/ame.md).

## Value

A data frame with one row per non-intercept covariate and columns:
`covariate`, `ext_ame`, `ext_ci_lo`, `ext_ci_hi`, `int_ame`,
`int_ci_lo`, `int_ci_hi`, `total_ame`, `total_ci_lo`, `total_ci_hi`,
`ext_share`, `int_share`, `sign_pattern`.

## See also

[`ame`](https://joonho112.github.io/hurdlebb/reference/ame.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ame_result <- ame(fit, cholesky = chol)
decomp <- ame_decomposition(ame_result)
print(decomp)
} # }
```
