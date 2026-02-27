# Conditional mean of the zero-truncated Beta-Binomial

Returns the conditional mean of the zero-truncated Beta-Binomial, \\E\[Y
\mid Y \> 0\] = n \mu / (1 - p_0)\\.

## Usage

``` r
compute_ztbb_mean(n, mu, kappa)
```

## Arguments

- n:

  Integer vector of trial sizes (\\n \ge 1\\).

- mu:

  Numeric vector of means in \\(\varepsilon, 1 - \varepsilon)\\.

- kappa:

  Numeric vector of concentrations (\\\> 0\\).

## Value

A numeric vector of conditional means.

## Details

Requires \\n \ge 1\\ and \\\mu \in (\varepsilon, 1 - \varepsilon)\\ so
that the zero-truncated distribution is well-defined (i.e., \\p_0 \<
1\\).

## References

Ghosal, S., Ghosh, S., and Moores, M. (2020). “Hierarchical
beta-binomial models for batch effects in cytometry data.” *Journal of
the Royal Statistical Society: Series A*, **183**(4), 1579–1601.

## See also

Other distributions:
[`compute_p0()`](https://joonho112.github.io/hurdlebb/reference/compute_p0.md),
[`dbetabinom()`](https://joonho112.github.io/hurdlebb/reference/dbetabinom.md),
[`dhurdle_betabinom()`](https://joonho112.github.io/hurdlebb/reference/dhurdle_betabinom.md),
[`dztbetabinom()`](https://joonho112.github.io/hurdlebb/reference/dztbetabinom.md),
[`hurdle_mean()`](https://joonho112.github.io/hurdlebb/reference/hurdle_mean.md),
[`hurdle_variance()`](https://joonho112.github.io/hurdlebb/reference/hurdle_variance.md),
[`pbetabinom()`](https://joonho112.github.io/hurdlebb/reference/pbetabinom.md),
[`rbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rbetabinom.md),
[`rhurdle_betabinom()`](https://joonho112.github.io/hurdlebb/reference/rhurdle_betabinom.md),
[`rztbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rztbetabinom.md)

## Examples

``` r
# Compare analytical mean with empirical mean
compute_ztbb_mean(10, mu = 0.3, kappa = 5)
#> [1] 3.495269

set.seed(1)
mean(rztbetabinom(10000, n = 10, mu = 0.3, kappa = 5))
#> [1] 3.4795
```
