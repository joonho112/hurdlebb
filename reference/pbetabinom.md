# Cumulative distribution function of the Beta-Binomial distribution

Computes the CDF \\P(Y \le q)\\ by summing the PMF from 0 to `floor(q)`,
element-wise.

## Usage

``` r
pbetabinom(q, n, mu, kappa, lower.tail = TRUE)
```

## Arguments

- q:

  Numeric vector of quantiles.

- n:

  Integer vector of trial sizes (\\n \ge 0\\).

- mu:

  Numeric vector of means, each in \\\[0, 1\]\\.

- kappa:

  Numeric vector of concentrations, each \\\> 0\\.

- lower.tail:

  Logical; if `TRUE` (default), returns \\P(Y \le q)\\; otherwise \\P(Y
  \> q)\\.

## Value

A numeric vector of probabilities.

## Details

Numerical stability is maintained by performing the summation in the log
domain using the log-sum-exp trick.

## References

Ghosal, S., Ghosh, S., and Moores, M. (2020). “Hierarchical
beta-binomial models for batch effects in cytometry data.” *Journal of
the Royal Statistical Society: Series A*, **183**(4), 1579–1601.

## See also

Other distributions:
[`compute_p0()`](https://joonho112.github.io/hurdlebb/reference/compute_p0.md),
[`compute_ztbb_mean()`](https://joonho112.github.io/hurdlebb/reference/compute_ztbb_mean.md),
[`dbetabinom()`](https://joonho112.github.io/hurdlebb/reference/dbetabinom.md),
[`dhurdle_betabinom()`](https://joonho112.github.io/hurdlebb/reference/dhurdle_betabinom.md),
[`dztbetabinom()`](https://joonho112.github.io/hurdlebb/reference/dztbetabinom.md),
[`hurdle_mean()`](https://joonho112.github.io/hurdlebb/reference/hurdle_mean.md),
[`hurdle_variance()`](https://joonho112.github.io/hurdlebb/reference/hurdle_variance.md),
[`rbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rbetabinom.md),
[`rhurdle_betabinom()`](https://joonho112.github.io/hurdlebb/reference/rhurdle_betabinom.md),
[`rztbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rztbetabinom.md)

## Examples

``` r
# CDF at each value
pbetabinom(0:5, n = 5, mu = 0.3, kappa = 10)
#> [1] 0.2307692 0.5454545 0.7972028 0.9370629 0.9895105 1.0000000

# Upper tail
pbetabinom(2, n = 5, mu = 0.3, kappa = 10, lower.tail = FALSE)
#> [1] 0.2027972
```
