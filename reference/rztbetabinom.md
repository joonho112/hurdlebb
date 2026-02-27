# Random generation from the Zero-Truncated Beta-Binomial

Generates random draws from the zero-truncated Beta-Binomial via
rejection sampling on the standard Beta-Binomial.

## Usage

``` r
rztbetabinom(nn, n, mu, kappa)
```

## Arguments

- nn:

  Number of draws to generate (positive integer).

- n:

  Integer vector of trial sizes (\\n \ge 1\\), recycled to length `nn`.

- mu:

  Numeric vector of means in \\(\varepsilon, 1 - \varepsilon)\\,
  recycled to length `nn`.

- kappa:

  Numeric vector of concentrations (\\\> 0\\), recycled to length `nn`.

## Value

An integer vector of length `nn`, each element in \\\\1, \ldots,
n_i\\\\.

## Details

The function is vectorised over `n`, `mu`, and `kappa` using a
parameter-grouping strategy. Unique parameter combinations are
identified, and rejection sampling is performed once per group with
smart batch sizing based on the analytical \\p_0\\.

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
[`pbetabinom()`](https://joonho112.github.io/hurdlebb/reference/pbetabinom.md),
[`rbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rbetabinom.md),
[`rhurdle_betabinom()`](https://joonho112.github.io/hurdlebb/reference/rhurdle_betabinom.md)

## Examples

``` r
set.seed(42)
x <- rztbetabinom(1000, n = 10, mu = 0.3, kappa = 5)
min(x)  # always >= 1
#> [1] 1
mean(x)
#> [1] 3.438

# Vectorised parameters
rztbetabinom(6, n = c(5, 10, 20), mu = c(0.2, 0.5, 0.8), kappa = c(3, 10, 50))
#> [1]  5  4 14  2  2 15
```
