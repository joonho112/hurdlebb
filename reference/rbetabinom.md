# Random generation from the Beta-Binomial distribution

Generates random draws via the hierarchical representation: \\p \sim
\text{Beta}(a, b)\\, \\Y \mid p \sim \text{Binomial}(n, p)\\.

## Usage

``` r
rbetabinom(nn, n, mu, kappa)
```

## Arguments

- nn:

  Number of draws to generate (positive integer).

- n:

  Integer vector of trial sizes (\\n \ge 0\\).

- mu:

  Numeric vector of means, each in \\\[0, 1\]\\.

- kappa:

  Numeric vector of concentrations, each \\\> 0\\.

## Value

An integer vector of length `nn`.

## Details

All parameter vectors are recycled to length `nn`. The generation is
fully vectorised (no loop).

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
[`rhurdle_betabinom()`](https://joonho112.github.io/hurdlebb/reference/rhurdle_betabinom.md),
[`rztbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rztbetabinom.md)

## Examples

``` r
set.seed(42)
x <- rbetabinom(1000, n = 10, mu = 0.3, kappa = 5)
hist(x, breaks = -0.5:10.5)

mean(x)  # approximately 3
#> [1] 2.965

# Vectorised parameters
rbetabinom(4, n = c(5, 10), mu = c(0.2, 0.8), kappa = c(3, 20))
#> [1] 1 7 0 5
```
