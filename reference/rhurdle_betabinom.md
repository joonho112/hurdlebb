# Random generation from the Hurdle Beta-Binomial

Generates random draws from the two-part hurdle model. First draws \\Z
\sim \text{Bernoulli}(q)\\, then for \\Z = 1\\ draws \\Y \sim
\text{ZT-BetaBin}(n, \mu, \kappa)\\.

## Usage

``` r
rhurdle_betabinom(nn, n, q, mu, kappa)
```

## Arguments

- nn:

  Number of draws to generate (positive integer).

- n:

  Integer vector of trial sizes (\\n \ge 1\\).

- q:

  Numeric vector of participation probabilities, each in \\\[0, 1\]\\.

- mu:

  Numeric vector of intensity means in \\(\varepsilon, 1 -
  \varepsilon)\\.

- kappa:

  Numeric vector of concentrations (\\\> 0\\).

## Value

An integer vector of length `nn`.

## Details

All parameter vectors (`n`, `q`, `mu`, `kappa`) are recycled to length
`nn`. The function is fully vectorised.

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
[`rztbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rztbetabinom.md)

## Examples

``` r
set.seed(42)
x <- rhurdle_betabinom(1000, n = 10, q = 0.7, mu = 0.3, kappa = 5)
mean(x == 0)  # approximately 0.3
#> [1] 0.293
hist(x, breaks = -0.5:10.5)

```
