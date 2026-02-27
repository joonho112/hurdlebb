# Compute P(Y = 0) under the Beta-Binomial distribution

Computes the probability of observing zero under a \\\text{BetaBin}(n,
\mu, \kappa)\\ distribution using the stable lgamma formula.

## Usage

``` r
compute_p0(n, mu, kappa)
```

## Arguments

- n:

  Integer vector of trial sizes (\\n \ge 0\\).

- mu:

  Numeric vector of means, each in \\\[0, 1\]\\.

- kappa:

  Numeric vector of concentrations, each \\\> 0\\.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\.

## Details

The formula is \$\$ p_0 = \frac{B(b + n,\\ a)}{B(b,\\ a)} =
\exp\bigl\[\Gamma(b + n) + \ln\Gamma(\kappa) - \ln\Gamma(b) -
\ln\Gamma(\kappa + n)\bigr\], \$\$ where \\a = \mu\kappa\\ and \\b =
(1 - \mu)\kappa\\.

**Boundary handling:**

- \\\mu = 0\\: Returns 1 (point mass at zero).

- \\\mu = 1\\: Returns 0 (point mass at \\n\\).

- \\n = 0\\: Returns 1 (trivially, no trials).

The result is clamped to \\\[0, 1\]\\ to guard against floating-point
drift.

## References

Ghosal, S., Ghosh, S., and Moores, M. (2020). “Hierarchical
beta-binomial models for batch effects in cytometry data.” *Journal of
the Royal Statistical Society: Series A*, **183**(4), 1579–1601.

## See also

Other distributions:
[`compute_ztbb_mean()`](https://joonho112.github.io/hurdlebb/reference/compute_ztbb_mean.md),
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
# Small n: compare with dbetabinom
compute_p0(5, mu = 0.3, kappa = 10)
#> [1] 0.2307692
dbetabinom(0, n = 5, mu = 0.3, kappa = 10)
#> [1] 0.2307692

# Boundary cases
compute_p0(5, mu = 0, kappa = 10)    # 1
#> [1] 1
compute_p0(5, mu = 1, kappa = 10)    # 0
#> [1] 0
compute_p0(0, mu = 0.3, kappa = 10)  # 1
#> [1] 1
```
