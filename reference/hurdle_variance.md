# Variance under the Hurdle Beta-Binomial

Computes the unconditional variance via the law of total variance: \$\$
\text{Var}\[Y\] = q \cdot V\_{\text{ZT}} + q (1 - q) \cdot
(E\_{\text{ZT}})^2, \$\$ where \\E\_{\text{ZT}} = n\mu / (1 - p_0)\\ is
the zero-truncated mean and \\V\_{\text{ZT}}\\ is the zero-truncated
variance.

## Usage

``` r
hurdle_variance(n, q, mu, kappa)
```

## Arguments

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

A numeric vector of variances (non-negative).

## Details

The zero-truncated variance is derived from the untruncated moments:
\$\$ V\_{\text{ZT}} = \frac{E\_{\text{BB}}\[Y^2\]}{1 - p_0} -
\left(\frac{E\_{\text{BB}}\[Y\]}{1 - p_0}\right)^2, \$\$ where \$\$
E\_{\text{BB}}\[Y^2\] = \text{Var}\_{\text{BB}} + (n\mu)^2 = n\mu(1 -
\mu) \frac{n + \kappa}{1 + \kappa} + n^2 \mu^2. \$\$

The result is clamped to be non-negative via `pmax(., 0)`.

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
[`pbetabinom()`](https://joonho112.github.io/hurdlebb/reference/pbetabinom.md),
[`rbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rbetabinom.md),
[`rhurdle_betabinom()`](https://joonho112.github.io/hurdlebb/reference/rhurdle_betabinom.md),
[`rztbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rztbetabinom.md)

## Examples

``` r
hurdle_variance(n = 10, q = 0.7, mu = 0.3, kappa = 5)
#> [1] 5.635486

# q = 0 always gives 0
hurdle_variance(n = 10, q = 0, mu = 0.3, kappa = 5)
#> [1] 0
```
