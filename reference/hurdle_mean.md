# Expected value under the Hurdle Beta-Binomial

Computes the unconditional mean \$\$ E\[Y\] = q \cdot \frac{n \mu}{1 -
p_0}, \$\$ where \\p_0 = P(Y = 0)\\ under the (untruncated)
Beta-Binomial.

## Usage

``` r
hurdle_mean(n, q, mu, kappa)
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

A numeric vector of expected values.

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
[`hurdle_variance()`](https://joonho112.github.io/hurdlebb/reference/hurdle_variance.md),
[`pbetabinom()`](https://joonho112.github.io/hurdlebb/reference/pbetabinom.md),
[`rbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rbetabinom.md),
[`rhurdle_betabinom()`](https://joonho112.github.io/hurdlebb/reference/rhurdle_betabinom.md),
[`rztbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rztbetabinom.md)

## Examples

``` r
hurdle_mean(n = 10, q = 0.7, mu = 0.3, kappa = 5)
#> [1] 2.446688

# q = 0 always gives 0
hurdle_mean(n = 10, q = 0, mu = 0.3, kappa = 5)
#> [1] 0
```
