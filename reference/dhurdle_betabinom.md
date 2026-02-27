# Probability mass function of the Hurdle Beta-Binomial

Computes the PMF of the two-part hurdle model: \$\$ P(Y = 0) = 1 - q,
\quad P(Y = y) = q \cdot f\_{\text{ZT}}(y \mid n, \mu, \kappa), \quad y
\in \\1, \ldots, n\\, \$\$ where \\f\_{\text{ZT}}\\ is the
zero-truncated Beta-Binomial PMF.

## Usage

``` r
dhurdle_betabinom(y, n, q, mu, kappa, log = FALSE)
```

## Arguments

- y:

  Integer vector of observed counts.

- n:

  Integer vector of trial sizes (\\n \ge 1\\).

- q:

  Numeric vector of participation probabilities, each in \\\[0, 1\]\\.

- mu:

  Numeric vector of intensity means, each in \\(\varepsilon, 1 -
  \varepsilon)\\.

- kappa:

  Numeric vector of concentrations (\\\> 0\\).

- log:

  Logical; if `TRUE`, return log-probabilities.

## Value

A numeric vector of (log-)probabilities.

## Details

The structural zero probability \\1 - q\\ is computed as `log1p(-q)` in
the log domain for numerical stability when \\q\\ is near 1.

## References

Ghosal, S., Ghosh, S., and Moores, M. (2020). “Hierarchical
beta-binomial models for batch effects in cytometry data.” *Journal of
the Royal Statistical Society: Series A*, **183**(4), 1579–1601.

## See also

Other distributions:
[`compute_p0()`](https://joonho112.github.io/hurdlebb/reference/compute_p0.md),
[`compute_ztbb_mean()`](https://joonho112.github.io/hurdlebb/reference/compute_ztbb_mean.md),
[`dbetabinom()`](https://joonho112.github.io/hurdlebb/reference/dbetabinom.md),
[`dztbetabinom()`](https://joonho112.github.io/hurdlebb/reference/dztbetabinom.md),
[`hurdle_mean()`](https://joonho112.github.io/hurdlebb/reference/hurdle_mean.md),
[`hurdle_variance()`](https://joonho112.github.io/hurdlebb/reference/hurdle_variance.md),
[`pbetabinom()`](https://joonho112.github.io/hurdlebb/reference/pbetabinom.md),
[`rbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rbetabinom.md),
[`rhurdle_betabinom()`](https://joonho112.github.io/hurdlebb/reference/rhurdle_betabinom.md),
[`rztbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rztbetabinom.md)

## Examples

``` r
# PMF over full support
dhurdle_betabinom(0:5, n = 5, q = 0.7, mu = 0.4, kappa = 8)
#> [1] 0.30000000 0.20157439 0.21708011 0.16600244 0.08872544 0.02661763

# Should sum to 1
sum(dhurdle_betabinom(0:5, n = 5, q = 0.7, mu = 0.4, kappa = 8))
#> [1] 1

# Log scale
dhurdle_betabinom(0:5, n = 5, q = 0.7, mu = 0.4, kappa = 8, log = TRUE)
#> [1] -1.203973 -1.601597 -1.527489 -1.795753 -2.422209 -3.626181
```
