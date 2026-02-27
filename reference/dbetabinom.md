# Probability mass function of the Beta-Binomial distribution

Computes the PMF of the Beta-Binomial distribution parameterised by mean
\\\mu \in \[0, 1\]\\ and concentration \\\kappa \> 0\\. The shape
parameters are \\a = \mu \kappa\\ and \\b = (1 - \mu) \kappa\\.

## Usage

``` r
dbetabinom(y, n, mu, kappa, log = FALSE)
```

## Arguments

- y:

  Integer vector of observed counts.

- n:

  Integer vector of trial sizes (\\n \ge 0\\).

- mu:

  Numeric vector of means, each in \\\[0, 1\]\\.

- kappa:

  Numeric vector of concentrations, each \\\> 0\\.

- log:

  Logical; if `TRUE`, return log-probabilities.

## Value

A numeric vector of (log-)probabilities, same length as the recycled
inputs.

## Details

The PMF is \$\$ f(y \mid n, \mu, \kappa) = \binom{n}{y} \frac{B(y + a,\\
n - y + b)}{B(a, b)}, \quad y \in \\0, 1, \ldots, n\\, \$\$ where \\a =
\mu \kappa\\ and \\b = (1 - \mu) \kappa\\.

**Boundary handling.** When \\\mu = 0\\ the Beta-Binomial degenerates to
a point mass at \\y = 0\\ (all probability on zero). When \\\mu = 1\\ it
degenerates to a point mass at \\y = n\\. These boundaries are handled
with explicit branches because the `lbeta` formula produces `NaN` when
\\a = 0\\ or \\b = 0\\.

All arguments are recycled to common length via `rep_len`.

## References

Ghosal, S., Ghosh, S., and Moores, M. (2020). “Hierarchical
beta-binomial models for batch effects in cytometry data.” *Journal of
the Royal Statistical Society: Series A*, **183**(4), 1579–1601.

## See also

Other distributions:
[`compute_p0()`](https://joonho112.github.io/hurdlebb/reference/compute_p0.md),
[`compute_ztbb_mean()`](https://joonho112.github.io/hurdlebb/reference/compute_ztbb_mean.md),
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
# Standard usage
dbetabinom(0:5, n = 5, mu = 0.3, kappa = 10)
#> [1] 0.23076923 0.31468531 0.25174825 0.13986014 0.05244755 0.01048951
dbetabinom(0:5, n = 5, mu = 0.3, kappa = 10, log = TRUE)
#> [1] -1.466337 -1.156182 -1.379326 -1.967112 -2.947942 -4.557380

# Boundary: mu = 0 gives point mass at y = 0
dbetabinom(0:3, n = 3, mu = 0, kappa = 5)
#> [1] 1 0 0 0

# Boundary: mu = 1 gives point mass at y = n
dbetabinom(0:3, n = 3, mu = 1, kappa = 5)
#> [1] 0 0 0 1
```
