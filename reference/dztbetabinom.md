# Probability mass function of the Zero-Truncated Beta-Binomial

Computes the PMF of the zero-truncated Beta-Binomial distribution for
\\y \in \\1, 2, \ldots, n\\\\.

## Usage

``` r
dztbetabinom(y, n, mu, kappa, log = FALSE)
```

## Arguments

- y:

  Integer vector of observed counts (\\y \ge 1\\).

- n:

  Integer vector of trial sizes (\\n \ge 1\\).

- mu:

  Numeric vector of means, each in \\(\varepsilon, 1 - \varepsilon)\\
  where \\\varepsilon\\ is machine epsilon.

- kappa:

  Numeric vector of concentrations, each \\\> 0\\.

- log:

  Logical; if `TRUE`, return log-probabilities.

## Value

A numeric vector of (log-)probabilities.

## Details

The zero-truncated PMF is \$\$ f\_{\text{ZT}}(y \mid n, \mu, \kappa) =
\frac{f\_{\text{BB}}(y \mid n, \mu, \kappa)}{1 - p_0}, \quad y \in \\1,
\ldots, n\\, \$\$ where \\p_0 = P(Y = 0)\\ under the untruncated
Beta-Binomial.

Numerical stability is ensured by computing in the log domain and using
`log1mexp()` for \\\log(1 - p_0)\\.

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
[`hurdle_mean()`](https://joonho112.github.io/hurdlebb/reference/hurdle_mean.md),
[`hurdle_variance()`](https://joonho112.github.io/hurdlebb/reference/hurdle_variance.md),
[`pbetabinom()`](https://joonho112.github.io/hurdlebb/reference/pbetabinom.md),
[`rbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rbetabinom.md),
[`rhurdle_betabinom()`](https://joonho112.github.io/hurdlebb/reference/rhurdle_betabinom.md),
[`rztbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rztbetabinom.md)

## Examples

``` r
# ZT-BB PMF (support starts at 1)
dztbetabinom(1:5, n = 5, mu = 0.3, kappa = 10)
#> [1] 0.40909091 0.32727273 0.18181818 0.06818182 0.01363636

# Sums to 1 over the support
sum(dztbetabinom(1:5, n = 5, mu = 0.3, kappa = 10))
#> [1] 1

# Log scale
dztbetabinom(1:5, n = 5, mu = 0.3, kappa = 10, log = TRUE)
#> [1] -0.8938179 -1.1169614 -1.7047481 -2.6855773 -4.2950153
```
