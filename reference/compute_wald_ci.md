# Wald Confidence Intervals from Sandwich Variance

The Wald confidence interval for parameter \\\theta_p\\ is \$\$
\mathrm{CI}\_{1-\alpha}(\theta_p) = \hat\theta_p \pm z\_{1-\alpha/2}\\
\sqrt{V\_{\mathrm{sand},pp}}, \$\$ where \\z\_{1-\alpha/2}\\ is the
standard normal quantile. Because the sandwich variance is
design-consistent under mild regularity conditions (Williams and
Savitsky, 2021, Theorem 3.2), Wald intervals achieve correct frequentist
coverage for the pseudo-true parameter.

## Usage

``` r
compute_wald_ci(theta_hat, V_sand, level = 0.95, param_labels = NULL)
```

## Arguments

- theta_hat:

  Numeric vector of length D containing point estimates (typically
  posterior means).

- V_sand:

  Numeric symmetric positive-(semi)definite matrix of dimension D by D,
  the sandwich variance. Only the diagonal entries are used for marginal
  Wald intervals.

- level:

  Numeric scalar in \\(0,1)\\. Confidence level. Default is `0.95`.

- param_labels:

  Optional character vector of length D giving parameter names. If
  `NULL` (default), names are taken from `names(theta_hat)` or generated
  as `param_1`, `param_2`, etc.

## Value

A data frame with one row per parameter and columns:

- `parameter`:

  Character: parameter label.

- `post_mean`:

  Numeric: point estimate.

- `se`:

  Numeric: sandwich standard error.

- `z_stat`:

  Numeric: Wald z-statistic.

- `p_value`:

  Numeric: two-sided p-value.

- `ci_lo`:

  Numeric: lower confidence limit.

- `ci_hi`:

  Numeric: upper confidence limit.

- `ci_width`:

  Numeric: interval width.

- `significant`:

  Logical: TRUE if p-value is below the nominal significance level.

## Details

Computes Wald confidence intervals using the sandwich variance
\\V\_{\mathrm{sand}}\\. This is a standalone function that does not
require a fitted model object or MCMC draws.

## Recommendation

Wald intervals are recommended as the primary inference device because:

1.  They depend only on the point estimate and \\V\_{\mathrm{sand}}\\,
    not on MCMC sampling variability.

2.  They are algebraically equivalent to the marginal quantiles of the
    Cholesky-corrected draws in large samples.

3.  They avoid the non-trivial Monte Carlo error in tail quantile
    estimation from finite MCMC samples.

## Significance testing

The Wald z-statistic is \\z_p = \hat\theta_p /
\sqrt{V\_{\mathrm{sand},pp}}\\, with two-sided p-value
\\2\\\Phi(-\|z_p\|)\\. A parameter is flagged significant when the
p-value falls below \\1 - \mathrm{level}\\.

## References

Williams, M. R. and Savitsky, T. D. (2021). Uncertainty estimation for
pseudo-Bayesian inference under complex sampling. *International
Statistical Review*, **89**(1), 72–107.
[doi:10.1111/insr.12376](https://doi.org/10.1111/insr.12376)

## See also

[`cholesky_correct`](https://joonho112.github.io/hurdlebb/reference/cholesky_correct.md)
for the full Cholesky recalibration,
[`sandwich_variance`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md)
for obtaining `V_sand`.

## Examples

``` r
# Minimal example with known values
theta_hat <- c(alpha_1 = -0.324, beta_1 = 0.090, log_kappa = 1.92)
V_sand <- diag(c(0.0027, 0.00035, 0.0041))
compute_wald_ci(theta_hat, V_sand, level = 0.95)
#>   parameter post_mean         se    z_stat       p_value       ci_lo      ci_hi
#> 1   alpha_1    -0.324 0.05196152 -6.235383  4.506744e-10 -0.42584272 -0.2221573
#> 2    beta_1     0.090 0.01870829  4.810702  1.504008e-06  0.05333243  0.1266676
#> 3 log_kappa     1.920 0.06403124 29.985362 1.522994e-197  1.79450107  2.0454989
#>     ci_width significant
#> 1 0.20368543        TRUE
#> 2 0.07333514        TRUE
#> 3 0.25099786        TRUE
```
