# Summary Method for hbb_fit Objects

Produces a comprehensive summary of a fitted hurdle Beta-Binomial model,
including fixed-effect estimates with standard errors and confidence
intervals, dispersion parameter summary, MCMC diagnostics, and
optionally design effect ratios from the sandwich variance.

## Usage

``` r
# S3 method for class 'hbb_fit'
summary(object, sandwich = NULL, level = 0.95, ...)
```

## Arguments

- object:

  An object of class `"hbb_fit"` returned by
  [`hbb`](https://joonho112.github.io/hurdlebb/reference/hbb.md).

- sandwich:

  An optional object of class `"hbb_sandwich"` returned by
  [`sandwich_variance`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md).
  If provided, produces sandwich-based Wald inference. If `NULL`
  (default), produces posterior-based inference.

- level:

  Numeric scalar in \\(0, 1)\\. Confidence/credible level. Default is
  `0.95`.

- ...:

  Currently unused; included for S3 method consistency.

## Value

An S3 object of class `"summary.hbb_fit"` containing:

- `fixed_effects`:

  Data frame with columns: `parameter`, `estimate`, `se`, `ci_lower`,
  `ci_upper`, `rhat`, `ess_bulk`. One row per fixed-effect parameter
  (\\D = 2P + 1\\).

- `dispersion`:

  Named list with `kappa_hat` (point estimate on the natural scale),
  `kappa_se` (delta-method SE or posterior SD), and `kappa_ci` (2-vector
  confidence limits on the \\\kappa\\ scale).

- `random_effects`:

  `NULL` for base/weighted models. For SVC models, a list with element
  `tau` (posterior mean of the random-effect scale parameters).

- `diagnostics`:

  Named list with `n_divergent`, `n_max_treedepth`, `ebfmi`, `max_rhat`,
  `min_ess_bulk`, `min_ess_tail`.

- `model_info`:

  Named list with `N`, `P`, `S` (or `NULL`), `zero_rate`, `model_type`,
  `formula`.

- `sandwich_used`:

  Logical: whether sandwich SEs were used.

- `DER`:

  Named numeric vector of Design Effect Ratios if `sandwich` was
  provided, `NULL` otherwise.

- `level`:

  The confidence level used.

- `call`:

  The original model call.

## Inference mode

When `sandwich` is supplied, standard errors and confidence intervals
are computed from the sandwich variance \\V\_{\mathrm{sand}}\\
(Wald-type intervals): \$\$ \mathrm{CI}\_{1-\alpha}(\theta_p) =
\hat\theta_p \pm z\_{(1+\mathrm{level})/2}\\
\sqrt{V\_{\mathrm{sand},pp}}. \$\$ When `sandwich` is `NULL`, standard
errors are posterior standard deviations and intervals are
quantile-based credible intervals from the MCMC draws.

## Dispersion

The dispersion parameter \\\kappa\\ controls overdispersion relative to
the Binomial. It is estimated on the log scale (`log_kappa`) and
back-transformed via the delta method: \$\$\mathrm{SE}(\hat\kappa) =
\hat\kappa \cdot \mathrm{SE}(\widehat{\log\kappa}).\$\$ The confidence
interval on the \\\kappa\\ scale is obtained by exponentiating the
interval for \\\log\kappa\\.

## References

Williams, M. R. and Savitsky, T. D. (2021). Uncertainty estimation for
pseudo-Bayesian inference under complex sampling. *International
Statistical Review*, **89**(1), 72–107.
[doi:10.1111/insr.12376](https://doi.org/10.1111/insr.12376)

## See also

[`print.summary.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/print.summary.hbb_fit.md)
for the print method,
[`coef.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/coef.hbb_fit.md),
[`vcov.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/vcov.hbb_fit.md),
[`sandwich_variance`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit  <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
sand <- sandwich_variance(fit)

# Summary with sandwich SEs (recommended for survey data)
s <- summary(fit, sandwich = sand, level = 0.95)
print(s)

# Summary with posterior SDs (for unweighted models)
s0 <- summary(fit)
print(s0)
} # }
```
