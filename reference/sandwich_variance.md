# Cluster-Robust Sandwich Variance Estimator

Computes the design-corrected sandwich variance-covariance matrix for
the fixed effects of a survey-weighted hurdle Beta-Binomial model. The
estimator uses the stratified cluster-robust formula

## Usage

``` r
sandwich_variance(fit)
```

## Arguments

- fit:

  An object of class `"hbb_fit"`, as returned by
  [`hbb()`](https://joonho112.github.io/hurdlebb/reference/hbb.md). Must
  be a weighted model (`model_type %in% c("weighted", "svc_weighted")`),
  since score generated quantities are only produced by survey-weighted
  Stan models.

## Value

An S3 object of class `"hbb_sandwich"` containing:

- `V_sand`:

  Numeric \\D \times D\\ sandwich variance matrix.

- `H_obs`:

  Numeric \\D \times D\\ observed Fisher information (block-diagonal:
  extensive \\P \times P\\ and intensive+kappa \\(P+1) \times (P+1)\\).

- `H_obs_inv`:

  Numeric \\D \times D\\ inverse of `H_obs`, representing the data-only
  variance (no design correction).

- `J_cluster`:

  Numeric \\D \times D\\ cluster-robust meat matrix.

- `Sigma_MCMC`:

  Numeric \\D \times D\\ posterior covariance of fixed-effect draws (for
  comparison, not used as bread).

- `scores`:

  Numeric \\N \times D\\ posterior mean score matrix.

- `DER`:

  Named numeric vector of length \\D\\. Design Effect Ratios:
  \\\mathrm{DER}\_p = V\_{\mathrm{sand}}\[p,p\] /
  H\_{\mathrm{obs}}^{-1}\[p,p\]\\.

- `param_labels`:

  Character vector of length \\D\\ with human-readable parameter labels.

- `D`:

  Integer. Total fixed-effect dimension.

- `N`:

  Integer. Number of observations.

- `P`:

  Integer. Number of covariates (including intercept).

- `model_type`:

  Character. The model type from `fit`.

- `survey_info`:

  List with survey design summaries: `n_strata`, `n_psu`, `df`,
  `n_singleton`.

- `matrix_diagnostics`:

  List of eigenvalue/condition diagnostics for `H_obs`, `J_cluster`,
  `V_sand`, and `Sigma_MCMC`.

- `nearPD_applied`:

  Logical. Whether nearPD correction was applied to V_sand.

- `call`:

  The matched call.

## Details

\$\$V\_{\mathrm{sand}} = H\_{\mathrm{obs}}^{-1}\\
J\_{\mathrm{cluster}}\\ H\_{\mathrm{obs}}^{-1}\$\$

where \\H\_{\mathrm{obs}}\\ is the block-diagonal observed Fisher
information (analytic for the extensive margin, empirical information
identity for the intensive margin plus dispersion), and
\\J\_{\mathrm{cluster}}\\ is the cluster-robust "meat" matrix that
accounts for within-PSU correlation and unequal survey weights.

## Why not Sigma_MCMC as bread

In hierarchical models with state random effects, the MCMC posterior
covariance of fixed effects absorbs prior and random-effect variance,
leading to \\\Sigma\_{\mathrm{MCMC}} \gg H\_{\mathrm{obs}}^{-1}\\. Using
\\\Sigma\_{\mathrm{MCMC}}\\ as bread yields astronomical Design Effect
Ratios (3000–26000 in the NSECE application). The explicit
\\H\_{\mathrm{obs}}\\ resolves this.

## Fixed-effect parameter vector

The \\D = 2P + 1\\ dimensional parameter vector is ordered as \$\$\theta
= (\alpha_1, \ldots, \alpha_P,\\ \beta_1, \ldots, \beta_P,\\
\log\kappa)\$\$ where \\\alpha\\ governs the extensive margin
(Bernoulli), \\\beta\\ governs the intensive margin (zero-truncated
Beta-Binomial), and \\\kappa\\ is the dispersion parameter.

## Design Effect Ratio

The DER for each parameter is \\\mathrm{DER}\_p =
V\_{\mathrm{sand}}\[p,p\] / H\_{\mathrm{obs}}^{-1}\[p,p\]\\. Expected
values are 1–5 for typical survey designs.

## References

Williams, M. R. and Savitsky, T. D. (2021). Uncertainty estimation for
pseudo-Bayesian inference under complex sampling. *International
Statistical Review*, **89**(1), 72–107.

## See also

[`compute_score_matrix()`](https://joonho112.github.io/hurdlebb/reference/compute_score_matrix.md),
[`compute_H_obs()`](https://joonho112.github.io/hurdlebb/reference/compute_H_obs.md),
[`compute_J_cluster()`](https://joonho112.github.io/hurdlebb/reference/compute_J_cluster.md),
[`compute_der()`](https://joonho112.github.io/hurdlebb/reference/compute_der.md)

Other sandwich:
[`compute_H_obs()`](https://joonho112.github.io/hurdlebb/reference/compute_H_obs.md),
[`compute_J_cluster()`](https://joonho112.github.io/hurdlebb/reference/compute_J_cluster.md),
[`compute_der()`](https://joonho112.github.io/hurdlebb/reference/compute_der.md),
[`compute_score_matrix()`](https://joonho112.github.io/hurdlebb/reference/compute_score_matrix.md),
[`print.hbb_sandwich()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_sandwich.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- hbb(
  y | trials(n_trial) ~ poverty + urban,
  data = my_data, weights = "weight",
  stratum = "vstratum", psu = "vpsu"
)
sand <- sandwich_variance(fit)
print(sand)
compute_der(sand)
} # }
```
