# Build a comparison table of naive, corrected, and Wald CIs

Build a comparison table of naive, corrected, and Wald CIs

## Usage

``` r
.build_comparison_table(
  theta_hat,
  theta_draws,
  theta_corrected,
  V_sand,
  Sigma_MCMC,
  H_obs_inv,
  param_labels,
  level
)
```

## Arguments

- theta_hat:

  Numeric vector length D: posterior means.

- theta_draws:

  M x D matrix: original MCMC draws.

- theta_corrected:

  M x D matrix: corrected draws.

- V_sand:

  D x D sandwich variance.

- Sigma_MCMC:

  D x D MCMC posterior covariance.

- H_obs_inv:

  D x D inverse observed Fisher information.

- param_labels:

  Character vector length D.

- level:

  Confidence level.

## Value

A data.frame with one row per parameter and 16 columns.
