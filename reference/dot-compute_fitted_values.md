# Compute fitted values from posterior means

Uses BLAS-optimised matrix multiplication to compute linear predictors
and inverse-logit transformations at the posterior mean of the
fixed-effect parameter vector. For SVC models, adds state-level random
effects via `.extract_delta_means()`.

## Usage

``` r
.compute_fitted_values(object)
```

## Arguments

- object:

  An `hbb_fit` object (already validated).

## Value

Named list with elements:

- `q_hat`:

  N-vector of participation probabilities.

- `mu_hat`:

  N-vector of conditional intensities.

- `kappa_hat`:

  Scalar dispersion parameter (\\\hat\kappa =
  \exp(\widehat{\log\kappa})\\, clamped to \\\[10^{-15}, 10^{15}\]\\).

- `theta_hat`:

  D-vector of posterior means.

## Edge cases handled

- Extreme `log_kappa`: `kappa_hat` is clamped via
  `pmin(exp(log_kappa), 1e15)` to prevent overflow.

- Non-finite posterior means: a warning is issued but computation
  proceeds (`plogis` handles +/-Inf gracefully).

- SVC delta extraction failure: a warning is issued and fitted values
  use fixed effects only.

- All `z = 0` (all structural zeros): `q_hat` will be small but not
  exactly zero (logistic never reaches 0); `mu_hat` is still computed
  but irrelevant.

- Draw matrix column mismatch: aborts with informative error.

- Insufficient MCMC draws (M \< 2): aborts with informative error.
