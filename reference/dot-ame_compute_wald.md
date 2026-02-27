# Compute Wald delta-method CIs for AME

Uses central-difference numerical gradient of the AME with respect to
the D-vector theta, combined with the sandwich variance V_sand, to
obtain standard errors via the delta method: SE(AME_k) = sqrt(grad_k'
V_sand grad_k).

## Usage

``` r
.ame_compute_wald(theta_hat, X, V_sand, P, D, level)
```

## Arguments

- theta_hat:

  Numeric D-vector: posterior mean.

- X:

  Numeric N x P matrix.

- V_sand:

  D x D sandwich variance matrix.

- P:

  Integer: number of covariates.

- D:

  Integer: total parameters (2P + 1).

- level:

  Confidence level.

## Value

Data frame with P rows and 13 columns: covariate, ext_point, ext_se,
ext_lo, ext_hi, int_point, int_se, int_lo, int_hi, total_point,
total_se, total_lo, total_hi.
