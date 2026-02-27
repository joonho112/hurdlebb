# Verify the Cholesky transform satisfies algebraic identities

Checks three properties of the affine correction: mean preservation,
finite-sample variance recovery, and the algebraic identity A Sigma A' =
V_sand.

## Usage

``` r
.verify_cholesky_transform(theta_corrected, theta_hat, A, Sigma_MCMC, V_sand)
```

## Arguments

- theta_corrected:

  Corrected M by D draw matrix.

- theta_hat:

  D-vector of posterior means.

- A:

  D by D transformation matrix.

- Sigma_MCMC:

  D by D MCMC covariance.

- V_sand:

  D by D sandwich variance.

## Value

Named list with numeric diagnostics and logical pass flags.
