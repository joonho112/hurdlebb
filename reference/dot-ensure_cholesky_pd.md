# Ensure a matrix is positive definite for Cholesky decomposition

Checks positive definiteness and applies corrections if needed. For
Sigma_MCMC, uses ridge regularisation (additive diagonal perturbation).
For V_sand, uses nearPD projection (Higham 2002).

## Usage

``` r
.ensure_cholesky_pd(mat, mat_name = "matrix", method = c("ridge", "nearpd"))
```

## Arguments

- mat:

  Symmetric matrix.

- mat_name:

  Character label for diagnostic messages.

- method:

  Character: `"ridge"` adds a diagonal ridge, `"nearpd"` projects via
  [`Matrix::nearPD`](https://rdrr.io/pkg/Matrix/man/nearPD.html).

## Value

A list with elements: `mat` (corrected matrix), `corrected` (logical),
`details` (character or NULL), `min_eig`, `max_eig`, `cond_number`.
