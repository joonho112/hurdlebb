# Safely compute lower-triangular Cholesky factor

Returns `t(chol(mat))`, i.e., the lower-triangular factor L such that
`mat = L %*% t(L)`. Wraps in tryCatch for informative error messages.

## Usage

``` r
.safe_chol_lower(mat, mat_name = "matrix")
```

## Arguments

- mat:

  A positive-definite matrix.

- mat_name:

  Name for error messages.

## Value

Lower-triangular Cholesky factor.
