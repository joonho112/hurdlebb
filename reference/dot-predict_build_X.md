# Build the design matrix for prediction

If `newdata` is `NULL`, returns the training design matrix
`object$hbb_data$X`. Otherwise, constructs a new design matrix by
extracting the covariates named in `object$formula$fixed`, applying the
stored centering and scaling transformations, and prepending an
intercept column.

## Usage

``` r
.predict_build_X(object, newdata)
```

## Arguments

- object:

  An `hbb_fit` object.

- newdata:

  A data frame or `NULL`.

## Value

Numeric matrix of dimension \\N\_{\mathrm{new}} \times P\\ (with
intercept column).

## Centering and scaling

The stored transformations `hbb_data$x_center` (named numeric vector of
training column means) and `hbb_data$x_scale` (named numeric vector of
training column SDs) are applied to newdata covariates to ensure the
design matrix is on the same scale as the training data. The
centering/scaling formula for covariate \\j\\ is: \$\$ x\_{ij}^\* =
\frac{x\_{ij} - \bar{x}\_j^{\mathrm{train}}} {s_j^{\mathrm{train}}},
\$\$ where \\\bar{x}\_j^{\mathrm{train}}\\ and \\s_j^{\mathrm{train}}\\
are the mean and SD from the training data.

## Transformation steps

1.  Extract columns named in `object$formula$fixed`.

2.  If training centred (`x_center` is not `NULL`), subtract the
    *training* column means via
    [`sweep()`](https://rdrr.io/r/base/sweep.html).

3.  If training scaled (`x_scale` is not `NULL`), divide by the
    *training* column SDs via
    [`sweep()`](https://rdrr.io/r/base/sweep.html).

4.  Prepend an intercept column of ones.

5.  Validate that `ncol(X_new) == P`.
