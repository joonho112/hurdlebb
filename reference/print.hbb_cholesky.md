# Print Method for Cholesky-Corrected Posterior Objects

Displays a structured summary of the Cholesky posterior recalibration,
including the transformation matrix diagonal (with shrinkage/inflation
labels), verification checks, DER summary, and an abbreviated comparison
table.

## Usage

``` r
# S3 method for class 'hbb_cholesky'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `"hbb_cholesky"` returned by
  [`cholesky_correct`](https://joonho112.github.io/hurdlebb/reference/cholesky_correct.md).

- digits:

  Integer; number of significant digits for numeric output. Default is
  4.

- ...:

  Additional arguments passed to
  [`format`](https://rdrr.io/r/base/format.html).

## Value

The object `x`, invisibly.

## See also

[`cholesky_correct`](https://joonho112.github.io/hurdlebb/reference/cholesky_correct.md)
