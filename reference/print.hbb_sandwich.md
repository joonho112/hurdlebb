# Print Method for hbb_sandwich Objects

Displays a compact summary of the sandwich variance estimator, including
the design effect ratios (DER), survey design information, and matrix
condition diagnostics. Wrapped in `tryCatch` so that printing never
fails even on malformed objects.

## Usage

``` r
# S3 method for class 'hbb_sandwich'
print(x, ...)
```

## Arguments

- x:

  An object of class `"hbb_sandwich"`.

- ...:

  Currently unused.

## Value

Invisibly returns `x`.

## See also

[`sandwich_variance()`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md),
[`compute_der()`](https://joonho112.github.io/hurdlebb/reference/compute_der.md)

Other sandwich:
[`compute_H_obs()`](https://joonho112.github.io/hurdlebb/reference/compute_H_obs.md),
[`compute_J_cluster()`](https://joonho112.github.io/hurdlebb/reference/compute_J_cluster.md),
[`compute_der()`](https://joonho112.github.io/hurdlebb/reference/compute_der.md),
[`compute_score_matrix()`](https://joonho112.github.io/hurdlebb/reference/compute_score_matrix.md),
[`sandwich_variance()`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md)
