# Compute the Cluster-Robust Meat Matrix

Aggregates weighted scores to the stratum–PSU level and computes the
Taylor linearisation variance estimator with finite population
correction (FPC) at the stratum level.

## Usage

``` r
compute_J_cluster(scores, hbb_data)
```

## Arguments

- scores:

  Numeric matrix of dimension \\N \times D\\. The posterior mean score
  matrix from
  [`compute_score_matrix()`](https://joonho112.github.io/hurdlebb/reference/compute_score_matrix.md).

- hbb_data:

  An S3 object of class `"hbb_data"` (from
  [`prepare_stan_data()`](https://joonho112.github.io/hurdlebb/reference/prepare_stan_data.md))
  or a list containing the fields `stratum_idx`, `psu_idx`, `w_tilde`,
  and `N`.

## Value

Numeric matrix of dimension \\D \times D\\. The cluster-robust meat
matrix (positive semi-definite).

## Details

### Algorithm

For each stratum \\h\\ with \\C_h\\ PSUs:

1.  Compute weighted score totals per PSU: \\s\_{hc} = \sum\_{i \in
    \mathrm{PSU}(h,c)} \tilde{w}\_i s_i\\

2.  Compute stratum mean: \\\bar{s}\_h = C_h^{-1} \sum_c s\_{hc}\\

3.  Center: \\\delta\_{hc} = s\_{hc} - \bar{s}\_h\\

4.  Accumulate with FPC: \\J_h = \frac{C_h}{C_h - 1} \sum_c \delta\_{hc}
    \delta\_{hc}^\top\\

Uses [`rowsum()`](https://rdrr.io/r/base/rowsum.html) for efficient
PSU-level aggregation. Singleton strata (\\C_h = 1\\) are skipped with a
warning.

## See also

[`sandwich_variance()`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md)

Other sandwich:
[`compute_H_obs()`](https://joonho112.github.io/hurdlebb/reference/compute_H_obs.md),
[`compute_der()`](https://joonho112.github.io/hurdlebb/reference/compute_der.md),
[`compute_score_matrix()`](https://joonho112.github.io/hurdlebb/reference/compute_score_matrix.md),
[`print.hbb_sandwich()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_sandwich.md),
[`sandwich_variance()`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md)
