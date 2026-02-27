# Compute Design Effect Ratios from a Sandwich Variance Object

Returns the Design Effect Ratio (DER) for each fixed-effect parameter:
\$\$\mathrm{DER}\_p = V\_{\mathrm{sand}}\[p,p\] /
H\_{\mathrm{obs}}^{-1}\[p,p\]\$\$

## Usage

``` r
compute_der(sandwich)
```

## Arguments

- sandwich:

  An object of class `"hbb_sandwich"` from
  [`sandwich_variance()`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md).

## Value

Named numeric vector of length \\D\\.

## Details

DER \> 1 means the survey design inflates variance beyond what a naive
(non-survey-corrected) analysis would produce. Expected values are
typically 1–5 for stratified cluster designs (Kish DEFF ~ 3–4).

The following warnings are issued:

- DER \< 0.5: survey design appears to reduce variance substantially,
  which is unusual and may signal a problem.

- DER \> 20: very large design effect, possibly indicating model
  misspecification or extreme weight variability.

## See also

[`sandwich_variance()`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md)

Other sandwich:
[`compute_H_obs()`](https://joonho112.github.io/hurdlebb/reference/compute_H_obs.md),
[`compute_J_cluster()`](https://joonho112.github.io/hurdlebb/reference/compute_J_cluster.md),
[`compute_score_matrix()`](https://joonho112.github.io/hurdlebb/reference/compute_score_matrix.md),
[`print.hbb_sandwich()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_sandwich.md),
[`sandwich_variance()`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md)
