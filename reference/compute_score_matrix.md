# Extract Posterior Mean Score Matrix from an hbb_fit

Extracts the posterior means of the Stan-generated score variables
(`score_ext`, `score_int`, `score_kappa`) and assembles them into an \\N
\times D\\ matrix, where \\D = 2P + 1\\: \\(\alpha_1, \ldots, \alpha_P,
\beta_1, \ldots, \beta_P, \log\kappa)\\.

## Usage

``` r
compute_score_matrix(fit)
```

## Arguments

- fit:

  An object of class `"hbb_fit"`. Must be a weighted model variant
  (`"weighted"` or `"svc_weighted"`).

## Value

Numeric matrix of dimension \\N \times D\\. Columns are ordered as
(score_ext columns 1 to P, score_int columns 1 to P, score_kappa).

## Details

Score variables are generated quantities in the weighted/svc_weighted
Stan models. CmdStanR names matrix variables as
`var[1,1], var[1,2], ..., var[1,P], var[2,1], ...` (first index varies
slowest, row-major ordering). Posterior means are computed across all
MCMC draws.

## See also

[`sandwich_variance()`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md)

Other sandwich:
[`compute_H_obs()`](https://joonho112.github.io/hurdlebb/reference/compute_H_obs.md),
[`compute_J_cluster()`](https://joonho112.github.io/hurdlebb/reference/compute_J_cluster.md),
[`compute_der()`](https://joonho112.github.io/hurdlebb/reference/compute_der.md),
[`print.hbb_sandwich()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_sandwich.md),
[`sandwich_variance()`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md)
