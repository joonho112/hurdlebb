# Compute the Block-Diagonal Observed Information Matrix

Computes the observed information matrix \\H\_{\mathrm{obs}}\\ for the
fixed effects of a hurdle Beta-Binomial model. Uses a block-diagonal
structure:

## Usage

``` r
compute_H_obs(fit, scores = NULL)
```

## Arguments

- fit:

  An object of class `"hbb_fit"`.

- scores:

  Optional \\N \times D\\ score matrix (as from
  [`compute_score_matrix()`](https://joonho112.github.io/hurdlebb/reference/compute_score_matrix.md)).
  If `NULL`, computed internally.

## Value

A list with components:

- `H_obs`:

  Numeric matrix (\\D \times D\\). Block-diagonal observed information.

- `H_obs_inv`:

  Numeric matrix (\\D \times D\\). Inverse of H_obs (the bread matrix).

- `H_ext`:

  Numeric matrix (\\P \times P\\). Extensive-margin Fisher information
  block.

- `H_int`:

  Numeric matrix (\\(P+1) \times (P+1)\\). Intensive + kappa block from
  information identity.

- `ridge_applied`:

  Logical. Whether ridge regularisation was needed.

## Details

- **Extensive margin** (alpha block, \\P \times P\\): analytic Fisher
  information from logistic regression.

- **Intensive + kappa block** (\\(P+1) \times (P+1)\\): empirical
  information identity using posterior mean scores.

## Extensive margin (H_ext)

For the base weighted model: \$\$H\_{\mathrm{ext}} = \sum\_{i=1}^N
\tilde{w}\_i\\ q_i (1-q_i)\\ X_i X_i^\top\$\$ where \\q_i =
\mathrm{logit}^{-1}(X_i^\top \hat{\alpha})\\.

For SVC models, the linear predictor includes state random effects:
\\q_i = \mathrm{logit}^{-1}(X_i^\top \hat{\alpha} + X_i^\top
\hat{\delta}^{\mathrm{ext}}\_{s\[i\]})\\.

## Intensive margin (H_int)

The zero-truncated Beta-Binomial Hessian is analytically complex, so we
use the information identity: \$\$H\_{\mathrm{int}} = \sum\_{i=1}^N
\tilde{w}\_i\\ s_i^{\mathrm{int}} (s_i^{\mathrm{int}})^\top\$\$ where
\\s_i^{\mathrm{int}} = (s\_{\beta,i}, s\_{\kappa,i})^\top\\ is the
\\(P+1)\\-dimensional intensive score at the posterior mean.

## See also

[`sandwich_variance()`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md),
[`compute_score_matrix()`](https://joonho112.github.io/hurdlebb/reference/compute_score_matrix.md)

Other sandwich:
[`compute_J_cluster()`](https://joonho112.github.io/hurdlebb/reference/compute_J_cluster.md),
[`compute_der()`](https://joonho112.github.io/hurdlebb/reference/compute_der.md),
[`compute_score_matrix()`](https://joonho112.github.io/hurdlebb/reference/compute_score_matrix.md),
[`print.hbb_sandwich()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_sandwich.md),
[`sandwich_variance()`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md)
