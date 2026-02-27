# Compute posterior summary for an AME draw matrix

Compute posterior summary for an AME draw matrix

## Usage

``` r
.ame_summarize(ame_mat, labels, level)
```

## Arguments

- ame_mat:

  M x P numeric matrix of AME draws.

- labels:

  Character vector of length P: covariate labels.

- level:

  Numeric scalar in (0,1): credible interval level.

## Value

Data frame with columns: covariate, post_mean, post_median, ci_lo,
ci_hi, post_sd, pr_positive.
