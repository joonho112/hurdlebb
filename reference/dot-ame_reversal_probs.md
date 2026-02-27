# Compute reversal probabilities for each non-intercept covariate

For each non-intercept covariate, computes the posterior probability
that the extensive and intensive AME components have opposing signs.

## Usage

``` r
.ame_reversal_probs(ext_draws, int_draws, labels, idx)
```

## Arguments

- ext_draws:

  M by P matrix of extensive AME draws.

- int_draws:

  M by P matrix of intensive AME draws.

- labels:

  Character P-vector of covariate labels.

- idx:

  Integer vector of non-intercept indices.

## Value

Named numeric vector of reversal probabilities.
