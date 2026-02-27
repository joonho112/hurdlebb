# Compute AME at a given theta vector

Compute AME at a given theta vector

## Usage

``` r
.ame_at_theta(theta, X, P)
```

## Arguments

- theta:

  Numeric vector of length D = 2P + 1.

- X:

  Numeric matrix N x P.

- P:

  Integer: number of covariates per margin.

## Value

Named list with elements: ext (P-vector), int (P-vector), total
(P-vector), mean_q (scalar), mean_mu (scalar).
