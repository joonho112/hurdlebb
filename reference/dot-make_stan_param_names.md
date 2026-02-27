# Generate CmdStanR parameter names

Produces
`c("alpha[1]", ..., "alpha[P]", "beta[1]", ..., "beta[P]", "log_kappa")`
matching CmdStanR naming conventions.

## Usage

``` r
.make_stan_param_names(P)
```

## Arguments

- P:

  Integer; number of covariates per margin.

## Value

Character vector of length 2P + 1.
