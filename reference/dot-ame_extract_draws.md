# Extract parameter draws for AME computation

If a cholesky object is provided, uses theta_corrected and theta_hat.
Otherwise, extracts draws from the CmdStanR fit.

## Usage

``` r
.ame_extract_draws(fit, cholesky, P)
```

## Arguments

- fit:

  An hbb_fit object.

- cholesky:

  Optional hbb_cholesky object.

- P:

  Integer: number of covariates per margin.

## Value

Named list: draws (M x D matrix), theta_hat (D-vector), M_total
(integer).
