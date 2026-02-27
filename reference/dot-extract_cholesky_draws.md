# Extract posterior draws from an hbb_fit for Cholesky correction

Validates that `fit$fit` exists and draws match expected dimensions.

## Usage

``` r
.extract_cholesky_draws(fit, param_names)
```

## Arguments

- fit:

  An hbb_fit object.

- param_names:

  Character vector of CmdStanR parameter names.

## Value

Numeric matrix of draws (M x D).
