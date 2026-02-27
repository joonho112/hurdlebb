# Fitted Values from an hbb_fit

Computes observation-level fitted values from the posterior means of the
fixed effects (and state random effects for SVC models).

## Usage

``` r
# S3 method for class 'hbb_fit'
fitted(object, type = c("response", "extensive", "intensive"), ...)
```

## Arguments

- object:

  An object of class `"hbb_fit"` returned by
  [`hbb`](https://joonho112.github.io/hurdlebb/reference/hbb.md).

- type:

  Character string: `"response"` (default), `"extensive"`, or
  `"intensive"`.

- ...:

  Currently unused; included for S3 method consistency.

## Value

Numeric vector of length \\N\\.

## Details

Three types are available:

- `"response"`:

  (Default.) The predicted enrollment proportion \\\hat{q}\_i \cdot
  \hat{\mu}\_i\\, i.e. the probability of a positive enrollment times
  the conditional enrollment share. All values lie in \\\[0, 1\]\\.

- `"extensive"`:

  The predicted participation probability \\\hat{q}\_i =
  \mathrm{logistic}(X_i' \hat\alpha)\\.

- `"intensive"`:

  The predicted conditional enrollment share \\\hat\mu_i =
  \mathrm{logistic}(X_i' \hat\beta)\\.

## SVC adjustment

For state-varying coefficient models
(`model_type %in% c("svc", "svc_weighted")`), the linear predictors
include state-level random effects: \$\$ \eta\_{\mathrm{ext},i} = X_i'
\hat\alpha + X_i' \hat\delta\_{\mathrm{ext}}\[s(i)\], \$\$ where
\\\hat\delta\\ is the posterior mean of the state random effects
extracted via `.extract_delta_means()`. If delta extraction fails, a
warning is issued and fitted values use fixed effects only.

## See also

[`residuals.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/residuals.hbb_fit.md),
[`coef.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/coef.hbb_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
yhat   <- fitted(fit)                      # q * mu (proportion)
q_hat  <- fitted(fit, type = "extensive")  # participation prob
mu_hat <- fitted(fit, type = "intensive")  # conditional share
} # }
```
