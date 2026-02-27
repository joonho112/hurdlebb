# Residuals from an hbb_fit

Computes residuals for the hurdle Beta-Binomial model. Two types are
supported:

## Usage

``` r
# S3 method for class 'hbb_fit'
residuals(object, type = c("response", "pearson"), ...)
```

## Arguments

- object:

  An object of class `"hbb_fit"` returned by
  [`hbb`](https://joonho112.github.io/hurdlebb/reference/hbb.md).

- type:

  Character string: `"response"` (default) or `"pearson"`.

- ...:

  Currently unused; included for S3 method consistency.

## Value

Numeric vector of length \\N\\.

## Details

- `"response"`:

  Raw residuals on the proportion scale: \\r_i = y_i / n_i - \hat{q}\_i
  \hat{\mu}\_i\\.

- `"pearson"`:

  Pearson residuals standardised by the model-implied standard deviation
  (count scale): \$\$ r_i^P = \frac{y_i - n_i \hat{q}\_i \hat{\mu}\_i}
  {\sqrt{\widehat{\mathrm{Var}}(Y_i)}}, \$\$ where
  \\\widehat{\mathrm{Var}}(Y_i)\\ is computed via
  [`hurdle_variance`](https://joonho112.github.io/hurdlebb/reference/hurdle_variance.md)
  using \\(n_i, \hat{q}\_i, \hat\mu_i, \hat\kappa)\\. Under correct
  model specification, Pearson residuals have approximately mean zero
  and variance one.

## Numerical stability (Pearson residuals)

Several guards are applied to prevent division by zero or NaN:

- \\\hat\mu_i\\ is clamped to \\(\varepsilon,\\ 1 - \varepsilon)\\ where
  \\\varepsilon = \sqrt{\mathtt{.Machine\\double.eps}}\\ (\\\approx 1.49
  \times 10^{-8}\\).

- \\\hat{q}\_i\\ is clamped to \\\[\varepsilon,\\ 1\]\\. (\\q = 0\\
  yields \\\mathrm{Var}\[Y\] = 0\\, making Pearson undefined.)

- \\\hat{\kappa}\\ is clamped to \\\le 10^{15}\\.

- The hurdle variance is floored at `.Machine$double.eps` to prevent
  division by zero.

- Non-finite residuals are replaced with `NA` and a warning is issued.

## See also

[`fitted.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/fitted.hbb_fit.md),
[`hurdle_variance`](https://joonho112.github.io/hurdlebb/reference/hurdle_variance.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
r_raw     <- residuals(fit)
r_pearson <- residuals(fit, type = "pearson")
hist(r_pearson, breaks = 50)
} # }
```
