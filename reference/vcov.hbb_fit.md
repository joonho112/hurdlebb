# Variance-Covariance Matrix of Fixed Effects

Returns the \\D \times D\\ variance-covariance matrix of the
fixed-effect parameter vector.

## Usage

``` r
# S3 method for class 'hbb_fit'
vcov(object, sandwich = NULL, ...)
```

## Arguments

- object:

  An object of class `"hbb_fit"` returned by
  [`hbb`](https://joonho112.github.io/hurdlebb/reference/hbb.md).

- sandwich:

  An optional object of class `"hbb_sandwich"` returned by
  [`sandwich_variance`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md),
  or `NULL` (default).

  If provided:

  :   The design-consistent sandwich variance \\V\_{\mathrm{sand}} =
      H\_{\mathrm{obs}}^{-1}\\ J\_{\mathrm{cluster}}\\
      H\_{\mathrm{obs}}^{-1}\\ is returned. This is the recommended
      choice for survey-weighted models.

  If `NULL`:

  :   The MCMC posterior covariance \\\Sigma\_{\mathrm{MCMC}} =
      \mathrm{Cov}(\theta^{(m)})\\ is returned. For hierarchical models
      this is prior-dominated and should be used only for unweighted
      base models.

- ...:

  Currently unused; included for S3 method consistency.

## Value

A named numeric matrix of dimension \\D \times D\\ where \\D = 2P + 1\\.
Row and column names are human-readable parameter labels (e.g.,
`alpha_intercept`, `beta_poverty`, `log_kappa`).

## Warning – prior domination

For hierarchical models with state-varying coefficients,
\\\Sigma\_{\mathrm{MCMC}}\\ absorbs prior and random-effect variance,
producing Prior Inflation ratios of 600–6500 for fixed effects (see the
sandwich variance documentation). Using \\\Sigma\_{\mathrm{MCMC}}\\ as a
variance estimate will substantially over-cover. The sandwich variance
is the appropriate measure of frequentist uncertainty.

## See also

[`coef.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/coef.hbb_fit.md),
[`sandwich_variance`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
vcov(fit)  # MCMC posterior covariance

# With sandwich variance (recommended for survey data)
sand <- sandwich_variance(fit)
vcov(fit, sandwich = sand)
} # }
```
