# Extract Fixed-Effect Coefficients from an hbb_fit

Returns the posterior mean of the fixed-effect parameter vector
\\\hat\theta = (\hat\alpha_1, \ldots, \hat\alpha_P,\\ \hat\beta_1,
\ldots, \hat\beta_P,\\ \widehat{\log\kappa})\\, optionally restricted to
a single margin.

## Usage

``` r
# S3 method for class 'hbb_fit'
coef(object, margin = c("both", "extensive", "intensive"), ...)
```

## Arguments

- object:

  An object of class `"hbb_fit"` returned by
  [`hbb`](https://joonho112.github.io/hurdlebb/reference/hbb.md).

- margin:

  Character string specifying which parameters to return:

  `"both"`

  :   (Default.) The full \\D = 2P + 1\\ parameter vector.

  `"extensive"`

  :   Only \\\hat\alpha\_{1:P}\\ (P-vector).

  `"intensive"`

  :   Only \\\hat\beta\_{1:P}\\ (P-vector).

- ...:

  Currently unused; included for S3 method consistency.

## Value

A named numeric vector of posterior means. Length depends on `margin`:
\\2P + 1\\ for `"both"`, \\P\\ for `"extensive"` or `"intensive"`.

## Mathematical note

The posterior mean is the Bayes estimator under squared-error loss. For
survey-weighted pseudo-posteriors, the posterior mean converges to the
pseudo-true parameter under standard regularity conditions (Williams and
Savitsky, 2021).

## See also

[`vcov.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/vcov.hbb_fit.md),
[`summary.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/summary.hbb_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
coef(fit)                       # full D-vector
coef(fit, margin = "extensive") # alpha only
coef(fit, margin = "intensive") # beta only
} # }
```
