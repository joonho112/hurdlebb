# Predictions from a Hurdle Beta-Binomial Model

Point predictions are computed at the posterior mean of the fixed
effects (and state random effects for SVC models): \$\$ \hat{y}\_i =
\hat{q}(x_i) \cdot \hat\mu(x_i), \$\$ where \\\hat{q}(x) =
\mathrm{logistic}(x'\hat\alpha)\\ and \\\hat\mu(x) =
\mathrm{logistic}(x'\hat\beta)\\.

When `interval = "credible"`, the full posterior distribution is
propagated through the link function to produce empirical credible
intervals. For each MCMC draw \\m\\: \$\$ \hat{y}\_i^{(m)} =
q^{(m)}(x_i) \cdot \mu^{(m)}(x_i), \$\$ and the \\(\alpha/2,\\
1-\alpha/2)\\ quantiles across draws give the interval bounds.

## Usage

``` r
# S3 method for class 'hbb_fit'
predict(
  object,
  newdata = NULL,
  type = c("response", "extensive", "intensive"),
  interval = c("none", "credible"),
  level = 0.95,
  ndraws = NULL,
  state = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `"hbb_fit"` returned by
  [`hbb`](https://joonho112.github.io/hurdlebb/reference/hbb.md).

- newdata:

  An optional data frame containing the covariates for prediction. If
  `NULL` (default), in-sample predictions are returned using the
  training design matrix.

- type:

  Character string: `"response"` (default), `"extensive"`, or
  `"intensive"`. See **Prediction types** above.

- interval:

  Character string: `"none"` (default) for point predictions only, or
  `"credible"` for posterior credible intervals.

- level:

  Numeric scalar in \\(0, 1)\\. Credible level when
  `interval = "credible"`. Default is `0.95`.

- ndraws:

  Optional positive integer. If specified, the posterior draws are
  thinned to `ndraws` before computing intervals. Ignored when
  `interval = "none"`. Default is `NULL` (use all draws).

- state:

  Character string or character vector specifying states for SVC
  prediction.

  Character of length 1:

  :   If `newdata` is provided, interpreted as a column name in
      `newdata`. Otherwise, interpreted as a single state label to apply
      to all observations.

  Character vector of length `nrow(newdata)`:

  :   Direct state labels for each observation.

  `NULL`:

  :   No state information. For SVC models, predictions use fixed
      effects only (with a warning).

- ...:

  Additional arguments (currently unused).

## Value

A data frame with one row per observation. Columns:

- `fit`:

  Numeric: point prediction (posterior mean).

- `lwr`:

  Numeric: lower credible bound. Present only when
  `interval = "credible"`.

- `upr`:

  Numeric: upper credible bound. Present only when
  `interval = "credible"`.

## Details

Generates point predictions and optional posterior credible intervals
from a fitted hurdle Beta-Binomial model. Supports both in-sample
prediction (using the training data) and out-of-sample prediction (using
`newdata`).

## Prediction types

Three types of predictions are available:

- `"response"`:

  (Default.) The composite enrollment proportion \\\hat{q}\_i \cdot
  \hat\mu_i \in \[0,1\]\\. This represents the unconditional expected
  proportion \\E\[Y_i/n_i\]\\.

- `"extensive"`:

  The participation probability \\\hat{q}\_i =
  \mathrm{logistic}(x_i'\hat\alpha) \in (0,1)\\.

- `"intensive"`:

  The conditional enrollment share \\\hat\mu_i =
  \mathrm{logistic}(x_i'\hat\beta) \in (0,1)\\.

## Newdata design matrix

When `newdata` is provided, the design matrix is constructed by:

1.  Extracting the covariate columns named in `object$formula$fixed`
    from `newdata`.

2.  Applying the same centering and scaling as the training data, using
    the stored `hbb_data$x_center` (column means) and `hbb_data$x_scale`
    (column SDs). This ensures that regression coefficients remain on
    the standardised scale.

3.  Prepending an intercept column of ones.

All covariate columns must be present in `newdata`. Missing columns
cause an informative error.

## State-varying coefficients (SVC)

For SVC models, out-of-sample prediction requires specifying `state`
(either a column name in `newdata` or a vector of state labels). The
state labels must match `hbb_data$group_levels` from the training data.
Unknown state labels cause an error.

If `state = NULL` for an SVC model, predictions use fixed effects only
(i.e., the population-average prediction without state random effects),
and a warning is issued.

## Credible interval theory

The posterior credible interval for prediction \\i\\ at level
\\1-\alpha\\ is \$\$ \bigl\[Q\_{\alpha/2}(\hat{y}\_i^{(1:M)}),\\
Q\_{1-\alpha/2}(\hat{y}\_i^{(1:M)})\bigr\], \$\$ where \\Q_p\\ denotes
the empirical \\p\\-quantile over \\M\\ posterior draws. This is a
*conditional* credible interval for the expected proportion \\E\[Y_i/n_i
\| \theta\]\\, not a predictive interval for \\Y_i\\ itself (which would
additionally account for sampling variability from the Beta-Binomial).

For computational efficiency, `ndraws` can be used to thin the draws
before propagation. With \\M = 4000\\ total draws, setting
`ndraws = 500` reduces computation 8-fold while typically preserving
interval accuracy.

## References

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., and
Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). Chapman and
Hall/CRC.

Ghosal, R., Ghosh, S. K., and Maiti, T. (2020). Two-part regression
models for longitudinal zero-inflated count data. *Journal of the Royal
Statistical Society: Series A*, **183**(4), 1603–1626.

## See also

[`fitted.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/fitted.hbb_fit.md)
for in-sample fitted values (simpler interface without newdata support),
[`residuals.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/residuals.hbb_fit.md)
for residuals,
[`summary.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/summary.hbb_fit.md)
for model summary.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)

# In-sample point predictions
pred <- predict(fit)

# Out-of-sample with credible intervals
new_df <- data.frame(poverty = c(0.1, 0.3, 0.5),
                     urban   = c(1, 0, 1))
pred_ci <- predict(fit, newdata = new_df,
                   interval = "credible", level = 0.95)

# Extensive-margin predictions only
pred_ext <- predict(fit, type = "extensive")

# SVC model with state assignment
pred_svc <- predict(fit_svc, newdata = new_df,
                    state = c("AL", "CA", "NY"),
                    interval = "credible")
} # }
```
