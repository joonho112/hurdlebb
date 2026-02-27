# Plot Diagnostics and Inferential Summaries for hbb_fit Objects

This method dispatches to four specialised internal helpers:

- `"coefficients"`:

  A forest plot of fixed-effect point estimates with confidence/credible
  intervals, faceted by model margin (extensive, intensive, dispersion).
  Uses
  [`summary.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/summary.hbb_fit.md)
  internally, so the same inference mode (Wald or posterior) applies.

- `"trace"`:

  MCMC trace plots showing the evolution of sampled parameter values
  across iterations, coloured by chain. Useful for diagnosing
  convergence and mixing.

- `"random_effects"`:

  Caterpillar plots of the posterior mean and credible interval for each
  state-level random effect \\\delta\_{s,k}\\, faceted by covariate.
  Only available for SVC models
  (`model_type %in% c("svc", "svc_weighted")`).

- `"residuals"`:

  A two-panel diagnostic combining residuals-versus-fitted and a normal
  QQ plot of Pearson residuals. Panels are arranged via patchwork.

## Usage

``` r
# S3 method for class 'hbb_fit'
plot(
  x,
  type = c("coefficients", "trace", "random_effects", "residuals"),
  sandwich = NULL,
  level = 0.95,
  margin = c("both", "extensive", "intensive"),
  pars = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `"hbb_fit"` returned by
  [`hbb`](https://joonho112.github.io/hurdlebb/reference/hbb.md).

- type:

  Character string: one of `"coefficients"` (default), `"trace"`,
  `"random_effects"`, or `"residuals"`.

- sandwich:

  An optional object of class `"hbb_sandwich"` returned by
  [`sandwich_variance`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md).
  Used only when `type = "coefficients"` to produce Wald-based
  intervals. If `NULL` (default), posterior-based intervals are shown.

- level:

  Numeric scalar in \\(0, 1)\\. Confidence/credible level for interval
  construction. Default is `0.95`.

- margin:

  Character string controlling which coefficients appear in the forest
  plot. One of `"both"` (default), `"extensive"`, or `"intensive"`.
  Ignored for other plot types.

- pars:

  Optional character vector of parameter names to include in the trace
  plot. If `NULL` (default), all \\D = 2P+1\\ fixed-effect parameters
  are plotted. Ignored for other plot types.

- ...:

  Additional arguments passed to internal helpers (currently unused).

## Value

A `ggplot` or `patchwork` object, returned invisibly. The plot is
printed as a side effect.

## Details

Produces diagnostic and summary plots for fitted Hurdle Beta-Binomial
models. Four plot types are available, selected via the `type` argument.

## ggplot2 requirement

This function requires ggplot2 (and patchwork for `type = "residuals"`).
Both are in `Suggests`; they are checked via
[`rlang::check_installed()`](https://rlang.r-lib.org/reference/is_installed.html)
and the user will be prompted to install them if absent.

## Inference mode (coefficients plot)

When `sandwich` is supplied, the forest plot shows Wald confidence
intervals derived from the sandwich variance: \$\$
\mathrm{CI}\_{1-\alpha}(\theta_p) = \hat\theta_p \pm
z\_{(1+\mathrm{level})/2}\\ \sqrt{V\_{\mathrm{sand},pp}}. \$\$ When
`sandwich = NULL`, quantile-based credible intervals from the MCMC
posterior are shown instead.

## Colour palette

The project palette is used throughout:

- `#4393C3` (blue) — extensive margin

- `#D6604D` (red) — intensive margin

- `#1B7837` (green) — reference lines, intervals

- `#762A83` (purple) — dispersion parameter

## References

Ghosal, R., Ghosh, S. K., and Maiti, T. (2020). Two-part regression
models for longitudinal zero-inflated count data. *Journal of the Royal
Statistical Society: Series A*, **183**(4), 1603–1626.

Williams, M. R. and Savitsky, T. D. (2021). Uncertainty estimation for
pseudo-Bayesian inference under complex sampling. *International
Statistical Review*, **89**(1), 72–107.
[doi:10.1111/insr.12376](https://doi.org/10.1111/insr.12376)

## See also

[`summary.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/summary.hbb_fit.md)
for the tabular summary reused by the coefficients plot,
[`residuals.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/residuals.hbb_fit.md)
for residual computation,
[`fitted.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/fitted.hbb_fit.md)
for fitted values,
[`ppc`](https://joonho112.github.io/hurdlebb/reference/ppc.md) and
[`plot.hbb_ppc`](https://joonho112.github.io/hurdlebb/reference/plot.hbb_ppc.md)
for posterior predictive checks.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)

# Forest plot with posterior intervals
plot(fit)

# Forest plot with sandwich (Wald) intervals
sand <- sandwich_variance(fit)
plot(fit, type = "coefficients", sandwich = sand)

# Trace plots for a subset of parameters
plot(fit, type = "trace", pars = c("alpha[1]", "beta[1]"))

# Residual diagnostics
plot(fit, type = "residuals")

# Random effects caterpillar (SVC models only)
plot(fit, type = "random_effects")
} # }
```
