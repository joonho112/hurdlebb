# Posterior Predictive Checks for Hurdle Beta-Binomial Models

The function evaluates two test statistics designed to probe the two
structural components of the hurdle Beta-Binomial model:

- Zero rate:

  \$\$T_0(\mathbf{y}) = \frac{1}{N} \sum\_{i=1}^{N} \mathbf{1}(y_i =
  0).\$\$ Targets the extensive margin: whether a center serves
  infant/toddler children at all. Under the fitted model,
  \\E\[T_0(\mathbf{y}^{\mathrm{rep}})\] = 1 - N^{-1}\sum_i q_i\\.

- IT share:

  \$\$T_s(\mathbf{y}) = \frac{1}{N^+} \sum\_{i:\\ y_i \> 0}
  \frac{y_i}{n_i},\$\$ where \\N^+ = \\\\i : y_i \> 0\\\\ is the count
  of participating providers. Targets the intensive margin: the
  enrollment share conditional on participation. Under the fitted model,
  \\E\[T_s(\mathbf{y}^{\mathrm{rep}}) \mid y^{\mathrm{rep}} \> 0\]
  \approx N^{-1}\sum_i \mu_i\\.

## Usage

``` r
ppc(
  fit,
  type = c("both", "zero_rate", "it_share"),
  method = c("stan", "simulate"),
  n_draws = NULL,
  level = 0.95
)
```

## Arguments

- fit:

  An object of class `"hbb_fit"` returned by
  [`hbb`](https://joonho112.github.io/hurdlebb/reference/hbb.md). Must
  contain a CmdStanMCMC fit with generated quantities `log_lik[N]` and
  `y_rep[N]`.

- type:

  Character string; which test statistics to compute. One of `"both"`
  (default), `"zero_rate"`, or `"it_share"`.

- method:

  Character string; how to obtain posterior predictive draws. One of
  `"stan"` (default, extracts `y_rep` from the Stan GQ block) or
  `"simulate"` (generates draws in R via
  [`rhurdle_betabinom`](https://joonho112.github.io/hurdlebb/reference/rhurdle_betabinom.md)).

- n_draws:

  Integer or `NULL`. Number of posterior draws to use. If `NULL`
  (default), all available draws are used. If specified, draws are
  subsampled by systematic thinning.

- level:

  Numeric scalar in \\(0, 1)\\. Coverage level for the posterior
  predictive interval. Default `0.95`.

## Value

An S3 object of class `"hbb_ppc"` containing:

- `observed`:

  Named list: `zero_rate` (scalar), `it_share` (scalar or `NA` if no
  positive observations), `N` (total observations), `N_pos` (count of
  positives).

- `predicted`:

  Named list with elements `zero_rate` and `it_share` (each `NULL` if
  not requested), where each non-null element is a list containing:
  `draws` (numeric vector of length \\M\\), `mean`, `median`, `ci`
  (named length-2 vector with entries `"lower"` and `"upper"`), `sd`.

- `coverage`:

  Named list with elements `zero_rate` and `it_share` (each `NULL` if
  not requested), where each non-null element contains `in_ci` (logical)
  and `ci` (named length-2 vector).

- `n_draws`:

  Integer; number of draws used.

- `n_total_draws`:

  Integer; total available draws before thinning.

- `type`:

  Character; the `type` argument used.

- `method`:

  Character; the `method` argument used.

- `level`:

  Numeric; the `level` argument used.

- `model_type`:

  Character; from `fit$model_type`.

## Details

Computes numerical posterior predictive checks (PPC) for a fitted
`hbb_fit` model by comparing observed test statistics to their posterior
predictive distributions.

## Theory and Motivation

A well-calibrated Bayesian model should produce posterior predictive
distributions \\p(\mathbf{y}^{\mathrm{rep}} \mid \mathbf{y})\\ that are
consistent with the observed data. Formally, the Bayesian p-value for a
test statistic \\T\\ is: \$\$ p_B =
\Pr\bigl(T(\mathbf{y}^{\mathrm{rep}}) \ge T(\mathbf{y}) \mid
\mathbf{y}\bigr), \$\$ which should be near \\0.5\\ for a well-specified
model and near \\0\\ or \\1\\ for systematic misspecification (Gelman et
al., 1996).

The zero rate \\T_0\\ and the IT share \\T_s\\ are chosen because they
directly correspond to the two structural parts of the hurdle model. A
defect in \\T_0\\ indicates misspecification of the Bernoulli
participation equation; a defect in \\T_s\\ indicates misspecification
of the zero-truncated Beta-Binomial intensity equation. Together they
provide a minimal but targeted diagnostic toolkit, following the
visualisation philosophy of Gabry et al. (2019).

## Coverage Criterion

A test statistic "passes" the PPC at level \\\alpha\\ if the observed
value falls within the \\(1-\alpha)\\ central posterior predictive
interval: \$\$ T(\mathbf{y}) \\\in\\
\bigl\[Q\_{\alpha/2}\bigl(T(\mathbf{y}^{\mathrm{rep}})\bigr),\\
Q\_{1-\alpha/2}\bigl(T(\mathbf{y}^{\mathrm{rep}})\bigr) \bigr\]. \$\$
Coverage is reported in the `coverage` element of the returned object.

## Extraction Methods

- `"stan"`:

  Extracts the posterior predictive replications
  \\\mathbf{y}^{\mathrm{rep}}\\ directly from the `y_rep[N]` array in
  the Stan generated quantities block. This is the preferred method: it
  uses the exact joint posterior and incurs no additional R-side
  simulation cost.

- `"simulate"`:

  Reconstructs \\\mathbf{y}^{\mathrm{rep}}\\ in R by composing posterior
  draws of \\(\alpha, \beta, \log\kappa)\\ with the hurdle Beta-Binomial
  PMF via
  [`rhurdle_betabinom`](https://joonho112.github.io/hurdlebb/reference/rhurdle_betabinom.md).
  Useful for cross-validating the Stan GQ block and for models where GQ
  draws are unavailable.

## References

Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A.
(2019). Visualisation in Bayesian workflow. *Journal of the Royal
Statistical Society: Series A*, **182**(2), 389–402.

Gelman, A., Meng, X.-L., and Stern, H. (1996). Posterior predictive
assessment of model fitness via realized discrepancies. *Statistica
Sinica*, **6**(4), 733–760.

Ghosal, R., Ghosh, S. K., and Maiti, T. (2020). Two-part regression
models for longitudinal zero-inflated count data. *Journal of the Royal
Statistical Society: Series A*, **183**(4), 1603–1626.

## See also

[`loo.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/loo.hbb_fit.md)
for LOO-CV model comparison,
[`hbb`](https://joonho112.github.io/hurdlebb/reference/hbb.md) for model
fitting,
[`rhurdle_betabinom`](https://joonho112.github.io/hurdlebb/reference/rhurdle_betabinom.md)
for the hurdle BetaBinomial sampler.

Other model-checking:
[`hbb_loo_compare()`](https://joonho112.github.io/hurdlebb/reference/hbb_loo_compare.md),
[`loo.hbb_fit()`](https://joonho112.github.io/hurdlebb/reference/loo.hbb_fit.md),
[`plot.hbb_ppc()`](https://joonho112.github.io/hurdlebb/reference/plot.hbb_ppc.md),
[`print.hbb_loo_compare()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_loo_compare.md),
[`print.hbb_ppc()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_ppc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit a model
fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data,
           weights = "weight")

# Full PPC (both statistics, Stan draws)
ppc_result <- ppc(fit)
print(ppc_result)
plot(ppc_result)

# Zero-rate only, using 500 draws
ppc_zr <- ppc(fit, type = "zero_rate", n_draws = 500)

# Cross-check via R simulation
ppc_sim <- ppc(fit, method = "simulate", n_draws = 200)
} # }
```
