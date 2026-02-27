# Leave-One-Out Cross-Validation for Hurdle Beta-Binomial Models

PSIS-LOO approximates the leave-one-out predictive density: \$\$ p(y_i
\mid \mathbf{y}\_{-i}) = \int p(y_i \mid \theta) \\ p(\theta \mid
\mathbf{y}\_{-i}) \\ d\theta, \$\$ using importance sampling with
weights proportional to \\1/p(y_i \mid \theta)\\. Raw importance weights
have heavy tails; PSIS stabilises them by fitting a generalised Pareto
distribution to the largest weights and replacing the upper tail.

The expected log predictive density (ELPD) is: \$\$
\widehat{\mathrm{elpd}}\_{\mathrm{LOO}} = \sum\_{i=1}^{N} \log p(y_i
\mid \mathbf{y}\_{-i}), \$\$ and the LOO information criterion is
\\\mathrm{LOOIC} = -2 \cdot \widehat{\mathrm{elpd}}\_{\mathrm{LOO}}\\.
Models with higher (less negative) ELPD predict held-out data better.

## Usage

``` r
# S3 method for class 'hbb_fit'
loo(x, ..., r_eff = TRUE, cores = 1L)
```

## Arguments

- x:

  An object of class `"hbb_fit"` returned by
  [`hbb`](https://joonho112.github.io/hurdlebb/reference/hbb.md). Must
  contain a CmdStanMCMC fit with the generated quantity `log_lik[N]`.

- ...:

  Additional arguments passed to
  [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html) (e.g.,
  `save_psis = TRUE`).

- r_eff:

  Logical; if `TRUE` (default), compute relative effective sample sizes
  to improve PSIS accuracy. If `FALSE`, sets all `r_eff = 1`.

- cores:

  Integer; number of cores for parallel PSIS computation. Default is
  `1`.

## Value

An object of class `"psis_loo"` as returned by
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html). Key
elements:

- `estimates`:

  3 x 2 matrix with rows `elpd_loo`, `p_loo`, `looic` and columns
  `Estimate`, `SE`.

- `diagnostics`:

  List with `pareto_k` (N-vector) and `n_eff` (N-vector).

- `pointwise`:

  N x 5 matrix of per-observation LOO summaries.

## Details

Computes Pareto-smoothed importance sampling LOO-CV (PSIS-LOO) for a
fitted `hbb_fit` object using the pointwise log-likelihood `log_lik[N]`
computed in the Stan generated quantities block.

## Pareto-k Diagnostic

The Pareto shape parameter \\k\\ for each observation provides a
diagnostic of LOO reliability:

- \\k \< 0.5\\: Good –LOO estimate is reliable.

- \\0.5 \le k \< 0.7\\: OK –LOO estimate is moderately reliable; the
  finite-moment condition is met.

- \\k \ge 0.7\\: Bad –LOO estimate may be unreliable; the importance
  weight distribution has infinite variance. Consider moment matching
  ([`loo::loo_moment_match`](https://mc-stan.org/loo/reference/loo_moment_match.html))
  or \\K\\-fold CV.

Pareto-k values are automatically reported via `cli` upon completion.

## Relative Effective Sample Size

When `r_eff = TRUE` (default), the relative effective sample size
\\r\_{\mathrm{eff},i} = \mathrm{ESS}\_i / M\\ is computed via
[`loo::relative_eff()`](https://mc-stan.org/loo/reference/relative_eff.html)
and passed to [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html)
to correct for within-chain autocorrelation. The correction uses the
chain membership vector constructed from CmdStanR metadata (draws are in
row-major chain order: chain 1 first, then chain 2, etc.).

## Survey-Weighted Models

For survey-weighted models the `log_lik[N]` array in the Stan GQ block
contains the *unweighted* nominal log-likelihoods \\\log p(y_i \mid
\theta)\\, not the pseudo-log-likelihoods. LOO therefore provides model
selection among alternative HBB specifications on a common likelihood
scale, but does not account for survey design weighting. For
design-based inference use the sandwich-corrected Wald confidence
intervals from
[`sandwich_variance`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md).

## References

Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*, **27**(5), 1413–1432.

Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A.
(2019). Visualisation in Bayesian workflow. *Journal of the Royal
Statistical Society: Series A*, **182**(2), 389–402.

## See also

[`hbb_loo_compare`](https://joonho112.github.io/hurdlebb/reference/hbb_loo_compare.md)
for multi-model LOO comparison,
[`ppc`](https://joonho112.github.io/hurdlebb/reference/ppc.md) for
posterior predictive checks,
[`loo::loo`](https://mc-stan.org/loo/reference/loo.html) for the
underlying PSIS-LOO engine.

Other model-checking:
[`hbb_loo_compare()`](https://joonho112.github.io/hurdlebb/reference/hbb_loo_compare.md),
[`plot.hbb_ppc()`](https://joonho112.github.io/hurdlebb/reference/plot.hbb_ppc.md),
[`ppc()`](https://joonho112.github.io/hurdlebb/reference/ppc.md),
[`print.hbb_loo_compare()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_loo_compare.md),
[`print.hbb_ppc()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_ppc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data)
loo_result <- loo(fit)
print(loo_result)
} # }
```
