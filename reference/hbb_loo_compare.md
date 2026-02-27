# Compare Hurdle Beta-Binomial Models via LOO-CV

Models are ranked from highest (best) to lowest (worst) ELPD. The ELPD
difference \\\Delta \mathrm{elpd}\_{ij} = \mathrm{elpd}\_i -
\mathrm{elpd}\_j\\ for the best model against each other model is
reported together with its standard error.

For consecutive pairs in the ranked list, the pairwise z-ratio is: \$\$
z\_{AB} = \frac{\Delta\\\mathrm{elpd}\_{AB}}
{\mathrm{SE}(\Delta\\\mathrm{elpd}\_{AB})}, \$\$ where \\\mathrm{SE}\\
is estimated from the pointwise ELPD differences via: \$\$
\mathrm{SE}(\Delta\\\mathrm{elpd}\_{AB}) = \sqrt{N} \cdot
\mathrm{SD}\bigl(\hat{\ell}\_i^{(A)} - \hat{\ell}\_i^{(B)}\bigr), \$\$
where \\\hat{\ell}\_i^{(A)}\\ is the per-observation LOO log-density for
model \\A\\. A magnitude \\\|z\| \> 2\\ is taken as substantial evidence
for the higher-ranked model (Vehtari et al., 2017).

## Usage

``` r
hbb_loo_compare(..., cores = 1L)
```

## Arguments

- ...:

  Two or more named `hbb_fit` objects, or a single named list of
  `hbb_fit` objects.

- cores:

  Integer; number of cores for parallel PSIS computation. Default `1`.

## Value

An S3 object of class `"hbb_loo_compare"` containing:

- `comparison`:

  Data frame from
  [`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html):
  models ranked best to worst by ELPD, with columns `model`, `elpd_loo`,
  `elpd_diff`, `se_diff`, `looic`.

- `loo_list`:

  Named list of `"psis_loo"` objects, one per model.

- `pairwise`:

  Data frame of pairwise z-ratios for consecutive model pairs (in
  ELPD-ranked order), with columns `comparison`, `delta_elpd`,
  `se_diff`, `z_ratio`.

- `pareto_k`:

  Named list of Pareto-k diagnostic vectors (one per model, length
  \\N\\).

- `model_names`:

  Character vector of model names (as supplied by the user).

- `n_models`:

  Integer: number of models compared.

## Details

Computes LOO-CV for each supplied `hbb_fit` object and calls
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
to rank models by expected log predictive density (ELPD). Also reports
pairwise z-statistics for consecutive model pairs.

## Input Format

Pass models as named arguments (e.g., `m0 = fit0, m1 = fit1, ...`) or as
a single named list. All models must share the same number of
observations \\N\\.

## References

Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*, **27**(5), 1413–1432.

## See also

[`loo.hbb_fit`](https://joonho112.github.io/hurdlebb/reference/loo.hbb_fit.md)
for single-model LOO-CV,
[`ppc`](https://joonho112.github.io/hurdlebb/reference/ppc.md) for
posterior predictive checks.

Other model-checking:
[`loo.hbb_fit()`](https://joonho112.github.io/hurdlebb/reference/loo.hbb_fit.md),
[`plot.hbb_ppc()`](https://joonho112.github.io/hurdlebb/reference/plot.hbb_ppc.md),
[`ppc()`](https://joonho112.github.io/hurdlebb/reference/ppc.md),
[`print.hbb_loo_compare()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_loo_compare.md),
[`print.hbb_ppc()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_ppc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit0 <- hbb(y | trials(n) ~ 1,              data = my_data)
fit1 <- hbb(y | trials(n) ~ poverty,        data = my_data)
fit2 <- hbb(y | trials(n) ~ poverty + urban, data = my_data)
comp <- hbb_loo_compare(m0 = fit0, m1 = fit1, m2 = fit2)
print(comp)
} # }
```
