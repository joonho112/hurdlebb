# Average Marginal Effects Decomposition for Hurdle Beta-Binomial Models

The hurdle Beta-Binomial model implies two channels through which
covariates affect the expected outcome. The extensive margin governs
whether a center serves infant/toddler children at all (\\q_i =
\mathrm{logistic}(X_i' \alpha)\\), while the intensive margin governs
the enrollment share conditional on participation (\\\mu_i =
\mathrm{logistic}(X_i' \beta)\\).

By the product rule, the marginal effect of covariate \\k\\ on \\E\[y_i
/ n_i\]\\ decomposes as: \$\$ \frac{\partial E\[y_i/n_i\]}{\partial
x\_{ik}} = \underbrace{\alpha_k \\ q_i(1-q_i) \\
\mu_i}\_{\text{extensive}} + \underbrace{\beta_k \\ \mu_i(1-\mu_i) \\
q_i}\_{\text{intensive}}. \$\$

The Average Marginal Effect averages this over all \\N\\ observations:
\$\$ \mathrm{AME}\_k = \frac{1}{N} \sum\_{i=1}^{N} \bigl\[\alpha_k \\
q_i(1-q_i) \\ \mu_i + \beta_k \\ \mu_i(1-\mu_i) \\ q_i\bigr\] =
\mathrm{AME}\_k^{\mathrm{ext}} + \mathrm{AME}\_k^{\mathrm{int}}. \$\$

## Usage

``` r
ame(fit, cholesky = NULL, sandwich = NULL, level = 0.95, n_draws = NULL)
```

## Arguments

- fit:

  An object of class `"hbb_fit"` returned by
  [`hbb`](https://joonho112.github.io/hurdlebb/reference/hbb.md). Must
  contain a CmdStanR fit with parameters `alpha[1:P]`, `beta[1:P]`,
  `log_kappa`, and an `hbb_data` element with design matrix `X`.

- cholesky:

  An optional object of class `"hbb_cholesky"` returned by
  [`cholesky_correct`](https://joonho112.github.io/hurdlebb/reference/cholesky_correct.md).
  If provided, uses the Cholesky-corrected draws
  (`cholesky$theta_corrected`) and posterior mean (`cholesky$theta_hat`)
  for AME computation. If `NULL` (default), uses raw MCMC draws from the
  fit object.

- sandwich:

  An optional object of class `"hbb_sandwich"` returned by
  [`sandwich_variance`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md).
  If provided, computes delta-method Wald confidence intervals for each
  AME component. If `NULL` (default), Wald results are omitted.

- level:

  Numeric scalar in \\(0, 1)\\. Confidence level for posterior credible
  intervals and Wald confidence intervals. Default is `0.95`.

- n_draws:

  Integer or `NULL`. Number of MCMC draws to use. If `NULL` (default),
  uses all available draws. If an integer, subsamples by systematic
  thinning to approximately `n_draws` draws.

## Value

An S3 object of class `"hbb_ame"` containing:

- `theta_draws_used`:

  Numeric matrix of dimension `M_use x D` containing the (possibly
  subsampled) draws.

- `theta_hat`:

  Numeric vector of length D: point estimate (posterior mean or Cholesky
  theta_hat).

- `ext_ame_draws`:

  Numeric matrix `M_use x P`: extensive AME draws for each covariate.

- `int_ame_draws`:

  Numeric matrix `M_use x P`: intensive AME draws for each covariate.

- `total_ame_draws`:

  Numeric matrix `M_use x P`: total AME draws for each covariate.

- `ext_summary`:

  Data frame with posterior summaries of extensive AME (7 columns:
  covariate, post_mean, post_median, ci_lo, ci_hi, post_sd,
  pr_positive).

- `int_summary`:

  Data frame with posterior summaries of intensive AME.

- `total_summary`:

  Data frame with posterior summaries of total AME.

- `decomp_table`:

  Data frame for non-intercept covariates with columns: covariate,
  ext_ame, ext_ci_lo, ext_ci_hi, int_ame, int_ci_lo, int_ci_hi,
  total_ame, total_ci_lo, total_ci_hi, ext_share, int_share,
  sign_pattern.

- `reversal_probs`:

  Named numeric vector: for each non-intercept covariate, the posterior
  probability of opposing signs between extensive and intensive
  components.

- `wald_summary`:

  Data frame of Wald-based AME inference, or `NULL` if `sandwich` was
  not provided.

- `ame_ext_hat`:

  Numeric vector of length P: extensive AME point estimates at
  `theta_hat`.

- `ame_int_hat`:

  Numeric vector of length P: intensive AME point estimates at
  `theta_hat`.

- `ame_total_hat`:

  Numeric vector of length P: total AME point estimates at `theta_hat`.

- `mean_q`:

  Numeric vector of length `M_use`: mean participation probability
  across observations per draw.

- `mean_mu`:

  Numeric vector of length `M_use`: mean conditional intensity across
  observations per draw.

- `N, P, D, M_total, M_use, level, cov_labels`:

  Dimensional and metadata scalars.

## Details

Computes the Average Marginal Effect (AME) of each covariate on the
expected IT enrollment share \\E\[y_i / n_i\] = q_i \mu_i\\, decomposed
into extensive-margin and intensive-margin contributions.

## Poverty Reversal

A covariate exhibits a "reversal" when its extensive and intensive AME
components have opposing signs. For the poverty variable in the NSECE
application: \\\alpha\_{\mathrm{poverty}} \< 0\\ (higher poverty reduces
participation) but \\\beta\_{\mathrm{poverty}} \> 0\\ (higher poverty
increases IT share among servers). The posterior probability of reversal
is \\\Pr(\mathrm{AME}\_k^{\mathrm{ext}} \< 0 \\\text{AND}\\
\mathrm{AME}\_k^{\mathrm{int}} \> 0)\\.

## Population-Average AME

This function computes the population-average AME (PA-AME) using global
(population-level) coefficients \\\alpha\\ and \\\beta\\ only, without
state random effects. This is consistent with the sandwich-corrected
inference framework and answers the question: "what is the average
marginal effect for a new (arbitrary) state?"

## Wald Inference

When `sandwich` is provided, the function also computes delta-method
Wald confidence intervals for each AME component. The numerical gradient
of \\\mathrm{AME}\_k(\theta)\\ with respect to \\\theta\\ is computed
via central differences (step size \\\epsilon = 10^{-5}\\), and the Wald
variance is: \$\$ \mathrm{Var}(\mathrm{AME}\_k) \approx \nabla\_\theta
\mathrm{AME}\_k(\hat\theta)^\top V\_{\mathrm{sand}}\\ \nabla\_\theta
\mathrm{AME}\_k(\hat\theta). \$\$

## Performance

The computation is vectorised over observations within each MCMC draw.
The key insight is that the N-vectors \\q'\\, \\\mu\\, \\\mu'\\, \\q\\
are independent of the covariate index \\k\\, so: \$\$
\mathrm{AME}\_k^{\mathrm{ext}} = \alpha_k \cdot \overline{q'(1-q) \mu},
\quad \mathrm{AME}\_k^{\mathrm{int}} = \beta_k \cdot
\overline{\mu'(1-\mu) q}, \$\$ reducing the per-draw cost from \\O(NP)\\
to \\O(N)\\.

## References

Ghosal, R., Ghosh, S. K., and Maiti, T. (2020). Two-part regression
models for longitudinal zero-inflated count data. *Journal of the Royal
Statistical Society: Series A*, **183**(4), 1603–1626.

## See also

[`ame_decomposition`](https://joonho112.github.io/hurdlebb/reference/ame_decomposition.md)
for extracting the decomposition table,
[`cholesky_correct`](https://joonho112.github.io/hurdlebb/reference/cholesky_correct.md)
for Cholesky-corrected draws,
[`sandwich_variance`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md)
for the sandwich variance.

## Examples

``` r
if (FALSE) { # \dontrun{
# After fitting:
fit <- hbb(y | trials(n_trial) ~ poverty + urban, data = my_data,
           weights = "weight")
sand <- sandwich_variance(fit)
chol <- cholesky_correct(fit, sand)

# Full AME decomposition with Wald comparison
ame_result <- ame(fit, cholesky = chol, sandwich = sand)
print(ame_result)

# Extract decomposition table
ame_decomposition(ame_result)
} # }
```
