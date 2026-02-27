# Cholesky Posterior Recalibration for Survey-Weighted Inference

In survey-weighted Bayesian models the pseudo-posterior covariance
\\\Sigma\_{\mathrm{MCMC}}\\ typically over-covers because the prior
contributes substantially to the curvature (Prior Inflation ratios of
600–6500 are common for fixed effects). The Cholesky correction
constructs an affine map \$\$ \theta^{\*(m)} = \hat\theta +
A\bigl(\theta^{(m)} - \hat\theta\bigr), \quad A = L\_{\mathrm{sand}}\\
L\_{\mathrm{MCMC}}^{-1}, \$\$ where \\L\_{\mathrm{sand}}\\ and
\\L\_{\mathrm{MCMC}}\\ are the lower-triangular Cholesky factors of the
sandwich variance \\V\_{\mathrm{sand}}\\ and the MCMC covariance
\\\Sigma\_{\mathrm{MCMC}}\\, respectively.

## Usage

``` r
cholesky_correct(fit, sandwich, level = 0.95)
```

## Arguments

- fit:

  An object of class `"hbb_fit"` returned by
  [`hbb`](https://joonho112.github.io/hurdlebb/reference/hbb.md) or
  similar model-fitting function. Must contain a CmdStanR fit accessible
  via `fit$fit$draws()`.

- sandwich:

  An object of class `"hbb_sandwich"` returned by
  [`sandwich_variance`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md).
  Must contain `V_sand`, `Sigma_MCMC`, `H_obs_inv`, `DER`, and
  `param_labels`.

- level:

  Numeric scalar in \\(0,1)\\. Confidence level for interval
  construction. Default is `0.95`.

## Value

An S3 object of class `"hbb_cholesky"` containing:

- `theta_corrected`:

  Numeric matrix of dimension M by D containing the Cholesky-corrected
  posterior draws.

- `theta_hat`:

  Numeric vector of length D: posterior means (column means of the
  original draws).

- `A`:

  The D by D affine transformation matrix \\A =
  L\_{\mathrm{sand}}\\L\_{\mathrm{MCMC}}^{-1}\\.

- `L_MCMC`:

  Lower-triangular Cholesky factor of the MCMC posterior covariance.

- `L_sand`:

  Lower-triangular Cholesky factor of the sandwich variance.

- `comparison_table`:

  Data frame with one row per parameter containing naive, corrected, and
  Wald confidence intervals, width ratios, DER, DER-vs-MCMC, prior
  inflation, and the square root of DER.

- `verification`:

  Named list of numerical verification checks: mean preservation error,
  relative variance error, A algebraic identity error, and logical
  pass/fail flags.

- `level`:

  The confidence level used.

- `D`:

  Integer: number of parameters (2P + 1).

- `M`:

  Integer: number of MCMC draws.

- `P`:

  Integer: number of covariates per margin.

- `pd_corrections`:

  List of PD correction details for each matrix (whether ridge or nearPD
  was applied, condition numbers).

## Details

Applies the affine Cholesky correction of Williams and Savitsky (2021,
Theorem 4.1) to MCMC posterior draws, replacing the prior-dominated
posterior covariance with the design-consistent sandwich variance while
preserving the posterior mean.

## Mathematical properties

The transformation satisfies three exact algebraic identities:

1.  **Mean preservation.** \\E\[\theta^\*\] = \hat\theta\\ because the
    affine map is centred at the posterior mean.

2.  **Covariance recovery.** \$\$ \mathrm{Cov}(\theta^\*) =
    A\\\Sigma\_{\mathrm{MCMC}}\\A^\top =
    L\_{\mathrm{sand}}\\L\_{\mathrm{MCMC}}^{-1}\\
    \Sigma\_{\mathrm{MCMC}}\\
    L\_{\mathrm{MCMC}}^{-\top}\\L\_{\mathrm{sand}}^\top =
    V\_{\mathrm{sand}}. \$\$ The middle cancellation uses
    \\\Sigma\_{\mathrm{MCMC}} =
    L\_{\mathrm{MCMC}}\\L\_{\mathrm{MCMC}}^\top\\.

3.  **Affine equivariance.** Quantile ordering is preserved under
    monotone marginal transformations, so credible intervals from the
    corrected draws are valid frequentist confidence intervals with
    correct nominal coverage.

## Prior domination and shrinkage

For parameters where the prior dominates the likelihood (Prior Inflation
\\\mathrm{PI}\_p = \Sigma\_{\mathrm{MCMC},pp} /
H\_{\mathrm{obs},pp}^{-1} \gg 1\\), the diagonal of \\A\\ is much less
than unity (typically 0.01–0.06). This *shrinks* the over-dispersed MCMC
draws towards \\\hat\theta\\, replacing prior-inflated credible
intervals with data-driven confidence intervals. This is the expected
and correct behaviour: the sandwich variance \\V\_{\mathrm{sand}}\\ is
smaller than both \\\Sigma\_{\mathrm{MCMC}}\\ and
\\H\_{\mathrm{obs}}^{-1}\\ when design effects are moderate.

For parameters where the prior contributes little (e.g., `log_kappa`
with \\\mathrm{PI} \approx 2.3\\), \\A\_{pp} \approx 1.06\\, producing
slight inflation consistent with design-effect adjustment.

## Design Effect Ratio

The DER is the key diagnostic for the correction magnitude: \$\$
\mathrm{DER}\_p =
\frac{V\_{\mathrm{sand},pp}}{H\_{\mathrm{obs},pp}^{-1}}. \$\$ Values in
the range 1–5 are typical for complex survey designs.

## Relationship to Wald inference

Wald confidence intervals \\\hat\theta_p \pm
z\_{1-\alpha/2}\\\sqrt{V\_{\mathrm{sand},pp}}\\ are algebraically
equivalent to the marginal quantiles of the corrected draws in large
samples. They are recommended as the primary reporting device because
they do not depend on MCMC sampling variability.

## References

Williams, M. R. and Savitsky, T. D. (2021). Uncertainty estimation for
pseudo-Bayesian inference under complex sampling. *International
Statistical Review*, **89**(1), 72–107.
[doi:10.1111/insr.12376](https://doi.org/10.1111/insr.12376)

Ghosal, R., Ghosh, S. K., and Maiti, T. (2020). Two-part regression
models for longitudinal zero-inflated count data. *Journal of the Royal
Statistical Society: Series A*, **183**(4), 1603–1626.

## See also

[`sandwich_variance`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md)
for the sandwich variance computation,
[`compute_wald_ci`](https://joonho112.github.io/hurdlebb/reference/compute_wald_ci.md)
for standalone Wald confidence intervals,
[`print.hbb_cholesky`](https://joonho112.github.io/hurdlebb/reference/print.hbb_cholesky.md)
for the print method.

## Examples

``` r
if (FALSE) { # \dontrun{
# After fitting and computing sandwich variance:
sand <- sandwich_variance(fit)
chol_obj <- cholesky_correct(fit, sand, level = 0.95)
print(chol_obj)

# Compare interval widths
chol_obj$comparison_table[, c("parameter", "naive_width",
                              "corrected_width", "wald_width")]
} # }
```
