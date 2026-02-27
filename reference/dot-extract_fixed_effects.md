# Build the fixed-effects summary table

Constructs a data frame with \\D\\ rows (one per fixed-effect parameter)
containing the posterior mean, standard error, confidence/credible
interval bounds, Rhat, and bulk ESS.

## Usage

``` r
.extract_fixed_effects(object, sandwich, level, param_labels)
```

## Arguments

- object:

  An `hbb_fit` object.

- sandwich:

  An optional `hbb_sandwich` object (or `NULL`).

- level:

  Numeric confidence level in (0, 1).

- param_labels:

  Character vector of length D.

## Value

Data frame with columns: parameter, estimate, se, ci_lower, ci_upper,
rhat, ess_bulk.

## Details

When `sandwich` is provided, SEs come from
\\\sqrt{\mathrm{diag}(V\_{\mathrm{sand}})}\\ and CIs are Wald intervals:
\\\hat\theta_p \pm z\_{(1+\mathrm{level})/2} \cdot \mathrm{SE}\_p\\.
Otherwise, SEs are posterior standard deviations and CIs are
quantile-based credible intervals from the MCMC draws.

Rhat and ESS are always extracted from the CmdStanR `$summary()` method,
wrapped in `tryCatch`.
