# Extract MCMC diagnostics from an hbb_fit

Safely extracts divergent transitions, max treedepth hits, E-BFMI, max
Rhat, and min ESS (bulk and tail) from the CmdStanMCMC object. All
extractions are wrapped in `tryCatch` so that this function never fails,
even if the fit object is corrupted or output files have been moved.

## Usage

``` r
.extract_diagnostics(object)
```

## Arguments

- object:

  An `hbb_fit` object.

## Value

Named list with elements: `n_divergent`, `n_max_treedepth`, `ebfmi`,
`max_rhat`, `min_ess_bulk`, `min_ess_tail`. Missing diagnostics are set
to `NA`.
