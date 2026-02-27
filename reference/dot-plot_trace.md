# MCMC trace plots for fixed-effect parameters

Extracts posterior draws from the CmdStanMCMC object and plots the
iteration-by-iteration trace for each selected parameter, coloured by
chain. Useful for visual assessment of convergence and mixing.

## Usage

``` r
.plot_trace(object, pars)
```

## Arguments

- object:

  An `hbb_fit` object.

- pars:

  Character vector of CmdStanR parameter names (e.g., `"alpha[1]"`,
  `"beta[2]"`) or `NULL` for all fixed effects.

## Value

A ggplot object.

## Details

Uses `format = "array"` to extract draws with chain structure preserved,
avoiding the need for chain-ID reconstruction from flat matrices.

## Chain-ID reconstruction

CmdStanR `$draws(format = "array")` returns an \\\mathrm{iter} \times
\mathrm{chains} \times \mathrm{params}\\ array, so chain identity is
directly available.
