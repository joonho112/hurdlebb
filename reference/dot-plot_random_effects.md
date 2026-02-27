# Caterpillar plot of state-varying random effects

For SVC models (`model_type %in% c("svc", "svc_weighted")`), extracts
the posterior summaries (mean and quantile-based intervals) of the state
random effects \\\delta\_{s,k}\\ and displays them as a caterpillar plot
faceted by covariate.

## Usage

``` r
.plot_random_effects(object, level)
```

## Arguments

- object:

  An `hbb_fit` object (must be SVC model).

- level:

  Numeric in (0,1) for credible interval width.

## Value

A ggplot object.

## Details

Each panel shows all \\S\\ states ordered by their posterior mean, with
horizontal bars spanning the `level` credible interval. The vertical
reference line at zero aids interpretation: states with intervals that
exclude zero exhibit statistically meaningful departures from the
fixed-effect mean.
