# Forest plot of fixed-effect estimates

Constructs a ggplot forest plot by calling
[`summary.hbb_fit()`](https://joonho112.github.io/hurdlebb/reference/summary.hbb_fit.md)
and extracting the `fixed_effects` data frame. Parameters are grouped
into three facets: Extensive (\\\alpha\\), Intensive (\\\beta\\), and
Dispersion (\\\log\kappa\\). A vertical reference line at zero
facilitates significance assessment.

## Usage

``` r
.plot_coefficients(object, sandwich, level, margin)
```

## Arguments

- object:

  An `hbb_fit` object.

- sandwich:

  An `hbb_sandwich` or NULL.

- level:

  Confidence level in (0,1).

- margin:

  One of "both", "extensive", "intensive".

## Value

A ggplot object.
