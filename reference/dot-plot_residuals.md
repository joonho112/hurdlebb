# Two-panel residual diagnostic plot

Constructs two diagnostic panels using existing S3 methods:

1.  **Residuals vs fitted**: Pearson residuals (vertical) against fitted
    values \\\hat{q}\_i \hat\mu_i\\ (horizontal). Under correct model
    specification, the residuals should scatter randomly around zero
    with no systematic pattern. A LOESS smoother (red) is overlaid to
    highlight trends.

2.  **Normal QQ plot**: Sorted Pearson residuals against theoretical
    standard normal quantiles \\\Phi^{-1}((i - 0.375)/(N + 0.25))\\.
    Departures from the diagonal indicate non-normality of the residual
    distribution.

## Usage

``` r
.plot_residuals(object)
```

## Arguments

- object:

  An `hbb_fit` object.

## Value

A patchwork object combining two ggplot panels.

## Details

The panels are combined via patchwork.

## Pearson residuals

The Pearson residual for observation \\i\\ is \$\$ r_i^P = \frac{y_i -
n_i \hat{q}\_i \hat\mu_i} {\sqrt{\widehat{\mathrm{Var}}(Y_i)}}, \$\$
where \\\widehat{\mathrm{Var}}(Y_i)\\ is the hurdle variance computed
from \\(n_i, \hat{q}\_i, \hat\mu_i, \hat\kappa)\\. See
[`hurdle_variance`](https://joonho112.github.io/hurdlebb/reference/hurdle_variance.md)
for the variance decomposition.
