# Compute point predictions at the posterior mean

Computes predictions using the posterior mean of the fixed-effect
parameter vector \\\hat\theta\\, optionally incorporating posterior mean
state random effects.

## Usage

``` r
.predict_compute_point(object, X_new, P, delta, type)
```

## Arguments

- object:

  An `hbb_fit` object.

- X_new:

  Numeric matrix N_new x P.

- P:

  Integer: number of covariates per margin.

- delta:

  List from `.predict_resolve_delta` or NULL.

- type:

  One of "response", "extensive", "intensive".

## Value

Numeric vector of length N_new.

## Details

The computation is vectorised via BLAS matrix multiplication: \$\$
\hat\eta\_{\mathrm{ext}} = X\_{\mathrm{new}} \hat\alpha +
X\_{\mathrm{new}} \hat\delta\_{\mathrm{ext}}\[s(\cdot)\], \$\$ and
similarly for the intensive margin, followed by the logistic inverse
link.
