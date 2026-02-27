# Compute posterior credible intervals via draw propagation

Propagates each MCMC draw through the inverse-link function to obtain
draw-level predictions, then computes empirical quantiles for the
credible interval.

## Usage

``` r
.predict_compute_interval(object, X_new, P, delta, type, level, ndraws)
```

## Arguments

- object:

  An `hbb_fit` object.

- X_new:

  Numeric matrix N_new x P.

- P:

  Integer: covariates per margin.

- delta:

  List from `.predict_resolve_delta` or NULL.

- type:

  One of "response", "extensive", "intensive".

- level:

  Numeric in (0,1).

- ndraws:

  Integer or NULL.

## Value

List with `lwr` (N_new-vector) and `upr` (N_new-vector).

## Algorithm

1.  Extract the M x D matrix of fixed-effect draws.

2.  If `ndraws` is specified and smaller than M, subsample M draws via
    deterministic evenly-spaced thinning.

3.  For each draw \\m\\, compute \\q_i^{(m)} = \mathrm{logistic}(x_i'
    \alpha^{(m)})\\ and \\\mu_i^{(m)} = \mathrm{logistic}(x_i'
    \beta^{(m)})\\, then form the prediction according to `type`.

4.  Compute \\(\alpha/2,\\ 1-\alpha/2)\\ quantiles across draws for each
    observation.

## Vectorisation

The draw-level computation is vectorised over observations using matrix
multiplication: for each draw \\m\\, the N_new-vector of predictions is
computed via \\X\_{\mathrm{new}} \alpha^{(m)}\\ etc., avoiding explicit
loops over observations.
