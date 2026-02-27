# Resolve state random effects for prediction

For non-SVC models or when `state` is `NULL`, returns `NULL`. For SVC
models with state information, maps state labels to integer indices via
`hbb_data$group_levels` and returns the posterior mean delta matrix
along with the state index vector.

## Usage

``` r
.predict_resolve_delta(object, state, N_new, newdata = NULL)
```

## Arguments

- object:

  An `hbb_fit` object.

- state:

  State specification (character vector, column name, or NULL).

- N_new:

  Number of prediction observations.

- newdata:

  Data frame or NULL.

## Value

`NULL` if no SVC adjustment, or a list with:

- `delta_hat`:

  S x K matrix of posterior mean deltas.

- `state_idx`:

  Integer vector of length N_new mapping observations to state indices.
