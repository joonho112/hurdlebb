# Prepare Stan Data for a Hurdle Beta-Binomial Model

Converts a formula specification and data frame into a validated list
suitable for passing to Stan. Constructs design matrices with optional
standardisation, builds group indices, normalises survey weights, and
computes all required dimensions.

## Usage

``` r
prepare_stan_data(
  formula,
  data,
  weights = NULL,
  stratum = NULL,
  psu = NULL,
  state_data = NULL,
  prior = NULL,
  center = TRUE,
  scale = TRUE
)
```

## Arguments

- formula:

  An object of class `"hbb_formula"` (from
  [`hbb_formula()`](https://joonho112.github.io/hurdlebb/reference/hbb_formula.md)),
  or a raw formula that will be parsed via
  [`hbb_formula()`](https://joonho112.github.io/hurdlebb/reference/hbb_formula.md).

- data:

  A data frame containing provider-level variables: the response,
  trials, fixed effects, and optionally the grouping variable.

- weights:

  Character string naming the column in `data` containing survey
  sampling weights, or `NULL` (default) for unweighted analysis.

- stratum:

  Character string naming the column in `data` containing sampling
  stratum identifiers, or `NULL`.

- psu:

  Character string naming the column in `data` containing primary
  sampling unit (PSU) identifiers, or `NULL`.

- state_data:

  A data frame of group-level (state-level) variables for policy
  moderators. Required if `formula$policy` is not `NULL`. Must contain a
  column matching the grouping variable for merging.

- prior:

  A prior specification list, or `NULL` (default) to use package
  defaults. (Prior specification API is under development.)

- center:

  Logical. If `TRUE` (default), subtract column means from non-intercept
  columns of X and V. The means are stored in the output for
  back-transformation.

- scale:

  Logical. If `TRUE` (default), divide non-intercept columns of X and V
  by their standard deviations after centering. The SDs are stored in
  the output for back-transformation.

## Value

An S3 object of class `"hbb_data"` (a list). See **Value** for full
field descriptions.

## Value fields

- `N`:

  Integer. Number of observations.

- `P`:

  Integer. Number of provider covariates including intercept.

- `S`:

  Integer. Number of groups (states). 1 if no grouping.

- `Q`:

  Integer. Number of policy covariates including intercept.

- `N_pos`:

  Integer. Number of positive observations (z == 1).

- `K`:

  Integer. Random effects dimension per group.

- `y`:

  Integer vector of length N. Response counts.

- `n_trial`:

  Integer vector of length N. Trial sizes.

- `z`:

  Integer vector of length N. Participation indicators.

- `X`:

  Numeric matrix N x P. Design matrix (column 1 = intercept).

- `state`:

  Integer vector of length N. Group indices 1..S.

- `V`:

  Numeric matrix S x Q, or `NULL` if no groups.

- `idx_pos`:

  Integer vector. Row indices where z == 1.

- `w_tilde`:

  Numeric vector of normalised weights, or `NULL`.

- `stratum_idx`:

  Integer vector of stratum indices, or `NULL`.

- `psu_idx`:

  Integer vector of PSU indices, or `NULL`.

- `n_strata`:

  Integer. Number of unique strata, or `NULL`.

- `n_psu`:

  Integer. Number of unique PSUs, or `NULL`.

- `x_center`:

  Named numeric vector of column means, or `NULL`.

- `x_scale`:

  Named numeric vector of column SDs, or `NULL`.

- `v_center`:

  Named numeric vector of V column means, or `NULL`.

- `v_scale`:

  Named numeric vector of V column SDs, or `NULL`.

- `formula`:

  The `hbb_formula` object.

- `prior`:

  The prior specification (or `NULL`).

- `group_levels`:

  Character vector of original group labels.

- `model_type`:

  Character. One of `"base"`, `"weighted"`, `"svc"`, `"svc_weighted"`.

## See also

[`hbb_formula()`](https://joonho112.github.io/hurdlebb/reference/hbb_formula.md),
[`validate_hbb_data()`](https://joonho112.github.io/hurdlebb/reference/validate_hbb_data.md)

Other data:
[`validate_hbb_data()`](https://joonho112.github.io/hurdlebb/reference/validate_hbb_data.md)

## Examples

``` r
# Minimal example with the small synthetic dataset
data(nsece_synth_small, package = "hurdlebb")
f <- hbb_formula(y | trials(n_trial) ~ poverty + urban)
d <- prepare_stan_data(f, nsece_synth_small)
d
#> Hurdle Beta-Binomial Data
#> -------------------------
#>   Observations (N) : 504 
#>   Covariates   (P) : 3 (intercept + 2 predictors) 
#>   Groups       (S) : 1 
#>   Policy vars  (Q) : 1 
#>   Positive obs     : 313 (62.1%) 
#>   Zero rate        : 0.379 
#>   RE dimension (K) : 0 
#>   Model type       : base 
#>   Survey weights   : no
#>   X standardised   : yes

# With grouping and weights
f2 <- hbb_formula(
  y | trials(n_trial) ~ poverty + urban + (1 | state_id)
)
d2 <- prepare_stan_data(f2, nsece_synth_small, weights = "weight")
d2
#> Hurdle Beta-Binomial Data
#> -------------------------
#>   Observations (N) : 504 
#>   Covariates   (P) : 3 (intercept + 2 predictors) 
#>   Groups       (S) : 51 
#>   Policy vars  (Q) : 1 
#>   Positive obs     : 313 (62.1%) 
#>   Zero rate        : 0.379 
#>   RE dimension (K) : 2 
#>   Model type       : weighted 
#>   Survey weights   : yes (range 0.05--14.92) 
#>   X standardised   : yes

# Full SVC model with policy moderators
data(nsece_state_policy, package = "hurdlebb")
f3 <- hbb_formula(
  y | trials(n_trial) ~ poverty + urban +
    (poverty + urban | state_id) +
    state_level(mr_pctile)
)
d3 <- prepare_stan_data(
  f3, nsece_synth_small,
  weights = "weight", state_data = nsece_state_policy
)
d3
#> Hurdle Beta-Binomial Data
#> -------------------------
#>   Observations (N) : 504 
#>   Covariates   (P) : 3 (intercept + 2 predictors) 
#>   Groups       (S) : 51 
#>   Policy vars  (Q) : 2 
#>   Positive obs     : 313 (62.1%) 
#>   Zero rate        : 0.379 
#>   RE dimension (K) : 6 
#>   Model type       : svc_weighted 
#>   Survey weights   : yes (range 0.05--14.92) 
#>   X standardised   : yes
```
