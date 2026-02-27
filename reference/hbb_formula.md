# Parse a Hurdle Beta-Binomial Model Formula

Converts a brms-style formula into a structured specification for a
two-part hurdle Beta-Binomial model. The fixed effects specified in the
formula enter **both** model margins (extensive and intensive) with
separate coefficient vectors \\\alpha\\ and \\\beta\\.

## Usage

``` r
hbb_formula(formula)
```

## Arguments

- formula:

  A two-sided formula of the form `y | trials(n) ~ predictors`. See
  **Details** for the full syntax.

## Value

An S3 object of class `"hbb_formula"` with components:

- `response`:

  Character. The response variable name.

- `trials`:

  Character. The trials variable name.

- `fixed`:

  Character vector of fixed effect term names, **excluding** the
  intercept (which is always implicit).

- `group`:

  Character or `NULL`. The grouping variable for random effects.

- `random`:

  Character vector of random effect terms, or `NULL` if no random
  effects. Contains `"1"` for intercept-only random effects, or variable
  names for state-varying coefficients.

- `policy`:

  Character vector of policy moderator term names, or `NULL` if none
  specified.

- `formula`:

  The original formula object.

- `svc`:

  Logical. `TRUE` if the model includes random slopes (state-varying
  coefficients), not just a random intercept.

## Details

### Formula syntax

The left-hand side **must** specify both the response variable and the
trials variable using the `| trials()` syntax:

    y | trials(n_trial) ~ ...

The right-hand side supports four types of terms:

1.  **Fixed effects**: Standard R formula terms (e.g., `x1 + x2`). An
    intercept is always included automatically. Use `y | trials(n) ~ 1`
    for an intercept-only model.

2.  **Random intercept**: `(1 | group)` specifies a random intercept
    grouped by `group`.

3.  **State-varying coefficients (SVC)**: `(x1 + x2 | group)` specifies
    random slopes for `x1` and `x2` (plus the implicit random
    intercept). All random slope variables must also appear as fixed
    effects.

4.  **Policy moderators**: `state_level(v1 + v2)` specifies cross-level
    interaction terms that predict the random effects. These require
    that random effects are also specified.

### Model mapping

The parsed formula maps to the model as follows:

- Fixed effects \\\to\\ the \\X\\ matrix (with auto-prepended intercept
  column)

- The grouping variable \\\to\\ the state index vector

- Policy moderators \\\to\\ the \\V\\ matrix (with auto-prepended
  intercept column)

- Both hurdle margins share the same \\X\\ design matrix

## See also

[`validate_hbb_formula()`](https://joonho112.github.io/hurdlebb/reference/validate_hbb_formula.md)
to check a formula against a dataset.

Other formula:
[`validate_hbb_formula()`](https://joonho112.github.io/hurdlebb/reference/validate_hbb_formula.md)

## Examples

``` r
# Intercept-only model
f1 <- hbb_formula(y | trials(n_trial) ~ 1)
f1
#> Hurdle Beta-Binomial Formula
#> ----------------------------
#>   Response : y 
#>   Trials   : n_trial 
#>   Fixed    : (intercept only)
#>   Random   : none
#>   Policy   : none
#>   SVC      : FALSE 

# Fixed effects only
f2 <- hbb_formula(y | trials(n_trial) ~ poverty + urban)
f2
#> Hurdle Beta-Binomial Formula
#> ----------------------------
#>   Response : y 
#>   Trials   : n_trial 
#>   Fixed    : intercept + poverty + urban 
#>   Random   : none
#>   Policy   : none
#>   SVC      : FALSE 

# Random intercept
f3 <- hbb_formula(y | trials(n_trial) ~ poverty + (1 | state_id))
f3
#> Hurdle Beta-Binomial Formula
#> ----------------------------
#>   Response : y 
#>   Trials   : n_trial 
#>   Fixed    : intercept + poverty 
#>   Random   : (1 | state_id )   [random intercept]
#>   Policy   : none
#>   SVC      : FALSE 

# SVC with policy moderators
f4 <- hbb_formula(
  y | trials(n_trial) ~ poverty + urban +
    (poverty + urban | state_id) +
    state_level(mr_pctile + tiered_reim)
)
f4
#> Hurdle Beta-Binomial Formula
#> ----------------------------
#>   Response : y 
#>   Trials   : n_trial 
#>   Fixed    : intercept + poverty + urban 
#>   Random   : (1 + poverty + urban | state_id )   [SVC]
#>   Policy   : state_level( mr_pctile + tiered_reim )
#>   SVC      : TRUE 
```
