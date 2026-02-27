# Specify Priors for a Hurdle Beta-Binomial Model

Creates a prior specification object that controls the prior
distributions used in the Stan model. All arguments are optional;
unspecified components receive the defaults from
[`default_prior()`](https://joonho112.github.io/hurdlebb/reference/default_prior.md).

## Usage

``` r
hbb_prior(
  alpha = NULL,
  beta = NULL,
  log_kappa = NULL,
  gamma = NULL,
  tau = NULL,
  lkj_eta = NULL
)
```

## Arguments

- alpha:

  Prior on the extensive-margin fixed effects \\\alpha\\ (logit scale).
  A named list with elements `dist` (character), `mean` (numeric
  scalar), and `sd` (positive numeric scalar). Currently only
  `dist = "normal"` is supported.

- beta:

  Prior on the intensive-margin fixed effects \\\beta\\ (logit scale).
  Same structure as `alpha`.

- log_kappa:

  Prior on the log-concentration parameter \\\log \kappa\\. Same
  structure as `alpha`.

- gamma:

  Prior on the policy moderator coefficients \\\Gamma\\. Same structure
  as `alpha`.

- tau:

  Prior on the random effect scale parameters \\\tau\\. Same structure
  as `alpha`. Note: in Stan, \\\tau\\ is declared with `lower = 0`, so a
  `Normal(mean, sd)` prior becomes a **half-normal** prior in practice.
  The `mean` component is typically set to 0 for half-normal priors.

- lkj_eta:

  Concentration parameter \\\eta \> 0\\ for the LKJ prior on the
  Cholesky factor of the correlation matrix \\L\_\Omega\\. Higher values
  shrink the correlation toward the identity matrix. A scalar \\\eta =
  1\\ gives a uniform prior over correlation matrices; \\\eta = 2\\
  (default) mildly favours the identity.

## Value

An S3 object of class `"hbb_prior"` with components:

- `alpha`:

  Named list with `dist`, `mean`, `sd`.

- `beta`:

  Named list with `dist`, `mean`, `sd`.

- `log_kappa`:

  Named list with `dist`, `mean`, `sd`.

- `gamma`:

  Named list with `dist`, `mean`, `sd`.

- `tau`:

  Named list with `dist`, `mean`, `sd`.

- `lkj_eta`:

  Positive numeric scalar.

## Details

### Extensibility

The `dist` field in each component currently accepts only `"normal"`.
The structure is designed to accommodate future distributions (e.g.,
Student-t) without breaking the API.

### Model-specific usage

Not all parameters are present in every model variant:

- **hbb_base** and **hbb_weighted**: use `alpha`, `beta`, `log_kappa`
  only.

- **hbb_svc** and **hbb_svc_weighted**: additionally use `gamma`, `tau`,
  and `lkj_eta`.

Unused prior components are silently ignored by the fitting function.

## See also

[`default_prior()`](https://joonho112.github.io/hurdlebb/reference/default_prior.md)
for the package defaults.

Other priors:
[`default_prior()`](https://joonho112.github.io/hurdlebb/reference/default_prior.md)

## Examples

``` r
# Default priors
default_prior()
#> Hurdle Beta-Binomial Prior Specification
#> ----------------------------------------
#>   alpha     ~ Normal(0, 2) 
#>   beta      ~ Normal(0, 2) 
#>   log_kappa ~ Normal(2, 1.5) 
#>   gamma     ~ Normal(0, 1) 
#>   tau       ~ HalfNormal(0, 1) 
#>   L_Omega   ~ LKJ(2) 

# Tighter prior on fixed effects
hbb_prior(alpha = list(dist = "normal", mean = 0, sd = 1))
#> Hurdle Beta-Binomial Prior Specification
#> ----------------------------------------
#>   alpha     ~ Normal(0, 1) 
#>   beta      ~ Normal(0, 2) 
#>   log_kappa ~ Normal(2, 1.5) 
#>   gamma     ~ Normal(0, 1) 
#>   tau       ~ HalfNormal(0, 1) 
#>   L_Omega   ~ LKJ(2) 

# Custom prior on dispersion and LKJ
hbb_prior(
  log_kappa = list(dist = "normal", mean = 3, sd = 1),
  lkj_eta   = 4
)
#> Hurdle Beta-Binomial Prior Specification
#> ----------------------------------------
#>   alpha     ~ Normal(0, 2) 
#>   beta      ~ Normal(0, 2) 
#>   log_kappa ~ Normal(3, 1) 
#>   gamma     ~ Normal(0, 1) 
#>   tau       ~ HalfNormal(0, 1) 
#>   L_Omega   ~ LKJ(4) 

# Override everything
p <- hbb_prior(
  alpha     = list(dist = "normal", mean = 0, sd = 1),
  beta      = list(dist = "normal", mean = 0, sd = 1),
  log_kappa = list(dist = "normal", mean = 2, sd = 1),
  gamma     = list(dist = "normal", mean = 0, sd = 0.5),
  tau       = list(dist = "normal", mean = 0, sd = 0.5),
  lkj_eta   = 5
)
print(p)
#> Hurdle Beta-Binomial Prior Specification
#> ----------------------------------------
#>   alpha     ~ Normal(0, 1) 
#>   beta      ~ Normal(0, 1) 
#>   log_kappa ~ Normal(2, 1) 
#>   gamma     ~ Normal(0, 0.5) 
#>   tau       ~ HalfNormal(0, 0.5) 
#>   L_Omega   ~ LKJ(5) 
```
