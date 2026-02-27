# Fit a Hurdle Beta-Binomial Model

Fits a two-part hurdle Beta-Binomial model using Hamiltonian Monte Carlo
(HMC) via CmdStan. The model comprises an extensive margin (Bernoulli
for whether a center serves infant/toddler children) and an intensive
margin (zero-truncated Beta-Binomial for the enrollment share among
servers). Both margins share the same design matrix with separate
coefficient vectors.

## Usage

``` r
hbb(
  formula,
  data,
  weights = NULL,
  stratum = NULL,
  psu = NULL,
  state_data = NULL,
  prior = default_prior(),
  chains = 4L,
  iter_warmup = 1000L,
  iter_sampling = 1000L,
  adapt_delta = 0.95,
  max_treedepth = 12L,
  seed = NULL,
  refresh = NULL,
  cpp_options = list(),
  ...
)
```

## Arguments

- formula:

  A two-sided formula of the form `y | trials(n) ~ predictors`, or an
  object of class
  [hbb_formula](https://joonho112.github.io/hurdlebb/reference/hbb_formula.md).
  See
  [`hbb_formula()`](https://joonho112.github.io/hurdlebb/reference/hbb_formula.md)
  for the full syntax including random effects and policy moderators.

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
  moderators. Required if the formula includes `state_level()` terms.
  Must contain a column matching the grouping variable for merging.

- prior:

  A prior specification object of class
  [hbb_prior](https://joonho112.github.io/hurdlebb/reference/hbb_prior.md),
  or `NULL` to use
  [`default_prior()`](https://joonho112.github.io/hurdlebb/reference/default_prior.md).
  See
  [`hbb_prior()`](https://joonho112.github.io/hurdlebb/reference/hbb_prior.md)
  for customisation.

- chains:

  Integer. Number of Markov chains. Default `4`.

- iter_warmup:

  Integer. Number of warmup iterations per chain. Default `1000`.

- iter_sampling:

  Integer. Number of post-warmup iterations per chain. Default `1000`.

- adapt_delta:

  Numeric in `(0, 1)`. Target average proposal acceptance probability
  during warmup adaptation. Higher values reduce divergent transitions
  at the cost of slower sampling. Default `0.95`.

- max_treedepth:

  Integer. Maximum tree depth for the NUTS sampler. Default `12`.

- seed:

  Integer or `NULL`. Random seed for reproducibility.

- refresh:

  Integer or `NULL`. How often to print progress (every `refresh`
  iterations). Set to `0` to suppress progress output. `NULL` uses the
  CmdStan default.

- cpp_options:

  Named list of C++ compilation options passed to
  [`cmdstanr::cmdstan_model()`](https://mc-stan.org/cmdstanr/reference/cmdstan_model.html).
  For example, `list(stan_threads = TRUE)` enables within-chain
  threading. Default [`list()`](https://rdrr.io/r/base/list.html).

- ...:

  Additional arguments passed to the CmdStanModel `$sample()` method.

## Value

An S3 object of class `"hbb_fit"`. See
[`is.hbb_fit()`](https://joonho112.github.io/hurdlebb/reference/is.hbb_fit.md)
and
[`print.hbb_fit()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_fit.md)
for methods. The object contains:

- `fit`:

  The `CmdStanMCMC` object from cmdstanr.

- `stan_data`:

  The list passed to Stan's `$sample()`.

- `hbb_data`:

  The full `hbb_data` object for downstream methods.

- `formula`:

  The `hbb_formula` object.

- `prior`:

  The `hbb_prior` object used.

- `model_type`:

  Character: one of `"base"`, `"weighted"`, `"svc"`, `"svc_weighted"`.

- `model_name`:

  Character: the Stan model file name (e.g., `"hbb_base"`).

- `call`:

  The matched call.

- `elapsed`:

  Numeric: wall-clock time in seconds.

## Details

### Model variants

The Stan model selected depends on the formula and whether survey
weights are provided:

|                    |                           |             |
|--------------------|---------------------------|-------------|
| **Model**          | **Formula pattern**       | **Weights** |
| `hbb_base`         | `y | trials(n) ~ x1 + x2` | No          |
| `hbb_weighted`     | same                      | Yes         |
| `hbb_svc`          | `... + (x1 | group)`      | No          |
| `hbb_svc_weighted` | `... + (x1 | group)`      | Yes         |

Note that random-intercept-only models (`(1 | state_id)`) currently use
the `base` or `weighted` variants, as the random intercept is absorbed
into the SVC framework only when random slopes are present.

### Workflow

`hbb()` performs the following steps:

1.  Parses the formula (if a raw formula rather than an `hbb_formula`).

2.  Validates the prior specification and MCMC arguments.

3.  Prepares Stan data via
    [`prepare_stan_data()`](https://joonho112.github.io/hurdlebb/reference/prepare_stan_data.md).

4.  Maps the prepared data and priors to the Stan data block.

5.  Compiles the appropriate Stan model (cached after first use).

6.  Runs MCMC sampling via CmdStan.

7.  Checks for MCMC pathologies (divergences, treedepth, E-BFMI).

8.  Returns an `hbb_fit` object for downstream analysis.

### Compilation

Stan models are compiled on first use and cached persistently. See
[`hbb_compile()`](https://joonho112.github.io/hurdlebb/reference/hbb_compile.md)
for details on the caching strategy and explicit pre-compilation.

## See also

[`hbb_formula()`](https://joonho112.github.io/hurdlebb/reference/hbb_formula.md),
[`prepare_stan_data()`](https://joonho112.github.io/hurdlebb/reference/prepare_stan_data.md),
[`hbb_prior()`](https://joonho112.github.io/hurdlebb/reference/hbb_prior.md),
[`hbb_compile()`](https://joonho112.github.io/hurdlebb/reference/hbb_compile.md),
[`is.hbb_fit()`](https://joonho112.github.io/hurdlebb/reference/is.hbb_fit.md)

Other fitting:
[`is.hbb_fit()`](https://joonho112.github.io/hurdlebb/reference/is.hbb_fit.md),
[`print.hbb_fit()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(nsece_synth_small, package = "hurdlebb")

# Base model (no weights, no random effects)
fit1 <- hbb(
  y | trials(n_trial) ~ poverty + urban,
  data = nsece_synth_small,
  chains = 2, iter_warmup = 500, iter_sampling = 500
)
print(fit1)

# Weighted model
fit2 <- hbb(
  y | trials(n_trial) ~ poverty + urban,
  data = nsece_synth_small,
  weights = "weight",
  chains = 2, iter_warmup = 500, iter_sampling = 500
)

# SVC model with policy moderators
data(nsece_state_policy, package = "hurdlebb")
fit3 <- hbb(
  y | trials(n_trial) ~ poverty + urban +
    (poverty + urban | state_id) +
    state_level(mr_pctile),
  data = nsece_synth_small,
  state_data = nsece_state_policy,
  weights = "weight",
  chains = 2, iter_warmup = 500, iter_sampling = 500
)
} # }
```
