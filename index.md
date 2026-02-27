# hurdlebb

**hurdlebb** fits Bayesian two-part hurdle Beta-Binomial models for
bounded discrete proportions with structural zeros. It is designed for
complex survey data where many units report zero counts and the rest
show overdispersed proportions — a pattern common in childcare
enrollment, health service utilization, and ecological count data.

The package uses [CmdStan](https://mc-stan.org/cmdstanr/) as its
computational backend and provides a complete post-estimation workflow
including survey-weighted inference, state-varying coefficients,
marginal effect decomposition, and model comparison.

## Key Features

- **Two-part hurdle framework** — Separates the decision to participate
  (extensive margin, Bernoulli) from the intensity of participation
  (intensive margin, zero-truncated Beta-Binomial).
- **Survey design support** — Pseudo-posterior inference with sampling
  weights, sandwich standard errors, and Cholesky posterior
  recalibration.
- **State-varying coefficients (SVC)** — Hierarchical random intercepts
  and slopes with optional cross-level policy moderators via
  `state_level()`.
- **Marginal effect decomposition** — Average marginal effects
  decomposed into extensive and intensive margin contributions, with
  reversal probability detection.
- **Model comparison** — PSIS-LOO cross-validation with multi-model
  comparison and Pareto-k diagnostics.
- **Posterior predictive checks** — Customizable PPC statistics for
  evaluating model fit.

## Installation

### 1. Install CmdStan (required)

hurdlebb requires [CmdStan](https://mc-stan.org/cmdstanr/) for MCMC
sampling. Install the R interface and the CmdStan toolchain:

``` r
# Install cmdstanr from the Stan R-universe
install.packages("cmdstanr", repos = c(
  "https://stan-dev.r-universe.dev",
  getOption("repos")
))

# Install the CmdStan backend
cmdstanr::install_cmdstan()
```

### 2. Install hurdlebb

``` r
# install.packages("remotes")
remotes::install_github("joonho112/hurdlebb")
```

## Quick Example

``` r
library(hurdlebb)

# Load synthetic NSECE data
data(nsece_synth_small)

# Fit a hurdle Beta-Binomial model
fit <- hbb(
  y | trials(n_trial) ~ poverty + urban,
  data   = nsece_synth_small,
  chains = 4,
  seed   = 42
)

# Model summary
summary(fit)

# Posterior predictive check
ppc_result <- ppc(fit)
plot(ppc_result)

# LOO cross-validation
loo_result <- loo(fit)
print(loo_result)
```

## The Model

The hurdle Beta-Binomial decomposes a bounded count outcome into two
parts:

**Part 1 — Extensive margin (participation):** A Bernoulli model for
whether the unit has a non-zero count. The participation probability *q*
is linked to covariates via a logit: logit(*q*) = *X*’*α*.

**Part 2 — Intensive margin (intensity):** A zero-truncated
Beta-Binomial for the count conditional on being positive. The expected
proportion *μ* is linked to covariates via a logit: logit(*μ*) =
*X*’*β*, with dispersion parameter *κ*.

Both margins share the same design matrix but have separate coefficient
vectors (*α* and *β*), making it straightforward to compare whether a
covariate pushes participation and intensity in the same or opposite
directions.

## Extended Workflow

Beyond basic model fitting, hurdlebb supports a full inferential
pipeline:

``` r
# Survey-weighted model
fit_w <- hbb(
  y | trials(n_trial) ~ poverty + urban,
  data    = nsece_synth_small,
  weights = "weight",
  stratum = "stratum",
  psu     = "psu",
  seed    = 42
)

# Sandwich variance and Cholesky recalibration
sand <- sandwich_variance(fit_w)
chol_fit <- cholesky_correct(fit_w, sand)

# Average marginal effects
ame_result <- ame(fit_w, sand)
ame_decomposition(ame_result)
```

## Vignettes

The package includes five vignettes that walk through the complete
analysis workflow:

| Vignette                                                                                                     | Topic                                                                |
|:-------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------|
| [Getting Started](https://joonho112.github.io/hurdlebb/articles/01-getting-started.md)                       | Data exploration, model fitting, diagnostics, PPC, and LOO           |
| [Survey Design](https://joonho112.github.io/hurdlebb/articles/02-survey-design.md)                           | Survey-weighted inference, sandwich variance, Cholesky recalibration |
| [State-Varying Coefficients](https://joonho112.github.io/hurdlebb/articles/03-state-varying-coefficients.md) | Random intercepts and slopes, shrinkage interpretation               |
| [Policy Moderators](https://joonho112.github.io/hurdlebb/articles/04-policy-moderators.md)                   | Cross-level interactions with `state_level()` terms                  |
| [Marginal Effects](https://joonho112.github.io/hurdlebb/articles/05-marginal-effects.md)                     | AME decomposition, extensive vs. intensive margins                   |

## Datasets

Three synthetic datasets are included, generated via Gaussian copula
from the 2019 National Survey of Early Care and Education (NSECE):

- **`nsece_synth`** — Full dataset (~6,785 providers, 51 states)
- **`nsece_synth_small`** — Stratified subsample (~500 providers) for
  quick demos
- **`nsece_state_policy`** — State-level policy variables (51 rows) for
  cross-level moderators

## Citation

If you use hurdlebb in your research, please cite:

    Lee, J. (2026). The poverty reversal in infant/toddler childcare: a Bayesian
    hurdle beta-binomial model with state-varying coefficients. arXiv preprint.

Or in R:

``` r
citation("hurdlebb")
```

## License

MIT © JoonHo Lee
