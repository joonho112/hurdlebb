# Package index

## Package

- [`hurdlebb`](https://joonho112.github.io/hurdlebb/reference/hurdlebb-package.md)
  [`hurdlebb-package`](https://joonho112.github.io/hurdlebb/reference/hurdlebb-package.md)
  : hurdlebb: Hurdle Beta-Binomial Models for Bounded Discrete
  Proportions

## Model Fitting

Main fitting function and formula specification

- [`hbb()`](https://joonho112.github.io/hurdlebb/reference/hbb.md) : Fit
  a Hurdle Beta-Binomial Model
- [`hbb_formula()`](https://joonho112.github.io/hurdlebb/reference/hbb_formula.md)
  : Parse a Hurdle Beta-Binomial Model Formula
- [`is.hbb_fit()`](https://joonho112.github.io/hurdlebb/reference/is.hbb_fit.md)
  : Test if an Object is an hbb_fit

## Data Preparation

Build and validate Stan-ready data structures

- [`prepare_stan_data()`](https://joonho112.github.io/hurdlebb/reference/prepare_stan_data.md)
  : Prepare Stan Data for a Hurdle Beta-Binomial Model
- [`validate_hbb_data()`](https://joonho112.github.io/hurdlebb/reference/validate_hbb_data.md)
  : Validate an hbb_data Object
- [`validate_hbb_formula()`](https://joonho112.github.io/hurdlebb/reference/validate_hbb_formula.md)
  : Validate an hbb_formula Against a Dataset

## Prior Specification

Customize or inspect model priors

- [`hbb_prior()`](https://joonho112.github.io/hurdlebb/reference/hbb_prior.md)
  : Specify Priors for a Hurdle Beta-Binomial Model
- [`default_prior()`](https://joonho112.github.io/hurdlebb/reference/default_prior.md)
  : Default Prior Specification for Hurdle Beta-Binomial Models

## Distributions

Beta-Binomial, zero-truncated, and hurdle distribution functions

- [`dbetabinom()`](https://joonho112.github.io/hurdlebb/reference/dbetabinom.md)
  : Probability mass function of the Beta-Binomial distribution
- [`pbetabinom()`](https://joonho112.github.io/hurdlebb/reference/pbetabinom.md)
  : Cumulative distribution function of the Beta-Binomial distribution
- [`rbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rbetabinom.md)
  : Random generation from the Beta-Binomial distribution
- [`dztbetabinom()`](https://joonho112.github.io/hurdlebb/reference/dztbetabinom.md)
  : Probability mass function of the Zero-Truncated Beta-Binomial
- [`rztbetabinom()`](https://joonho112.github.io/hurdlebb/reference/rztbetabinom.md)
  : Random generation from the Zero-Truncated Beta-Binomial
- [`dhurdle_betabinom()`](https://joonho112.github.io/hurdlebb/reference/dhurdle_betabinom.md)
  : Probability mass function of the Hurdle Beta-Binomial
- [`rhurdle_betabinom()`](https://joonho112.github.io/hurdlebb/reference/rhurdle_betabinom.md)
  : Random generation from the Hurdle Beta-Binomial
- [`compute_p0()`](https://joonho112.github.io/hurdlebb/reference/compute_p0.md)
  : Compute P(Y = 0) under the Beta-Binomial distribution
- [`compute_ztbb_mean()`](https://joonho112.github.io/hurdlebb/reference/compute_ztbb_mean.md)
  : Conditional mean of the zero-truncated Beta-Binomial
- [`hurdle_mean()`](https://joonho112.github.io/hurdlebb/reference/hurdle_mean.md)
  : Expected value under the Hurdle Beta-Binomial
- [`hurdle_variance()`](https://joonho112.github.io/hurdlebb/reference/hurdle_variance.md)
  : Variance under the Hurdle Beta-Binomial

## Design-Based Inference

Sandwich variance and Cholesky posterior recalibration for survey data

- [`sandwich_variance()`](https://joonho112.github.io/hurdlebb/reference/sandwich_variance.md)
  : Cluster-Robust Sandwich Variance Estimator
- [`cholesky_correct()`](https://joonho112.github.io/hurdlebb/reference/cholesky_correct.md)
  : Cholesky Posterior Recalibration for Survey-Weighted Inference
- [`compute_wald_ci()`](https://joonho112.github.io/hurdlebb/reference/compute_wald_ci.md)
  : Wald Confidence Intervals from Sandwich Variance

## Marginal Effects

Average marginal effect decomposition into extensive and intensive
margins

- [`ame()`](https://joonho112.github.io/hurdlebb/reference/ame.md) :
  Average Marginal Effects Decomposition for Hurdle Beta-Binomial Models
- [`ame_decomposition()`](https://joonho112.github.io/hurdlebb/reference/ame_decomposition.md)
  : Extract the AME Decomposition Table

## Model Comparison & Diagnostics

LOO cross-validation and posterior predictive checks

- [`loo(`*`<hbb_fit>`*`)`](https://joonho112.github.io/hurdlebb/reference/loo.hbb_fit.md)
  : Leave-One-Out Cross-Validation for Hurdle Beta-Binomial Models
- [`hbb_loo_compare()`](https://joonho112.github.io/hurdlebb/reference/hbb_loo_compare.md)
  : Compare Hurdle Beta-Binomial Models via LOO-CV
- [`ppc()`](https://joonho112.github.io/hurdlebb/reference/ppc.md) :
  Posterior Predictive Checks for Hurdle Beta-Binomial Models

## S3 Methods

Standard generic methods for hbb_fit objects

- [`coef(`*`<hbb_fit>`*`)`](https://joonho112.github.io/hurdlebb/reference/coef.hbb_fit.md)
  : Extract Fixed-Effect Coefficients from an hbb_fit
- [`vcov(`*`<hbb_fit>`*`)`](https://joonho112.github.io/hurdlebb/reference/vcov.hbb_fit.md)
  : Variance-Covariance Matrix of Fixed Effects
- [`fitted(`*`<hbb_fit>`*`)`](https://joonho112.github.io/hurdlebb/reference/fitted.hbb_fit.md)
  : Fitted Values from an hbb_fit
- [`residuals(`*`<hbb_fit>`*`)`](https://joonho112.github.io/hurdlebb/reference/residuals.hbb_fit.md)
  : Residuals from an hbb_fit
- [`predict(`*`<hbb_fit>`*`)`](https://joonho112.github.io/hurdlebb/reference/predict.hbb_fit.md)
  : Predictions from a Hurdle Beta-Binomial Model
- [`summary(`*`<hbb_fit>`*`)`](https://joonho112.github.io/hurdlebb/reference/summary.hbb_fit.md)
  : Summary Method for hbb_fit Objects
- [`nobs(`*`<hbb_fit>`*`)`](https://joonho112.github.io/hurdlebb/reference/nobs.hbb_fit.md)
  : Number of Observations in an hbb_fit
- [`plot(`*`<hbb_fit>`*`)`](https://joonho112.github.io/hurdlebb/reference/plot.hbb_fit.md)
  : Plot Diagnostics and Inferential Summaries for hbb_fit Objects
- [`plot(`*`<hbb_ppc>`*`)`](https://joonho112.github.io/hurdlebb/reference/plot.hbb_ppc.md)
  : Plot Method for hbb_ppc Objects
- [`print(`*`<hbb_fit>`*`)`](https://joonho112.github.io/hurdlebb/reference/print.hbb_fit.md)
  : Print Method for hbb_fit Objects
- [`print(`*`<hbb_sandwich>`*`)`](https://joonho112.github.io/hurdlebb/reference/print.hbb_sandwich.md)
  : Print Method for hbb_sandwich Objects
- [`print(`*`<hbb_cholesky>`*`)`](https://joonho112.github.io/hurdlebb/reference/print.hbb_cholesky.md)
  : Print Method for Cholesky-Corrected Posterior Objects
- [`print(`*`<hbb_ame>`*`)`](https://joonho112.github.io/hurdlebb/reference/print.hbb_ame.md)
  : Print Method for hbb_ame Objects
- [`print(`*`<hbb_loo_compare>`*`)`](https://joonho112.github.io/hurdlebb/reference/print.hbb_loo_compare.md)
  : Print Method for hbb_loo_compare Objects
- [`print(`*`<hbb_ppc>`*`)`](https://joonho112.github.io/hurdlebb/reference/print.hbb_ppc.md)
  : Print Method for hbb_ppc Objects
- [`print(`*`<summary.hbb_fit>`*`)`](https://joonho112.github.io/hurdlebb/reference/print.summary.hbb_fit.md)
  : Print Method for summary.hbb_fit Objects

## Stan Model Management

Compile and cache Stan models

- [`hbb_compile()`](https://joonho112.github.io/hurdlebb/reference/hbb_compile.md)
  : Compile Stan Models for hurdlebb
- [`hbb_clear_cache()`](https://joonho112.github.io/hurdlebb/reference/hbb_clear_cache.md)
  : Clear Compiled Stan Model Cache

## Sandwich Internals

Lower-level components for sandwich variance computation

- [`compute_score_matrix()`](https://joonho112.github.io/hurdlebb/reference/compute_score_matrix.md)
  : Extract Posterior Mean Score Matrix from an hbb_fit
- [`compute_H_obs()`](https://joonho112.github.io/hurdlebb/reference/compute_H_obs.md)
  : Compute the Block-Diagonal Observed Information Matrix
- [`compute_J_cluster()`](https://joonho112.github.io/hurdlebb/reference/compute_J_cluster.md)
  : Compute the Cluster-Robust Meat Matrix
- [`compute_der()`](https://joonho112.github.io/hurdlebb/reference/compute_der.md)
  : Compute Design Effect Ratios from a Sandwich Variance Object

## Datasets

Built-in synthetic datasets based on NSECE 2019

- [`nsece_synth`](https://joonho112.github.io/hurdlebb/reference/nsece_synth.md)
  : Synthetic NSECE 2019 Center-Based Provider Data
- [`nsece_synth_small`](https://joonho112.github.io/hurdlebb/reference/nsece_synth_small.md)
  : Small Synthetic NSECE 2019 Dataset
- [`nsece_state_policy`](https://joonho112.github.io/hurdlebb/reference/nsece_state_policy.md)
  : State-Level Childcare Policy Variables
