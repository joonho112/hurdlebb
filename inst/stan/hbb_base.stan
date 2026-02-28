// =============================================================================
// hbb_base.stan — Pooled Hurdle Beta-Binomial (no hierarchy, no weights)
// =============================================================================
// Package: hurdlebb
// Source:  Consolidated from M0 (hbb_m0.stan)
//
// Two-part hurdle model for bounded discrete proportions (y / n_trial):
//   Part 1 (Extensive margin): Bernoulli — does provider participate?
//     P(z_i = 1) = q_i = logit^{-1}(X_i alpha)
//   Part 2 (Intensive margin):  Zero-truncated Beta-Binomial — share among participants
//     Y_i | z_i=1 ~ ZT-BetaBin(n_i, mu_i, kappa)
//     mu_i = logit^{-1}(X_i beta),  kappa = exp(log_kappa)
//     a_i = mu_i * kappa,  b_i = (1 - mu_i) * kappa
//
// Hurdle log-likelihood:
//   z_i = 0:  ell_i = log(1 - q_i)
//   z_i = 1:  ell_i = log(q_i) + log f_BB(y_i | n_i, a_i, b_i) - log(1 - p0_i)
//   where p0_i = P(Y=0 | BB) = B(b_i + n_i, kappa) / B(b_i, kappa)
//
// Configurable priors via data block (no recompilation needed):
//   alpha ~ Normal(prior_alpha_mean, prior_alpha_sd)
//   beta  ~ Normal(prior_beta_mean, prior_beta_sd)
//   log_kappa ~ Normal(prior_kappa_mean, prior_kappa_sd)
//
// Generated quantities:
//   log_lik[N] — pointwise log-likelihood for LOO-CV
//   y_rep[N]   — posterior predictive replications (rejection sampling for ZT-BB)
// =============================================================================

data {
  int<lower=1> N;                    // number of observations
  int<lower=1> P;                    // number of covariates (including intercept)
  array[N] int<lower=0> y;           // outcome count (e.g., IT enrollment)
  array[N] int<lower=1> n_trial;     // number of trials (e.g., total 0-5 enrollment)
  array[N] int<lower=0, upper=1> z;  // participation indicator: z[i] = 1 iff y[i] > 0
  matrix[N, P] X;                    // design matrix (col 1 = intercept)

  // --- Configurable prior hyperparameters ---
  real prior_alpha_mean;              // mean for alpha ~ Normal(., .)
  real prior_beta_mean;               // mean for beta  ~ Normal(., .)
  real<lower=0> prior_alpha_sd;      // SD for alpha ~ Normal(., .)
  real<lower=0> prior_beta_sd;       // SD for beta  ~ Normal(., .)
  real prior_kappa_mean;             // mean for log_kappa ~ Normal(., .)
  real<lower=0> prior_kappa_sd;      // SD for log_kappa ~ Normal(., .)
}

parameters {
  vector[P] alpha;       // extensive-margin coefficients (logit scale)
  vector[P] beta;        // intensive-margin coefficients (logit scale)
  real log_kappa;        // log concentration parameter
}

transformed parameters {
  real<lower=0> kappa = exp(log_kappa);
}

model {
  // --- Priors (configurable via data) ---
  alpha ~ normal(prior_alpha_mean, prior_alpha_sd);
  beta ~ normal(prior_beta_mean, prior_beta_sd);
  log_kappa ~ normal(prior_kappa_mean, prior_kappa_sd);

  // --- Precompute linear predictors (vectorized) ---
  //   eta_ext[i] = X[i] * alpha
  //   eta_int[i] = X[i] * beta
  vector[N] eta_ext = X * alpha;
  vector[N] eta_int = X * beta;

  // --- Hurdle likelihood ---
  // TODO: consider reduce_sum for within-chain parallelism
  for (i in 1:N) {
    real q_i = inv_logit(eta_ext[i]);       // P(participate)
    real mu_i = inv_logit(eta_int[i]);      // conditional mean share
    real a_i = mu_i * kappa;                // BB shape alpha = mu * kappa
    real b_i = (1 - mu_i) * kappa;          // BB shape beta  = (1-mu) * kappa

    if (z[i] == 0) {
      // --- Structural zero: not participating ---
      target += log1m(q_i);

    } else {
      // --- Positive count: log(q_i) + log_fBB - log(1 - p0) ---

      // log P(Y=0 | BetaBin) via lgamma for numerical stability:
      //   p0 = B(b+n, kappa) / B(b, kappa)
      //      = Gamma(b+n) * Gamma(kappa) / [Gamma(b) * Gamma(kappa+n)]
      real log_p0_i = lgamma(b_i + n_trial[i]) + lgamma(kappa)
                    - lgamma(b_i) - lgamma(kappa + n_trial[i]);

      // log(1 - p0): log1m_exp computes log(1 - exp(x)) for x <= 0
      real log_1mp0_i = log1m_exp(log_p0_i);

      // Beta-Binomial log-PMF:
      //   log f_BB = lchoose(n, y) + lbeta(y+a, n-y+b) - lbeta(a, b)
      real log_fBB_i = lchoose(n_trial[i], y[i])
                     + lbeta(y[i] + a_i, n_trial[i] - y[i] + b_i)
                     - lbeta(a_i, b_i);

      // Hurdle contribution
      target += log(q_i) + log_fBB_i - log_1mp0_i;
    }
  }
}

generated quantities {
  // --- Pointwise log-likelihood (for LOO-CV) ---
  vector[N] log_lik;

  // --- Posterior predictive replications ---
  array[N] int<lower=0> y_rep;

  {
    vector[N] eta_ext = X * alpha;
    vector[N] eta_int = X * beta;

    for (i in 1:N) {
      real q_i = inv_logit(eta_ext[i]);
      real mu_i = inv_logit(eta_int[i]);
      real a_i = mu_i * kappa;
      real b_i = (1 - mu_i) * kappa;

      if (z[i] == 0) {
        log_lik[i] = log1m(q_i);
      } else {
        real log_p0_i = lgamma(b_i + n_trial[i]) + lgamma(kappa)
                      - lgamma(b_i) - lgamma(kappa + n_trial[i]);
        real log_1mp0_i = log1m_exp(log_p0_i);
        real log_fBB_i = lchoose(n_trial[i], y[i])
                       + lbeta(y[i] + a_i, n_trial[i] - y[i] + b_i)
                       - lbeta(a_i, b_i);
        log_lik[i] = log(q_i) + log_fBB_i - log_1mp0_i;
      }

      // --- Posterior predictive replication ---
      // Step 1: Draw participation from Bernoulli(q_i)
      int z_rep_i = bernoulli_rng(q_i);

      if (z_rep_i == 0) {
        y_rep[i] = 0;
      } else {
        // Step 2: Rejection sampling for zero-truncated Beta-Binomial
        //   Draw from BetaBin(n, a, b) via beta_rng + binomial_rng; reject y=0
        int draw = 0;
        int max_iter = 1000;   // safety cap to prevent infinite loops
        int iter = 0;
        while (draw == 0 && iter < max_iter) {
          real p_draw = beta_rng(a_i, b_i);
          draw = binomial_rng(n_trial[i], p_draw);
          iter += 1;
        }
        // If rejection sampling exhausted (extremely unlikely), set to 1
        y_rep[i] = (draw == 0) ? 1 : draw;
      }
    }
  }
}
