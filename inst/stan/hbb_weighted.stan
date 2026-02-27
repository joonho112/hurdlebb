// =============================================================================
// hbb_weighted.stan — Pooled Hurdle Beta-Binomial + Survey Weights + Scores
// =============================================================================
// Package: hurdlebb
// Source:  New combination — pooled (no random effects) with survey weights
//          and score vectors for sandwich variance estimation
//
// Two-part hurdle model for bounded discrete proportions (y / n_trial):
//   Part 1 (Extensive margin): Bernoulli — does provider participate?
//     P(z_i = 1) = q_i = logit^{-1}(X_i alpha)
//   Part 2 (Intensive margin):  Zero-truncated Beta-Binomial — share among participants
//     Y_i | z_i=1 ~ ZT-BetaBin(n_i, mu_i, kappa)
//     mu_i = logit^{-1}(X_i beta),  kappa = exp(log_kappa)
//     a_i = mu_i * kappa,  b_i = (1 - mu_i) * kappa
//
// KEY DIFFERENCE FROM hbb_base.stan:
//   hbb_base:     target += ell_i                     (unweighted)
//   hbb_weighted: target += w_tilde[i] * ell_i        (weighted pseudo-posterior)
//
// Weights:
//   w_tilde[i] = w[i] * N / sum(w)   (normalized to sum to N)
//   This ensures the pseudo-posterior is on the same "effective sample size"
//   scale as the unweighted posterior (Savitsky & Toth, 2016).
//
// Configurable priors via data block (no recompilation needed):
//   alpha ~ Normal(0, prior_alpha_sd)
//   beta  ~ Normal(0, prior_beta_sd)
//   log_kappa ~ Normal(prior_kappa_mean, prior_kappa_sd)
//
// Generated quantities:
//   log_lik[N]        — UNWEIGHTED pointwise log-likelihood for LOO-CV
//   y_rep[N]          — posterior predictive replications
//   score_ext[N, P]   — extensive-margin score: d ell_i / d alpha
//   score_int[N, P]   — intensive-margin score: d ell_i / d beta
//   score_kappa[N]    — dispersion score: d ell_i / d log_kappa
//
// NOTE: Score vectors are UNWEIGHTED individual scores.
//   Weight multiplication and cluster aggregation happen in R for the
//   sandwich variance estimator: V_sand = H^{-1} J_cluster H^{-1}
//
// Score derivations (pooled model — no delta terms):
//   eta_ext_i = X_i * alpha   (no random effects)
//   eta_int_i = X_i * beta    (no random effects)
//
//   Extensive:
//     score_ext[i] = (z_i - q_i) * X_i   (Bernoulli chain rule)
//
//   Intensive (z_i = 1 only; zero for z_i = 0):
//     S_BB_mu = kappa * [psi(y+a) - psi(n-y+b) - psi(a) + psi(b)]
//     Lambda_i = kappa * [psi(b+n) - psi(b)]
//     trunc_corr_mu = p0 * Lambda / (1 - p0)
//     score_mu = S_BB_mu - trunc_corr_mu              (MINUS sign)
//     score_int[i] = score_mu * mu(1-mu) * X_i
//
//   Dispersion (z_i = 1 only; zero for z_i = 0):
//     S_BB_kappa = kappa * [mu(psi(y+a)-psi(a)) + (1-mu)(psi(n-y+b)-psi(b))
//                           + psi(kappa) - psi(n+kappa)]
//     trunc_sum = sum_{j=1}^{n-1} j / [(b+j)(kappa+j)]
//     trunc_corr_kappa = p0 * mu * kappa * trunc_sum / (1 - p0)
//     score_kappa = S_BB_kappa - trunc_corr_kappa     (MINUS sign)
// =============================================================================

data {
  int<lower=1> N;                    // number of observations
  int<lower=1> P;                    // number of covariates (including intercept)
  array[N] int<lower=0> y;           // outcome count
  array[N] int<lower=1> n_trial;     // number of trials
  array[N] int<lower=0, upper=1> z;  // participation indicator: z[i] = 1 iff y[i] > 0
  matrix[N, P] X;                    // design matrix (col 1 = intercept)

  // --- Survey weights ---
  vector<lower=0>[N] w_tilde;       // normalized survey weights: sum(w_tilde) = N

  // --- Configurable prior hyperparameters ---
  real<lower=0> prior_alpha_sd;      // SD for alpha ~ Normal(0, .)
  real<lower=0> prior_beta_sd;       // SD for beta  ~ Normal(0, .)
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
  alpha ~ normal(0, prior_alpha_sd);
  beta ~ normal(0, prior_beta_sd);
  log_kappa ~ normal(prior_kappa_mean, prior_kappa_sd);

  // --- Precompute linear predictors (vectorized, pooled — no delta) ---
  //   eta_ext[i] = X[i] * alpha
  //   eta_int[i] = X[i] * beta
  vector[N] eta_ext = X * alpha;
  vector[N] eta_int = X * beta;

  // --- WEIGHTED Likelihood (pseudo-posterior) ---
  //   target += w_tilde[i] * ell_i
  // TODO: consider reduce_sum for within-chain parallelism
  for (i in 1:N) {
    real q_i = inv_logit(eta_ext[i]);       // P(participate)
    real mu_i = inv_logit(eta_int[i]);      // conditional mean share
    real a_i = mu_i * kappa;                // BB shape alpha
    real b_i = (1 - mu_i) * kappa;          // BB shape beta

    if (z[i] == 0) {
      // --- Structural zero: weighted ---
      target += w_tilde[i] * log1m(q_i);

    } else {
      // --- Positive count: weighted ---
      real log_p0_i = lgamma(b_i + n_trial[i]) + lgamma(kappa)
                    - lgamma(b_i) - lgamma(kappa + n_trial[i]);
      real log_1mp0_i = log1m_exp(log_p0_i);
      real log_fBB_i = lchoose(n_trial[i], y[i])
                     + lbeta(y[i] + a_i, n_trial[i] - y[i] + b_i)
                     - lbeta(a_i, b_i);

      target += w_tilde[i] * (log(q_i) + log_fBB_i - log_1mp0_i);
    }
  }
}

generated quantities {
  // --- Pointwise log-likelihood (UNWEIGHTED for LOO-CV) ---
  vector[N] log_lik;

  // --- Posterior predictive replications ---
  array[N] int<lower=0> y_rep;

  // --- Score vectors for sandwich variance (UNWEIGHTED individual scores) ---
  matrix[N, P] score_ext;           // d ell_i / d alpha
  matrix[N, P] score_int;           // d ell_i / d beta
  vector[N] score_kappa;            // d ell_i / d log_kappa

  {
    // Precompute linear predictors (pooled — no delta)
    vector[N] eta_ext = X * alpha;
    vector[N] eta_int = X * beta;

    for (i in 1:N) {
      real q_i = inv_logit(eta_ext[i]);
      real mu_i = inv_logit(eta_int[i]);
      real a_i = mu_i * kappa;
      real b_i = (1 - mu_i) * kappa;

      // =========================================================
      // Extensive-margin score: s_i^ext = (z_i - q_i) * X_i
      //   Chain rule: d ell_i / d alpha_p = (d ell_i / d eta_ext) * X_ip
      //   For Bernoulli hurdle: d ell_i / d eta_ext = z_i - q_i
      //   (applies to both z=0 and z=1 cases)
      // =========================================================
      score_ext[i] = (z[i] - q_i) * to_row_vector(X[i]);

      if (z[i] == 0) {
        // --- Structural zero ---
        log_lik[i] = log1m(q_i);
        score_int[i] = rep_row_vector(0, P);
        score_kappa[i] = 0;

      } else {
        // --- Positive count ---
        real log_p0_i = lgamma(b_i + n_trial[i]) + lgamma(kappa)
                      - lgamma(b_i) - lgamma(kappa + n_trial[i]);
        real log_1mp0_i = log1m_exp(log_p0_i);
        real log_fBB_i = lchoose(n_trial[i], y[i])
                       + lbeta(y[i] + a_i, n_trial[i] - y[i] + b_i)
                       - lbeta(a_i, b_i);
        real p0_i = exp(log_p0_i);

        log_lik[i] = log(q_i) + log_fBB_i - log_1mp0_i;

        // =========================================================
        // Intensive-margin score: s_i^int (CORRECTED MINUS sign)
        //   S_BB_mu = kappa [psi(y+a) - psi(n-y+b) - psi(a) + psi(b)]
        //   Lambda_i = kappa [psi(b+n) - psi(b)]
        //   trunc_corr = p0 * Lambda / (1-p0)
        //   score_mu = S_BB_mu - trunc_corr    <-- MINUS sign
        //   score_int = score_mu * mu(1-mu) * X_i
        //
        // Chain rule: d ell_i / d beta_p = (d ell_i / d mu) * mu(1-mu) * X_ip
        //   where mu(1-mu) is the logistic derivative d mu / d eta_int
        // =========================================================
        real S_BB_mu = kappa * (digamma(y[i] + a_i)
                               - digamma(n_trial[i] - y[i] + b_i)
                               - digamma(a_i) + digamma(b_i));
        real Lambda_i = kappa * (digamma(b_i + n_trial[i]) - digamma(b_i));
        real trunc_corr_mu = p0_i * Lambda_i / (1 - p0_i);
        real score_mu_i = S_BB_mu - trunc_corr_mu;  // MINUS sign

        score_int[i] = score_mu_i * mu_i * (1 - mu_i) * to_row_vector(X[i]);

        // =========================================================
        // Dispersion score: s_i^kappa (w.r.t. log_kappa, MINUS sign)
        //   S_BB_kappa = kappa[mu(psi(y+a)-psi(a)) + (1-mu)(psi(n-y+b)-psi(b))
        //                + psi(kappa) - psi(n+kappa)]
        //   trunc_corr_kappa = p0*mu*kappa/(1-p0) * sum_{j=1}^{n-1} j/[(b+j)(kappa+j)]
        //   score_kappa = S_BB_kappa - trunc_corr_kappa
        //
        // Note: the kappa multiplier at front accounts for d/d(log kappa)
        //   via chain rule: d/d(log kappa) = kappa * d/d(kappa)
        // =========================================================
        real S_BB_kappa = kappa * (
          mu_i * (digamma(y[i] + a_i) - digamma(a_i))
          + (1 - mu_i) * (digamma(n_trial[i] - y[i] + b_i) - digamma(b_i))
          + digamma(kappa) - digamma(n_trial[i] + kappa)
        );

        real trunc_sum = 0;
        for (j in 1:(n_trial[i] - 1)) {
          trunc_sum += 1.0 * j / ((b_i + j) * (kappa + j));
        }
        real trunc_corr_kappa = p0_i * mu_i * kappa * trunc_sum / (1 - p0_i);

        score_kappa[i] = S_BB_kappa - trunc_corr_kappa;  // MINUS sign
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
