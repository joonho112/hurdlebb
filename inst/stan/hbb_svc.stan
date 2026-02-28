// =============================================================================
// hbb_svc.stan — State-Varying Coefficients + Optional Policy Moderators (unweighted)
// =============================================================================
// Package: hurdlebb
// Source:  Consolidated from M3b (hbb_m3b.stan)
//
// Two-part hurdle model for bounded discrete proportions (y / n_trial):
//   Part 1 (Extensive margin): Bernoulli — does provider participate?
//   Part 2 (Intensive margin):  Zero-truncated Beta-Binomial — share among participants
//
// Hierarchical structure:
//   State-varying coefficients (SVC) with policy moderators that EXPLAIN
//   variation in the random effects across groups.
//
//   delta[s] = Gamma * v_state[s]' + epsilon[s]
//
//   where:
//     delta[s]       = K-vector of group-specific deviations (K = 2P)
//                      head(delta[s], P) = extensive-margin deviations
//                      tail(delta[s], P) = intensive-margin deviations
//     Gamma          = K x Q matrix of policy moderator coefficients
//     v_state[s]     = Q-vector of group-level covariates (col 1 = intercept)
//     Gamma * v[s]'  = policy-explained component of group variation
//     epsilon[s]     ~ N(0, Sigma_epsilon) = residual group variation
//     Sigma_epsilon  = diag(tau) * Omega * diag(tau)
//
//   Non-centered parameterization (NCP) for residuals:
//     z_eps[s] ~ N(0, I_K)   iid for s = 1,...,S
//     epsilon[s] = L_Sigma * z_eps[s]
//     L_Sigma = diag_pre_multiply(tau, L_Omega)
//
// Linear predictors:
//   eta_ext[i] = X[i] * (alpha + head(delta[state[i]], P))
//   eta_int[i] = X[i] * (beta  + tail(delta[state[i]], P))
//
// Identifiability note:
//   v_state column 1 = all 1s (group intercept). Gamma[,1] absorbs the mean
//   of delta across groups. The residual epsilon[s] is truly mean-zero.
//   alpha and beta remain identifiable as individual-level fixed effects.
//
// Configurable priors via data block (no recompilation needed):
//   alpha     ~ Normal(prior_alpha_mean, prior_alpha_sd)
//   beta      ~ Normal(prior_beta_mean, prior_beta_sd)
//   log_kappa ~ Normal(prior_kappa_mean, prior_kappa_sd)
//   vec(Gamma) ~ Normal(0, prior_gamma_sd)
//   tau       ~ Normal+(0, prior_tau_sd)
//   L_Omega   ~ LKJ_cholesky(prior_lkj_eta)
//
// Generated quantities:
//   log_lik[N]       — pointwise log-likelihood for LOO-CV
//   y_rep[N]         — posterior predictive replications
//   Omega[K,K]       — residual correlation matrix (from L_Omega)
// =============================================================================

data {
  int<lower=1> N;                    // number of observations
  int<lower=1> P;                    // number of individual-level covariates (incl. intercept)
  int<lower=1> S;                    // number of groups (e.g., states)
  int<lower=1> Q;                    // number of group-level covariates (incl. intercept)
  array[N] int<lower=0> y;           // outcome count
  array[N] int<lower=1> n_trial;     // number of trials
  array[N] int<lower=0, upper=1> z;  // participation indicator: z[i] = 1 iff y[i] > 0
  matrix[N, P] X;                    // individual-level design matrix (col 1 = intercept)
  array[N] int<lower=1, upper=S> state;  // group index for each observation
  matrix[S, Q] v_state;             // group-level design matrix (col 1 = intercept)

  // --- Configurable prior hyperparameters ---
  real prior_alpha_mean;              // mean for alpha ~ Normal(., .)
  real prior_beta_mean;               // mean for beta  ~ Normal(., .)
  real<lower=0> prior_alpha_sd;      // SD for alpha ~ Normal(., .)
  real<lower=0> prior_beta_sd;       // SD for beta  ~ Normal(., .)
  real prior_kappa_mean;             // mean for log_kappa ~ Normal(., .)
  real<lower=0> prior_kappa_sd;      // SD for log_kappa ~ Normal(., .)
  real<lower=0> prior_gamma_sd;      // SD for vec(Gamma) ~ Normal(0, .)
  real<lower=0> prior_tau_sd;        // SD for tau ~ Normal+(0, .)
  real<lower=0> prior_lkj_eta;       // LKJ concentration for L_Omega
}

transformed data {
  int K = 2 * P;  // joint random effect dimension (extensive + intensive)
}

parameters {
  vector[P] alpha;                   // extensive-margin fixed effects (logit scale)
  vector[P] beta;                    // intensive-margin fixed effects (logit scale)
  real log_kappa;                    // log concentration parameter

  // --- Policy moderator coefficients ---
  // Gamma[k, q] = effect of group covariate q on random effect dimension k
  //   Rows 1:P     -> extensive margin
  //   Rows (P+1):K -> intensive margin
  matrix[K, Q] Gamma;

  // --- Residual random effects (NCP) ---
  vector<lower=0>[K] tau;            // residual scales: tau[1:P] = ext, tau[(P+1):K] = int
  cholesky_factor_corr[K] L_Omega;   // Cholesky factor of K x K residual correlation matrix
  array[S] vector[K] z_eps;          // raw NCP variates (K-dimensional per group)
}

transformed parameters {
  real<lower=0> kappa = exp(log_kappa);

  // --- Recover state-varying coefficients ---
  // delta[s] = Gamma * v_state[s]' + epsilon[s]
  //   Gamma * v[s]' = policy-explained component (K-vector)
  //   epsilon[s]    = residual group variation (via NCP)
  array[S] vector[K] delta;
  {
    // Form L_Sigma ONCE outside the state loop
    matrix[K, K] L_Sigma = diag_pre_multiply(tau, L_Omega);
    for (s in 1:S) {
      // epsilon[s] = L_Sigma * z_eps[s]  (residual, mean-zero)
      vector[K] epsilon_s = L_Sigma * z_eps[s];
      // delta[s] = policy-explained + residual
      // v_state[s] is row_vector[Q]; transpose to column vector for Gamma multiply
      delta[s] = Gamma * v_state[s]' + epsilon_s;
    }
  }
}

model {
  // --- Priors: fixed effects (configurable via data) ---
  alpha ~ normal(prior_alpha_mean, prior_alpha_sd);
  beta ~ normal(prior_beta_mean, prior_beta_sd);
  log_kappa ~ normal(prior_kappa_mean, prior_kappa_sd);

  // --- Prior: policy moderator coefficients ---
  to_vector(Gamma) ~ normal(0, prior_gamma_sd);

  // --- Priors: residual random effects ---
  tau ~ normal(0, prior_tau_sd);           // half-normal (lower=0 enforced by constraint)
  L_Omega ~ lkj_corr_cholesky(prior_lkj_eta);  // LKJ prior on K x K residual correlation
  for (s in 1:S)
    z_eps[s] ~ std_normal();               // iid N(0, I_K) for NCP

  // --- Precompute fixed-effect linear predictors (vectorized) ---
  vector[N] eta_ext_fixed = X * alpha;
  vector[N] eta_int_fixed = X * beta;

  // --- Hurdle likelihood ---
  // TODO: consider reduce_sum for within-chain parallelism
  for (i in 1:N) {
    // Extract margin-specific deviations from joint delta vector
    //   head(delta[s], P) = delta_ext (first P elements)
    //   tail(delta[s], P) = delta_int (last P elements)
    real eta_ext_i = eta_ext_fixed[i] + X[i] * head(delta[state[i]], P);
    real eta_int_i = eta_int_fixed[i] + X[i] * tail(delta[state[i]], P);

    real q_i = inv_logit(eta_ext_i);       // P(participate)
    real mu_i = inv_logit(eta_int_i);      // conditional mean share
    real a_i = mu_i * kappa;               // BB shape alpha
    real b_i = (1 - mu_i) * kappa;         // BB shape beta

    if (z[i] == 0) {
      // --- Structural zero: not participating ---
      target += log1m(q_i);

    } else {
      // --- Positive count: log(q_i) + log_fBB - log(1 - p0) ---

      // log P(Y=0 | BetaBin) via lgamma
      real log_p0_i = lgamma(b_i + n_trial[i]) + lgamma(kappa)
                    - lgamma(b_i) - lgamma(kappa + n_trial[i]);

      // log(1 - p0): log1m_exp for numerical stability
      real log_1mp0_i = log1m_exp(log_p0_i);

      // Beta-Binomial log-PMF
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

  // --- FULL K x K residual correlation matrix ---
  // Omega has block structure:
  //   Omega[1:P, 1:P]           = within-extensive residual correlations
  //   Omega[(P+1):K, (P+1):K]   = within-intensive residual correlations
  //   Omega[1:P, (P+1):K]       = cross-margin residual correlations
  // NOTE: This is the RESIDUAL correlation (after removing policy effects).
  corr_matrix[K] Omega = multiply_lower_tri_self_transpose(L_Omega);

  {
    vector[N] eta_ext_fixed = X * alpha;
    vector[N] eta_int_fixed = X * beta;

    for (i in 1:N) {
      // Extract margin-specific deviations from joint delta vector
      real eta_ext_i = eta_ext_fixed[i] + X[i] * head(delta[state[i]], P);
      real eta_int_i = eta_int_fixed[i] + X[i] * tail(delta[state[i]], P);

      real q_i = inv_logit(eta_ext_i);
      real mu_i = inv_logit(eta_int_i);
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
