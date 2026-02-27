// =============================================================================
// hbb_svc_weighted.stan — SVC + Survey Weights + Scores + Policy Moderators
// =============================================================================
// Package: hurdlebb
// Source:  Consolidated from M3b-W (hbb_m3b_weighted.stan)
//
// Two-part hurdle model for bounded discrete proportions (y / n_trial):
//   Part 1 (Extensive margin): Bernoulli — does provider participate?
//   Part 2 (Intensive margin):  Zero-truncated Beta-Binomial — share among participants
//
// Hierarchical structure:
//   State-varying coefficients (SVC) with policy moderators.
//   delta[s] = Gamma * v_state[s]' + epsilon[s]
//     head(delta[s], P) = extensive-margin deviations
//     tail(delta[s], P) = intensive-margin deviations
//     Sigma_epsilon = diag(tau) * Omega * diag(tau)
//     NCP: epsilon[s] = L_Sigma * z_eps[s],  z_eps[s] ~ N(0, I_K)
//
// Linear predictors (include delta terms):
//   eta_ext[i] = X[i] * (alpha + head(delta[state[i]], P))
//   eta_int[i] = X[i] * (beta  + tail(delta[state[i]], P))
//
// KEY DIFFERENCE FROM hbb_svc.stan:
//   hbb_svc:          target += ell_i                  (unweighted)
//   hbb_svc_weighted: target += w_tilde[i] * ell_i     (weighted pseudo-posterior)
//
// Weights:
//   w_tilde[i] = w[i] * N / sum(w)   (normalized to sum to N)
//   Pseudo-posterior on same "effective sample size" scale (Savitsky & Toth, 2016).
//
// Configurable priors via data block (no recompilation needed):
//   alpha     ~ Normal(0, prior_alpha_sd)
//   beta      ~ Normal(0, prior_beta_sd)
//   log_kappa ~ Normal(prior_kappa_mean, prior_kappa_sd)
//   vec(Gamma) ~ Normal(0, prior_gamma_sd)
//   tau       ~ Normal+(0, prior_tau_sd)
//   L_Omega   ~ LKJ_cholesky(prior_lkj_eta)
//
// Generated quantities:
//   log_lik[N]        — UNWEIGHTED pointwise log-likelihood for LOO-CV
//   y_rep[N]          — posterior predictive replications
//   Omega[K,K]        — residual correlation matrix
//   score_ext[N,P]    — extensive-margin score: d ell_i / d alpha
//   score_int[N,P]    — intensive-margin score: d ell_i / d beta
//   score_kappa[N]    — dispersion score: d ell_i / d log_kappa
//
// NOTE: Score vectors are UNWEIGHTED individual scores.
//   Weight multiplication and cluster aggregation happen in R for the
//   sandwich variance estimator: V_sand = H^{-1} J_cluster H^{-1}
//
// Score derivations (SVC model — eta includes delta[state[i]] terms):
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
  int<lower=1> P;                    // number of individual-level covariates (incl. intercept)
  int<lower=1> S;                    // number of groups (e.g., states)
  int<lower=1> Q;                    // number of group-level covariates (incl. intercept)
  array[N] int<lower=0> y;           // outcome count
  array[N] int<lower=1> n_trial;     // number of trials
  array[N] int<lower=0, upper=1> z;  // participation indicator: z[i] = 1 iff y[i] > 0
  matrix[N, P] X;                    // individual-level design matrix (col 1 = intercept)
  array[N] int<lower=1, upper=S> state;  // group index for each observation
  matrix[S, Q] v_state;             // group-level design matrix (col 1 = intercept)

  // --- Survey weights ---
  vector<lower=0>[N] w_tilde;       // normalized survey weights: sum(w_tilde) = N

  // --- Configurable prior hyperparameters ---
  real<lower=0> prior_alpha_sd;      // SD for alpha ~ Normal(0, .)
  real<lower=0> prior_beta_sd;       // SD for beta  ~ Normal(0, .)
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
  cholesky_factor_corr[K] L_Omega;   // Cholesky factor of K x K residual correlation
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
      delta[s] = Gamma * v_state[s]' + epsilon_s;
    }
  }
}

model {
  // --- Priors: fixed effects (configurable via data) ---
  alpha ~ normal(0, prior_alpha_sd);
  beta ~ normal(0, prior_beta_sd);
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

  // --- WEIGHTED Likelihood (pseudo-posterior) ---
  //   target += w_tilde[i] * ell_i
  // TODO: consider reduce_sum for within-chain parallelism
  for (i in 1:N) {
    // Linear predictors include delta[state[i]] terms
    real eta_ext_i = eta_ext_fixed[i] + X[i] * head(delta[state[i]], P);
    real eta_int_i = eta_int_fixed[i] + X[i] * tail(delta[state[i]], P);

    real q_i = inv_logit(eta_ext_i);       // P(participate)
    real mu_i = inv_logit(eta_int_i);      // conditional mean share
    real a_i = mu_i * kappa;               // BB shape alpha
    real b_i = (1 - mu_i) * kappa;         // BB shape beta

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

  // --- FULL K x K residual correlation matrix ---
  // Omega has block structure:
  //   Omega[1:P, 1:P]           = within-extensive residual correlations
  //   Omega[(P+1):K, (P+1):K]   = within-intensive residual correlations
  //   Omega[1:P, (P+1):K]       = cross-margin residual correlations
  // NOTE: This is the RESIDUAL correlation (after removing policy effects).
  corr_matrix[K] Omega = multiply_lower_tri_self_transpose(L_Omega);

  // --- Score vectors for sandwich variance (UNWEIGHTED individual scores) ---
  // s_i = d log f(y_i | theta) / d theta
  // Weight multiplication and cluster aggregation happen in R.
  matrix[N, P] score_ext;           // d ell_i / d alpha
  matrix[N, P] score_int;           // d ell_i / d beta
  vector[N] score_kappa;            // d ell_i / d log_kappa

  {
    vector[N] eta_ext_fixed = X * alpha;
    vector[N] eta_int_fixed = X * beta;

    for (i in 1:N) {
      // Linear predictors include delta[state[i]] terms
      real eta_ext_i = eta_ext_fixed[i] + X[i] * head(delta[state[i]], P);
      real eta_int_i = eta_int_fixed[i] + X[i] * tail(delta[state[i]], P);

      real q_i = inv_logit(eta_ext_i);
      real mu_i = inv_logit(eta_int_i);
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
