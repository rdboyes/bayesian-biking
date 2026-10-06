// Hierarchical lagged kinetic model of heart rate driven by power, across rides.
// Within ride r, latent HR relaxes toward
//   a_ref[r] + Bt[r] * elapsed + b[r] * (P - P_ref)
// with time constant tau[r]; residuals follow an AR(1) process.
// Ride-level (a_ref, log b, log tau, Bt) are correlated draws around population
// means, with a bike offset on a_ref and log b. Rides are split across threads
// with reduce_sum.
//
// Sampling notes: every parameter is declared on unit scale and shifted/scaled
// in transformed parameters, and the baseline is anchored at P_ref rather than
// 0 W (HR at 0 W is an extrapolation and strongly anti-correlated with b).
// Ride noise levels are centred (every ride pins its own noise down tightly,
// and a non-centred form mixed badly). Ride effects stay non-centred: short
// rides say little about their own drift or tau, and centring them caused
// divergences.
functions {
  real partial_sum(array[] int rides, int start, int end,
                   array[] int first, array[] int len,
                   vector dt, vector P, vector hr, vector elapsed, real P_ref,
                   vector a_ref, vector b, vector tau, vector Bt, vector mu0,
                   vector sigma, real phi) {
    real lp = 0;
    real sd_scale = 1 / sqrt(1 - square(phi));
    for (i in 1:size(rides)) {
      int r = rides[i];
      int s = first[r];
      int n = len[r];
      vector[n] e;
      real mu = mu0[r];
      e[1] = hr[s] - mu;
      for (j in 2:n) {
        int t = s + j - 1;
        mu += (1 - exp(-dt[t] / tau[r]))
              * (a_ref[r] + Bt[r] * elapsed[t] + b[r] * (P[t] - P_ref) - mu);
        e[j] = hr[t] - mu;
      }
      lp += normal_lpdf(e[1] | 0, sigma[r] * sd_scale)
            + normal_lpdf(e[2:n] | phi * e[1:(n - 1)], sigma[r]);
    }
    return lp;
  }
}
data {
  int<lower=1> R;                    // rides
  int<lower=1> N;                    // total observations
  array[R] int<lower=1> first;       // index of each ride's first observation
  array[R] int<lower=2> len;         // observations per ride
  vector<lower=0>[N] dt;             // seconds since previous observation
  vector<lower=0>[N] P;              // estimated power, in units of 100 W
  vector[N] hr;                      // heart rate, bpm
  vector<lower=0>[N] elapsed;        // hours since ride start
  vector<lower=0, upper=1>[R] trek;  // 1 if ridden on the Trek 520
  real<lower=0> P_ref;               // reference power for a_ref, units of 100 W
  int<lower=1> grainsize;
}
transformed data {
  array[R] int ride_idx;
  vector[R] hr_first;
  // Prior centres and scales for (a_ref, log b, log tau, Bt)
  vector[4] pop_loc = [150, log(25), log(60), 0]';
  vector[4] pop_scale = [20, 0.5, 0.7, 10]';
  vector[4] sd_scale = [15, 0.5, 0.7, 5]';
  vector[2] trek_scale = [10, 0.3]';
  for (r in 1:R) {
    ride_idx[r] = r;
    hr_first[r] = hr[first[r]];
  }
}
parameters {
  vector[4] pop_raw;                 // population means, unit scale
  vector<lower=0>[4] pop_sd_raw;     // between-ride sds, unit scale
  cholesky_factor_corr[4] L;         // between-ride correlation
  matrix[4, R] z;                    // ride effects (non-centred)
  vector[2] trek_raw;                // Trek offset on a_ref and log b
  vector[R] mu0_raw;                 // latent HR at each ride's start
  real log_sigma_raw;
  real<lower=0> sd_log_sigma;
  vector[R] log_sigma;               // ride noise levels (centred)
  real<lower=-1, upper=1> phi;
}
transformed parameters {
  vector[4] pop_mu = pop_loc + pop_scale .* pop_raw;
  vector[4] pop_sd = sd_scale .* pop_sd_raw;
  vector[2] beta_trek = trek_scale .* trek_raw;
  real mu_log_sigma = log(3) + 0.5 * log_sigma_raw;
  matrix[R, 4] theta = rep_matrix(pop_mu', R)
                       + (diag_pre_multiply(pop_sd, L) * z)';
  vector[R] a_ref = theta[:, 1] + beta_trek[1] * trek;   // HR at P_ref, bpm
  vector[R] b = exp(theta[:, 2] + beta_trek[2] * trek);  // bpm per 100 W
  vector[R] tau = exp(theta[:, 3]);                      // s
  vector[R] Bt = theta[:, 4];                            // bpm per hour
  vector[R] a = a_ref - b * P_ref;                       // HR at 0 W, bpm
  vector[R] mu0 = hr_first + 15 * mu0_raw;
  vector[R] sigma = exp(log_sigma);
}
model {
  pop_raw ~ std_normal();
  pop_sd_raw ~ std_normal();
  L ~ lkj_corr_cholesky(2);
  to_vector(z) ~ std_normal();
  trek_raw ~ std_normal();
  mu0_raw ~ std_normal();
  log_sigma_raw ~ std_normal();
  sd_log_sigma ~ normal(0, 0.5);
  log_sigma ~ normal(mu_log_sigma, sd_log_sigma);
  phi ~ normal(0.5, 0.3);

  target += reduce_sum(partial_sum, ride_idx, grainsize,
                       first, len, dt, P, hr, elapsed, P_ref,
                       a_ref, b, tau, Bt, mu0, sigma, phi);
}
generated quantities {
  matrix[4, 4] Omega = multiply_lower_tri_self_transpose(L);
  real pop_b = exp(pop_mu[2]);       // typical bpm per 100 W (Miele)
  real pop_tau = exp(pop_mu[3]);     // typical time constant, s
}
