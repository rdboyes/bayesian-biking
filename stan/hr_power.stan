// Lagged kinetic model of heart rate driven by power.
// Latent HR relaxes toward a steady state a + Bt * elapsed + b * P with time
// constant tau; Bt captures slow upward drift in baseline HR over the ride.
// Observation errors follow an AR(1) process.
data {
  int<lower=2> N;
  vector<lower=0>[N] dt;       // seconds since previous observation
  vector<lower=0>[N] P;        // estimated power, in units of 100 W
  vector[N] hr;                // heart rate, bpm
  vector<lower=0>[N] elapsed;  // hours since ride start
}
parameters {
  real a;                          // steady-state HR at 0 W, ride start
  real<lower=0> b;                 // bpm per 100 W
  real<lower=0> tau;               // response time constant, s
  real mu0;                        // latent HR at first observation
  real<lower=0> sigma;             // innovation sd
  real<lower=-1, upper=1> phi;     // AR(1) coefficient of residuals
  real Bt;                         // drift in baseline HR, bpm per hour
}
transformed parameters {
  vector[N] mu;
  mu[1] = mu0;
  for (t in 2:N) {
    real target_hr = a + Bt * elapsed[t] + b * P[t];
    mu[t] = mu[t - 1] + (1 - exp(-dt[t] / tau)) * (target_hr - mu[t - 1]);
  }
}
model {
  vector[N] e = hr - mu;
  real sd1 = sigma / sqrt(1 - square(phi));

  a ~ normal(90, 20);
  b ~ normal(30, 15);
  tau ~ lognormal(log(45), 0.5);
  mu0 ~ normal(hr[1], 10);
  sigma ~ normal(0, 5);
  phi ~ normal(0.5, 0.3);
  Bt ~ normal(0, 10);

  e[1] ~ normal(0, sd1);
  e[2:N] ~ normal(phi * e[1:(N - 1)], sigma);
}
generated quantities {
  vector[N] hr_rep;
  vector[N] log_lik;
  {
    vector[N] e = hr - mu;
    real sd1 = sigma / sqrt(1 - square(phi));
    real e_rep = normal_rng(0, sd1);
    hr_rep[1] = mu[1] + e_rep;
    log_lik[1] = normal_lpdf(e[1] | 0, sd1);
    for (t in 2:N) {
      e_rep = normal_rng(phi * e_rep, sigma);
      hr_rep[t] = mu[t] + e_rep;
      log_lik[t] = normal_lpdf(e[t] | phi * e[t - 1], sigma);
    }
  }
}
