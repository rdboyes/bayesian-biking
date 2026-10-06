library(cmdstanr)
library(posterior)
library(dplyr)
library(ggplot2)

# Fit the hierarchical model across all rides from 03_read_all.R.
# Usage: Rscript 04_fit_all.R [iter_warmup] [iter_sampling]

args <- as.integer(commandArgs(trailingOnly = TRUE))
iter_warmup   <- if (length(args) >= 1) args[1] else 1000
iter_sampling <- if (length(args) >= 2) args[2] else 1000
chains <- 4
threads_per_chain <- 2

d <- readRDS(file.path("output", "all_rides.rds"))
rides <- d$info |>
  filter(keep) |>
  arrange(date) |>
  mutate(r = row_number())
bins <- d$bins |>
  inner_join(select(rides, ride_id, r), by = "ride_id") |>
  arrange(r, t)

idx <- bins |>
  mutate(row = row_number()) |>
  group_by(r) |>
  summarise(first = min(row), len = n(), .groups = "drop")

stan_data <- list(
  R = nrow(rides),
  N = nrow(bins),
  first = idx$first,
  len = idx$len,
  dt = bins$dt,
  P = bins$P / 100,
  hr = bins$hr,
  elapsed = bins$t / 3600,
  trek = as.numeric(rides$bike == "Trek 520"),
  P_ref = 1.5,  # 150 W
  grainsize = 1
)

# Start near the prior centres (all unit-scale parameters near 0) with a little
# jitter between chains, rather than Stan's default uniform(-2, 2) inits.
init_fn <- function() {
  list(
    pop_raw = rnorm(4, 0, 0.2),
    pop_sd_raw = runif(4, 0.5, 1),
    L = diag(4),
    z = matrix(rnorm(4 * stan_data$R, 0, 0.1), 4, stan_data$R),
    trek_raw = c(0, 0),
    mu0_raw = rep(0, stan_data$R),
    log_sigma_raw = 0,
    sd_log_sigma = 0.2,
    log_sigma = rep(log(3), stan_data$R),
    phi = 0.5
  )
}

mod <- cmdstan_model("stan/hr_power_hier.stan",
                     cpp_options = list(stan_threads = TRUE))
fit <- mod$sample(
  data = stan_data, init = init_fn, seed = 1,
  chains = chains, parallel_chains = chains,
  threads_per_chain = threads_per_chain,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
  refresh = max(1, (iter_warmup + iter_sampling) %/% 10)
)

pop_pars <- c("pop_mu", "pop_sd", "beta_trek", "pop_b", "pop_tau",
              "mu_log_sigma", "sd_log_sigma", "phi", "Omega[1,2]", "Omega[1,3]",
              "Omega[1,4]", "Omega[2,3]", "Omega[2,4]", "Omega[3,4]")
print(fit$summary(pop_pars), n = 50, width = 120)
print(fit$diagnostic_summary())

ride_sum <- fit$summary(c("a", "a_ref", "b", "tau", "Bt", "sigma"))
cat(sprintf("Max R-hat: %.3f | min bulk ESS: %.0f\n",
            max(ride_sum$rhat, na.rm = TRUE), min(ride_sum$ess_bulk, na.rm = TRUE)))

fit$save_object(file.path("output", "fit_all.rds"))

# ---- ride-level estimates over time -------------------------------------------
ride_est <- ride_sum |>
  mutate(param = sub("\\[.*", "", variable),
         r = as.integer(sub(".*\\[(\\d+)\\]", "\\1", variable))) |>
  filter(param %in% c("a", "b", "tau", "Bt")) |>
  inner_join(rides, by = "r") |>
  mutate(param = recode(param,
                        a   = "a: HR at 0 W (bpm)",
                        b   = "b: bpm per 100 W",
                        tau = "tau: time constant (s)",
                        Bt  = "Bt: drift (bpm/hour)"))
write.csv(select(ride_est, ride_id, date, bike, param, mean, median, q5, q95),
          file.path("output", "fit_all_rides.csv"), row.names = FALSE)

p_time <- ggplot(ride_est, aes(date, median, colour = bike)) +
  geom_linerange(aes(ymin = q5, ymax = q95), alpha = 0.4) +
  geom_point(size = 1) +
  facet_wrap(~param, scales = "free_y", ncol = 1) +
  labs(x = NULL, y = NULL, colour = NULL,
       subtitle = "Per-ride posterior median and 90% interval") +
  theme_minimal() +
  theme(legend.position = "top")
ggsave(file.path("output", "fit_all_over_time.png"), p_time,
       width = 9, height = 9, dpi = 120)

# ---- observed vs fitted latent HR for a sample of rides ------------------------
med <- function(v) ride_sum$median[match(sprintf("%s[%d]", v, seq_len(nrow(rides))), ride_sum$variable)]
mu0_med <- fit$summary("mu0", "median")$median
a_m <- med("a"); b_m <- med("b"); tau_m <- med("tau"); Bt_m <- med("Bt")

latent <- function(r) {
  s <- idx$first[r]; n <- idx$len[r]; rows <- s:(s + n - 1)
  mu <- numeric(n); mu[1] <- mu0_med[r]
  for (j in 2:n) {
    t <- rows[j]
    target <- a_m[r] + Bt_m[r] * stan_data$elapsed[t] + b_m[r] * stan_data$P[t]
    mu[j] <- mu[j - 1] + (1 - exp(-stan_data$dt[t] / tau_m[r])) * (target - mu[j - 1])
  }
  data.frame(r = r, t = bins$t[rows], hr = bins$hr[rows], mu = mu)
}
show <- unique(round(seq(1, nrow(rides), length.out = 6)))
traj <- bind_rows(lapply(show, latent)) |>
  inner_join(select(rides, r, date, bike), by = "r") |>
  mutate(label = paste(date, bike))
p_traj <- ggplot(traj, aes(t / 60)) +
  geom_line(aes(y = hr), linewidth = 0.3) +
  geom_line(aes(y = mu), colour = "steelblue", linewidth = 0.5) +
  facet_wrap(~label, scales = "free_x", ncol = 2) +
  labs(x = "Elapsed time (min)", y = "Heart rate (bpm)",
       subtitle = "black: observed; blue: latent HR at posterior medians") +
  theme_minimal()
ggsave(file.path("output", "fit_all_trajectories.png"), p_traj,
       width = 10, height = 8, dpi = 120)
