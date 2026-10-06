library(cmdstanr)
library(posterior)
library(dplyr)
library(ggplot2)

ride_id <- "5835787118"
d <- readRDS(file.path("output", paste0("ride_", ride_id, ".rds")))
bins <- d$bins

stan_data <- list(
  N  = nrow(bins),
  dt = bins$dt,
  P  = bins$P / 100,
  hr = bins$hr
)

mod <- cmdstan_model("stan/hr_power.stan")
fit <- mod$sample(data = stan_data, chains = 4, parallel_chains = 4,
                  seed = 1, refresh = 500)

pars <- c("a", "b", "tau", "mu0", "sigma", "phi")
print(fit$summary(pars), width = 120)
print(fit$diagnostic_summary())
fit$save_object(file.path("output", paste0("fit_", ride_id, ".rds")))

# ---- (a) posteriors of the key parameters ------------------------------------
dr <- as_draws_df(fit$draws(c("a", "b", "tau")))
post <- bind_rows(
  data.frame(param = "a: HR at 0 W (bpm)",       value = dr$a),
  data.frame(param = "b: bpm per 100 W",         value = dr$b),
  data.frame(param = "tau: time constant (s)",   value = dr$tau)
)
p_post <- ggplot(post, aes(value)) +
  geom_density(fill = "grey70", colour = NA) +
  facet_wrap(~param, scales = "free") +
  labs(x = NULL, y = NULL) +
  theme_minimal()
ggsave(file.path("output", paste0("fit_", ride_id, "_posteriors.png")), p_post,
       width = 9, height = 3, dpi = 120)

# ---- (b) observed HR vs latent trajectory and posterior predictive -----------
mu_mean <- fit$summary("mu", "mean")$mean
rep_q <- fit$summary("hr_rep", ~quantile(.x, c(0.05, 0.95)))
traj <- bins |>
  mutate(mu = mu_mean, lo = rep_q$`5%`, hi = rep_q$`95%`)
p_traj <- ggplot(traj, aes(t / 60)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = hr), linewidth = 0.3) +
  geom_line(aes(y = mu), colour = "steelblue", linewidth = 0.5) +
  labs(x = "Elapsed time (min)", y = "Heart rate (bpm)",
       subtitle = "black: observed; blue: latent HR (mean) with 90% predictive band") +
  theme_minimal()
ggsave(file.path("output", paste0("fit_", ride_id, "_trajectory.png")), p_traj,
       width = 9, height = 4, dpi = 120)

# ---- (c) steady-state HR vs power ---------------------------------------------
P_grid <- seq(0, quantile(bins$P, 0.99), length.out = 50)
ss <- sapply(P_grid, function(p) dr$a + dr$b * p / 100)
ss_df <- data.frame(
  P   = P_grid,
  mid = apply(ss, 2, median),
  lo  = apply(ss, 2, quantile, 0.05),
  hi  = apply(ss, 2, quantile, 0.95)
)
p_ss <- ggplot(ss_df, aes(P)) +
  geom_point(data = bins, aes(P, hr), alpha = 0.15, size = 0.8) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = mid), colour = "steelblue") +
  labs(x = "Estimated power (W, 10 s bins)", y = "Heart rate (bpm)",
       subtitle = "Steady-state HR = a + b * P (median, 90% interval); points are raw bins") +
  theme_minimal()
ggsave(file.path("output", paste0("fit_", ride_id, "_steady_state.png")), p_ss,
       width = 6, height = 4, dpi = 120)
