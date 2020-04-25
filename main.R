for(f in list.files("scripts", full.names = T)) source(f)

p_observed = c(empirical_timing_dist, rep(0, nrow(dat_diff)-length(empirical_timing_dist)))

modeled_state = "massachusetts"
col_idx = which(colnames(dat_diff) == modeled_state)

stan_cases = dat_diff[, col_idx]
stan_tests = tests[, col_idx]
stan_timesteps = nrow(stan_tests)
stan_states = ncol(stan_tests)

stan_data = list(
  timesteps = stan_timesteps,
  states = stan_states,
  cases = stan_cases,
  tests = stan_tests,
  cum_p_observed = cum_p_observed,
  p_observed = p_observed
)

mod = stan_model("stan_models/rt-v7-reporting-error.stan")

fit = sampling(
  mod, 
  stan_data,
  chains = 2, 
  cores = 2,
  iter = 1000
)

post = rstan::extract(fit)

stan_trace(fit, pars = c("tau", "step_size"))

post$smoothed_cases %>% apply(c(2,3), mean) %>% plot
