for(f in list.files("scripts", full.names = T)) source(f)

p_observed = c(empirical_timing_dist, rep(0, nrow(dat_diff)-length(empirical_timing_dist)))

cases_to_onsets = matrix(0, nrow = nrow(dat_diff), ncol = nrow(dat_diff))
for(i in 1:nrow(dat_diff)) {
  cases_to_onsets[i, i:nrow(dat_diff)]  = p_observed[1:(nrow(dat_diff)-i+1)]
}

cases_to_onsets[1:4, 1:4]
cases_to_onsets[50:54, 50:54]

modeled_state = c("massachusetts", "michigan")
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
  iter = 2000
)

post = rstan::extract(fit)

stan_trace(fit, pars = c("tau"))

post$smoothed_cases %>% apply(c(2,3), mean) %>% .[, 1] %>%  plot
post$smoothed_cases %>% apply(c(2,3), mean) %>% .[, 2] %>%  plot
post$smoothed_onsets %>% apply(c(2,3), mean) %>% .[, 2] %>%  plot
post$log_smoothed_onsets %>% apply(c(2,3), mean) %>% .[, 2] %>%  plot


post$theta %>% apply(c(2,3), mean) %>% .[, 1] %>% plot
post$theta %>% apply(c(2,3), mean) %>% .[, 1] %>% diff %>% plot
post$rt %>% apply(c(2,3), mean) %>% .[, 1] %>% plot
post$rt %>% apply(c(2,3), mean) %>% .[, 2] %>% plot

