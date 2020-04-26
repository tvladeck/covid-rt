for(f in list.files("scripts", full.names = T)) source(f)

p_observed = c(empirical_timing_dist, rep(0, nrow(dat_diff)-length(empirical_timing_dist)))

cases_to_onsets = matrix(0, nrow = nrow(dat_diff), ncol = nrow(dat_diff))
for(i in 1:nrow(dat_diff)) {
  cases_to_onsets[i, i:nrow(dat_diff)]  = p_observed[1:(nrow(dat_diff)-i+1)]
}

cases_to_onsets[1:4, 1:4]
cases_to_onsets[50:54, 50:54]

onsets_to_cases = solve(cases_to_onsets)

onsets_to_cases[1:4, 1:4]
onsets_to_cases[51:54, 51:54]

onsets = 1:54
onsets_to_cases %*% onsets

onsets = rnorm(54, 0, .1)
onsets_to_cases %*% onsets

impl_cases_1 = onsets_to_cases %*% (post$smoothed_onsets %>% apply(c(2,3), mean) %>% .[, 1])


modeled_state = c("massachusetts", "michigan")
col_idx = which(colnames(dat_diff) == modeled_state)

col_idx = c(23, 24, 34)

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

mod = stan_model("stan_models/rt-v8-reporting-error.stan")

fit = sampling(
  mod, 
  stan_data,
  chains = 2, 
  cores = 2,
  iter = 2000
)

post = rstan::extract(fit)

stan_trace(fit, pars = c("tau"))
stan_trace(fit, pars = c("theta_steps[2,3]"))

post$lambda %>% apply(c(2,3), mean) %>% .[, 1] %>%  ts.plot()
post$lambda %>% apply(c(2,3), mean) %>% .[, 2] %>%  ts.plot()
post$lambda %>% apply(c(2,3), mean) %>% .[, 3] %>%  ts.plot()
post$lambda_steps %>% apply(c(2,3), mean) %>% .[, 3] %>%  ts.plot()

stan_data$cases[, 1] %>% ts.plot()
stan_data$tests[, 1] %>% ts.plot()
post$inferred_onsets %>% apply(c(2,3), mean) %>% .[, 1] %>%  ts.plot()

stan_data$cases[, 2] %>% ts.plot()
post$inferred_onsets %>% apply(c(2,3), mean) %>% .[, 2] %>%  ts.plot()

stan_data$cases[, 3] %>% ts.plot()
stan_data$tests[, 3] %>% ts.plot()
post$smoothed_cases %>% apply(c(2,3), mean) %>% .[, 3] %>%  ts.plot()

post$smoothed_onsets %>% apply(c(2,3), mean) %>% .[, 1] %>%  ts.plot()
post$smoothed_onsets %>% apply(c(2,3), mean) %>% .[, 2] %>%  ts.plot()
post$smoothed_onsets %>% apply(c(2,3), mean) %>% .[, 3] %>%  ts.plot()

impl_theta_1 = 
  log(
    (post$smoothed_onsets %>% apply(c(2,3), mean) %>% .[-1, 1] / cum_p_observed[-1]) / 
    (post$smoothed_onsets %>% apply(c(2,3), mean) %>% .[-54, 1] / cum_p_observed[-54])
  )

ts.plot(impl_theta_1 * 4 + 1)


impl_theta_2 = 
  log(
    (post$smoothed_onsets %>% apply(c(2,3), mean) %>% .[-1, 2] / cum_p_observed[-1]) / 
      (post$smoothed_onsets %>% apply(c(2,3), mean) %>% .[-54, 2] / cum_p_observed[-54])
  )

ts.plot(impl_theta_2 * 4 + 1)


post$theta %>% apply(c(2,3), mean) %>% .[, 1] %>% ts.plot()
post$rt %>% apply(c(2,3), mean) %>% .[, 1] %>% plot
post$rt %>% apply(c(2,3), mean) %>% .[, 2] %>% plot

plot_rt_from_posterior(post, stan_cases, date_vector)
