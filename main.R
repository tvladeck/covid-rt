for(f in list.files("scripts", full.names = T)) source(f)

p_observed = c(empirical_timing_dist, rep(0, nrow(dat_diff)-length(empirical_timing_dist)))

timesteps = nrow(dat_diff)
onset_to_cases = matrix(0, nrow = timesteps, ncol = timesteps)
for(i in 1:nrow(dat_diff)) {
  onset_to_cases[i:timesteps, i]  = p_observed[1:(timesteps-i+1)]
}

cases_to_onsets = solve(onset_to_cases)

onset_to_cases[1:4, 1:4]
cases_to_onsets[1:4, 1:4]

obs_cases = dat_diff$michigan
obs_cases %>% ts.plot()

ll_onsets = function(.onsets, tau) {
  .obs_onsets = .onsets * cum_p_observed
  cases = onset_to_cases %*% .obs_onsets
  ll = sum(dnbinom(obs_cases, mu = cases, size = tau, log = T))
  return(ll)
}

find_onsets = function(par) {
  .onsets = par[1:length(obs_cases)]
  tau = par[1+length(obs_cases)]
  ll_onsets(.onsets, tau)
}

found_onsets = 
  optim(rep(1, 55), find_onsets, method = "L-BFGS-B", lower = rep(0.1, 55),
      control = list(fnscale = -1, parscale = rep(10, 55)))

found_cases = (onset_to_cases %*% (found_onsets$par[1:54] * cum_p_observed)) 
ts.plot(cbind(found_cases, obs_cases, found_onsets$par[1:54]), col = c("red", "blue", "green"))


ll_thetasteps_smoothdir = function(initial_onsetload, tau, theta_steps) {
  thetas = cumsum(theta_steps)
  log_delta_onset = rep(0, stan_timesteps)
  log_delta_onset[1] = log(initial_onsetload)
  log_delta_onset[2:stan_timesteps] = thetas
  onsets = exp(cumsum(log_delta_onset)) 
  observed_onsets = onsets * cum_p_observed
  cases = onset_to_cases %*% observed_onsets
  ll_thetasteps = sum(dnorm(theta_steps, mean = 0, sd = tau[1], log = T))
  ll_cases = sum(dnbinom(obs_cases, mu = cases, size = tau[2], log = T))
  ll = ll_thetasteps + ll_cases
  return(ll)
}

find_onsets_smoothdir = function(par) {
  initial_onsetload = par[1]
  tau = rep(0, 2)
  tau[1] = par[2]
  tau[2] = par[3]
  theta_steps = par[4:56]
  ll_thetasteps_smoothdir(initial_onsetload, tau, theta_steps)
}




found_onsets = 
  optim(c(rep(1, 3), rep(0, 53)), find_onsets_smoothdir, method = "L-BFGS-B", 
        lower = c(rep(0.01, 3), rep(-Inf, 53)),
        control = list(fnscale = -1, parscale = rep(1, 56)))

par = found_onsets$par

implied_onsets = function(par) {
  initial_onsetload = par[1]
  tau = rep(0, 2)
  tau[1] = par[2]
  tau[2] = par[3]
  theta_steps = par[4:56]
  
  thetas = cumsum(theta_steps)
  log_delta_onset = rep(0, stan_timesteps)
  log_delta_onset[1] = log(initial_onsetload)
  log_delta_onset[2:stan_timesteps] = thetas
  onsets = exp(cumsum(log_delta_onset)) 
  observed_onsets = onsets * cum_p_observed
  cases = onset_to_cases %*% observed_onsets
  return(
    list(
      onsets = onsets,
      observed_onsets = observed_onsets,
      cases = cases,
      thetas = thetas
    )
  )
}





modeled_state = c("massachusetts", "michigan", "new_york")
col_idx = unlist(map(modeled_state, ~ which(colnames(dat_diff) == .x)))

stan_cases = dat_diff[, col_idx]
stan_tests = tests[, col_idx]
stan_timesteps = nrow(stan_tests)
stan_states = ncol(stan_tests)
date_vector = dat_diff$date

stan_data = list(
  timesteps = stan_timesteps,
  states = stan_states,
  cases = stan_cases,
  tests = stan_tests,
  cum_p_observed = cum_p_observed,
  p_observed = p_observed
)

mod = stan_model("stan_models/rt-v8-reporting-error.stan")

initfn = function(x) {
  list(
    theta_steps = matrix(0, ncol = stan_states, nrow = stan_timesteps-1)
  )
}

fit = optimizing(
  mod, 
  stan_data,
  iter = 1000,
  init = initfn,
  as_vector = F,
  verbose = T
)

fit = sampling(
  mod, 
  stan_data,
  chains = 2, 
  cores = 2,
  iter = 1000,
  init = initfn
)

post = rstan::extract(fit)

stan_trace(fit, pars = c("tau"))
stan_trace(fit, pars = c("theta_steps[2,3]", "initial_onsetload"))

post$lambda %>% apply(c(2,3), mean) %>% .[, 1] %>%  ts.plot()
post$lambda %>% apply(c(2,3), mean) %>% .[, 2] %>%  ts.plot()
post$lambda %>% apply(c(2,3), mean) %>% .[, 3] %>%  ts.plot()
post$lambda_steps %>% apply(c(2,3), mean) %>% .[, 3] %>%  ts.plot()

stan_data$cases[, 1] %>% ts.plot()
stan_data$tests[, 1] %>% ts.plot()

stan_data$cases[, 2] %>% ts.plot()
stan_data$tests[, 2] %>% ts.plot()

stan_data$cases[, 3] %>% ts.plot()
stan_data$tests[, 3] %>% ts.plot()

plot_par_from_posterior("scaled_inferred_onsets", post, stan_data$cases,date_vector, 1)
plot_par_from_posterior("inferred_onsets", post, stan_data$cases,date_vector, 1)
plot_par_from_posterior("rt", post, stan_data$cases,date_vector, 1)
plot_par_from_posterior("lambda", post, stan_data$cases,date_vector, 1)

log(post$inferred_theta) %>% apply(c(2,3), mean) %>% .[, 1] %>%  ts.plot()

post$rt %>% apply(c(2,3), mean) %>% .[, 2] %>%  ts.plot()
post$rt %>% apply(c(2,3), mean) %>% .[, 3] %>%  ts.plot()

impl_theta_1 = 
  log(
    (post$inferred_onsets %>% apply(c(2,3), mean) %>% .[-1, 1] / cum_p_observed[-1]) / 
    (post$inferred_onsets %>% apply(c(2,3), mean) %>% .[-54, 1] / cum_p_observed[-54])
  )

ts.plot(impl_theta_1 * 4 + 1)


impl_theta_2 = 
  log(
    (post$inferred_onsets %>% apply(c(2,3), mean) %>% .[-1, 2] / cum_p_observed[-1]) / 
      (post$inferred_onsets %>% apply(c(2,3), mean) %>% .[-54, 2] / cum_p_observed[-54])
  )

ts.plot(impl_theta_2 * 4 + 1)

impl_theta_3 = 
  log(
    (post$inferred_onsets %>% apply(c(2,3), mean) %>% .[-1, 3] / cum_p_observed[-1]) / 
      (post$inferred_onsets %>% apply(c(2,3), mean) %>% .[-54, 3] / cum_p_observed[-54])
  )

ts.plot(impl_theta_3)


post$theta %>% apply(c(2,3), mean) %>% .[, 1] %>% ts.plot()
post$rt %>% apply(c(2,3), mean) %>% .[, 1] %>% plot
post$rt %>% apply(c(2,3), mean) %>% .[, 2] %>% plot

plot_rt_from_posterior(post, stan_cases, date_vector)
