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
