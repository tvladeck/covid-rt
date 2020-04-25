for(f in list.files("scripts", full.names = T)) source(f)



p_observed = c(empirical_timing_dist, rep(0, nrow(dat_diff)-length(empirical_timing_dist)))

onsets_to_cases = matrix(0, nrow = nrow(dat_diff), ncol = nrow(dat_diff))
for(i in 1:nrow(dat_diff)) {
  onsets_to_cases[i:nrow(dat_diff), i] = p_observed[1:(nrow(dat_diff)-(i-1))]
}

cases_to_onsets = solve(onsets_to_cases)

observed_cases = dat_diff$california

ll_inferred_onsets = function(.inferred_onsets) {
  
  implied_cases = onsets_to_cases %*% .inferred_onsets
  
  error = sum((dat_diff$california-implied_cases)^2)
  
  return(error)
}

ml_inferred_onsets = optim(rep(0, 53), ll_inferred_onsets,
                            control = list(parscale = rep(100, 53)),
                            lower = rep(0, 53), method = "L-BFGS-B")


implied_cases_from_scaled_inferred_onsets = function(.inferred_onsets) {
  
  implied_cases = onsets_to_cases %*% .inferred_onsets
  
  return(implied_cases)
  
}

ts.plot(
  cbind(
    implied_cases_from_scaled_inferred_onsets(ml_inferred_onsets$par),
    dat_diff$california
  ),
  col = c("red", "blue")
)




stan_mod_generative = stan_model("stan_models/rt-v6-generative.stan")


data_generative = list(
  cases = dat_diff[, 6],
  states = 1,
  timesteps = nrow(dat_diff),
  p_observed = p_observed,
  cum_p_observed = cum_p_observed
)

fit_generative = sampling(
  stan_mod_generative,
  data_generative,
  chains = 2,
  iter = 5000,
  cores = 2
)

post_generative = rstan::extract(fit_generative)

ts.plot(post_generative$expected_cases %>% apply(c(2,3), mean), ylab = "expected_cases")
ts.plot(post_generative$onsets %>% apply(c(2,3), mean), ylab = "expected_onsets")

stan_mod_generative_2 = stan_model("stan_models/rt-v6-generative-v2.stan")

fit_generative_2 = sampling(
  stan_mod_generative_2,
  data_generative,
  chains = 2,
  iter = 5000,
  cores = 2
)

stan_trace(fit_generative_2, pars = c("step_size", "theta[1,1]", "initial_seed"))
print(fit_generative_2, pars = c("step_size", "tau"))
stan_trace(fit_generative_2, pars = c("step_size", "tau"))

post_generative_2 = rstan::extract(fit_generative_2)

ts.plot(post_generative_2$expected_cases %>% apply(c(2,3), median), ylab = "expected_cases")
ts.plot(post_generative_2$onsets %>% apply(c(2,3), median), ylab = "expected_onsets")
ts.plot(post_generative_2$theta %>% apply(c(2,3), median), ylab = "theta")
ts.plot(post_generative_2$rt %>% apply(c(2,3), median), ylab = "rt")

onsets = post_generative_2$onsets %>% apply(c(2,3), median)
theta_ratio = log(onsets[-53]/onsets[-1])
theta_sample =  post_generative_2$theta %>% apply(c(2,3), mean)


log_lik_generative_2 = loo::extract_log_lik(fit_generative_2, merge_chains = F)
reff_generative_2 = loo::relative_eff(exp(log_lik_generative_2))


loo_generative_2 <- loo::loo(log_lik_generative_2, r_eff = reff_generative_2, cores = 2) 

log_lik_generative = loo::extract_log_lik(fit_generative, merge_chains = F)
reff_generative = loo::relative_eff(exp(log_lik_generative))


loo_generative <- loo::loo(log_lik_generative, r_eff = reff_generative, cores = 2) 




stan_mod_base = stan_model("stan_models/rt-v3.stan")
stan_mod_shutdown = stan_model("stan_models/rt-v4-plus-shutdown.stan")

ITER = 10e3
WARMUP = 8e3

fit_base = sampling(
  stan_mod_base, 
  stan_data,
  chains = 2,
  cores = 2,
  iter = ITER,
  warmup = WARMUP
)

fit_shutdown = sampling(
  stan_mod_shutdown, 
  stan_data,
  chains = 2,
  cores = 2,
  iter = ITER,
  warmup = WARMUP
)

post_base = rstan::extract(fit_base)
post_shutdown = rstan::extract(fit_shutdown)

log_lik_base = loo::extract_log_lik(fit_base, merge_chains = F)
log_lik_shutdown = loo::extract_log_lik(fit_shutdown, merge_chains = F)
reff_base = loo::relative_eff(exp(log_lik_base))
reff_shutdown = loo::relative_eff(exp(log_lik_shutdown))

# reff_base[which(is.na(reff_base))] = mean(reff_base, na.rm = T)
# reff_shutdown[which(is.na(reff_shutdown))] = mean(reff_shutdown, na.rm = T)


loo_base <- loo(log_lik_base, r_eff = reff_base, cores = 2)
loo_shutdown <- loo(log_lik_shutdown, r_eff = reff_shutdown, cores = 2)

loo_base
loo_shutdown

loo::pareto_k_ids(loo_base, threshold = 0.5)

a = loo::loo_compare(loo_base, loo_shutdown)

wts <- loo::loo_model_weights(list(loo_base, loo_shutdown), method = "stacking")

print(a)

stan_trace(fit_shutdown, "shutdown_impact_on_rt")
print(fit_shutdown, "shutdown_impact_on_rt")

plot_rt_from_posterior(post_shutdown, stan_data$cases, date_vector)



print(fit, "step_size")
print(fit, "shutdown_impact_on_rt")
print(fit, "tau")
stan_trace(fit, "rho")
stan_trace(fit, "shutdown_impact_on_rt")
stan_trace(fit, "tau")
stan_trace(fit, "step_size")
stan_trace(fit, "theta[20,1]")
print(fit, c("step_size", "theta[20,1]"))


step_size = rstan::extract(fit, permuted =F)[, 1, 2] 
a = get_sampler_params(fit)
energy = map(a, ~ .x[2001:3000, "energy__"]) %>% 
  reduce(c)

plot(cbind(step_size, energy))

plot_par_from_posterior("rt", post, stan_data$cases, date_vector, 1)

newrt = summarize_rt_from_posterior(post, stan_data$cases, date_vector)

plot_rt_from_posterior(post, stan_data$cases, date_vector)

plot_state_rt_from_posterior("california", post, stan_data$cases, date_vector)

print(fit, "mean_step_size")
stan_trace(fit, "step_size")
stan_trace(fit, "rt[11,10]")

stan_trace(fit, "shutdown_impact_on_rt[12]")
print(fit, "shutdown_impact_on_rt")