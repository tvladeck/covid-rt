for(f in list.files("scripts", full.names = T)) source(f)

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

mod = stan_model("stan_models/rt-v9-just-reporting-error.stan")

fit = sampling(
  mod, 
  stan_data,
  chains = 2, 
  cores = 2,
  iter = 1000
)


log_lik <- loo::extract_log_lik(fit, merge_chains = FALSE)
r_eff <- loo::relative_eff(exp(log_lik))
loo_ll <- loo::loo(log_lik, r_eff = r_eff, cores = 2)

# fine - we 
log_lik[loo::pareto_k_ids(loo_ll, threshold = 1)]
log_lik[loo::pareto_k_ids(loo_ll, threshold = .7)]

print(loo_ll)

post = rstan::extract(fit)

modeled_fit = post$fitted_cases %>% apply(c(2,3), mean) 

plot_modeled_vs_actual = function(modeled, actual) {
  
  modeled = modeled/mean(modeled) * mean(actual)
  
  cbind.data.frame(modeled, actual) %>% 
    mutate(idx = 1:nrow(.)) %>% 
    gather(-idx, key = series, value = estimate) %>% 
    ggplot() + 
    aes(x = idx, y = estimate, color = series) + 
    geom_line()
  
  
}


cbind(, ) %>% ts %>% plot

cbind(modeled_fit[, 2], stan_data$cases[, 2]) %>% ts %>% plot

cbind(modeled_fit[, 3], stan_data$cases[, 3]) %>% ts %>% plot

