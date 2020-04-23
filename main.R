for(f in list.files("scripts", full.names = T)) source(f)

# stan_mod = stan_model("stan_models/impact-of-shutdown.stan")
# stan_mod = stan_model("stan_models/rt.stan")
stan_mod = stan_model("stan_models/rt-v3.stan")

fit = sampling(
  stan_mod, 
  stan_data,
  chains = 2,
  cores = 2,
  iter = 3000,
  warmup = 2000
)

post = rstan::extract(fit)



stan_trace(fit, "step_size")
print(fit, "mean_step_size")


step_size = rstan::extract(fit, permuted =F)[, 1, 2] 
a = get_sampler_params(fit)
energy = map(a, ~ .x[4001:8000, "energy__"]) %>% 
  reduce(c)

plot(cbind(step_size, energy))


newrt = summarize_rt_from_posterior(post, stan_data$cases, date_vector)

plot_rt_from_posterior(post, stan_data$cases, date_vector)

plot_state_rt_from_posterior("california", post, stan_data$cases, date_vector)

print(fit, "mean_step_size")
stan_trace(fit, "step_size")
stan_trace(fit, "rt[11,10]")

stan_trace(fit, "shutdown_impact_on_rt[12]")
print(fit, "shutdown_impact_on_rt")