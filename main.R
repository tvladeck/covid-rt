for(f in list.files("scripts", full.names = T)) source(f)

# stan_mod = stan_model("stan_models/impact-of-shutdown.stan")
# stan_mod = stan_model("stan_models/rt.stan")
stan_mod = stan_model("stan_models/rt-v3.stan")

fit = sampling(
  stan_mod, 
  stan_data,
  chains = 2,
  cores = 2,
  iter = 1000
)

post = rstan::extract(fit)



newrt = summarize_rt_from_posterior(post, stan_data$cases, date_vector)

plot_rt_from_posterior(post, stan_data$cases, date_vector)


stan_trace(fit, "shutdown_impact_on_rt[12]")
print(fit, "shutdown_impact_on_rt")