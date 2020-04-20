for(f in list.files("scripts", full.names = T)) source(f)

stan_mod = stan_model("rt.stan")

fit = sampling(
  stan_mod, 
  stan_data,
  chains = 2,
  cores = 2,
  iter = 4000
)

post = rstan::extract(fit)

summarize_rt_from_posterior(post)

plot_rt_from_posterior(post)





