for(f in list.files("scripts", full.names = T)) source(f)

stan_mod = stan_model("stan_models/case-curve-base.stan")

fit = sampling(
  stan_mod,
  stan_data,
  cores = 2,
  chains = 2,
  iter = 2000
)

stan_trace(fit, "fitted_cases[10, 1]")

post = rstan::extract(fit)

plot_par_from_posterior("fitted_cases", post, stan_data$cases_raw, date_vector, 1)






