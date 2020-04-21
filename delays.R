# remotes::install_github("epiforecasts/NCoVUtils", dependencies = TRUE)



stan_mod = stan_model("stan_models/time-to-test.stan")

idk = sampling(stan_mod, list(n = nrow(delay), timing = as.numeric(delay$delta)), cores = 2)

pars = rstan::extract(idk)

empirical_timing_dist = 
  delay$delta %>% 
  ceiling %>% 
  table %>% 
  prop.table 

empirical_timing_dist %>% cumsum

empirical_timing_dist_cap = empirical_timing_dist[1:18]

ser = rev(c(rep(0, 17), dat_multivar_with_shutdowns$alabama))
apply_1d_filter(empirical_timing_dist_cap, ser) %>% sum

apply_1d_filter_rev_pad(empirical_timing_dist_cap, dat_multivar_with_shutdowns$alabama) %>% plot
  
