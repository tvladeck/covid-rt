library(foreach)

windows = floor(nrow(dat_diff)/7)

theta_mean = 0.3352583
theta_sd = 0.3032382

rolling_sim_mass = 
  foreach(i = 1:(nrow(dat_diff)-7), .combine = rbind.data.frame) %do% {
  
  if(i > 1) {
    theta_mean = mean(post$theta[, 6, 1])
    theta_sd   = sd(post$theta[, 6, 1])
  }  
    
  shutdown_grid = 
    shutdown_dates %>% 
    map(~ convert_shutdown_dates_to_date_vector(.x, dat_diff, 7)) %>% 
    reduce(cbind.data.frame) %>% 
    setNames(names(shutdown_dates)) %>% 
    select(massachusetts)
  
  dat_recompiled = 
    dat_diff %>% 
    select(-date) %>% 
    map_df(~ round(apply_1d_filter_rev_pad(empirical_timing_dist, .x))) %>% 
    # map(~ c(rep(0, 19), .x)) %>%
    # map_df(~ rollsum(.x, 20))
    identity()
  
  
  
  shutdown_grid_capped = shutdown_grid %>% 
    slice(1:nrow(dat_recompiled))
  
  dat_recompiled_with_shutdowns = 
    dat_recompiled %>% 
    select(colnames(shutdown_grid_capped)) %>% 
    slice(i:(i+6)) %>%
    identity()
  
  date_vector = dat_diff %>% 
    slice(i:(i+6)) %>% 
    pull(date)
  
  cum_p_observed = convert_filter_to_cumsum(empirical_timing_dist, nrow(dat_recompiled_with_shutdowns))
  
  stan_data = list(
    timesteps = nrow(dat_recompiled_with_shutdowns),
    states = ncol(dat_recompiled_with_shutdowns),
    cases = dat_recompiled_with_shutdowns,
    shutdowns = shutdown_grid_capped,
    cum_p_observed = cum_p_observed,
    window = 20,
    simulation_steps = 15,
    theta_mean = theta_mean,
    theta_sd = theta_sd
  )
  
  fit = sampling(
    stan_mod, 
    stan_data,
    chains = 2,
    cores = 2,
    iter = 3000,
    warmup = 2000
  )
  
  post = rstan::extract(fit)
  
  newrt = summarize_rt_from_posterior(post, stan_data$cases, date_vector)
  
  newrt
}


rolling_sim_mass %>% 
  group_by(date,state) %>% 
  summarize_at(vars(mean:lower), mean) %>% 
  gather(-date, -state, key = series, value = rt) %>% 
  arrange(date) %>% 
  ggplot() + 
    aes(x = date, y = rt, color = series, lty = series) + 
    geom_point(size = .1) + 
    geom_line() + 
    scale_color_manual("", values = c("mean" = "red", "lower" = "grey", "upper" = "grey")) + 
    theme_bw() + 
    scale_linetype_manual("", values = c("mean" = 1, "lower" = 2, "upper" = 2)) + 
    theme(legend.position = "none") + 
    labs(y = to_title_case("rt"), x = "") + 
    geom_hline(yintercept = 1) + 
    facet_wrap(~ state)

