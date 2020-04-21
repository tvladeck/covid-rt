dat = read_csv(url)
linelist = NCoVUtils::get_linelist()



delay = 
  linelist %>% 
  select(date_onset_symptoms, date_confirmation) %>% 
  na.omit %>% 
  mutate(delta = difftime(date_confirmation, date_onset_symptoms, units = "days")) %>% 
  filter(delta > 0)

empirical_timing_dist = 
  delay$delta %>% 
  ceiling %>% 
  table %>% 
  prop.table 

dat_diff = 
  dat %>% 
  select(date, state, cases) %>% 
  spread(state, cases) %>% 
  setNames(to_snake_case(colnames(.))) %>% 
  filter(date > lubridate::ymd("2020-03-01")) %>% 
  mutate_at(vars(-date), ~ ifelse(is.na(.x), 0, .x)) %>% 
  mutate_at(vars(-date), ~ c(NA, diff(.x))) %>% 
  .[-1, ]

shutdown_grid = 
  shutdown_dates %>% 
  map(~ convert_shutdown_dates_to_date_vector(.x, dat_diff)) %>% 
  reduce(cbind.data.frame) %>% 
  setNames(names(shutdown_dates))

dat_recompiled = 
  dat_diff %>% 
  select(-date) %>% 
  map_df(~ round(apply_1d_filter_rev_pad(empirical_timing_dist, .x))) %>% 
  map(~ c(rep(0, 19), .x)) %>%
  map_df(~ rollsum(.x, 20))
  # identity()

date_vector = dat_diff$date[1:nrow(dat_recompiled)]

shutdown_grid_capped = shutdown_grid[1:nrow(dat_recompiled), ]

dat_recompiled_with_shutdowns = 
  dat_recompiled %>% 
  select(colnames(shutdown_grid_capped))

cum_p_observed = convert_filter_to_cumsum(empirical_timing_dist, nrow(dat_recompiled_with_shutdowns))

stan_data = list(
  timesteps = nrow(dat_recompiled_with_shutdowns),
  states = ncol(dat_recompiled_with_shutdowns),
  cases = dat_recompiled_with_shutdowns,
  shutdowns = shutdown_grid_capped,
  cum_p_observed = cum_p_observed,
  window = 20,
  simulation_steps = 15
)

