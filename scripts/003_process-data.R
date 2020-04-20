dat = read_csv(url)


dat_multivar = 
  dat %>% 
  filter(state %in% STATES) %>% 
  select(date, state, cases) %>% 
  spread(state, cases) %>% 
  setNames(to_snake_case(colnames(.))) %>% 
  filter(date > lubridate::ymd("2020-03-01")) %>% 
  mutate_at(vars(-date), ~ ifelse(is.na(.x), 0, .x)) %>% 
  mutate_at(vars(-date), function(x) {
    diff(x) %>% 
      {. + 1} %>% 
      {c(rep(0, WINDOW), .)} %>% 
      rollsum(., WINDOW) 
  }) %>% 
  .[-1, ]


shutdown_grid = 
  shutdown_dates %>% 
  map(~ convert_shutdown_dates_to_date_vector(.x)) %>% 
  reduce(cbind.data.frame) %>% 
  setNames(names(shutdown_dates))


dat_multivar_with_shutdowns = 
  dat_multivar %>% 
  select(names(shutdown_grid))


stan_data = list(
  timesteps = nrow(dat_multivar_with_shutdowns),
  states = ncol(dat_multivar_with_shutdowns),
  cases = dat_multivar_with_shutdowns,
  shutdowns = shutdown_grid
)