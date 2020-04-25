dat = read_csv(url)
linelist = NCoVUtils::get_linelist()

state_abbrev = 
  read_csv("https://raw.githubusercontent.com/jasonong/List-of-US-States/master/states.csv") %>% 
  mutate_at(vars(State), to_snake_case) %>% 
  setNames(tolower(colnames(.)))


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
  mutate(date = lubridate::ymd(date)) %>% 
  rename(cases = positive) %>% 
  rename(abbreviation = state) %>% 
  left_join(state_abbrev, by = "abbreviation") %>% 
  select(date, state, cases) %>% 
  filter(!is.na(state)) %>% 
  spread(state, cases) %>% 
  setNames(to_snake_case(colnames(.))) %>% 
  filter(date > lubridate::ymd("2020-03-01")) %>% 
  mutate_at(vars(-date), ~ ifelse(is.na(.x), 0, .x)) %>% 
  mutate_at(vars(-date), ~ c(NA, diff(.x))) %>% 
  .[-1, ]

shutdown_grid = 
  shutdown_dates %>% 
  map(~ convert_shutdown_dates_to_date_vector(.x, dat_diff, 3)) %>% 
  reduce(cbind.data.frame) %>% 
  setNames(names(shutdown_dates)) %>% 
  select(1, 3:20) %>% 
  select(-hawaii)

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
  identity()

date_vector = dat_diff$date[1:nrow(dat_recompiled_with_shutdowns)]
  
cum_p_observed = convert_filter_to_cumsum(empirical_timing_dist, nrow(dat_recompiled_with_shutdowns))

stan_data = list(
  timesteps = nrow(dat_recompiled_with_shutdowns),
  states = ncol(dat_recompiled_with_shutdowns),
  cases = dat_recompiled_with_shutdowns,
  shutdowns = shutdown_grid_capped,
  cum_p_observed = cum_p_observed,
  window = 20,
  simulation_steps = 15,
  init_theta_mean = 0.3352583,
  init_theta_sd = 0.3032382
)

