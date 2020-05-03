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

cases = 
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

tests = 
  dat %>% 
  mutate(date = lubridate::ymd(date)) %>% 
  rename(cases =  totalTestResults) %>% 
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

shutdowns = 
  shutdown_dates %>% 
  map(~ convert_shutdown_dates_to_date_vector(.x, cases, 3)) %>% 
  reduce(cbind.data.frame) %>% 
  setNames(names(shutdown_dates)) 

active_states = colnames(shutdowns)

cases_convolved = 
  cases %>% 
  select(-date) %>% 
  map_df(~ round(apply_1d_filter_rev_pad(empirical_timing_dist, .x))) %>% 
  identity()

deaths = 
  dat %>% 
  rename(abbreviation = state) %>% 
  left_join(state_abbrev, by = "abbreviation") %>% 
  select(state, date, deathIncrease) %>% 
  na.omit %>% 
  mutate(death_date = ymd(date)) %>% 
  select(-date) %>% 
  spread(state, deathIncrease) %>% 
  mutate_all(~ ifelse(is.na(.x), 0, .x))

date_vector = cases$date
  
cum_p_observed = convert_filter_to_cumsum(empirical_timing_dist, nrow(cases_convolved))

tests = 
  tests %>% 
  select(colnames(dat_recompiled_with_shutdowns)) %>% 
  mutate_all( ~ ifelse(.x < 0, 0, .x))

stan_data = list(
  timesteps = nrow(tests),
  states = length(active_states),
  tests = select(tests, active_states),
  cases_raw = select(cases, active_states),
  cases_convolved = select(cases_convolved, active_states),
  deaths = select(deaths, active_states),
  shutdowns = shutdowns,
  cum_p_observed = cum_p_observed,
  max_scaling_factor = 10,
  init_theta_mean = 0.3352583,
  init_theta_sd = 0.3032382
)



