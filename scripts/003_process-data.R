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

shutdown_grid = 
  shutdown_dates %>% 
  map(~ convert_shutdown_dates_to_date_vector(.x, dat_diff, 3)) %>% 
  reduce(cbind.data.frame) %>% 
  setNames(names(shutdown_dates)) %>% 
  select(1, 3:10) %>% 
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

tests_dat = 
  tests %>% select(colnames(dat_recompiled_with_shutdowns)) %>% 
  mutate_all( ~ ifelse(.x < 0, 0, .x))

stan_data = list(
  timesteps = nrow(dat_recompiled_with_shutdowns),
  states = ncol(dat_recompiled_with_shutdowns),
  tests = tests_dat,
  cases = dat_recompiled_with_shutdowns,
  shutdowns = shutdown_grid_capped,
  cum_p_observed = cum_p_observed,
  window = 20,
  simulation_steps = 15,
  init_theta_mean = 0.3352583,
  init_theta_sd = 0.3032382
)




deaths_dat = dat %>% 
  select(state, date, deathIncrease) %>% 
  mutate(death_date = ymd(date)) %>% 
  select(-date)

cases_dat = dat %>% 
  select(state, date, positiveIncrease) %>% 
  mutate(case_date = ymd(date)) %>% 
  select(-date)

all_dates = expand.grid(
  case_date = unique(cases_dat$case_date),
  death_date = unique(deaths_dat$death_date)
) %>% 
  filter(case_date < death_date)

cases_dat_expanded = all_dates %>% 
  left_join(cases_dat)

deaths_dat_expanded = deaths_dat %>% 
  left_join(cases_dat_expanded) %>% 
  mutate(days_ago = as.integer(difftime(death_date, case_date, units = "days"))) 

fin_dat1 = 
  deaths_dat_expanded %>% 
  select(state, deathIncrease, death_date, positiveIncrease, days_ago) %>% 
  spread(days_ago, positiveIncrease)

colnames(fin_dat1)[4:ncol(fin_dat1)] = str_c("d_", colnames(fin_dat1)[4:ncol(fin_dat1)]) 

fin_dat2 = 
  fin_dat1 %>% 
  mutate_at(vars(deathIncrease, starts_with("d_")), ~ ifelse(is.na(.x), 0, .x))

all_zeros = rowSums(fin_dat2[, c(2, 4:ncol(fin_dat2))]) == 0

fin_dat3 = fin_dat2[!all_zeros, ] %>% 
  filter(deathIncrease >= 0)


fin_dat_simple = fin_dat3[, c(2, 4:70)]
summary(lm(deathIncrease ~ -1 + ., fin_dat_simple))
coef(lm(deathIncrease ~ -1 + ., fin_dat_simple)) %>% plot
summary(glm(deathIncrease ~ -1 + ., fin_dat_simple, family = "poisson"))
coef(glm(deathIncrease ~ -1 + ., fin_dat_simple, family = "poisson")) %>% exp

library(glmnet)
