library(tidyverse)
library(KFAS)
library(zoo)
library(snakecase)
library(tictoc)
library(lubridate)
library(rstan)

url = 'https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv'
dat = read_csv(url)

WINDOW = 20
SERIAL_INTERVAL = 4
GAMMA = 1 / SERIAL_INTERVAL
STATES = dat$state %>% unique
DIM = length(STATES)

stan_mod = stan_model("impact-of-shutdown.stan")


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


shutdown_dates = 
  list(
    alabama = ymd(20200404),
    alaska = ymd(20200328),
    california = ymd(20200319),
    colorado = ymd(20200326),
    connecticut = ymd(20200323),
    delaware = ymd(20200324),
    district_of_columbia = ymd(20200401),
    florida = ymd(20200403),
    georgia = ymd(20200403),
    hawaii = ymd(20200325),
    idaho = ymd(20200325),
    illinois = ymd(20200321),
    indiana = ymd(20200324),
    kansas = ymd(20200330),
    kentucky = ymd(20200326),
    louisiana = ymd(20200323),
    maine = ymd(20200402),
    maryland = ymd(20200330),
    massachusetts = ymd(20200324),
    michigan = ymd(20200324),
    minnesota =ymd(20200331),
    mississippi = ymd(20200403),
    missouri = ymd(20200406),
    montana = ymd(20200328),
    nevada = ymd(20200401),
    new_hampshire = ymd(20200327),
    new_jersey = ymd(20200321),
    new_mexico = ymd(20200322)
  )

convert_shutdown_dates_to_date_vector = function(date) {
  
  case_when(
    dat_multivar$date < date ~ 0,
    dat_multivar$date > date + days(7) ~ 1,
    TRUE ~ as.numeric(difftime(dat_multivar$date, date, units = "days")) * 1/7
  )
  
}

shutdown_grid = shutdown_dates %>% 
  map(~ convert_shutdown_dates_to_date_vector(.x)) %>% 
  reduce(cbind.data.frame) %>% 
  setNames(names(shutdown_dates))


dat_multivar_with_shutdowns = dat_multivar %>% 
  select(names(shutdown_grid))

stan_mod = stan_model("impact-of-shutdown.stan")

fit = sampling(
  stan_mod, 
  list(
    timesteps = nrow(dat_multivar_with_shutdowns),
    states = ncol(dat_multivar_with_shutdowns),
    cases = dat_multivar_with_shutdowns,
    shutdowns = shutdown_grid
  ),
  chains = 2,
  cores = 2,
  iter = 4000
)

post = rstan::extract(fit)
stan_trace(fit, "shutdown_impact_on_rt")
stan_dens(fit, "shutdown_impact_on_rt")
print(fit, "shutdown_impact_on_rt")


itp1 = as.matrix(dat_multivar[-1, 2:ncol(dat_multivar)])
it = as.matrix(dat_multivar[-nrow(dat_multivar), 2:ncol(dat_multivar)])

observation_matrix = model.matrix(
  ~ 1 + f,
  data = data.frame(f = factor(1:DIM)),
  contrasts = list(f = "contr.sum")
)

mod_multivar = SSModel(
  itp1 ~ -1 + SSMcustom(
    Z = observation_matrix,
    T = diag(DIM),
    R = diag(DIM),
    a1 = rep(0, DIM),
    P1 = diag(c(0, rep(1, DIM-1))),
    P1inf = diag(c(1, rep(0, DIM-1))),
    Q = diag(DIM),
    n = nrow(itp1)
  ),
  u = it,
  distribution = "poisson"
)



update_fn = function(pars, mod) {
  QQ = diag(exp(pars[1:DIM]))
  mod$Q[, , 1] = QQ
  
  mod
}

# mod_multivar_fit = fitSSM(mod_multivar, rep(-5, DIM+1), update_fn, method = "BFGS")
# saveRDS(mod_multivar_fit, file = "mod_multivar_fit.rds")
# takes about 13 minutes to fit
mod_multivar_fit = readRDS("mod_multivar_fit.rds")

mod_multivar_filtered = KFS(mod_multivar_fit$model, c("state", "mean"), c("state", "mean"))

Z_augment = mod_multivar_filtered$model$Z[, ,1] %>% 
  rbind("average" = c(1, rep(0, DIM-1)))

theta = map(
  1:nrow(it), 
  ~ t(Z_augment %*% mod_multivar_filtered$alphahat[.x, ])
) %>% 
  reduce(rbind) %>% 
  as.data.frame 

compute_se_from_indices = function(var_index) {
  z_index = which(Z_augment[var_index, ] != 0)
  
  series = mod_multivar_filtered$P[z_index, z_index, ]
  if(is.null(dim(series))) return(sqrt(series))
  
  sqrt(apply(series, 3, sum))
}

ses = 1:(DIM+1) %>% 
  map(~ compute_se_from_indices(.x)) %>% 
  cbind.data.frame %>% 
  setNames(colnames(theta)) %>% 
  .[-1, ]

theta_upper = theta + 1.96 * ses
theta_lower = theta - 1.96 * ses

rt = theta/GAMMA + 1
rt_lower = theta_lower/GAMMA + 1
rt_upper = theta_upper/GAMMA + 1


rts_by_state = 
  list(rt, rt_lower, rt_upper) %>% 
  map2(
    c("mean", "lower", "upper"),
    ~ mutate(.x, idx = 1:n()) %>% 
        gather(-idx, key = state, value = !!.y)
  ) %>% 
  reduce(~ left_join(.x, .y, by = c("idx", "state")))

rts_by_state %>% 
  ggplot() + 
  aes(x = idx, y = mean, ymax = upper, ymin = lower) + 
  geom_line(color = "grey") + 
  geom_ribbon(alpha = 0.5) + 
  facet_wrap(~ state) + 
  geom_hline(yintercept = 1) +
  coord_cartesian(ylim = c(-1, 5))
