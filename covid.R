library(tidyverse)
library(KFAS)
library(bsts)
library(rstan)

#### dat ####

url = 'https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv'
dat = read_csv(url)



#### test ####

WINDOW = 20
SERIAL_INTERVAL = 4
GAMMA = 1 / SERIAL_INTERVAL

STATE = "New York"

series = 
  dat %>% 
  filter(state == STATE) %>%
  filter(cases>0) %>% 
  pull(cases) %>% 
  diff %>% 
  {. + 1} %>% 
  {c(rep(0, WINDOW-1), .)} %>% 
  rollsum(., WINDOW) 


stan_mod = stan_model("under-report.stan")

idfk = sampling(stan_mod, list(report_window = 20, timesteps = length(series), reported_cases = series),
                iter = 10000, chains = 2, cores = 2)

traceplot(idfk, "report_probs[7]")

dates = dat %>% 
  filter(state == STATE) %>%
  filter(cases>0) %>% 
  pull(date) %>% 
  .[c(-1,-2)]

ktm1 = series[-length(series)]
kt = series[-1]

mod = SSModel(
  kt ~ 1, u = ktm1,
  distribution = "poisson",
)

mod$Q[1,1,1] = NA

mod_fit = fitSSM(mod, c(1,1))

fitted(mod_fit$model) %>% plot

mod_fit_filtered = KFS(mod_fit$model, c("state"), c("state"))

ir  = tibble(
  mean_estimate = mod_fit_filtered$a[, 1],
  upper = mean_estimate + 1.96 * sqrt(mod_fit_filtered$P[1,1,]),
  lower = mean_estimate - 1.96 * sqrt(mod_fit_filtered$P[1,1,])
)[-1, ]


rt = (ir / GAMMA + 1) %>% 
  mutate(date = dates)

rt %>% 
  filter(date > lubridate::ymd("20200301")) %>% 
  ggplot() + 
  aes(x = date, y = mean_estimate, ymin = lower, ymax = upper) + 
  geom_line() + 
  geom_ribbon(alpha = .5) + 
  geom_hline(yintercept = 1)









