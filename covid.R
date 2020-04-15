library(tidyverse)
library(KFAS)
library(bsts)

#### dat ####

url = 'https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv'
dat = read_csv(url)


#### Before ####

N_LAGS = 7

x1 = rep(0, N_LAGS)

sum(x1)

plot(x1, type = "l", main = "R0 by duration", ylab = "", xlab = "day")

Px1 = rep(1/10*N_LAGS, N_LAGS)

dat_state = 
  dat %>% 
  filter(state == "New York") %>% 
  mutate(cases = c(NA, diff(cases))) %>% 
  {
    tmp <- .
    for(i in 1:N_LAGS) {
      var = str_c("lag_cases_", i)
      tmp <- mutate(tmp, !!var := lag(cases, i))
    }
    tmp
  } %>% 
  na.omit %>% 
  filter(cumsum(cases)>0) 


fm = as.formula(str_c("~ ", str_c("lag_cases_", 1:N_LAGS, collapse = "+")))


ss_spec = 
  SSModel(
    cases ~ -1 + SSMregression(fm, dat_state, Q = diag(NA, nrow = 1), remove.intercept = T) + 
      SSMseasonal(7),
    data = dat_state,
    distribution = "gaussian"
  )


ss_spec_fit = fitSSM(ss_spec, rep(0, N_LAGS+1))

ss_spec_filtered = KFS(ss_spec_fit$model, c("state", "mean"), c("state", "mean"))

variances = rep(0, nrow(dat_state))
for(i in 1:nrow(dat_state)) variances[i] = sum(ss_spec_filtered$V[1:N_LAGS, 1:N_LAGS, i])
ses = sqrt(variances)

tibble(
  mean = rowSums(ss_spec_filtered$alphahat[, 1:N_LAGS]),
) %>% mutate(
  lower = mean - 1.96 * ses,
  upper = mean + 1.96 * ses,
  date = dat_state$date
) %>% 
  ggplot() + 
  aes(x = date, ymin = lower, ymax = upper, y = mean) + 
  geom_line() + 
  geom_ribbon(alpha = .5)


ts.plot(cbind(ss_spec_filtered$muhat, dat_state$cases), col = c("red", "blue"))

imp <- simulateSSM(ss_spec_fit$model, antithetics = TRUE)

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









