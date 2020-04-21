for(f in list.files("scripts", full.names = T)) source(f)


a = stan_data$cases %>% 
  mutate_all(~ c(NA, diff(.x)/cum_p_observed[-nrow(stan_data$cases)])) %>% 
  .[-1, ] %>% 
  mutate_all(~ c(NA, .x[-1]/.x[-(nrow(stan_data$cases)-1)]))


stan_mod = stan_model("stan_models/rt-v2.stan")

dat_frame = dat_diff[, 2:10]

stan_data = list(
  timesteps = nrow(dat_frame),
  states = ncol(dat_frame),
  cases = dat_frame,
  n = nrow(delay), 
  timing = as.numeric(delay$delta),
  MAX_COUNT = dat_frame %>% max * 100
)


fit = sampling(
  stan_mod, 
  stan_data,
  chains = 2,
  cores = 2,
  iter = 4000
)

post=rstan::extract(fit)

post$expected_counts


# stan_mod = stan_model("stan_models/impact-of-shutdown.stan")
# stan_mod = stan_model("stan_models/rt.stan")

fit = sampling(
  stan_mod, 
  stan_data,
  chains = 2,
  cores = 2,
  iter = 4000
)

post = rstan::extract(fit)

stan_trace(fit, "shutdown_impact_on_rt[12]")
print(fit, "shutdown_impact_on_rt")

newrt = summarize_rt_from_posterior(post, stan_data$cases, date_vector)

plot_rt_from_posterior(post, stan_data$cases, date_vector)

state_abbrev = 
  read_csv("https://raw.githubusercontent.com/jasonong/List-of-US-States/master/states.csv") %>% 
  mutate_at(vars(State), to_snake_case)

existingrt = 
  read_csv("https://d14wlfuexuxgcm.cloudfront.net/covid/rt.csv") %>% 
  left_join(state_abbrev, by = c("state"="Abbreviation"))


comparison = 
  existingrt %>% 
  left_join(newrt, by = c("date" = "date", "State" = "state")) %>% 
  select(-State) %>% 
  na.omit

plot_state_comparisons = function() {
  comparison %>% 
    # gather(-state, -date, key = series, value = rt) %>% 
    # filter(state == !!state) %>% 
    ggplot() + 
    aes(x = date) + 
    geom_line(aes(y = ML), col = "blue") + 
    geom_line(aes(y = mean), col = "red") + 
    geom_ribbon(aes(ymin = Low_50, ymax = High_50), fill = "lightblue", alpha = .75) + 
    geom_ribbon(aes(ymin = Low_90, ymax = High_90), fill = "lightblue", alpha = .75) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "pink", alpha = .75) +
    # labs(y = "Rt", subtitle = state, x = "") +
    labs(y = "Rt",  x = "") +
    theme_bw() + 
    geom_hline(yintercept = 1) + 
    facet_wrap(~ state)+ 
    coord_cartesian(ylim = c(-.75, 4))

}
