
# does our filter work? ---------------------------------------------------

scaled_california_cases = stan_data$cases$california/cum_p_observed

fil = empirical_timing_dist

cbind(
  apply_1d_filter_rev_pad(fil, dat_diff$california),
  dat_diff$california
) %>% 
  ts.plot(col = c("red", "blue"))

cbind(
  filtered = apply_1d_filter_rev_pad(fil, dat_diff$california)[-1],
  ratio = cum_p_observed[-1]*(apply_1d_filter_rev_pad(fil, dat_diff$california) / cum_p_observed)[-48]
)


# how does this compare vs. previous estimates? ---------------------------



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