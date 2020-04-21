convert_shutdown_dates_to_date_vector = function(date) {
  
  case_when(
    dat_multivar$date < date ~ 0,
    dat_multivar$date >= date ~ 1
  )
  
}

summarize_rt_from_posterior = function(post, orig_dat, date_vec) {
  
  mean_mtx  = post$rt %>% apply(c(2,3), mean)
  upper_mtx = post$rt %>% apply(c(2,3), function(x) quantile(x, .95))
  lower_mtx = post$rt %>% apply(c(2,3), function(x) quantile(x, .05))
  
  r = list(mean_mtx, upper_mtx, lower_mtx) %>% 
    map2(c("mean", "upper", "lower"), 
         ~ .x %>% 
           as.data.frame %>% 
           setNames(colnames(orig_dat)) %>% 
           mutate(date = date_vector[-1]) %>% 
           gather(-date, key = state, value = !!.y)) %>% 
    reduce(~ left_join(.x, .y, by = c("date", "state")))
  
  return(r)
  
}

plot_rt_from_posterior = function(...) {
  s = 
    summarize_rt_from_posterior(...) %>% 
    gather(-date, -state, key = series, value = rt)
  
  
  
  ggplot(s) + 
    aes(x = date, y = rt, color = series, lty = series) + 
    geom_point(size = .1) + 
    geom_line() + 
    scale_color_manual("", values = c("mean" = "red", "lower" = "grey", "upper" = "grey")) + 
    theme_bw() + 
    scale_linetype_manual("", values = c("mean" = 1, "lower" = 2, "upper" = 2)) + 
    theme(legend.position = "none") + 
    labs(y = "Rt", x = "") + 
    facet_wrap(~ state)
  
}

check_pp = function(post, dat_grid) {
  
  exp_theta = post$theta %>% exp
  mu = exp_theta
  
  idk = map(1:dim(exp_theta)[1], function(i) {
    exp_theta[i, , ] * dat_grid[1:dim(exp_theta)[2], ]
  })
  
  idk2 = map(idk, ~ apply(.x, 2, function(x) {rpois(length(x), x)}))
  
  
  idk3 = map(idk2, ~ as.data.frame(.x) %>% mutate(idx = 1:nrow(.)) %>% 
               gather(-idx, key = state, value = pred_cases))
  
  idk4 = reduce(idk3[1:100], rbind.data.frame)
  for(i in 2:dim(mu)[1]/100) {
    
    tic(str_c("trying ", i))
    idk4 = rbind.data.frame(
      idk4, 
      reduce(idk3[(100*(i-1)+1):(100*i)], rbind.data.frame)
    )
    toc()
    
  }
  
  idk5 = idk4 %>% 
    group_by(idx, state) %>% 
    summarize(
      pred_cases_05 = quantile(pred_cases, .05),
      pred_cases_50 = median(pred_cases),
      pred_cases_95 = quantile(pred_cases, .95)
    )
  
  idk6 = idk5 %>% 
    left_join(
      dat_grid %>% 
        mutate(idx = 1:nrow(.)) %>% 
        gather(-idx, key=state, value=actual), 
      by = c("idx", "state")
    ) 
  
  return(idk6)
}


apply_1d_filter = function(fil, ser) {
  
  lser = length(ser)
  lfil = length(fil)
  
  newser = rep(0, lser)
  
  for(i in 1:lser) {
    idx = i:min(i+lfil-1, lser)
    newser[idx] = newser[idx] + ser[i] * fil[1:length(idx)]
  }
  
  return(newser)
  
}

apply_1d_filter_rev_pad = function(fil, ser, min_fil_len = 15, cap = T) {
  
  lfil = length(fil)
  lser = length(ser)
  
  
  ser2 = rev(c(rep(0, min_fil_len-1), ser))
  
  if(cap) lfin = lser-min_fil_len else lfin = lser2
  
  r = rev(apply_1d_filter(fil, ser2))[1:lfin]
  
  return(r)
  
}
