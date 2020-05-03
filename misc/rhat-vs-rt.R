par(mfrow=c(1,3))

state_idx = 2

ex_rt = post$rt %>% apply(c(2,3), mean) %>% .[, state_idx]

scaling = 1.1

rhat_delta = function(x) {
  d = c(x[1]*scaling, x)
  l = length(d)
  sum((ex_rt - d[2:l] - log(d[2:l]/d[1:(l-1)])*4)**2)
}

produce_rhat = function(x) {
  d = c(x[1]*scaling, x)
  l = length(d)
  return(d[2:l] + log(d[2:l]/d[1:(l-1)])*4)
}

find_rt = optim(par = ex_rt+rnorm(length(ex_rt), sd = 0.001), fn = rhat_delta, lower=rep(0.01, length(ex_rt)), method = "L-BFGS-B")
find_rt$value
find_rt$par

mt = str_c(colnames(stan_data$cases)[state_idx], "; K = ", scaling)

cbind(find_rt$par, ex_rt) %>% 
  ts.plot(col = c("red", "blue"), ylim = c(.8, 1.6), main = mt, xlim=c(20,60))

legend("topright", c("adjusted", "unadjusted"), col = c("red", "blue"), lty = 1)

# rt_0.9 = find_rt$par
# rt_1 = find_rt$par
# rt_1_1 = find_rt$par

cbind(rt_0.9, rt_1, rt_1_1) %>% 
  ts.plot(col = c("red", "green", "blue"))
legend("topright", c("K=.9", "K=1", "K=1.1"), col = c("red", "green", "blue"), lty = 1)

par(mfrow=c(1,1))