
vec = stan_data$cases$california
ret = rep(0, 1000)
for(i in 1:length(vec)) {
  retstar = rep(0, 100)
  n = vec[i]
  for(j in 1:n) {
    d = sample(delay$delta, 1)
    retstar[d] = retstar[d] + 1
  }
  ret[i:(i+99)] = ret[i:(i+99)] + retstar
}

sim_data = rpois(52, ret[1:52])
sim_data = rnbinom(52, mu = ret[1:52], size = 10)

ts.plot(cbind(sim_data, dat_diff$california), col = c("red", "blue"))

ts.plot(cbind(
  apply_1d_filter_rev_pad(empirical_timing_dist, sim_data),
  vec
), col = c("red", "blue")
)
