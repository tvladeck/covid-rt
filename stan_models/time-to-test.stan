data {
  int n;
  real timing[n];
}
parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
}
model {
  timing ~ gamma(alpha, beta);
}
generated quantities {
  real time_draw;
  time_draw = gamma_rng(alpha, beta);
}

