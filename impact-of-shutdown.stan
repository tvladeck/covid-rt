data {

  int timesteps; 
  int states;
  int cases[timesteps, states];
  matrix[timesteps, states] shutdowns;
  
} 

parameters {
  
  real<lower=0> serial_interval;
  real<lower=0> sigma;
  matrix[timesteps-1, states] theta_steps;
  
}

transformed parameters {
  matrix[timesteps-1, states] theta;
  
  for(s in 1:states) {
    theta[, s] = cumulative_sum(theta_steps[, s]);
  }
}

model {
  
  serial_interval ~ gamma(2, .5);
  
  sigma ~ normal(0, 0.2)T[0, ];
  
  theta_steps[1, ] ~ normal(0, 1);
  to_vector(theta_steps[2:(timesteps-1), ]) ~ normal(0, sigma);
  
  for(t in 2:timesteps) {
    for(s in 1:states) {
      cases[t, s] ~ poisson(cases[t-1, s] * exp(theta[t-1, s]));
    } 
  }
  
}

