data {

  int timesteps; 
  int states;
  int cases[timesteps, states];
  matrix[timesteps, states] shutdowns;
  
} 

parameters {
  
  real<lower=0> serial_interval;
  real<lower=0> sigma;
  matrix[timesteps-1, states] intcpt_steps;
  
  real impact_of_shutdown;
}

transformed parameters {
  matrix[timesteps-1, states] theta;
  matrix[timesteps-1, states] intcpt;
  real<lower=0> gam;
  matrix[timesteps-1, states] rt;
  
  gam = 1/serial_interval;
  
  for(s in 1:states) {
    intcpt[, s] = cumulative_sum(intcpt_steps[, s]);
  }
  
  theta = intcpt + impact_of_shutdown * shutdowns[2:timesteps, ];
  rt = theta/gam + 1;
  
}

model {
  
  impact_of_shutdown ~ normal(0, .5);
  
  serial_interval ~ gamma(2, .5);
  
  sigma ~ normal(0, 0.2)T[0, ];
  intcpt_steps[1, ] ~ normal(0, 1);
  to_vector(intcpt_steps[2:(timesteps-1), ]) ~ normal(0, sigma);
  
  for(t in 2:timesteps) {
    for(s in 1:states) {
      cases[t, s] ~ poisson(cases[t-1, s] * exp(theta[t-1, s]));
    } 
  }
  
}

