data {
  int report_window;
  int timesteps;
  int reported_cases[timesteps];
}

parameters {
  real theta[timesteps];
  real<lower=0> sigma;
  real beginning_cases; 
  real<lower=0, upper=1> report_prob;
}

transformed parameters {
  vector[timesteps] expected_cases;
  
  expected_cases[1] = beginning_cases * exp(theta[1]);
  
  for(t in 2:timesteps) {
    expected_cases[t] = expected_cases[t-1] * exp(theta[t]);
  }
}

model {
  
  
  sigma ~ normal(0, .1)T[0, ];
  
  theta[1] ~ normal(0, 1);
  for(t in 2:timesteps) {
    theta[t] ~ normal(theta[t-1], sigma);
  }
  
  
  
  reported_cases ~ poisson(report_prob*expected_cases);
    
  
}
