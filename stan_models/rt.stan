data {

  int timesteps; 
  int states;
  int cases[timesteps, states];
  
} 

parameters {
  real<lower=0> serial_interval;
  real<lower=0> step_size;
  matrix[timesteps-1, states] theta_steps;
}

transformed parameters {
  matrix[timesteps-1, states] theta;
  real<lower=0> gam;
  matrix[timesteps-1, states] rt;

  for(s in 1:states) {
    theta[, s] = cumulative_sum(theta_steps[, s]);
  }
  
  gam = 1/serial_interval;
  rt = theta/gam + 1;

}

model {
  
  serial_interval ~ gamma(6, 1.5);
  
  step_size ~ normal(0, 0.2)T[0, ];
  theta_steps[1, ] ~ normal(0, 1);
  to_vector(theta_steps[2:(timesteps-1), ]) ~ normal(0, step_size);
  
  for(t in 2:timesteps) {
    for(s in 1:states) {
      real mu;
      mu = fmax(cases[t-1, s] * exp(theta[t-1, s]), 0.1); // guard against log(0) errors with init case
      cases[t, s] ~ poisson(mu);
    } 
  }
  
}

generated quantities {
  // matrix[timesteps-1, states] pred_cases; 
  // 
  // for(t in 1:(timesteps-1)) {
  //   for(s in 1:states) {
  //     pred_cases[t, s] = poisson_rng(cases[t, s] * exp(theta[t, s]));
  //   } 
  // }
  
}

