data {

  int timesteps; 
  int states;
  int cases[timesteps, states];
  int window;
  vector[timesteps] cum_p_observed;
  
} 

transformed data {
  matrix[timesteps, states] scaled_cases;
  int cum_cases[timesteps, states];
  matrix[timesteps, states] cum_scaled_cases;
  
  for(s in 1:states) {
    scaled_cases[, s] = to_vector(cases[, s]) ./ cum_p_observed;    
    for(t in 1:timesteps) {
      int idx_low = max(1, t-window+1);
      // leaving this here as a note to myself
      // we can do the math either way, with cumulatives or diffs
      // but with diffs we only have to worry about scaling _on that day_
      // not about integrating the scaling vector across the vector of diffs 
      // to get the cumulative amount
      cum_cases[t, s] = sum(cases[idx_low:t, s]);
      cum_scaled_cases[t, s] = sum(scaled_cases[idx_low:t, s]);
    }
  }
  
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
      real expected_cases = cum_p_observed[t] * scaled_cases[t-1, s] * exp(theta[t-1, s]);
      mu = fmax(expected_cases, 0.1); // guard against log(0) errors with init case
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

