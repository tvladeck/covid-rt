data {
  int timesteps; 
  int states;
  int cases_raw[timesteps, states]; 
  int tests[timesteps, states]; 
  real max_scaling_factor;
} 

parameters {
  matrix[timesteps, states] lambda_steps;
  vector<lower=0>[3] tau;
}

transformed parameters {
  matrix[timesteps, states] lambda;
  matrix[timesteps, states] fitted_cases;
  
  for(s in 1:states){
    lambda[, s] = cumulative_sum(lambda_steps[, s]);
    fitted_cases[, s] = exp(lambda[, s]);
  }

}

model {
  tau ~ gamma(.5, .5);
  
  lambda_steps[1, ] ~ normal(0, tau[1]);
  
  to_vector(lambda_steps[2:(timesteps), ]) ~ normal(0, tau[2]);
  
  for(s in 1:states) {
    for(t in 1:timesteps) {
      if(tests[t, s] > 0) {
        real scaled_tests = tests[t, s] / max(tests[, s]);
        real mutest = fmax(scaled_tests, .1);
        cases_raw[t, s] ~ neg_binomial_2(mutest .* exp(lambda[t, s]), tau[3]);
      }
    }
  }

} 

generated quantities {
  vector[timesteps * states] log_lik;
  matrix[timesteps, states] log_lik_matrix;
  
  for(s in 1:states) {
    for(t in 1:timesteps) {
      if(tests[t, s] > 0) {
        real scaled_tests = tests[t, s] / max(tests[, s]);
        real mutest = fmax(scaled_tests, .1);
        log_lik_matrix[t, s] = neg_binomial_2_lpmf(cases_raw[t, s] | mutest .* exp(lambda[t, s]), tau[3]);
      } else {
        log_lik_matrix[t, s] = 0;
      }
    }
  }
  
  log_lik = to_vector(log_lik_matrix);
  
}



