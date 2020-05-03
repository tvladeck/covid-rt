data {
  int timesteps; 
  int states;
  int cases[timesteps, states]; // raw
  int tests[timesteps, states]; // raw
  int deaths[timesteps, states]; // raw
} 

parameters {
  matrix[timesteps, states] lambda_steps;
  matrix[timesteps, states] rho_steps;
  vector<lower=0>[6] tau;
}

transformed parameters {
  matrix[timesteps, states] lambda;
  matrix[timesteps, states] rho;
  matrix[timesteps, states] expected_cases;
  matrix[timesteps, states] fitted_cases;
  matrix[timesteps, states] fitted_tests;
  matrix[timesteps, states] expected_deaths;
  vector[timesteps] death_weights;
  
  for(s in 1:states){
    lambda[, s]          = cumulative_sum(lambda_steps[, s]);
    rho[, s]             = cumulative_sum(rho_steps[, s]);
    fitted_tests[, s]    = exp(rho[, s]);
    expected_cases[, s]  = exp(lambda[, s]) .* to_vector(tests[, s]);
    fitted_cases[, s]    = exp(lambda[, s]) .* fitted_tests[, s];
    for(t in 1:timesteps) {
      
    }
  }

}

model {
  tau ~ gamma(.5, .5);
  
  lambda_steps[1, ] ~ normal(0, tau[1]);
  rho_steps[1, ] ~ normal(0, tau[2]);
  
  to_vector(lambda_steps[2:(timesteps), ]) ~ normal(0, tau[3]);
  to_vector(rho_steps[2:(timesteps), ]) ~ normal(0, tau[4]);
  
  for(s in 1:states) {
    for(t in 1:timesteps) {
      tests[t, s] ~ neg_binomial_2(fitted_tests[t, s], tau[5]);
      if(expected_cases[t, s] > 0) cases[t, s] ~ neg_binomial_2(expected_cases[t, s], tau[6]);
    }
  }

} 

generated quantities {
  vector[timesteps * states] log_lik;
  matrix[timesteps, states] log_lik_matrix;
  
  for(s in 1:states) {
    for(t in 1:timesteps) {
      log_lik_matrix[t, s] = neg_binomial_2_lpmf(cases[t, s] | fitted_cases[t, s], tau[5]);
    }
  }
  
  log_lik = to_vector(log_lik_matrix);
  
}



