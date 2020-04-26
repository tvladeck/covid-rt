data {
  int timesteps; 
  int states;
  int cases[timesteps, states]; // raw
  int tests[timesteps, states]; // raw
  vector[timesteps] p_observed; // make this estimated
  vector[timesteps] cum_p_observed; // make this estimated given timing block
} 
transformed data {
  matrix[timesteps, timesteps] cases_to_onsets = rep_matrix(0, timesteps, timesteps);
  matrix[timesteps, timesteps] onset_to_cases = rep_matrix(0, timesteps, timesteps);
  
  for(i in 1:timesteps) {
    cases_to_onsets[i, i:timesteps] = p_observed[1:(timesteps-(i-1))]';
  }
  
  onset_to_cases = inverse(cases_to_onsets);
}
parameters {
  matrix[timesteps-1, states] theta_steps;
  vector<lower=0>[states] initial_onsetload;
  real<lower=0> serial_interval;
  vector<lower=0>[2] tau;
}
transformed parameters {
  matrix[timesteps-1, states] theta;
  matrix[timesteps, states] inferred_onsets;
  matrix[timesteps, states] lambda;
  real<lower=0> gam;
  
  for(s in 1:states) {
    theta[, s] = cumulative_sum(theta_steps[, s]);
    inferred_onsets[1, s] = initial_onsetload[s];
    for(t in 1:(timesteps-1)) {
      inferred_onsets[t+1, s] = inferred_onsets[t, s] * exp(theta[t, s]);
    }
    lambda[, s] = onset_to_cases * inferred_onsets[, s];
    
    for(t in 1:timesteps) {
      lambda[t, s] = fmax(lambda[t, s], 0.01);
    }
  }
  
  gam = 1/serial_interval;

}

model {
  tau ~ gamma(.5, .5);
  
  initial_onsetload ~ gamma(1, .1);
  
  theta_steps[1, ] ~ normal(0, 1);
  
  to_vector(theta_steps[2:(timesteps-1), ]) ~ normal(0, tau[1]);
  
  // https://epiforecasts.io/covid/
  // an uncertain serial interval with a mean of 4.7 days (95% CrI: 3.7, 6.0) 
  // and a standard deviation of 2.9 days (95% CrI: 1.9, 4.9) [7] .
  serial_interval ~ gamma(2.63, 0.56);
  
  for(s in 1:states) {
    for(t in 1:timesteps) {
      if(tests[t, s] > 0) {
        real scaled_tests = tests[t, s] / max(tests[, s]);
        real mutest = fmax(scaled_tests, .1);
        cases[t, s] ~ neg_binomial_2(mutest .* lambda[t, s], tau[2]);
      }
    }
  }

}

generated quantities {

}
