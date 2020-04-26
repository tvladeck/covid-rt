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
  
  for(i in 1:timesteps) {
    cases_to_onsets[i, i:timesteps] = p_observed[1:(timesteps-(i-1))]';
  }
}

parameters {
  real<lower=0> serial_interval;
  matrix[timesteps, states] lambda_steps;
  vector<lower=0>[5] tau;
}



transformed parameters {
  matrix[timesteps-1, states] theta;
  matrix[timesteps, states] lambda;
  real<lower=0> gam;
  matrix[timesteps-1, states] rt;
  matrix[timesteps, states]   smoothed_cases;
  matrix[timesteps, states]   smoothed_onsets;
  matrix[timesteps, states]   log_smoothed_onsets;
  matrix[timesteps, states]   upscaled_log_smoothed_onsets;
  matrix[timesteps-1, states] expected_log_smoothed_onsets;
  matrix[timesteps-1, states] theta_steps;

  for(s in 1:states) {
    lambda[, s] = cumulative_sum(lambda_steps[, s]);
  }
  
  smoothed_cases = exp(lambda);
  
  for(s in 1:states) {
    smoothed_onsets[, s] = cases_to_onsets * smoothed_cases[, s];
  }
  
  log_smoothed_onsets = log(smoothed_onsets);
  
  for(s in 1:states) {
    upscaled_log_smoothed_onsets[, s] = log_smoothed_onsets[, s] - log(cum_p_observed);
    theta[, s] = upscaled_log_smoothed_onsets[2:timesteps, s] - upscaled_log_smoothed_onsets[1:(timesteps-1), s];
    theta_steps[2:(timesteps-1), s] = theta[2:timesteps-1, s] - theta[1:(timesteps-2), s];
    theta_steps[1, s] = theta[1, s];
  }

  gam = 1/serial_interval;
  rt = theta/gam + 1;


}

model {
  tau ~ gamma(.5, .5);
  
  
  // https://epiforecasts.io/covid/
  // an uncertain serial interval with a mean of 4.7 days (95% CrI: 3.7, 6.0) 
  // and a standard deviation of 2.9 days (95% CrI: 1.9, 4.9) [7] .
  serial_interval ~ gamma(2.626635, 0.5588585);
  
  theta_steps[1, ] ~ normal(0, 1);
  lambda_steps[1, ] ~ normal(0, tau[1]);
  
  to_vector(theta_steps[2:(timesteps-1), ]) ~ normal(0, tau[2]);
  to_vector(lambda_steps[2:timesteps, ]) ~ normal(0, tau[3]);
  
  
  for(s in 1:states) {
    for(t in 1:timesteps) {
      if(tests[t, s] > 0) {
        real scaled_tests = tests[t, s] / max(tests[, s]);
        real mutest = fmax(scaled_tests, .1);
        cases[t, s] ~ neg_binomial_2(mutest .* smoothed_cases[t, s], tau[4]);
      }
    }
  }

  
}

generated quantities {
  // vector[(timesteps-1)*states] log_lik;
  // matrix[timesteps-1, states] loglikm;
  // 
  // for(t in 1:(timesteps-1)) {
  //   for(s in 1:states) {
  //     real mu;
  //     mu = fmax(expected_onsets_today[t, s], 0.1); 
  //     loglikm[t, s] = poisson_lpmf(cases[t+1, s] | mu);
  //   } 
  // }
  // 
  // log_lik = to_vector(loglikm);
    
  
}

