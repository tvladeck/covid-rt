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
  real<lower=0> step_size;
  matrix[timesteps-1, states] theta_steps;
  matrix[timesteps, states] lambda_steps;
  vector<lower=0>[4] tau;
}



transformed parameters {
  matrix[timesteps-1, states] theta;
  matrix[timesteps, states] lambda;
  real<lower=0> gam;
  matrix[timesteps-1, states] rt;
  matrix[timesteps, states]   smoothed_cases;
  matrix[timesteps, states]   smoothed_onsets;
  matrix[timesteps, states]   log_smoothed_onsets;
  matrix[timesteps-1, states] inferred_onsets_yesterday;
  matrix[timesteps-1, states] expected_onsets_today;

  for(s in 1:states) {
    theta[, s] = cumulative_sum(theta_steps[, s]);
    lambda[, s] = cumulative_sum(lambda_steps[, s]);
  }
  
  smoothed_cases = exp(lambda);
  
  for(s in 1:states) {
    smoothed_onsets[, s] = cases_to_onsets * smoothed_cases[, s];
  }
  
  log_smoothed_onsets = log(smoothed_onsets);
  
  for(s in 1:states) {
    inferred_onsets_yesterday[, s] = to_vector(smoothed_onsets[1:(timesteps-1), s]) ./ cum_p_observed[1:(timesteps-1)];
    expected_onsets_today[, s] = cum_p_observed[2:timesteps] .* inferred_onsets_yesterday[, s] .* exp(theta[, s]);
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
  
  step_size ~ normal(0, tau[1])T[0, ];
  
  theta_steps[1, ] ~ normal(0, 1);
  lambda_steps[1, ] ~ normal(0, tau[2]);
  
  to_vector(theta_steps[2:(timesteps-1), ]) ~ normal(0, step_size);
  to_vector(lambda_steps[2:(timesteps-1), ]) ~ normal(0, step_size);
  
  
  for(s in 1:states) {
    for(t in 1:timesteps) {
      if(tests[t, s] > 0) {
        real scaled_tests = tests[t, s] / max(tests[, s]);
        real mutest = fmax(scaled_tests, .1);
        cases[t, s] ~ neg_binomial_2(mutest .* smoothed_cases[t, s], tau[3]);
      }
    }
  }
  
  for(t in 1:(timesteps-1)) {
    for(s in 1:states) {
      real mu;
      mu = log(fmax(expected_onsets_today[t, s], 0.1));
      // smoothed_onsets[t+1, s] ~ poisson(mu);
      log_smoothed_onsets[t+1, s] ~ normal(mu, tau[4]);
      // include correction for log absolute determinant of transformation
      // https://mc-stan.org/docs/2_21/stan-users-guide/changes-of-variables.html
      target += -log_smoothed_onsets[t+1, s];
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

