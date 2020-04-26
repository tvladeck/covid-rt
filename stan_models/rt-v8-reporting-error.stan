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
  matrix[timesteps, states] onsets;
  matrix[timesteps, states] inferred_onsets;
  matrix[timesteps, states] scaled_inferred_onsets;
  matrix[timesteps-1, states] theta;
  matrix[timesteps-1, states] inferred_theta;
  matrix[timesteps-1, states] rt;
  matrix[timesteps, states] lambda;
  real<lower=0> gam;
  
  for(s in 1:states) {
    vector[timesteps] log_delta_onset; 
    log_delta_onset[1] = initial_onsetload[s];
    theta[, s] = cumulative_sum(theta_steps[, s]);
    log_delta_onset[2:timesteps] = theta[, s];
    onsets[, s] = exp(cumulative_sum(log_delta_onset));
    lambda[, s] = onset_to_cases * onsets[, s];
    inferred_onsets[, s] = cases_to_onsets * lambda[, s];
    scaled_inferred_onsets[, s] = inferred_onsets[, s] ./ cum_p_observed;
    inferred_theta[, s] = log(scaled_inferred_onsets[2:timesteps, s] ./ scaled_inferred_onsets[1:(timesteps-1), s]);
  }
  
  gam = 1/serial_interval;
  rt = inferred_theta/gam + 1;

}

model {
  tau ~ gamma(.5, .5);
  
  initial_onsetload ~ gamma(1, .1);
  
  to_vector(theta_steps[1:(timesteps-1), ]) ~ normal(0, tau[1]);
  
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

