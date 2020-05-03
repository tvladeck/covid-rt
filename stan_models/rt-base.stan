

data {

  int timesteps; 
  int states;
  int cases_convolved[timesteps, states]; 
  vector[timesteps] cum_p_observed; 
  real init_theta_mean;
  real init_theta_sd;
  
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
  matrix[timesteps-1, states] inferred_cases_yesterday;
  matrix[timesteps-1, states] expected_cases_today;

  for(s in 1:states) {
    theta[, s] = cumulative_sum(theta_steps[, s]);
    inferred_cases_yesterday[, s] = to_vector(cases_convolved[1:(timesteps-1), s]) ./ cum_p_observed[1:(timesteps-1)];
    expected_cases_today[, s] = cum_p_observed[2:timesteps] .* inferred_cases_yesterday[, s] .* exp(theta[, s]);
  }
  
  gam = 1/serial_interval;
  rt = theta/gam + 1;

}

model {
  
  // https://epiforecasts.io/covid/
  // an uncertain serial interval with a mean of 4.7 days (95% CrI: 3.7, 6.0) 
  // and a standard deviation of 2.9 days (95% CrI: 1.9, 4.9) [7] .
  serial_interval ~ gamma(2.626635, 0.5588585);
  
  step_size ~ normal(0, 0.2)T[0, ];
  
  
  theta_steps[1, ] ~ normal(init_theta_mean, init_theta_sd);
  to_vector(theta_steps[2:(timesteps-1), ]) ~ normal(0, step_size);
  
  for(t in 1:(timesteps-1)) {
    for(s in 1:states) {
      real mu;
      mu = fmax(expected_cases_today[t, s], 0.1); 
      cases_convolved[t+1, s] ~ poisson(mu);
    } 
  }
  
}

generated quantities {
  vector[(timesteps-1)*states] log_lik;
  matrix[timesteps-1, states] loglikm;
  
  for(t in 1:(timesteps-1)) {
    for(s in 1:states) {
      real mu;
      mu = fmax(expected_cases_today[t, s], 0.1); 
      loglikm[t, s] = poisson_lpmf(cases_convolved[t+1, s] | mu);
    } 
  }
  
  log_lik = to_vector(loglikm);
    
}

