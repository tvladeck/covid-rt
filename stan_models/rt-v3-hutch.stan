// to do
// - bring the timing distribution into stan with the line list
// - do the convolution in stan as well
//   - this will make the dependent var a parameter
//   - https://www.youtube.com/watch?v=KOIudAB6vJ0
//   - or can i just round it s.t. i get an int? 

data {

  int timesteps; 
  int states;
  int cases[timesteps, states]; // diffed & convolved // eventually move this to xformed parameters
  vector[timesteps] cum_p_observed; // make this estimated given timing block
  real init_theta_mean;
  real init_theta_sd;
  // real gamma_shape;
  // real gamma_rate;
  
} 

parameters {
  real<lower=0> serial_interval;
  real<lower=0> step_size;
  matrix[timesteps-1, states] rt_steps;
}

transformed parameters {
  matrix[timesteps-1, states] theta;
  real<lower=0> gam;
  matrix[timesteps-1, states] rt;
  matrix[timesteps-1, states] inferred_cases_yesterday;
  matrix[timesteps-1, states] expected_cases_today;


  gam = 1/serial_interval;
  
  for(s in 1:states) {
    vector[timesteps-1] rt_ratio;
    rt[, s] = cumulative_sum(rt_steps[, s]);
    theta[, s] = gam * (rt[, s] - 1);
    rt_ratio[1] = 1;
    rt_ratio[2:(timesteps-1)] = rt[2:(timesteps-1), s] ./ rt[1:(timesteps-2), s];
    inferred_cases_yesterday[, s] = to_vector(cases[1:(timesteps-1), s]) ./ cum_p_observed[1:(timesteps-1)];
    expected_cases_today[, s] = cum_p_observed[2:timesteps] .* rt_ratio .* inferred_cases_yesterday[, s] .* exp(theta[, s]);
  }

}

model {
  
  // https://epiforecasts.io/covid/
  // an uncertain serial interval with a mean of 4.7 days (95% CrI: 3.7, 6.0) 
  // and a standard deviation of 2.9 days (95% CrI: 1.9, 4.9) [7] .
  serial_interval ~ gamma(2.626635, 0.5588585);
  
  step_size ~ normal(0, 0.2)T[0, ];
  
  rt_steps[1, ] ~ normal(2, 2);
  to_vector(rt_steps[2:(timesteps-1), ]) ~ normal(0, step_size);
  
  for(t in 1:(timesteps-1)) {
    for(s in 1:states) {
      real mu;
      mu = fmax(expected_cases_today[t, s], 0.1); 
      cases[t+1, s] ~ poisson(mu);
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
      loglikm[t, s] = poisson_lpmf(cases[t+1, s] | mu);
    } 
  }
  
  log_lik = to_vector(loglikm);
    
  
}

