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
  matrix[timesteps, states] shutdowns;
  
} 

parameters {
  real<lower=0> serial_interval;
  
  real<lower=0> step_size;
  
  matrix[timesteps-1, states] intcpt_steps;
  
  real mean_shutdown_effect;
  vector[states] state_shutdown_effects_std;
  real<lower=0> tau;
}

transformed parameters {
  real<lower=0> gam;
  matrix[timesteps-1, states] rt;

  matrix[timesteps-1, states] theta;
  matrix[timesteps-1, states] intcpt;
  
  vector[states] shutdown_impact_on_rt;
  vector[states] state_shutdown_effects;
  
  matrix[timesteps-1, states] inferred_cases_yesterday;
  matrix[timesteps-1, states] expected_cases_today;

  
  state_shutdown_effects = mean_shutdown_effect + state_shutdown_effects_std * tau;

  for(s in 1:states) {
    intcpt[, s] = cumulative_sum(intcpt_steps[, s]);
    theta[, s] = intcpt[, s] + state_shutdown_effects[s] * shutdowns[2:timesteps, s];
    inferred_cases_yesterday[, s] = to_vector(cases[1:(timesteps-1), s]) ./ cum_p_observed[1:(timesteps-1)];
    expected_cases_today[, s] = cum_p_observed[2:timesteps] .* inferred_cases_yesterday[, s] .* exp(theta[, s]);
  }
  
  gam = 1/serial_interval;
  rt = theta/gam + 1;
  shutdown_impact_on_rt = state_shutdown_effects/gam;

}

model {
  
  tau ~ normal(0, .5)T[0, ];
  state_shutdown_effects_std ~ std_normal();
  mean_shutdown_effect ~ normal(0, 1);
  
  serial_interval ~ gamma(6, 1.5);
  
  step_size ~ normal(0, 0.05);
  
  intcpt_steps[1, ] ~ normal(2, 2);
  to_vector(intcpt_steps[2:(timesteps-1), ]) ~ normal(0, step_size);
  
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

