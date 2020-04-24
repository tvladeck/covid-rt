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
  
  real log_step_size;
  
  matrix[timesteps-1, states] intcpt_steps;
  
  real mean_shutdown_effect;
  vector[states] state_shutdown_effects_std;
  real<lower=0> tau;
}

transformed parameters {
  real<lower=0> gam;
  matrix[timesteps-1, states] rt;
  
  real<lower=0> step_size = exp(log_step_size);
  matrix[timesteps-1, states] theta;
  matrix[timesteps-1, states] intcpt;
  
  vector[states] shutdown_impact_on_rt;
  vector[states] state_shutdown_effects;
  
  state_shutdown_effects = mean_shutdown_effect + state_shutdown_effects_std * tau;

  for(s in 1:states) {
    intcpt[, s] = cumulative_sum(intcpt_steps[, s]);
    theta[, s] = intcpt[, s] + state_shutdown_effects[s] * shutdowns[2:timesteps, s];
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
  
  log_step_size ~ normal(0, 10);
  intcpt_steps[1, ] ~ normal(0, 1);
  to_vector(intcpt_steps[2:(timesteps-1), ]) ~ normal(0, step_size);
  
  for(t in 2:timesteps) {
    for(s in 1:states) {
      real mu;
      // we scale up yesterday's cases by the probability that we have observed it
      real inferred_cases_yesterday = (cases[t-1, s]/cum_p_observed[t-1]);
      // we scale down today's cases by the fact that we have not observed many of the later cases
      real expected_cases_today = cum_p_observed[t] * inferred_cases_yesterday * exp(theta[t-1, s]);
      // here guard against log(0) errors with init case
      // this 0.1 can be thought of as a seeding probability
      mu = fmax(expected_cases_today, 0.1); 
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

