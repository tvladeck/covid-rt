data {

  int timesteps; 
  int states;
  int cases[timesteps, states];
  matrix[timesteps, states] shutdowns;
  int simulation_steps;
  
} 

parameters {
  real<lower=0> serial_interval;
  real<lower=0> step_size;
  real<lower=0> tau;
  matrix[timesteps-1, states] intcpt_steps;
  real mean_shutdown_effect;
  vector[states] state_shutdown_effects_std;
  
}

transformed parameters {
  matrix[timesteps-1, states] theta;
  matrix[timesteps-1, states] intcpt;
  real<lower=0> gam;
  matrix[timesteps-1, states] rt;
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
  
  serial_interval ~ gamma(6, 1.5);
  
  step_size ~ normal(0, 0.2)T[0, ];
  intcpt_steps[1, ] ~ normal(0, 1);
  to_vector(intcpt_steps[2:(timesteps-1), ]) ~ normal(0, step_size);
  
  for(t in 2:timesteps) {
    for(s in 1:states) {
      real mu;
      mu = fmax(cases[t-1, s] * exp(theta[t-1, s]), 0.1); // guard against log(0) errors with init case
      cases[t, s] ~ poisson(mu);
    } 
  }
  
}

generated quantities {
  matrix[simulation_steps, states] theta_forward;

  for(s in 1:states) {
    theta_forward[, s] = theta[timesteps-1, s] + 
      to_vector(cumulative_sum(normal_rng(rep_vector(0, simulation_steps), step_size)));
  }
  
}

