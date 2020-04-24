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
  matrix[timesteps-1, states] theta_steps;
}

transformed parameters {
  matrix[timesteps-1, states] theta;
  real<lower=0> gam;
  matrix[timesteps-1, states] rt;

  for(s in 1:states) {
    theta[, s] = cumulative_sum(theta_steps[, s]);
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

