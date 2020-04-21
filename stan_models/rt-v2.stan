functions {
  
  vector reverse(vector x) {
    int n = num_elements(x);
    vector[n] r;
    
    for(i in 1:n) {
      r[n-i+1] = x[i];
    }
    
    return(r);
  }
  
}

data {
  
  int timesteps;
  int states;
  int cases[timesteps, states]; // DIFFED
  
  int n;
  real timing[n];
  
  int MAX_COUNT; 
}

parameters {
  
  real<lower=0> serial_interval;
  real<lower=0> step_size;
  matrix[timesteps-1, states] theta_steps;
  
  real<lower=0> alpha; // gamma distribution parameter
  real<lower=0> beta; // gamma distribution parameter
  
  vector<lower=0>[states] initial_count; // 
  
}

transformed parameters {
  matrix[timesteps, states] unobserved_counts;
  matrix[timesteps-1, states] theta;
  vector[timesteps] prob_observed_that_day;
  matrix[timesteps, states] expected_counts;
  
  // given a case we observe today, with what probability did it land X days ago
  // we use this to build up our expected cases from a vector of 
  for(t in 1:timesteps) {
    prob_observed_that_day[t] = gamma_cdf(t, alpha, beta) - gamma_cdf(t-1, alpha, beta);
  }
  
  for(s in 1:states) {
    theta[, s] = cumulative_sum(theta_steps[, s]);
  }
  
  for(s in 1:states) {
    unobserved_counts[1, s] = initial_count[s];
    for(t in 2:timesteps) {
      unobserved_counts[t, s] = fmax(unobserved_counts[t-1, s] * exp(theta[t-1,s]), MAX_COUNT);
    }
  }
  
  
  for(s in 1:states) {
    for(t in 1:timesteps) {
      vector[t] rev_unobs_counts = reverse(unobserved_counts[1:t, s]);
      expected_counts[t, s] = dot_product(prob_observed_that_day[1:t], rev_unobs_counts);
    }
    
  }
  
} 

model {
  
  
  // this handles the offset between onset and reporting
  timing ~ gamma(alpha, beta);
  
  // initial_count // not sure
  initial_count ~ gamma(0.1, 0.1);
  
  // 
  serial_interval ~ gamma(6, 1.5);
  
  // our random walk
  step_size ~ normal(0, 0.2)T[0, ];
  theta_steps[1, ] ~ normal(0, 1);
  to_vector(theta_steps[2:(timesteps-1), ]) ~ normal(0, step_size);
  
  
  // actual versus expected cases
  for(s in 1:states) cases[, s] ~ poisson(expected_counts[, s]);
  
}
