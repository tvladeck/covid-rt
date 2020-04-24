// https://mc-stan.org/docs/2_23/stan-users-guide/fit-gp-section.html

data {

  int timesteps; 
  int states;
  int cases[timesteps, states]; // diffed & convolved // eventually move this to xformed parameters
  vector[timesteps] cum_p_observed; // make this estimated given timing block
  
} 

transformed data {
  int N = timesteps-1;
  real x[N];
  real delta = 1e-9;
  
  
  for(i in 1:N) x[i] = 1;
}

parameters {
  matrix[timesteps-1, states] theta;
  real<lower=0> serial_interval;
  
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  vector[N] eta;
  
  
}

transformed parameters {
  real<lower=0> gam;
  matrix[timesteps-1, states] rt;
  
  gam = 1/serial_interval;
  rt = theta/gam + 1;
}

model {
  
  vector[N] f;
  {
    matrix[N, N] L_K;
    matrix[N, N] K = cov_exp_quad(x, alpha, rho);

    // diagonal elements
    for (n in 1:N)
      K[n, n] = K[n, n] + delta;

    L_K = cholesky_decompose(K);
    f = L_K * eta;
  }

  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma ~ std_normal();
  eta ~ std_normal();
  
  for(s in 1:states) theta[, s] ~ normal(f, sigma);
  
  serial_interval ~ gamma(6, 1.5);
  
  for(t in 2:timesteps) {
    for(s in 1:states) {
      real muhat;
      // we scale up yesterday's cases by the probability that we have observed it
      real inferred_cases_yesterday = (cases[t-1, s]/cum_p_observed[t-1]);
      // we scale down today's cases by the fact that we have not observed many of the later cases
      real expected_cases_today = cum_p_observed[t] * inferred_cases_yesterday * exp(theta[t-1, s]);
      // here guard against log(0) errors with init case
      // this 0.1 can be thought of as a seeding probability
      muhat = fmax(expected_cases_today, 0.1); 
      cases[t, s] ~ poisson(muhat);
    } 
  }
  
}

