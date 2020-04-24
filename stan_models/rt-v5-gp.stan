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
  
} 

transformed data {
  int N = timesteps-1;
  real x[N];
  vector[N] mu = rep_vector(0, N);
  
  for(i in 1:N) x[i] = 1;
}

parameters {
  matrix[timesteps-1, states] theta;
  real<lower=0> serial_interval;
  
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

transformed parameters {
  real<lower=0> gam;
  matrix[timesteps-1, states] rt;
  
  gam = 1/serial_interval;
  rt = theta/gam + 1;
}

model {
  
  matrix[N, N] L_K;
  matrix[N, N] K = cov_exp_quad(x, alpha, rho);
  real sq_sigma = square(sigma);

  // diagonal elements
  for (n in 1:N)
    K[n, n] = K[n, n] + sq_sigma;

  L_K = cholesky_decompose(K);

  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma ~ std_normal();
  
  for(s in 1:states) theta[, s] ~ multi_normal_cholesky(mu, L_K);
  
  
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

generated quantities {
  // matrix[timesteps-1, states] pred_cases; 
  // 
  // for(t in 1:(timesteps-1)) {
  //   for(s in 1:states) {
  //     pred_cases[t, s] = poisson_rng(cases[t, s] * exp(theta[t, s]));
  //   } 
  // }
  
}

