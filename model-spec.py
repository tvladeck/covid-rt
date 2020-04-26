with pm.Model() as model:
  
  tau = pm.Gamma('tau', .5, .5)
  lam = pm.GaussianRandomWalk('lam', sd = tau, shape = len(cases))
  
  rho = pm.Gamma('rho', .5, .5)
  
  expected_cases = pm.Deterministic(exp(lam) * tests)
  
  filtered_expected_cases = expected_cases[non_zero_days]
  filtered_cases          = cases[non_zero_days]
  
  likelihood = pm.NegativeBinomial('cases', mu = filtered_expected_cases, 
                                    alpha = rho, observed = filtered_cases)
  
  what_we_want_to_use = exp(lam)

  
  
