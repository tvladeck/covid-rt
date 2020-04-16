library("lme4", quietly = TRUE)
y_split <- split(sleepstudy["Reaction"], sleepstudy["Subject"])
p <- length(y_split)
y <- matrix(unlist(y_split), ncol = p,
            dimnames = list(NULL, paste("Subject", names(y_split))))

dataf <- split(sleepstudy, sleepstudy["Subject"])

P1 <- as.matrix(.bdiag(replicate(p, matrix(NA, 2, 2), simplify = FALSE)))

model_lmm <- SSModel(
  y ~ -1 +
    SSMregression(
      rep(list(~ Days), p), 
      type = "common", 
      data = dataf,
      remove.intercept = FALSE
    ) +
    SSMregression(
      rep(list(~ Days), p), 
      data = dataf,
      remove.intercept = FALSE, P1 = P1
    ),
    H = diag(NA, p)
  )

update_lmm <- function(pars, model) {
  P1 <- diag(exp(pars[1:2]))
  P1[1, 2] <- pars[3]
  P1 <- crossprod(P1)
  model["P1", states = 3:38] <-
    as.matrix(.bdiag(replicate(p, P1, simplify = FALSE)))
  model["H"] <- diag(exp(pars[4]), p)
  model
}

fit_lmm <- fitSSM(model_lmm, c(1, 1, 1, 5), update_lmm, method = "BFGS")

idk = KFS(fit_lmm$model, "state")



model_lmm2 <- SSModel(
  y ~ -1 +
    SSMregression(
      ~ 1, 
      type = "common", 
      remove.intercept = FALSE
    ) +
    SSMregression(
      ~ 1, 
      remove.intercept = FALSE, P1 = P1
    ),
  H = diag(NA, p)
)

url = 'https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv'
dat = read_csv(url)

dat_multivar = dat %>% 
  filter(state %in% c("New York", "Washington", "Louisiana")) %>% 
  select(date, state, cases) %>% 
  spread(state, cases) %>% 
  setNames(c("date", "ny", "wa", "la")) %>% 
  filter(date > lubridate::ymd("2020-03-01")) %>% 
  mutate_at(vars(ny, wa, la), ~ ifelse(is.na(.x), 0, .x)) %>% 
  mutate_at(vars(ny, wa, la), function(x) {
    diff(x) %>% 
      {. + 1} %>% 
      {c(rep(0, WINDOW), .)} %>% 
      rollsum(., WINDOW) 
  }) %>% 
  .[-1, ]



itp1 = as.matrix(dat_multivar[-1, 2:4])
it = as.matrix(dat_multivar[-nrow(dat_multivar), 2:4])
dummy_data = data.frame(intcpt = rep(1, nrow(it)))

mod_multivar = SSModel(
  itp1 ~ -1 + 
    SSMregression(
      ~ intcpt, data = dummy_data, remove.intercept = T, type = "common", Q = NA
    ) + 
    SSMregression(
      ~ intcpt, data = dummy_data, remove.intercept = T, Q = diag(NA, 3)
    ),
  u = it,
  distribution = "poisson"
)

diag(mod_multivar$P1)    = c(0, rep(1, 3))
diag(mod_multivar$P1inf) = c(1, rep(0, 3))

mod_multivar_fn = function(pars, mod) {
  
  QQ = diag(exp(pars[1:4]))
  # QQ[1, 2:4] = pars[5:7]
  # QQ[2:4, 1] = pars[5:7]
  
  mod$Q[, , 1] = QQ
  
  mod

}

mod_multivar_fit = fitSSM(mod_multivar, rep(-5, 4), mod_multivar_fn, method = "BFGS")

mod_multivar_filter = KFS(mod_multivar_fit$model, c("state", "mean"), c("state", "mean"))

cbind(fitted(mod_multivar_filter)[, 3], itp1[, 3]) %>% diff %>%  ts.plot()


idk = MASS::mvrnorm(n = 1000, mu = c(0, 0), Sigma = mod_multivar_filter$P[1:2, 1:2, 44])
idk %>% rowSums %>% var
mod_multivar_filter$P[1:2, 1:2, 44] %>% sum



