library(rstan)

model <- stan_model("seir.stan")

dt <- readRDS("SIRsample.rds")
head(dt)
plot(dt$onsets)

data <- list(
  cases = dt$onsets,
  n_days = length(dt$onsets),
  t0 = 0,
  tswitch = 10,
  N = 1e5,
  use_likelihood = 1 # lets us explore the priors
)
data$ts <- seq(1, data$n_days, by = 1)

options(mc.cores = 4)
fit_nuts <- sampling(model,
                     data = data,
                     chains = 4,
                     seed = 0)
