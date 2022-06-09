library(rstan)
library(msde)
cores <- parallel::detectCores()-1
options(mc.cores = cores)

data.names <- c("x")
param.names <- c("r", "K", "sigma")
logisticmod <- sde.make.model(ModelFile = "logistic.h",
                              data.names = data.names,
                              param.names = param.names)

N <- 50
r <- 0.25
K <- 100
sigma <- 0.1
dT <- 1
x0 <- 1

tseq <- seq(1,N,by=dT)
theta0 <- c(.25,100,.5)
logisticsim <- sde.sim(model = logisticmod, x0 = x0, theta = theta0,
                       nobs = N-1, # N-1 steps forward
                       dt = dT,
                       dt.sim = dT/100) # internal observation time

# plot data
Xobs <- c(x0, logisticsim$data) # include first observation
clrs <- c("orange","darkgreen")
plot(tseq,Xobs,pch=19,col=clrs[1])

stancode <- '
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real<lower=0> r;
  real<lower=1> K;
  real<lower=0> sigma;
}
model {
  real mu;
  real stddev;

  r ~ exponential(0.1); // Wide prior
  K ~ exponential(0.001); // Wide prior
  sigma ~ exponential(0.1); // Wide prior
  for (i in 2:N) {
    mu = y[i-1]+r*y[i-1]*(1-y[i-1]/K);
    stddev = sqrt(y[i-1]*sigma);
    y[i] ~ normal(mu,stddev);
  }
}'

mod <- stan_model(model_code = stancode) #, verbose = TRUE)
fit <- sampling(mod, data = list(y=Xobs,N=N))
# fit <- sampling(mod, data = list(y=y,N=N))

print(fit)
plot(fit)

params = extract(fit)

par(mfrow=c(1,3))
hist(params$r,main="",xlab=expression(r),col='skyblue',border='skyblue')
abline(v=r)
hist(params$K,main="",xlab=expression(K),col='skyblue',border='skyblue')
abline(v=K)
hist(params$sigma,main="",xlab=expression(sigma),col='skyblue',border='skyblue',xlim=c(0,1))
abline(v=sigma)
