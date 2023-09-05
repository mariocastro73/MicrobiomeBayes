library(forecast)
N <- 100
y <- arima.sim(n=N,list(ar=c(.5)))

library("rstan")
cores <- parallel::detectCores()-1
options(mc.cores = cores)


stancode <- '
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  for (n in 2:N) {
    y[n] ~ normal(alpha + beta * y[n-1], sigma);
  }
}'

mod <- stan_model(model_code = stancode, verbose = TRUE)
# fit <- sampling(mod, data = list(y=y,n=n))
fit <- sampling(mod, 
  data = list(y=y,N=N),    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  refresh = 1000          # show progress every 'refresh' iterations
  )


print(fit)
plot(fit)
params = extract(fit)


par(mfrow=c(2,3))
ts.plot(params$alpha,xlab="Iterations",ylab=expression(alpha))
ts.plot(params$beta,xlab="Iterations",ylab=expression(beta))
ts.plot(params$sigma,xlab="Iterations",ylab=expression(sigma))

hist(params$alpha,main="",xlab=expression(alpha),col='skyblue',border='skyblue')
abline(v=0)
hist(params$beta,main="",xlab=expression(beta),col='skyblue',border='skyblue')
abline(v=0.5)
hist(params$sigma,main="",xlab=expression(sigma),col='skyblue',border='skyblue')          
abline(v=1)