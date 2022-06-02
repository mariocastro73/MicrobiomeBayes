# One the sdeModel class is created as the C++ level, it is compiled in R using the following commands:

require(msde)
set.seed(1234)
## Loading required package: msde
# put lotvolModel.h in the working directory
data.names <- c("H", "L")
param.names <- c("alpha", "beta", "gamma")
lvmod <- sde.make.model(ModelFile = "LV.h",
                        data.names = data.names,
                        param.names = param.names)

# simulation parameters
theta0 <- c(alpha = .5, beta = .0025, gamma = .3) # true parameter values
# x0 <- c(H = 150, L = 79) # initial SDE values
x0 <- c(H = 200, L = 150) # initial SDE values
N <- 50 # number of observations
dT <- 1 # time between observations (years)

# simulate data
lvsim <- sde.sim(model = lvmod, x0 = x0, theta = theta0,
                 nobs = N-1, # N-1 steps forward
                 dt = dT,
                 dt.sim = dT/100) # internal observation time

# plot data
Xobs <- rbind(c(x0), lvsim$data) # include first observation
tseq <- (1:N-1)*dT # observation times
clrs <- c("orange", "darkgreen")
# par(mfrow=c(1,1))
# 
# par(mar = c(4, 4, 1, 0)+.1)
# plot(x = 0, type = "n", xlim = range(tseq), ylim = range(Xobs),
#      xlab = "Time (years)", ylab = "Population")
# lines(tseq, Xobs[,"H"], type = "o", pch = 16, col = clrs[1])
# lines(tseq, Xobs[,"L"], type = "o", pch = 16, col = clrs[2])
# legend("topleft", legend = c("Hare", "Lynx"), fill = clrs)


# initialize the posterior sampler
init <- sde.init(model = lvmod, x = Xobs, dt = dT,
                 m = 1, theta = c(.1, .1, .1))

nsamples <- 2e4
burn <- 2e3
lvpost <- sde.post(model = lvmod, init = init,
                   hyper = NULL, #prior specification
                   nsamples = nsamples, burn = burn)

# posterior histograms
tnames <- expression(alpha, beta, gamma)
par(mfrow = c(2,2))
for(ii in 1:lvmod$nparams) {
  hist(lvpost$params[,ii], breaks = 25, freq = FALSE,
       xlab = tnames[ii],
       main = parse(text = paste0("p[1](", tnames[ii], "*\" | \"*bold(Y))")))
  # superimpose true parameter value
  abline(v = theta0[ii], lwd = 4, lty = 2)
}

# plot data
par(mar = c(4, 4, 1, 0)+.1)
plot(x = 0, type = "n", xlim = range(tseq), ylim = c(0,max(2*Xobs)),
     xlab = "Time (years)", ylab = "Population")


for(i in 1:100) {
  theta0 <- lvpost$params[i,]
  lvsim <- sde.sim(model = lvmod, x0 = x0, theta = theta0,
                   nobs = N-1, # N-1 steps forward
                   dt = dT,
                   dt.sim = dT/100) # internal observation time
  
  # plot data
  Xsamp <- rbind(c(x0), lvsim$data) # include first observation
  lines(tseq, Xsamp[,"H"],lwd=0.5,lty=2, col = clrs[1])
  lines(tseq, Xsamp[,"L"],lwd=0.5,lty=2, col = clrs[2])
}

lines(tseq, Xobs[,"H"], type = "o", pch = 24, col =1, bg= clrs[1])
lines(tseq, Xobs[,"L"], type = "o", pch = 24, col =1, bg= clrs[2])
legend("topleft", legend = c("Hare", "Lynx"), fill = clrs)
