# Mario 01/06/2022
#install.packages(c('msde','GillespieSSA2','ggplot2'))
library(msde)
library(GillespieSSA2)



set.seed(123)
source('auxfunctions.R')
# Build model
# put lotvolModel.h in the working directory
data.names <- c("x")
param.names <- c("r", "K", "sigma")
logisticmod <- sde.make.model(ModelFile = "logistic.h",
                        data.names = data.names,
                        param.names = param.names)

x0 <- 10
N <- 50 # number of observations
dT <- 1
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


# initialize the posterior sampler
init <- sde.init(model = logisticmod, x = as.matrix(Xobs), dt = dT,
                 m = 1, theta = c(1, 10, .1))

nsamples <- 2e4
burn <- 2e3
logisticpost <- sde.post(model = logisticmod, init = init,
                   hyper = NULL, #prior specification
                   nsamples = nsamples, burn = burn)

# posterior histograms
tnames <- expression(r, K, sigma)
par(mfrow = c(2,2))
for(ii in 1:logisticmod$nparams) {
  hist(logisticpost$params[,ii], breaks = 25, freq = FALSE,
       xlab = tnames[ii],
       main = parse(text = paste0("p[1](", tnames[ii], "*\" | \"*bold(Y))")),
       xlim=c(min(logisticpost$params[,ii])/1.2,max(logisticpost$params[,ii])*1.2))
  # superimpose true parameter value
  abline(v = theta0[ii], lwd = 4, lty = 2,col='orange')
}

# plot data
# par(mfrow=c(1,1))
par(mar = c(4, 4, 1, 0)+.1)
clrs <- c("orange","darkgreen")
plot(x = 0, type = "n", xlim = range(tseq), ylim = c(0,max(Xobs)*1.2),
     xlab = "Time (years)", ylab = "Population")


x0 <- Xobs[1]

# Funciones auxiliares (para pintar cosucas)
colq <- function(X,q) apply(X, 2, function(X) as.numeric(quantile(X,q)))

plotPostPred <- function(theta0,lwd=0.1,lty=2, col = 'orange') {
  
  logisticsim <- sde.sim(model = logisticmod, x0 = x0, theta = theta0,
                         nobs = N-1, # N-1 steps forward
                         dt = dT,
                         dt.sim = dT/100) # internal observation time
  
  # plot data
  Xsamp <- c(x0, logisticsim$data) # include first observation
  lines(tseq, Xsamp,lwd=lwd,lty=lty,col=col)
  return(Xsamp)
}
####

temp <- c()
for(i in 1:500) {
  theta0 <- logisticpost$params[i,]
  temp <- rbind(temp,plotPostPred(theta0))
}

lines(tseq, Xobs, type = "p", pch = 24, col =1, bg= clrs[1])

legend("bottomright", legend = c("x", "posterior"), pch=c(19,-1),lty=c(-1,2),col='orange')


lines(tseq,colq(temp,.5))
lines(tseq,colq(temp,.025),col='gray')
lines(tseq,colq(temp,.975),col='gray')

