# Mario 29/05/2022
#install.packages(c('msde','GillespieSSA2','ggplot2'))
library(msde)
library(GillespieSSA2)
set.seed(123)

source('auxfunctions.R')
# Markov-Chain simulation using library GillespieSSA2
sim_name <- "Lotka Predator-Prey model"
params <- c(c1 = .4, c2 = .005, c3 = .2)
N <- 50 # number of observations
final_time <- N
initial_state <- c(Y1 = 100, Y2 = 75)

reactions <- list(
  reaction("c1 * Y1", c(Y1 = +1)),
  reaction("c2 * Y1 * Y2", c(Y1 = -1, Y2 = +1)),
  reaction("c3 * Y2", c(Y2 = -1))
)

out <- ssa(
  initial_state = initial_state,
  reactions = reactions,
  params = params,
  final_time = final_time,
  method = ssa_exact(),
  # method = ssa_etl(.25),
  sim_name = sim_name
) 
plot_ssa(out)

# Build model
# put lotvolModel.h in the working directory
data.names <- c("H", "L")
param.names <- c("alpha", "beta", "gamma")
lvmod <- sde.make.model(ModelFile = "LV.h",
                        data.names = data.names,
                        param.names = param.names)


#####

Xobs <- out$state[seq(1,length(out$time),length=N),]
tseq <- out$time[seq(1,length(out$time),length=N)]
colnames(Xobs) <- c("H","L")

dT <- 1
# initialize the posterior sampler
init <- sde.init(model = lvmod, x = Xobs, dt = dT,
                 m = 1, theta = c(.1, .01, .1))

nsamples <- 2e4
burn <- 2e4
lvpost <- sde.post(model = lvmod, init = init,
                   hyper = NULL, #prior specification
                   nsamples = nsamples, burn = burn)

# posterior histograms
tnames <- expression(alpha, beta, gamma)
par(mfrow = c(2,3))
theta0 <- params  # true parameter values
for(ii in 1:lvmod$nparams) {
  hist(lvpost$params[,ii], breaks = 25, freq = FALSE,
       xlab = tnames[ii],
       main = parse(text = paste0("p[1](", tnames[ii], "*\" | \"*bold(Y))")),
       xlim=c(min(lvpost$params[,ii])/1.1,max(lvpost$params[,ii])*1.1))
  # superimpose true parameter value
  abline(v = theta0[ii], lwd = 4, lty = 2,col='orange')
}

# plot data
# par(mfrow=c(1,1))
par(mar = c(4, 4, 1, 0)+.1)
clrs <- c("orange","darkgreen")
plot(x = 0, type = "n", xlim = range(tseq), ylim = c(0,max(Xobs)*2),
     xlab = "Time (years)", ylab = "Population")
lines(out$time, out$state[,1], type = "p", pch = 19, col = clrs[1],cex=.5)
lines(out$time, out$state[,2], type = "p", pch = 19, col = clrs[2],cex=.5)

x0 <- Xobs[1,]

for(i in 1:100) {
  theta0 <- lvpost$params[i,]
  lvsim <- sde.sim(model = lvmod, x0 = x0, theta = theta0,
                   nobs = N-1, # N-1 steps forward
                   dt = dT,
                   dt.sim = dT/100) # internal observation time
  
  # plot data
  Xsamp <- rbind(c(x0), lvsim$data) # include first observation
  lines(tseq, Xsamp[,"H"],lwd=0.1,lty=2, col = clrs[1])
  lines(tseq, Xsamp[,"L"],lwd=0.1,lty=2, col = clrs[2])
}

legend("topleft", legend = c("Hare", "Lynx"), fill = clrs)


plot(x = 0, type = "n", xlim=c(0,max(Xobs[,1])*1.5),ylim =c(0,max(Xobs[,2])*1.5),
     xlab = "Hare", ylab = "Lynx")

x0 <- Xobs[1,]

for(i in 1:100) {
  theta0 <- lvpost$params[i,]
  lvsim <- sde.sim(model = lvmod, x0 = x0, theta = theta0,
                   nobs = N-1, # N-1 steps forward
                   dt = dT,
                   dt.sim = dT/100) # internal observation time
  
  # plot data
  Xsamp <- rbind(c(x0), lvsim$data) # include first observation
  lines(Xsamp[,"H"],Xsamp[,"L"],lwd=0.2,lty=2, col = 'darkgray')
  
}
lines(out$state[,1],out$state[,2], type = "p", pch = 19, col = clrs[1],cex=.5)
