library(deSolve)
library(msde)

set.seed(1234) # semilla
# Modelo deSolve
population <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dx1 <- x1*(r1+b11*x1+b12*x2+b13*x3)
    dx2 <- x2*(r2+b21*x1+b22*x2+b23*x3)
    dx3 <- x3*(r3+b31*x1+b32*x2+b33*x3)
    
    list(c(dx1,dx2,dx3))
  })
}
df <- read.csv('fourpop.csv',header=F)[,-2]
colnames(df) <- c("x1","x2","x3")
matplot(df,type='l')
t <- 1:length(df[,1])
df <- cbind(t,df)
par(mar=c(5,5,2,2))
matplot(t,df[,-1], type = "o",pch=19,cex=1.2,lwd=2,cex.axis=1.5,cex.lab=2,xlab="Time",ylab="Population") # Plot con ruido
legend('topleft',legend=c("x1","x2","x3"),col=1:3,lty=1:3,lwd=2,pch=19,cex=1.5)

# Creamos modelo "msde"
data.names <- c("x1", "x2", "x3")
param.names <- c(names(params),"sigma")
LK3 <- sde.make.model(ModelFile = "3species-safe.h", # La madre del cordero
                      data.names = data.names,
                      param.names = param.names)

# dT <- .25 # intervalo de muestreo
dT <- 1
m <- 1
# 
# Xobs <- as.matrix(df)
Xobs <- as.matrix(df[seq(1,length(t),by=m),-1])
Xobs[Xobs<0] <- 0
init <- sde.init(model = LK3, x = Xobs, dt = dT,
                 theta=rep(.1,13)) 

# Ajuste bayesiano (usando msde)
nsamples <- 5e5
burn <- 2e3

LK3.posterior <- sde.post(model = LK3, init = init,
                          hyper = NULL, #prior plano
                          nsamples = nsamples, burn = burn)
# tnames <- expression(theta[1], theta[2], theta[3], 
#                      theta[4], theta[5], theta[6],
#                      theta[7], theta[8], theta[9],
#                      theta[10], theta[11], theta[12],
#                      sigma)
tnames <- expression(r[1] ,r[2] ,r[3],
                     b[11],b[21],b[31], 
                     b[12],b[22],b[32],
                     b[13],b[23],b[33],
                     sigma)
theta0 <- c(params,sigma=sigma)

x11("",18,10)
par(mfrow=c(3,5))
par(mar=c(5.5,5.5,3,5.5))
col <- "#ffa50088"
cex.size <- 2
for(ii in 1:LK3$nparams) {
  hist(LK3.posterior$params[,ii],  freq = FALSE,border=col,col=col,cex.main=cex.size,cex=cex.size,cex.lab=cex.size, cex.axis=cex.size,#breaks = 25,
       xlab = tnames[ii],xlim=c(min(theta0[ii],0,LK3.posterior$params[,ii]),max(theta0[ii],0,LK3.posterior$params[,ii])),
       main = parse(text = paste0("p(", tnames[ii], "*\" | \"*bold(Data))")))
  # superimpose true parameter value
  abline(v = theta0[ii], lwd = 2, lty = 2,col=1)
  abline(v = 0, lwd = 2, lty = 3,col='red')
}
dev.copy2pdf(file="posteriors-fourpop.pdf")

colq <- function(X,q) apply(X, 2, function(X) as.numeric(quantile(X,q))) # Función auxiliar (quantiles)
compara <- as.data.frame(cbind(real=theta0,median=colq(LK3.posterior$params,.5),
                               lowq=colq(LK3.posterior$params,.025),hiq=colq(LK3.posterior$params,.975)))
print(compara)
# for(i in 1:12) cat(tnames[i],params[i],":",length(which(LK3.posterior$params[,i]*sign(params[i])<0))/length(LK3.posterior$params[,i])*100,"\n")

for(i in 1:12) cat(sprintf("%s & %+.4f & %+.4f& %.1f%% & %.1f%%\\\\\n",tnames[i],params[i],median(LK3.posterior$params[,i]),
                           (params[i]-median(LK3.posterior$params[,i]))/params[i]*100,length(which(LK3.posterior$params[,i]*sign(params[i])<0))/length(LK3.posterior$params[,i])*100))


# Pintamos unas cuantas trayectorias
x11("",8,6)
par(mfrow=c(1,1))
par(mar=c(5,5,3,5))
with(simulacion,matplot(t, df[,-1], type='n',lty=1:3,cex=1.4,cex.lab=1.5,cex.axis=1.5, 
                        xlab = "time", ylab = "Populations",ylim=c(0,1)))

matplot(t[seq(1,length(t),by=m)],Xobs, type = "p",add=TRUE,pch=19,lwd=0) # Plot con ruido

ntraj <- 5
rango <- seq(nsamples-ntraj,nsamples)
for(i in rango) {
  traj <- data.frame(ode(y = init.cond, times = t, 
                         func = population, parms = LK3.posterior$params[i,1:12]))
  with(traj,matplot(t, cbind(x1,x2,x3), type='l',lwd=.1,lty=3,add=TRUE))
  
}
dev.copy2pdf(file='trajectories_thick_fourpop.pdf')

