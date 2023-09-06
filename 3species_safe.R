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

init.cond <- c(x1=10,x2=50,x3=10) # Cond. inciales
params <- c(r1 = 0.5,     r2 = 0.019,    r3 = 0.9, 
            b11= -0.02,      b21 = -0.001,       b31 = 0.0015, 
            b12= 0.0089,  b22=-0.05,     b32 = 0.025,
            b13 = 0.0049, b23 = 0.0057, b33= -0.02) # Parametros

tmax=15 # tiempo de simulación
dt <- 0.25
t <- seq(0,tmax,by=dt) 

# Simulamos la ODE
simulacion <- data.frame(ode(y = init.cond, times = t, func = population, parms = params))
par(mfrow=c(1,1))
with(simulacion,matplot(time, cbind(x1,x2,x3), type='l',lwd=1.5,lty=1:3, 
                        xlab = "time", ylab = "Populations",ylim=c(0,1.3*max(x3))))

# Añadimos ruido (gaussiano blanco)
sigma <- .5
df <- as.data.frame(cbind(simulacion$time,sapply(simulacion[,2:4],function(x) x+ rnorm(length(t),0,sigma))))
# df <- cbind(simulacion$time,sapply(simulacion[,2:4],function(x) rlnorm(length(t),log(x),sigma)))
df[df<0] <- 0 # Quitamos los valores negativos
colnames(df) <- colnames(simulacion)

matplot(df[,1], df[,-1], type = "p",add=TRUE,pch=19,cex=.7) # Plot con ruido
legend('topright',legend=c("x1","x2","x3"),col=1:3,lty=1:3,lwd=2,pch=19)

# Creamos modelo "msde"
data.names <- c("x1", "x2", "x3")
param.names <- c(names(params),"sigma")
LK3 <- sde.make.model(ModelFile = "3species-safe.h", # La madre del cordero
                      data.names = data.names,
                      param.names = param.names)

# dT <- .25 # intervalo de muestreo
dT <- 1
m <- 4
# Xobs <- as.matrix(df[,-1])
Xobs <- as.matrix(df[seq(1,length(t),by=m),-1])
init <- sde.init(model = LK3, x = Xobs, dt = dT,
                 # m=1,theta=rep(.1,13)) # No nos mojamos, todos los "priors" iguales
                 m=m,theta=rep(.1,13)) # No nos mojamos, todos los "priors" iguales

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
dev.copy2pdf(file="posteriors.pdf")

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
with(simulacion,matplot(time, cbind(x1,x2,x3), type='n',lty=1:3,cex=1.4,cex.lab=1.5,cex.axis=1.5, 
                        xlab = "time", ylab = "Populations",ylim=c(0,1.3*max(x3))))

matplot(t[seq(1,length(t),by=m)],Xobs, type = "p",add=TRUE,pch=19,lwd=0) # Plot con ruido

ntraj <- 500
rango <- seq(nsamples-ntraj,nsamples)
for(i in rango) {
  traj <- data.frame(ode(y = init.cond, times = t, 
                         func = population, parms = LK3.posterior$params[i,1:12]))
  with(traj,matplot(time, cbind(x1,x2,x3), type='l',lwd=.1,lty=3,add=TRUE))
  
}
dev.copy2pdf(file='trajectories.pdf')

