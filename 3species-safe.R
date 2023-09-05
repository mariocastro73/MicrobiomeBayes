library(deSolve)
library(msde)

set.seed(123) # semilla
# Modelo deSolve
population <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dx1 <- x1*(theta1+theta4*x1+theta7*x2+theta10*x3)
    dx2 <- x2*(theta2+theta5*x1+theta8*x2+theta11*x3)
    dx3 <- x3*(theta3+theta6*x1+theta9*x2+theta12*x3)
    
    list(c(dx1,dx2,dx3))
  })
}

init.cond <- c(x1=10,x2=10,x3=10) # Cond. inciales
params <- c(theta1 = 0.5,     theta2 = 0.2,    theta3 = 0.9, 
            theta4=-0.4,      theta5 = 0,       theta6 = 0.003, 
            theta7= 0.005,  theta8=-0.05,     theta9 = 0.002,
            theta10 = 0.005, theta11 = 0.006, theta12=-0.02) # Parametros

tmax=15 # tiempo de simulación
t <- seq(0,tmax,by=.25) 

# Simulamos la ODE
simulacion <- data.frame(ode(y = init.cond, times = t, func = population, parms = params))

# A pintar
par(mfrow=c(1,1))
with(simulacion,matplot(time, cbind(x1,x2,x3), type='l',lwd=1.5,lty=1:3, 
                        xlab = "time", ylab = "Populations",ylim=c(0,1.3*max(x3))))

# Añadimos ruido (gaussiano blanco)
sigma <- 1
df <- as.data.frame(cbind(simulacion$time,sapply(simulacion[,2:4],function(x) x+ rnorm(length(t),0,sigma))))
df[df<0] <- 0 # Quitamos los valores negativos
colnames(df) <- colnames(simulacion)

matplot(df[,1], df[,-1], type = "p",add=TRUE,pch=19,cex=.7) # Plot con ruido
legend('topright',legend=c("x1","x2","x3"),col=1:3,lty=1:3,lwd=2,pch=19)
# rm(.Random.seed)

# Creamos modelo "msde"
data.names <- c("x1", "x2", "x3")
param.names <- c(names(params),"sigma")
# LK3 <- sde.make.model(ModelFile = "3species-safe.h", # La madre del cordero
#                       data.names = data.names,
#                       param.names = param.names)

# dT <- .25 # intervalo de muestreo
dT <- 1
# Xobs <- as.matrix(df[,-1])
Xobs <- as.matrix(df[seq(1,61,by=4),-1])
init <- sde.init(model = LK3, x = Xobs, dt = dT,
                 # m=1,theta=rep(.1,13)) # No nos mojamos, todos los "priors" iguales
                 m=4,theta=rep(.01,13)) # No nos mojamos, todos los "priors" iguales

# Ajuste bayesiano (usando msde)
nsamples <- 2e5
burn <- 1e4

LK3.posterior <- sde.post(model = LK3, init = init,
                          hyper = NULL, #prior plano
                          nsamples = nsamples, burn = burn)
tnames <- expression(theta[1], theta[2], theta[3], 
                     theta[4], theta[5], theta[6],
                     theta[7], theta[8], theta[9],
                     theta[10], theta[11], theta[12],
                     sigma)
theta0 <- params

# x11()
par(mfrow=c(3,5))

for(ii in 1:LK3$nparams) {
  plot(density(LK3.posterior$params[,ii]), breaks = 25, freq = FALSE,
       xlab = tnames[ii],
       main = parse(text = paste0("p[1](", tnames[ii], "*\" | \"*bold(Data))"))) 
  # hist(LK3.posterior$params[,ii], breaks = 25, freq = FALSE,
       # xlab = tnames[ii],
       # main = parse(text = paste0("p[1](", tnames[ii], "*\" | \"*bold(Data))")))
  # superimpose true parameter value
  abline(v = theta0[ii], lwd = 4, lty = 2,col='orange')
  abline(v = 0, lwd = 2, lty = 3,col='red')
}

for(i in 1:length(theta0)) {
  print(mean(100*sapply(LK3.posterior$params[,i],function(x) x*sign(theta0[i])<0)))
}

colq <- function(X,q) apply(X, 2, function(X) as.numeric(quantile(X,q))) # Función auxiliar (quantiles)
compara <- as.data.frame(cbind(real=c(theta0,sigma=sigma),median=colq(LK3.posterior$params,.5),
                               lowq=colq(LK3.posterior$params,.025),hiq=colq(LK3.posterior$params,.975)))
print(compara)

# Pintamos unas cuantas trayectorias
par(mfrow=c(1,1))

with(simulacion,matplot(time, cbind(x1,x2,x3), type='n',lwd=1.5,lty=1:3, 
                        xlab = "time", ylab = "Populations",ylim=c(0,1.3*max(x3))))

matplot(t[seq(1,61,by=4)],Xobs, type = "p",add=TRUE,pch=19,cex=.7) # Plot con ruido

ntraj <- 1000
rango <- seq(nsamples-ntraj,nsamples)
for(i in rango) {
  traj <- data.frame(ode(y = init.cond, times = t, 
                         func = population, parms = LK3.posterior$params[i,1:12]))
  with(traj,matplot(time, cbind(x1,x2,x3), type='l',lwd=.1,lty=3,add=TRUE))
  
}

