library(deSolve)
library(msde)

set.seed(123) # semilla
# Modelo deSolve
population <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dx1 <- x1*(r1+b11*x1+b12*x2)
    dx2 <- x2*(r2+b21*x1+b22*x2)
    
    list(c(dx1,dx2))
  })
}

init.cond <- c(x1=10,x2=14) # Cond. inciales
params <- c(r1 = 0.25, r2 = 0.19,
            b11=-0.07, b21=-0.03,
            b12= 0.2, b22= -0.0015) # Parametros

tmax=10 # tiempo de simulación
t <- seq(0,tmax,length=50) 

# Simulamos la ODE
simulacion <- data.frame(ode(y = init.cond, times = t, func = population, parms = params))

# A pintar
par(mfrow=c(1,1))
with(simulacion,matplot(time, cbind(x1,x2), type='l',lwd=1.5,lty=1:3, 
                        xlab = "time", ylab = "Populations",ylim=c(0,1.3*max(x1,x2))))

# Añadimos ruido (gaussiano blanco)
sigma <- 1
df <- as.data.frame(cbind(simulacion$time,sapply(simulacion[,2:3],function(x) x+ rnorm(length(t),0,sigma))))
df[df<0] <- 0 # Quitamos los valores negativos
colnames(df) <- colnames(simulacion)

matplot(df[,1], df[,-1], type = "p",add=TRUE,pch="+",cex=.7) # Plot con ruido
legend('topleft',legend=c("x1","x2"),col=1:2,lty=1:2,lwd=2,pch="+")
# rm(.Random.seed)

# Creamos modelo "msde"
data.names <- c("x1", "x2")
param.names <- c("r1","r2","b11","b21","b12","b22","sigma")
LK2 <- sde.make.model(ModelFile = "2species.h", # La madre del cordero
                        data.names = data.names,
                        param.names = param.names)

# dT <- .25 # intervalo de muestreo
dT <- 1
# Xobs <- as.matrix(df[,-1])
Xobs <- as.matrix(df[seq(1,length(simulacion$time),by=dT),-1])
matplot(t[seq(1,length(simulacion$time),by=dT)],Xobs, type = "p",add=TRUE,pch=0,cex=.7) # Plot con ruido

init <- sde.init(model = LK2, x = Xobs, dt = dT,
                 # m=1,theta=rep(.1,13)) # No nos mojamos, todos los "priors" iguales
                 m=4,theta=rep(.1,7)) # No nos mojamos, todos los "priors" iguales

# Ajuste bayesiano (usando msde)
nsamples <- 2e4
burn <- 2e3


LK2.posterior <- sde.post(model = LK2, init = init,
                   hyper = NULL, #prior plano
                   nsamples = nsamples, burn = burn)
tnames <- expression(r[1], r[2], 
                     b[11], b[21], 
                     b[12], b[22], 
                     sigma)
theta0 <- c(params,sigma=sigma)

par(mfrow=c(3,3))

for(ii in 1:LK2$nparams) {
  hist(LK2.posterior$params[,ii], breaks = 25, freq = FALSE,
       xlab = tnames[ii],
       main = parse(text = paste0("p[1](", tnames[ii], "*\" | \"*bold(Data))")))
  # superimpose true parameter value
  abline(v = theta0[ii], lwd = 4, lty = 2,col='orange')
}

colq <- function(X,q) apply(X, 2, function(X) as.numeric(quantile(X,q))) # Función auxiliar (quantiles)
compara <- as.data.frame(cbind(real=theta0,median=colq(LK2.posterior$params,.5),
                    lowq=colq(LK2.posterior$params,.025),hiq=colq(LK2.posterior$params,.975)))
# View(compara)
print(compara)

# Pintamos unas cuantas trayectorias
par(mfrow=c(1,1))
with(simulacion,matplot(time, cbind(x1,x2), type='l',lwd=1.5,lty=1:3, 
                        xlab = "time", ylab = "Populations",ylim=c(0,1.3*max(x1,x2))))

matplot(t[seq(1,length(simulacion$time),by=dT)],Xobs, type = "p",add=TRUE,pch=19,cex=.7) # Plot con ruido

ntraj <- 100
rango <- seq(nsamples-ntraj,nsamples)
for(i in rango) {
    traj <- data.frame(ode(y = init.cond, times = t, 
                           func = population, parms = LK2.posterior$params[i,1:7]))
    with(traj,matplot(time, cbind(x1,x2), type='l',lwd=.1,lty=3,add=TRUE))
    
}
