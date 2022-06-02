library(deSolve)
library(msde)

set.seed(123) # semilla
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

init.cond <- c(x1=10,x2=10,x3=10) # Cond. inciales
params <- c(r1 = 0.5, r2 = 0.19, r3 = 0.9, 
            b11=-0.4, b22=-0.05, b33=-0.02, 
            b12= 0.00489, b21 = 0, b13 = 0.0049,
            b31 = 0.0025, b23 = 0.0057, b32 = 0.0025) # Parametros

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
param.names <- c("r1","r2","r3","b11","b12","b13","b21","b22","b23","b31","b32","b33","sigma")
LK3 <- sde.make.model(ModelFile = "3species.h", # La madre del cordero
                        data.names = data.names,
                        param.names = param.names)

# dT <- .25 # intervalo de muestreo
dT <- 1
# Xobs <- as.matrix(df[,-1])
Xobs <- as.matrix(df[seq(1,61,by=4),-1])
init <- sde.init(model = LK3, x = Xobs, dt = dT,
                 # m=1,theta=rep(.1,13)) # No nos mojamos, todos los "priors" iguales
                 m=4,theta=rep(.1,13)) # No nos mojamos, todos los "priors" iguales

# Ajuste bayesiano (usando msde)
nsamples <- 2e4
burn <- 2e3


LK3.posterior <- sde.post(model = LK3, init = init,
                   hyper = NULL, #prior plano
                   nsamples = nsamples, burn = burn)
tnames <- expression(r[1], r[2], r[3], 
                     b[11], b[22], b[33], 
                     b[12], b[21], b[13], 
                     b[31], b[23], b[32],
                     sigma)
theta0 <- params

par(mfrow=c(3,5))

for(ii in 1:LK3$nparams) {
  hist(LK3.posterior$params[,ii], breaks = 25, freq = FALSE,
       xlab = tnames[ii],
       main = parse(text = paste0("p[1](", tnames[ii], "*\" | \"*bold(Data))")))
  # superimpose true parameter value
  abline(v = theta0[ii], lwd = 4, lty = 2,col='orange')
}

colq <- function(X,q) apply(X, 2, function(X) as.numeric(quantile(X,q))) # Función auxiliar (quantiles)
compara <- as.data.frame(cbind(real=c(theta0,sigma=sigma),median=colq(LK3.posterior$params,.5),
                    lowq=colq(LK3.posterior$params,.025),hiq=colq(LK3.posterior$params,.975)))
# View(compara)
print(compara)

par(mfrow=c(1,1))
k <- sort(compara$real,index.return=TRUE)$ix
with(compara[k,],{
     plot(real,median);
     # abline(0,1);
     abline(h=0);
     abline(v=0);
     lines(real,lowq);
     lines(real,hiq)}
)

# Pintamos unas cuantas trayectorias
par(mfrow=c(1,1))
with(simulacion,matplot(time, cbind(x1,x2,x3), type='l',lwd=1.5,lty=1:3, 
                        xlab = "time", ylab = "Populations",ylim=c(0,1.3*max(x3))))

matplot(df[,1], df[,-1], type = "p",add=TRUE,pch=19,cex=.7) # Plot con ruido

ntraj <- 1000
rango <- seq(nsamples-ntraj,nsamples)
for(i in rango) {
    traj <- data.frame(ode(y = init.cond, times = t, 
                           func = population, parms = LK3.posterior$params[i,1:12]))
    with(traj,matplot(time, cbind(x1,x2,x3), type='l',lwd=.1,lty=3,add=TRUE))
    
}
