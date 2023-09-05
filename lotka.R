library(deSolve)
library(msde)

set.seed(12) # semilla
# Modelo deSolve
lotka <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dx0 <- x0*(theta0+theta3*x0+theta6*x1+theta9*x2)
    dx1 <- x1*(theta1+theta4*x0+theta7*x1+theta10*x2)
    dx2 <- x2*(theta2+theta5*x0+theta8*x1+theta11*x2)
    
    list(c(dx0,dx1,dx2))
    
  })
}

init.cond <- c(x0=10,x1=50,x2=10) # Cond. inciales
params <- c(theta0 = 0.5,     theta1 = 0.2,    theta2 = 0.9, 
            theta3= -0.4,      theta4 = 0,       theta5 = 0.003, 
            theta6= 0.005,  theta7=-0.05,     theta8 = 0.002,
            theta9 = 0.005, theta10 = 0.006, theta11= -0.02) 
# x = np.array([0.5,0.2,0.9,-0.4,0.,0.003,0.005,-0.05,0.002,0.005,0.006,-0.02])


tmax=15 # tiempo de simulación
dt <- 0.25
t <- seq(0,tmax,by=dt) 

# Simulamos la ODE
simulacion <- data.frame(ode(y = init.cond, times = t, func = lotka, parms = params))
par(mfrow=c(1,1))
with(simulacion,matplot(time, cbind(x0,x1,x2), add=F,type='l',lwd=1.5,lty=1:3, 
                        xlab = "time", ylab = "lotkas",ylim=c(0,1.1*max(x0,x1,x2))))



################################3
params2 = c(theta0=-1.99258654,theta1=0.16074929,theta2=0.95160099,
          theta3=-5.14265934,theta4=-0.07376196,theta5=0.10222699,
          theta6=0.88975381,theta7=-0.03621809,theta8=-0.01650715,
          theta9=0.07567895,theta10=0.0071044,theta11=-0.02146429)

simulacion2 <- data.frame(ode(y = init.cond, times = t, func = lotka, parms = params2))
# par(mfrow=c(1,1))
with(simulacion2,matplot(time, cbind(x0,x1,x2), add=TRUE,type='p',lwd=1.5,lty=1:3, 
                        xlab = "time", ylab = "lotkas",ylim=c(0,1.1*max(x0,x1,x2))))


matplot(simulacion$time,simulacion2[,-1]-simulacion[,-1],type='o')

RSS <- function(parameters) {
  names(parameters) <- c("theta0","theta1","theta2",
                         "theta3","theta4","theta5",
                         "theta6","theta7","theta8",
                         "theta9","theta10","theta11")
  simulacion2 <- ode(y = init.cond, times = t, func = lotka, parms = parameters)
  sum(colSums((simulacion[,-1]-simulacion2[,-1])^2))
}
RSS(params2)
Opt <- optim(params2*1.1, RSS, method = "L-BFGS-B", lower = rep(-1,12), upper = rep(1,12)) # optimize with some sensible conditions


params3 <- c(theta0=-1.44830118e+02,theta1=-1.95469140e+00,theta2=3.41887195e+01,
             theta3=-2.62791543e+02,theta4=-3.75997031e+00,theta5=6.01871605e+01,
             theta6=4.89074479e+01,theta7=6.52921379e-01,theta8=-1.12133794e+01,
             theta9=4.10153448e+00,theta10=6.57906691e-02,theta11=-9.58888522e-01)

simulacion3 <- data.frame(ode(y = init.cond, times = t, func = lotka, parms = params3))
# par(mfrow=c(1,1))
with(simulacion2,matplot(time, cbind(x0,x1,x2), add=TRUE,type='p',lwd=1.5,lty=1:3, 
                         xlab = "time", ylab = "lotkas",ylim=c(0,1.1*max(x0,x1,x2))))


lotkaRG <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dx1 <- x1*(theta1+theta7*x1+theta10*x2)
    dx2 <- x2*(theta2+theta8*x1+theta11*x2)
    
    list(c(dx1,dx2))
    
  })
}

init.cond <- c(x1=50,x2=10) # Cond. inciales
paramsRG <- c(theta1=0.1175097730560013, 
               theta2=1.218270185731395, theta7=-0.046836877839611, 
              theta8=-0.0121041123752068, 
               theta10=0.00710671110732864, theta11=-0.023951387079958433)


simulacionRG <- data.frame(ode(y = init.cond, times = t, func = lotkaRG, parms = paramsRG))

with(simulacionRG,matplot(time, cbind(-0.551122 + 0.186107*x1 + 0.0156076*x2,x1,x2), add=T,type='p',lwd=1.5,lty=1:3, 
                         xlab = "time", ylab = "lotkas",ylim=c(0,1.1*max(x1,x2))))

