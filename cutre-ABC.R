library(deSolve)



population <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dx1 <- x1*(theta0 + theta2 * x1 + theta4 * x2)
    dx2 <-  x2*(theta1 + theta3 * x1 + theta5 * x2)
    list(c(dx1,dx2))
  })
}

params <- c(theta1=2.5/45, theta2=2/45)
params <- c(theta0=1, theta1=-1, theta2=0,theta3=0.02,theta4=-0.01,theta5=0)
init.cond <- c(x1=20,x2=20)

tmax=15 # tiempo de simulación
t <- seq(0,tmax,length=15) 

# Simulamos la ODE
sim.ode <- function(params) {
  return(data.frame(ode(y = init.cond, times = t, func = population, parms = params)))
}

data<- sim.ode(params)
out <- data
out$x1 <- rlnorm(length(data$time),log(data$x1),.2)
out$x2 <- rlnorm(length(data$time),log(data$x2),.2)
matplot(data$time,data[,-1],type='l')
matplot(out$time,out[,-1],pch=19,add=TRUE)

loss <- function(data,params,tol=1e-2) {
  return(mean(colMeans((data[,-1]-out[,-1])^2/data[,-1]^2)))
}
print(loss(data,out))
