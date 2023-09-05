#install.packages("remotes")
# remotes::install_github("jpahle/CoRC")
library(CoRC)

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


# CoRC::loadModel('lotka.antimony___main_sbml.cps')
# CoRC::loadModel('lotka-basic.cps')
res <- runTimeCourse(15)$result

tmax=15 # tiempo de simulación
t <- seq(0,tmax,by=.25) 

output <- list()
matplot(res$Time,res[,-1],type='l')
for(i in 1:10)  {
  print(i)
  out <- runParameterEstimation()
  if(out$main$root_mean_square<0.56) {
    output[[i]] <- out
    params <- out$parameters$value
    names(params) <- c("theta10","theta11","theta12","theta1","theta2","theta3",
                       "theta4","theta5","theta6","theta7","theta8","theta9")
    simulacion <- data.frame(ode(y = init.cond, times = t, func = population, parms = params))
    with(simulacion,matplot(time, cbind(x1,x2,x3), type='p',pch="+",add=TRUE,lwd=1.5,lty=1:3, 
                            xlab = "time", ylab = "Populations",ylim=c(0,1.3*max(x3))))
    
  }
  print(out$main$root_mean_square)
}
# out <- runParameterEstimation(randomize_start_values = TRUE,method = "GeneticAlgorithm",
                              # calculate_statistics = TRUE)



# saveRDS(output,'batch.rda')
