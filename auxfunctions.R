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
