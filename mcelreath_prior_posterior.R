# prior - likelihood conflict
# install.packages(c("coda","mvtnorm","devtools","loo","dagitty","shape"))
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# devtools::install_github("rmcelreath/rethinking")
library(cmdstanr)
library(rethinking)

yobs <- 0

mtt <- ulam(
  alist(
    y ~ dstudent(2,mu,1),
    mu ~ dstudent(2,10,1)
  ), data=list(y=yobs) , chains=4 , iter=2000 )

mnn <- ulam(
  alist(
    y ~ dnorm(mu,1),
    mu ~ dnorm(10,1)
  ), data=list(y=yobs) , chains=4 , iter=2000)
mtn <- ulam(
  alist(
    y ~ dstudent(2,mu,1),
    mu ~ dnorm(10,1)
  ), data=list(y=yobs) , chains=4 , iter=2000)
mnt <- ulam(
  alist(
    y ~ dnorm(mu,1),
    mu ~ dstudent(2,10,1)
  ), data=list(y=yobs) , chains=4 , iter=2000)

# plot
par(mfrow=c(2,2),cex=1.05)
ymax <- 0.53
xlwd <- 1.5
postcol <- 2
xadj <- 0.8

p <- extract.samples(mnn)
dens(p$mu, xlim=c(-5,15), ylim=c(0,ymax), lwd=xlwd+1 , col=postcol, xlab="" ,adj=xadj )
#mtext("normal prior, normal likelihood")
curve( dnorm(yobs,x,1) , add=TRUE , lty=1 , lwd=xlwd ) # lik
curve( dnorm(x,10,1) , add=TRUE , lty=2 , lwd=xlwd ) # prior

text(0,0.42,"likelihood")
text(10,0.42,"prior")

p <- extract.samples(mtt)
dens(p$mu , xlim=c(-5,15) , ylim=c(0,ymax) , lwd=xlwd+1 , col=postcol , xlab="" ,adj=xadj )
#mtext("t prior, t likelihood")
curve( dstudent(yobs,2,x,1) , add=TRUE , lty=1 , lwd=xlwd ) # lik
curve( dstudent(x,2,10,1) , add=TRUE , lty=2 , lwd=xlwd ) # prior

text(0,0.42,"likelihood")
text(10,0.42,"prior")

p <- extract.samples(mnt)
dens(p$mu, xlim=c(-5,15), ylim=c(0,ymax), lwd=xlwd+1 , col=postcol, xlab="" ,adj=xadj )
#mtext("t prior, normal likelihood")
curve( dnorm(yobs,x,1) , add=TRUE , lty=1 , lwd=xlwd) # lik
curve( dstudent(x,2,10,1) , add=TRUE , lty=2 , lwd=xlwd) # prior

text(0,0.42,"likelihood")
text(10,0.42,"prior")

p <- extract.samples(mtn)
dens(p$mu, xlim=c(-5,15), ylim=c(0,ymax), lwd=xlwd+1 , col=postcol, xlab="" ,adj=xadj )
#mtext("normal prior, t likelihood")
curve( dstudent(yobs,2,x,1) , add=TRUE , lty=1 , lwd=xlwd) # lik
curve( dnorm(x,10,1) , add=TRUE , lty=2 , lwd=xlwd) # prior

text(0,0.42,"likelihood")
text(10,0.42,"prior")