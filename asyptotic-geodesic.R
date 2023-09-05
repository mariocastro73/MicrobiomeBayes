library(readr)
geoxs <- read_csv("~/Documentos/Dropbox/projects/galeano/MBAM/MBAM-fork/geoxs.csv",
col_names = FALSE)
library(readr)
geots <- read_csv("~/Documentos/Dropbox/projects/galeano/MBAM/MBAM-fork/geots.csv",
col_names = FALSE)

matplot(geots[geots>1],geoxs[geots>1,],type='l')

xs <- as.data.frame(geoxs)
x11()
par(mfrow=c(6,11))
for(i in 1:11) {
  for(j in (i+1):12) {
    plot(abs(xs[geots>1,i]),abs(xs[geots>1,j]),type='l',xlab=sprintf("x%s",i),ylab=sprintf("x%s",j),log='xy')
  }
}

par(mfrow=c(2,2))
i=1
j=4
plot(xs[geots>1,i],xs[geots>1,j],type='l',xlab=sprintf("x%s",i-1),ylab=sprintf("x%s",j-1))
abline(fit <- lm(xs[geots>1,j]~xs[geots>1,i]-1),col=2)
title(sprintf("slope=%s",coef(fit)))

i=1
j=7
plot(xs[geots>1,i],xs[geots>1,j],type='l',xlab=sprintf("x%s",i-1),ylab=sprintf("x%s",j-1))
abline(fit <- lm(xs[geots>1,j]~xs[geots>1,i]-1),col=2)
title(sprintf("slope=%s",coef(fit)))

i=4
j=7
plot(xs[geots>1,i],xs[geots>1,j],type='l',xlab=sprintf("x%s",i-1),ylab=sprintf("x%s",j-1))
abline(fit <- lm(xs[geots>1,j]~xs[geots>1,i]-1),col=2)
title(sprintf("slope=%s",coef(fit)))

i=4
j=10
plot(xs[geots>1,i],xs[geots>1,j],type='l',xlab=sprintf("x%s",i-1),ylab=sprintf("x%s",j-1))
abline(fit <- lm(xs[geots>1,j]~xs[geots>1,i]-1),col=2)
title(sprintf("slope=%s",coef(fit)))

par(mfrow=c(1,1))
plot(geots[geots>0],xs[geots>0,4],type='l',xlab=expression(tau),ylab=expression(x[3]))
