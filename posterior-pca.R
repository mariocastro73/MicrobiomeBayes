library(FactoMineR)
library(Factoshiny)
aux <- LK3.posterior$params[-(1:10000),]
aux <- LK3.posterior$params
res.shiny <- PCAshiny(aux)

