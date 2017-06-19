## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
#devtools::install_github("dpique/oncoMix1")
library(oncoMix)

## ------------------------------------------------------------------------
data(dfTumor, dfNormal, package="oncoMix")


dfNml = as.data.frame(matrix(data = rgamma(n = 150, shape = 2, rate = 2), nrow = 15, ncol = 10))
rownames(dfNml) = paste0("patient.n", 1:nrow(dfNml))
colnames(dfNml) = paste0("gene", 1:ncol(dfNml))
dfTumor = as.data.frame(matrix(data = rgamma(n = 150, shape = 4, rate = 3), nrow = 15, ncol = 10))
rownames(dfTumor) = paste0("patient.t", 1:nrow(dfTumor))
colnames(dfTumor) = paste0("gene", 1:ncol(dfTumor))
mmParams = mixModelParams(dfNml, dfTumor)



## ------------------------------------------------------------------------
#oncoMix::

## ----cars----------------------------------------------------------------
summary(cars)

## ----pressure, echo=FALSE------------------------------------------------
plot(pressure)

