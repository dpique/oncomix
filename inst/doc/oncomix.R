## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = F, message = F)

## ------------------------------------------------------------------------
#devtools::install_github("dpique/oncomix")
library(oncomix)

## ------------------------------------------------------------------------
library(ggplot2)
oncomix::oncoMixIdeal()


## ------------------------------------------------------------------------
data(dfNmlIsof, dfTumorIsof, package="oncomix")

dim(dfNmlIsof)
dfNmlIsof[1:5, 1:5] #take a look at the matrix of mRNA expression data from adjacent normal samples

dim(dfTumorIsof)
dfTumorIsof[1:5, 1:5] #take a look at the matrix of mRNA expression data from tumors

mmParams = oncomix::mixModelParams(dfNmlIsof, dfTumorIsof) #fits the mixture models, will take a few mins
head(mmParams)

## ---- fig.width = 7, fig.height= 6.5-------------------------------------
scatterMixPlot(mmParams)

## ---- fig.width = 7, fig.height= 6.5-------------------------------------
scatterMixPlot(mmParams, selIndThresh = .99)

## ------------------------------------------------------------------------
library(ggplot2)
qplot(mmParams[,"SI"]) + theme_classic() + xlab("Selectivity Index")


## ------------------------------------------------------------------------
mmParams.df = as.data.frame(mmParams)
topGeneQuant = oncomix::topGeneQuants(mmParams.df, deltMu2Thresh = 90, deltMu1Thresh = 10, siThresh = .99)
print(topGeneQuant)
topGeneTbl = oncomix::topGeneTable(mmParams.df, N = 10) #want the top 10 isoforms based on a custom score
print(topGeneTbl)


## ------------------------------------------------------------------------
isof = "uc002jxc.2"
plotGeneHist(mmParams, dfNmlIsof, dfTumorIsof, isof)


## ------------------------------------------------------------------------
sessionInfo()

