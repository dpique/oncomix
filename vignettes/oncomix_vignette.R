## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
#devtools::install_github("dpique/oncomix")
#setwd("\\\\data.einstein.yu.edu/home/dpique/dpLabNotebook/bimodality_brca_tcga")
#library(devtools)
#install_github("dpique/oncomix")
#install.packages("oncomix")
#vignette("oncomix_vignette")

#remove.packages("oncomix")
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

mmParams = mixModelParams(dfNmlIsof, dfTumorIsof) #will take a few mins
head(mmParams)

## ---- fig.width = 7, fig.height= 6.5-------------------------------------
scatterMixPlot(mmParams)

