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

## ------------------------------------------------------------------------
topGeneQuant = oncomix::topGeneQuants(mmParams, deltMu2Thresh = 99, deltMu1Thresh = 10, siThresh = .99)
print(topGeneQuant)

## ------------------------------------------------------------------------
mmParams.top10 = mmParams[1:10,]
print(mmParams.top10)

## ------------------------------------------------------------------------
isof = "uc002jxc.2"
plotGeneHist(mmParams, dfNmlIsof, dfTumorIsof, isof)

## ---- fig.width = 7, fig.height= 6.5-------------------------------------
scatterMixPlot(mmParams)

## ---- fig.width = 7, fig.height= 6.5-------------------------------------
scatterMixPlot(mmParams, selIndThresh = .99)

## ------------------------------------------------------------------------
scatterMixPlot(mmParams = mmParams, gene_labels = rownames(mmParams.top10))

## ------------------------------------------------------------------------
#install.packages("RMySQL")
library(RMySQL)

#read in a table of known human oncogenes from the ONGene database
ongene = read.table("http://ongene.bioinfo-minzhao.org/ongene_human.txt", header = T, sep = "\t")

#send a sql query to UCSC to map the human oncogenes to ucsc isoform ids
ucsc_genome <- dbConnect(MySQL(), user="genome", 
                         host="genome-mysql.cse.ucsc.edu", db = 'hg19')
createGeneQuery = function(name){ #name is a character vector
  p1 = paste(name, collapse = '\',\'')
  p2 = paste('(\'',p1, '\')',sep ="")
  return(p2)
}
gene_q = createGeneQuery(ongene$OncogeneName)
query_res = dbGetQuery(ucsc_genome, paste0("SELECT kgID, geneSymbol FROM kgXref WHERE geneSymbol IN ", gene_q, " ;"))
dbDisconnect(ucsc_genome)

#Merge the query_res & mmParams dataframes
query_res$kgID.s = substr(query_res$kgID, 1,8)
mmParams$kgID.s = substr(rownames(mmParams), 1,8)
mmParams$kgID = rownames(mmParams)
mmParams.m = merge(mmParams, query_res, by = "kgID.s", all.x = T)
rownames(mmParams.m) = mmParams.m$kgID.x

# Show the top 5 isoforms with the highest score 
#in our dataset that map to known oncogenes
mmParams.m = mmParams.m[with(mmParams.m, order(-score)), ]
mmParams.m.s = subset(mmParams.m, !is.na(geneSymbol))[1:5,]
print(mmParams.m.s)

scatterMixPlot(mmParams = mmParams, gene_labels = rownames(mmParams.m.s))

## ------------------------------------------------------------------------
library(ggplot2)
library(RColorBrewer)
col = brewer.pal(3, "Dark2")
ggplot(mmParams.m, aes(x = score, y = ..density.., fill=is.na(geneSymbol))) +
  geom_histogram(data=subset(mmParams.m, is.na(geneSymbol)), fill = col[2], alpha = 0.5)+ 
  geom_histogram(data=subset(mmParams.m, !is.na(geneSymbol)), fill = col[3], alpha = 0.5)+ theme_classic() + xlab("OncoMix Score")+ theme_classic() 

## ------------------------------------------------------------------------
sessionInfo()

