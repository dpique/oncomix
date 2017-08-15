#' Generate the parameters for two 2-component Gaussian mixture models with equal variances
#'
#' This function allows you to generate the parameters for two 2-component Gaussian mixture model
#' with equal variances from 2 matrices of data with a priori labels (eg tumor vs normal.) This application
#' was originally intended for matrices of gene expression data treated with 2 conditions.
#'
#' @param dfNml A dataframe of normal data containing patients as rows and genes as columns
#' @param dfTumor A dataframe of tumor data containing patients as rows and genes as columns
#' @keywords oncoMix, mixture-model, two-component
#' @return Returns a list of length 2, each element of which contains the parameters
#' for the Normal or the Tumor data in a 3 x n matrix, where n is the number of patient samples
#' @importFrom mclust Mclust mclustBIC
#' @export
#' @examples
#' dfNml = as.data.frame(matrix(data = rgamma(n = 150, shape = 2, rate = 2), nrow = 15, ncol = 10))
#' rownames(dfNml) = paste0("patient.n", 1:nrow(dfNml))
#' colnames(dfNml) = paste0("gene", 1:ncol(dfNml))
#'
#' dfTumor = as.data.frame(matrix(data = rgamma(n = 150, shape = 4, rate = 3), nrow = 15, ncol = 10))
#' rownames(dfTumor) = paste0("patient.t", 1:nrow(dfTumor))
#' colnames(dfTumor) = paste0("gene", 1:ncol(dfTumor))
#'
#' mmParams = mixModelParams(dfNml, dfTumor)

mixModelParams = function(dfNml, dfTumor) {
  params_normal <- apply(dfNml, 2, function(x) {
    y <- mclust::Mclust(data = x, G=2, modelNames = "E", verbose = F)
    z <- c(n.mu = y$parameters$mean, n.var = y$parameters$variance$sigmasq, n.pi.1 = y$parameters$pro[1])
    #n.class = y$classification
    return(z) })

  params_tumor <- apply(dfTumor, 2, function(x) {
    y <- mclust::Mclust(data = x, G=2, modelNames = "E", verbose = F)
    z <- c(t.mu = y$parameters$mean, t.var = y$parameters$variance$sigmasq, t.pi.1 = y$parameters$pro[1])
    #t.class = y$classification
    return(z) })

  params = rbind(params_normal, params_tumor)

  deltaMu2 = params["t.mu.2",] - params["n.mu.2",]
  deltaMu1 = params["t.mu.1",] - params["n.mu.1",]
  boundaryTumor = (params["t.mu.2",] + params["t.mu.1",]) / 2

  si_calc = function(vectNml, boundTumor){ #calculate the selectivity index
    x = sum(boundTumor > vectNml) / length(vectNml)
    return(x)
  }

  SI = sapply(1:ncol(dfNml), function(i) si_calc(dfNml[,i], boundaryTumor[i]))

  params = rbind(params, deltaMu2, deltaMu1, SI)
  return(t(params))
}


#' Calculate the selectivity Index from getMixModelParams
#'
#' This function allows you to generate the parameters for two 2-component mixture models
#' with equal variances
#'
#' @param mmParams The output from the getMixModelParams function. Will utilize the deltaMu2
#' and deltaMu1 rows
#' @param dfNml The normal dataframe. Will be used to calculate the SI.
#' @keywords oncoMix, visualization, two-component
#' @return Returns a vector, equal in length to the number of genes. 1 SI per gene.
#' @export
#' @examples
#' selectivityIndex(mmParams, dfNml)
#' @seealso \code{\link{mixModelParams}}

selectivityIndex <- function(mmParams, dfNml){
  boundaryTumor <- (mmParams[,"t.mu.2"] + mmParams[,"t.mu.1"]) / 2 #this is the boundary between the
  #tumor samples (classified into hi and low expression)
  selInd = colSums(dfNml < boundaryTumor) / nrow(dfNml)
  #How many nml samples are below this threshold?
  if(colnames(dfNml) == rownames(boundaryTumor)){
    names(selInd) = rownames(mmParams)
    return(selInd)
  } else {
    return("Error! colnames from dfNml and rownames from mmParams do not match")
  }
}


#' Plot the selectivity Index from getMixModelParams
#'
#' This function allows you to generate the parameters for two 2-component mixture models
#' with equal variances
#'
#' @param mmParams The output from the mixModelParams function. Will utilize the deltaMu2
#' and deltaMu1 rows
#' @param dfNml The normal dataframe. Will be used to calculate the SI.
#' @keywords oncoMix, visualization, two-component
#' @return Returns a vector, equal in length to the number of genes. 1 SI per gene.
#' @export
#' @examples
#' selectivityIndex(mmParams, dfNml)
#' @seealso \code{\link{mixModelParams}}

plotSelectivityIndex <- function(mmParams, dfNml){
  boundaryTumor <- (mmParams[,"t.mu.2"] + mmParams[,"t.mu.1"]) / 2 #this is the boundary between the
  #tumor samples (classified into hi and low expression)
  selInd = colSums(dfNml < boundaryTumor) / nrow(dfNml)
  #How many nml samples are below this threshold?
  if(colnames(dfNml) == rownames(boundaryTumor)){
    names(selInd) = rownames(mmParams)
    return(selInd)
  } else {
    return("Error! colnames from dfNml and rownames from mmParams do not match")
  }
}


#' Plot a histogram of gene expression values from tumor and adjacent normal tissue.
#'
#' This function allows you to plot a histogram of gene expression values
#' from tumor and adjacent normal tissue with the option of including the
#' best fitting Gaussian curve.
#'
#' @param mmParams The output from the getMixModelParams function.
#' @param dfNml The normal dataframe.
#' @param dfTumor The tumor dataframe.
#' @param isof The gene isoform to visualize
#' @keywords oncoMix, visualization, Gaussian, two-component
#' @return Returns a histogram of the gene expression values from the two groups.
#' @export
#' @examples
#' plotGeneHist(mmParams, dfNml, dfTumor, isof)
#' @seealso \code{\link{mixModelParams}}

plotGeneHist <- function(mmParams, dfNml, dfTumor, isof){
  tidy_df <- as.data.frame(cbind(as.numeric(c(dfTumor[,isof], dfNml[,isof])),
                  as.factor(c(rep("tumor",nrow(dfTumor)), rep("normal",nrow(dfNml))))),
                  stringsAsFactors = F)
  colnames(tidy_df) <- c("expr", "type")

  p1 = ggplot(tidy_df, aes(x=expr, color=as.factor(type), fill=as.factor(type),
                      group=as.factor(type))) +
    theme_classic() +
    geom_histogram(data=subset(tidy_df,type == 1),fill = "#F8766D",
                   alpha = 0.2, aes(y=..density..)) +
    geom_histogram(data=subset(tidy_df,type == 2),fill = "#00BFC4",
                   alpha = 0.2, aes(y=..density..)) +
    geom_rug(alpha=0.3, show.legend=FALSE) +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),legend.position = "none",
          plot.title = element_text(size = 12),
          axis.text=element_text(size=8), axis.title=element_text(size=8)) +
    ggtitle(paste0(isof, " : SI = ",
                   round(mmParams[isof,"SI"],4)))+
    xlab(expression(Log[2] *"(TPM Reads)")) +
    stat_function(fun = dnorm, colour = "#F8766D",
                  args = list(mean = mmParams[isof,"n.mu.2"],
                              sd=sqrt(mmParams[isof,"n.var"]))) +
    stat_function(fun = dnorm, colour = "#F8766D",
                  args = list(mean = mmParams[isof,"n.mu.1"],
                              sd=sqrt(mmParams[isof,"n.var"]))) +
    stat_function(fun = dnorm, colour = "#00BFC4",
                  args = list(mean = mmParams[isof,"t.mu.2"],
                              sd=sqrt(mmParams[isof,"t.var"]))) +
    stat_function(fun = dnorm, colour = "#00BFC4",
                  args = list(mean = mmParams[isof,"t.mu.1"],
                              sd=sqrt(mmParams[isof,"t.var"])))
  print(p1)
}


#' Generate a scatter plot with the output from mixModelParams
#'
#' This function allows you to generate the parameters for two 2-component mixture models
#' with equal variances
#'
#' @param mmParams The output from the mixModelParams function. Will utilize the deltaMu2
#' and deltaMu1 rows
#' @keywords oncoMix, visualization, two-component
#' @return Returns a ggplot object that you can plot
#' @import ggplot2 ggrepel
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @examples
#' scatterMixPlot(mmParams)
#' @seealso \code{\link{mixModelParams}}

scatterMixPlot <- function(mmParams, selIndThresh = 1, gene_labels = NULL){
  mmParams = as.data.frame(mmParams)
  one_over_alpha = diff(range(mmParams$deltaMu2))
  alpha1 = 1/one_over_alpha

  quants = c(0.01, 0.10, 0.50, 0.90, 0.99) #add in the quantiles
  colors_red=RColorBrewer::brewer.pal(n=length(quants), name="Reds")

  deltaMu2Quant <- quantile(mmParams[,"deltaMu2"], quants)
  deltaMu1Quant <- quantile(1/(abs(mmParams[,"deltaMu1"]) + alpha1), quants)

  x = ggplot(data = as.data.frame(mmParams), aes(x = deltaMu2, y = 1/(abs(deltaMu1) + alpha1))) +
    theme_classic() +
    geom_hline(yintercept = deltaMu1Quant, col = colors_red, size = c(1,1,1,1,1)) +
    geom_vline(xintercept = deltaMu2Quant, col = colors_red, size = c(1,1,1,1,1)) +
    geom_point(alpha= 0.5) +
    xlab(expression(paste(Delta, mu[2]))) +
    ylab(expression(paste(frac(1, paste(Delta, mu[1], " + ", alpha)))))
  #print(x)

  if(selIndThresh < 1){
    mmParams.si = mmParams[mmParams$SI > selIndThresh,]
    x = x + geom_point(data = as.data.frame(mmParams.si),
                       aes(x = deltaMu2, y = 1/(abs(deltaMu1)+alpha1)),
                       size = 10, alpha=0.1,
                       col=colors_red[length(colors_red)],
                       fill=colors_red[length(colors_red)]) +
      ggtitle(bquote(Distribution~of~Mixture~Model~Parameters*","~alpha~"="~.(round(alpha1,2))*", SI >"~.(selIndThresh)))
  } else if(!is.null(gene_labels)){
    mmParams.si = mmParams[gene_labels,]
    mmParams.si$gene_labels = gene_labels
    x = x + geom_point(data = as.data.frame(mmParams.si),
                       aes(x = deltaMu2, y = 1/(abs(deltaMu1)+alpha1)),
                       size = 10, alpha=0.1,
                       col=colors_red[length(colors_red)],
                       fill=colors_red[length(colors_red)]) +
      geom_text_repel(data = mmParams.si, aes(x = deltaMu2, y = 1/(abs(deltaMu1)+alpha1)),
                label = rownames(mmParams.si)) +
      ggtitle(bquote(Distribution~of~Mixture~Model~Parameters*","~alpha~"="~.(round(alpha1,2))))

  } else{
    x = x + ggtitle(bquote(Distribution~of~Mixture~Model~Parameters*","~alpha~"="~.(round(alpha1,2))))
  }
  return(x)
}

#' Identify genes that meet pre-specified quantiles
#'
#' This function allows you to subset genes that are above pre-specified quantiles
#' and that most closely resemble the distribution of oncogenes.
#'
#' @param mmParams The output from the mixModelParams function.
#' @param deltMu2Thresh The percentile threshold for the deltaMu2 statistic.
#' All genes exceeding this percentile threshold will be selected.
#' @param deltMu1Thresh The percentile threshold for the deltaMu1 statistic.
#' All genes exceeding this percentile threshold will be selected.
#' @param siThresh The threshold for the selectivity index statistic (between 0 - 1).
#' All genes exceeding this threshold will be selected.
#' @keywords subsetting
#' @return Returns a dataframe containing all genes meeting the prespecified thresholds.
#' @export
#' @examples
#' topGeneQuants(mmParams)
#' @seealso \code{\link{mixModelParams}}

topGeneQuants = function(mmParams, deltMu2Thresh = 90, deltMu1Thresh = 10, siThresh = .99){
  mmParams.df = as.data.frame(mmParams)
  deltaMu2Quant = quantile(mmParams.df$deltaMu2, deltMu2Thresh*.01)
  deltaMu1Quant = quantile(abs(mmParams.df$deltaMu1), {100-deltMu1Thresh}*.01)

  mmParams.df = as.data.frame(mmParams)
  tfVect = mmParams.df$deltaMu2 > deltaMu2Quant & abs(mmParams.df$deltaMu1) < deltaMu1Quant & mmParams.df$SI > siThresh
  mmParams.df.quantSubset = mmParams.df[tfVect,]
  return(mmParams.df.quantSubset)
}



#' Identify the top N genes based on a score
#'
#' This function allows you to identify the top N genes that most
#' closely match the distributional profiles of a theoretical oncogene.
#'
#' @param mmParams The output from the mixModelParams function.
#' @param N An integer specified by the user to indicate how many genes they would like to return.
#' @keywords subsetting
#' @return Returns a dataframe containing the top N genes (where N is specified by the user)
#' computed using a custom score that incorporates differences between component means and the variance.
#' @export
#' @examples
#' topGeneQuants(mmParams, N = 10)
#' @seealso \code{\link{mixModelParams}}


#returns a list of genes
topGeneTable = function(mmParams, N=nrow(mmParams)){
  #returns the top N genes based on a score
  mmParams.df = as.data.frame(mmParams)
  mmParams.df$score = mmParams.df$SI*{(mmParams.df$deltaMu2 -  mmParams.df$deltaMu1) - (mmParams.df$n.var + mmParams.df$t.var)}
  mmParams.df.s = mmParams.df[with(mmParams.df, order(-score)), ]
  mmParams.df.s.subset = mmParams.df.s[1:N,]
  return(mmParams.df.s.subset)
}
