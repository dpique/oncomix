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
    y <- mclust::Mclust(data = x, G=2, modelNames = "E")
    z <- c(n.mu = y$parameters$mean, n.var = y$parameters$variance$sigmasq, n.pi.1 = y$parameters$pro[1])
    #n.class = y$classification
    return(z) })

  params_tumor <- apply(dfTumor, 2, function(x) {
    y <- mclust::Mclust(data = x, G=2, modelNames = "E")
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


#' Generate a scatter plot with the output from mixModelParams
#'
#' This function allows you to generate the parameters for two 2-component mixture models 
#' with equal variances
#' 
#' @param mixModelParams The output from the getMixModelParams function. Will utilize the deltaMu2
#' and deltaMu1 rows
#' @keywords oncoMix, visualization, two-component
#' @return Returns a ggplot object that you can plot
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @examples
#' scatterMixPlot(mixModelParams)
#' @seealso \code{\link{mixModelParams}} 

scatterMixPlot <- function(mixModelParams, selIndThresh = 1){
  mixModelParams = as.data.frame(mixModelParams)
  one_over_alpha = diff(range(mixModelParams$deltaMu2))
  alpha1 = 1/one_over_alpha
  
  quants = c(0.01, 0.10, 0.50, 0.90, 0.99) #add in the quantiles
  colors_red=RColorBrewer::brewer.pal(n=length(quants), name="Reds")
  
  deltaMu2Quant <- quantile(mixModelParams[,"deltaMu2"], quants)
  deltaMu1Quant <- quantile(1/(abs(mixModelParams[,"deltaMu1"]) + alpha1), quants)

  x = ggplot(data = as.data.frame(mixModelParams), aes(x = deltaMu2, y = 1/(abs(deltaMu1) + alpha1))) + 
    theme_classic() + 
    geom_hline(yintercept = deltaMu1Quant, col = colors_red, size = c(1,1,1,1,1)) + 
    geom_vline(xintercept = deltaMu2Quant, col = colors_red, size = c(1,1,1,1,1)) +
    geom_point(alpha= 0.5) + 
    xlab(expression(paste(Delta, mu[2]))) + 
    ylab(expression(paste(frac(1, paste(Delta, mu[1], " + ", alpha))))) 
  #print(x)

  if(selIndThresh < 1){
    mixModelParams.si = mixModelParams[mixModelParams$SI > selIndThresh,]
    x = x + geom_point(data = as.data.frame(mixModelParams.si), 
                       aes(x = deltaMu2, y = 1/(abs(deltaMu1)+alpha1)), 
                       size = 10, alpha=0.1, 
                       col=colors_red[length(colors_red)], 
                       fill=colors_red[length(colors_red)]) +
      ggtitle(bquote(Distribution~of~Mixture~Model~Parameters*","~alpha~"="~.(round(alpha1,2))*", SI >"~.(selIndThresh)))
  } else {
    x = x + ggtitle(bquote(Distribution~of~Mixture~Model~Parameters*","~alpha~"="~.(round(alpha1,2))))
  }
  return(x)
}



