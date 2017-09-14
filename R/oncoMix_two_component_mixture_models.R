#' Generate the parameters for two 2-component Gaussian mixture models
#' with equal variances
#'
#' This function allows you to generate the parameters for two 2-component
#' Gaussian mixture model with equal variances from 2 matrices of data with
#' a priori labels (eg tumor vs normal.) This application was originally
#' intended for matrices of gene expression data treated with 2 conditions.
#'
#' @param dfNml A dataframe of normal data containing patients as rows and
#' genes as columns
#' @param dfTumor A dataframe of tumor data containing patients as rows and
#' genes as columns
#' @keywords oncoMix, mixture-model, two-component
#' @return Returns a list of length 2, each element of which contains the
#' parameters for the Normal or the Tumor data in a 3 x n matrix, where n
#' is the number of patient samples
#' @importFrom mclust Mclust mclustBIC
#' @export
#' @examples
#' dfNml <- as.data.frame(matrix(data=rgamma(n=150, shape=2, rate=2),
#' nrow=15, ncol=10))
#' rownames(dfNml) <- paste0("patientN", 1:nrow(dfNml))
#' colnames(dfNml) <- paste0("gene", 1:ncol(dfNml))
#'
#' dfTumor <- as.data.frame(matrix(data=rgamma(n=150, shape=4, rate=3),
#' nrow=15, ncol=10))
#' rownames(dfTumor) <- paste0("patientT", 1:nrow(dfTumor))
#' colnames(dfTumor) <- paste0("gene", 1:ncol(dfTumor))
#'
#' mmParams <- mixModelParams(dfNml, dfTumor)

mixModelParams <- function(dfNml, dfTumor) {
    paramsNormal <- apply(dfNml, 2, function(x) {
        y <- mclust::Mclust(data=x, G=2, modelNames="E", verbose=FALSE)
        z <- c(nMu=y$parameters$mean,
            nVar=y$parameters$variance$sigmasq,
            nPi1=y$parameters$pro[1])
        return(z)})

    paramsTumor <- apply(dfTumor, 2, function(x) {
        y <- mclust::Mclust(data=x, G=2, modelNames="E", verbose=FALSE)
        z <- c(tMu=y$parameters$mean,
            tVar=y$parameters$variance$sigmasq,
            tPi1=y$parameters$pro[1])
        return(z)})

    params <- rbind(paramsNormal, paramsTumor)
    rownames(params)[c(1:2, 5:6)] = c("nMu1","nMu2","tMu1", "tMu2")
    deltaMu2 <- params["tMu2",]-params["nMu2",]
    deltaMu1 <- params["tMu1",]-params["nMu1",]
    boundaryTumor <- (params["tMu2",] + params["tMu1",]) / 2

    siCalc <- function(vectNml, boundTumor){ #calculate the sel. idx
        x <- sum(boundTumor > vectNml) / length(vectNml)
        return(x)
    }
    SI <- sapply(1:ncol(dfNml),
            function(i) siCalc(dfNml[,i], boundaryTumor[i]))

    params <- rbind(params, deltaMu2, deltaMu1, SI)
    mmParamsDf <- data.frame(t(params))
    score <- NA
    mmParamsDf$score <- mmParamsDf$SI*{
        (mmParamsDf$deltaMu2-mmParamsDf$deltaMu1)-
        (mmParamsDf$nVar+mmParamsDf$tVar)}
    mmParamsDfS <- mmParamsDf[with(mmParamsDf, order(-score)), ]
    return(mmParamsDfS)
}


#' Plot a histogram of gene expression values from
#' tumor and adjacent normal tissue.
#'
#' This function allows you to plot a histogram of gene expression values from
#' tumor and adjacent normal tissue with the option of including the best
#' fitting Gaussian curve.
#'
#' @param mmParams The output from the getMixModelParams function.
#' @param dfNml The normal dataframe.
#' @param dfTumor The tumor dataframe.
#' @param isof The gene isoform to visualize
#' @keywords oncoMix, visualization, Gaussian, two-component
#' @return Returns a histogram of the gene expression
#' values from the two groups.
#' @export
#' @examples
#' dfNml <- as.data.frame(matrix(data=rgamma(n=150, shape=2, rate=2),
#' nrow=15, ncol=10))
#' rownames(dfNml) <- paste0("patientN", 1:nrow(dfNml))
#' colnames(dfNml) <- paste0("gene", 1:ncol(dfNml))
#'
#' dfTumor <- as.data.frame(matrix(data=rgamma(n=150, shape=4, rate=3),
#' nrow=15, ncol=10))
#' rownames(dfTumor) <- paste0("patientT", 1:nrow(dfTumor))
#' colnames(dfTumor) <- paste0("gene", 1:ncol(dfTumor))
#'
#' mmParams <- mixModelParams(dfNml, dfTumor)
#' isof <- rownames(mmParams)[1]
#' plotGeneHist(mmParams, dfNml, dfTumor, isof)
#' @seealso \code{\link{mixModelParams}}

plotGeneHist <- function(mmParams, dfNml, dfTumor, isof){
    tidyDf <- as.data.frame(cbind(as.numeric(c(dfTumor[,isof], dfNml[,isof])),
        as.factor(c(rep("tumor",nrow(dfTumor)), rep("normal",nrow(dfNml))))),
        stringsAsFactors=FALSE)
    colnames(tidyDf) <- c("expr", "type")
    expr <- type <- ..density.. <- NULL # Setting the variables to NULL first
    p1 <- ggplot(tidyDf, aes(x=expr, color=as.factor(type),
        fill=as.factor(type),
        group=as.factor(type))) +
    theme_classic() +
    geom_histogram(data=subset(tidyDf,type == 1),fill="#F8766D",
        alpha=0.2, aes(y=..density..)) +
    geom_histogram(data=subset(tidyDf,type == 2),fill="#00BFC4",
        alpha=0.2, aes(y=..density..)) +
    geom_rug(alpha=0.3, show.legend=FALSE) +
    theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="none",
        plot.title=element_text(size=12),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8)) +
    ggtitle(paste0(isof, " : SI = ",
        round(mmParams[isof,"SI"],4))) +
    xlab(expression(Log[2] *"(TPM Reads)")) +
    stat_function(fun="dnorm", colour="#F8766D",
        args=list(mean=mmParams[isof,"nMu2"],
        sd=sqrt(mmParams[isof,"nVar"]))) +
    stat_function(fun="dnorm", colour="#F8766D",
        args=list(mean=mmParams[isof,"nMu1"],
        sd=sqrt(mmParams[isof,"nVar"]))) +
    stat_function(fun="dnorm", colour="#00BFC4",
        args=list(mean=mmParams[isof,"tMu2"],
        sd=sqrt(mmParams[isof,"tVar"]))) +
    stat_function(fun="dnorm", colour="#00BFC4",
        args=list(mean=mmParams[isof,"tMu1"],
        sd=sqrt(mmParams[isof,"tVar"])))
    print(p1)
}


#' Generate a scatter plot with the output from mixModelParams
#'
#' This function allows you to generate the parameters for two 2-component
#' mixture models with equal variances
#'
#' @param mmParams The output from the mixModelParams function. Will utilize
#' the deltaMu2 and deltaMu1 rows
#' @param selIndThresh This is the selectivity index threshold to use. All
#' genes with SI values above this threshold will be highlighted in purple.
#' @param geneLabels A character vector of gene names used to label the
#' genes with that name on the scatter plot.
#' @keywords oncoMix, visualization, two-component
#' @return Returns a ggplot scatter object that can be plotted
#' @import ggplot2 ggrepel stats
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @examples
#' dfNml <- as.data.frame(matrix(data=rgamma(n=150, shape=2, rate=2),
#' nrow=15, ncol=10))
#' rownames(dfNml) <- paste0("patientN", 1:nrow(dfNml))
#' colnames(dfNml) <- paste0("gene", 1:ncol(dfNml))
#'
#' dfTumor <- as.data.frame(matrix(data=rgamma(n=150, shape=4, rate=3),
#' nrow=15, ncol=10))
#' rownames(dfTumor) <- paste0("patientT", 1:nrow(dfTumor))
#' colnames(dfTumor) <- paste0("gene", 1:ncol(dfTumor))
#'
#' mmParams <- mixModelParams(dfNml, dfTumor)
#' scatterMixPlot(mmParams)
#' @seealso \code{\link{mixModelParams}}

scatterMixPlot <- function(mmParams, selIndThresh=1, geneLabels=NULL){
    mmParams <- as.data.frame(mmParams)
    oneOverAlpha <- diff(range(mmParams$deltaMu2))
    alpha1 <- 1/oneOverAlpha

    quants <- c(0.01, 0.10, 0.50, 0.90, 0.99) #add in the quantiles
    colorsRed <- RColorBrewer::brewer.pal(n=length(quants), name="Reds")

    deltaMu2Quant <- stats::quantile(mmParams[,"deltaMu2"], quants)
    deltaMu1Quant <- stats::quantile(1/(abs(mmParams[,"deltaMu1"]) + alpha1),
                        quants)
    deltaMu2 <- deltaMu1 <- NULL
    x <- ggplot(data=as.data.frame(mmParams),
        aes(x=deltaMu2, y=1/(abs(deltaMu1) + alpha1))) +
        theme_classic() +
        geom_hline(yintercept=deltaMu1Quant,
            col=colorsRed, size=c(1,1,1,1,1)) +
        geom_vline(xintercept=deltaMu2Quant,
            col=colorsRed, size=c(1,1,1,1,1)) +
        geom_point(alpha=0.5) +
        xlab(expression(paste(Delta, mu[2]))) +
        ylab(expression(paste(frac(1, paste(Delta, mu[1], " + ", alpha)))))
    if(selIndThresh < 1){
        mmParamsSi <- mmParams[mmParams$SI > selIndThresh,]
        x <- x + geom_point(data=as.data.frame(mmParamsSi),
            aes(x=deltaMu2, y=1/(abs(deltaMu1)+alpha1)),
            size=10, alpha=0.1,
            col=colorsRed[length(colorsRed)],
            fill=colorsRed[length(colorsRed)]) +
        ggtitle(bquote(Distribution~of~Mixture~Model~
            Parameters*","~alpha~"="~.(round(alpha1,2))*", SI >"~
            .(selIndThresh)))
    } else if(!is.null(geneLabels)){
        mmParamsSi <- mmParams[geneLabels,]
        mmParamsSi$geneLabels <- geneLabels
        x <- x + geom_point(data=as.data.frame(mmParamsSi),
            aes(x=deltaMu2, y=1/(abs(deltaMu1)+alpha1)),
            size=10, alpha=0.1,
            col=colorsRed[length(colorsRed)],
            fill=colorsRed[length(colorsRed)]) +
        ggrepel::geom_text_repel(data=mmParamsSi, aes(x=deltaMu2,
            y=1/(abs(deltaMu1)+alpha1)),
            label=rownames(mmParamsSi)) +
        ggtitle(bquote(Distribution~of~Mixture~Model~
            Parameters*","~alpha~"="~.(round(alpha1,2))))
    } else {
        x <- x + ggtitle(bquote(Distribution~of~Mixture~Model~
            Parameters*","~alpha~"="~.(round(alpha1,2))))
    }
    return(x)
}

#' Identify genes that meet pre-specified quantiles
#'
#' This function allows you to subset genes that are above pre-specified
#' quantiles and that most closely resemble the distribution of oncogenes.
#'
#' @param mmParams The output from the mixModelParams function.
#' @param deltMu2Thr The percentile threshold for the deltaMu2 statistic.
#' All genes exceeding this percentile threshold will be selected.
#' @param deltMu1Thr The percentile threshold for the deltaMu1 statistic.
#' All genes exceeding this percentile threshold will be selected.
#' @param siThr The threshold for the selectivity index statistic
#' (between 0-1). All genes exceeding this threshold will be selected.
#' @keywords subsetting
#' @return Returns a dataframe containing all genes meeting the prespecified
#' thresholds.
#' @import ggplot2 ggrepel stats
#' @export
#' @examples
#' dfNml <- as.data.frame(matrix(data=rgamma(n=150, shape=2, rate=2),
#' nrow=15, ncol=10))
#' rownames(dfNml) <- paste0("patient_n", 1:nrow(dfNml))
#' colnames(dfNml) <- paste0("gene", 1:ncol(dfNml))
#'
#' dfTumor <- as.data.frame(matrix(data=rgamma(n=150, shape=4, rate=3),
#' nrow=15, ncol=10))
#' rownames(dfTumor) <- paste0("patient_t", 1:nrow(dfTumor))
#' colnames(dfTumor) <- paste0("gene", 1:ncol(dfTumor))
#'
#' mmParams <- mixModelParams(dfNml, dfTumor)
#' topGeneQuants(mmParams)
#' @seealso \code{\link{mixModelParams}}

topGeneQuants <- function(mmParams, deltMu2Thr=90, deltMu1Thr=10, siThr=.99){
    mmParamsDf <- as.data.frame(mmParams)
    deltaMu2Quant <- stats::quantile(mmParamsDf$deltaMu2, deltMu2Thr*.01)
    deltaMu1Quant <- stats::quantile(abs(mmParamsDf$deltaMu1),
        {100-deltMu1Thr}*.01)
    mmParamsDf <- as.data.frame(mmParams)
    tfVect <- mmParamsDf$deltaMu2 > deltaMu2Quant &
        abs(mmParamsDf$deltaMu1) < deltaMu1Quant & mmParamsDf$SI > siThr
    mmParamsDfQuantSubset <- mmParamsDf[tfVect,]
    return(mmParamsDfQuantSubset)
}


#' Human Breast Cancer RNA-sequencing data from TCGA - Tumor Tissue
#'
#' @name dfTumorIsof
#' @docType data
#' @author Daniel Pique \email{daniel.pique@med.einstein.yu.edu}
#' @references \url{https://gdc.cancer.gov/}
#' @keywords RNA expression, Breast Tumor

NULL


#' Human Breast Cancer RNA-sequencing data from TCGA - Adj. Normal Tissue
#'
#' @name dfNmlIsof
#' @docType data
#' @author Daniel Pique \email{daniel.pique@med.einstein.yu.edu}
#' @references \url{https://gdc.cancer.gov/}
#' @keywords RNA expression, Breast Tissue

NULL
