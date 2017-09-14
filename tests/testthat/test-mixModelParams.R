library(oncomix)
context("output for mixModelParams")
test_that("mixModelParams takes in and returns a DataFrame", {
    dfNml <- as.data.frame(matrix(data=rgamma(n=150, shape=2, rate=2),
    nrow=15, ncol=10))
    rownames(dfNml) <- paste0("patientN", 1:nrow(dfNml))
    colnames(dfNml) <- paste0("gene", 1:ncol(dfNml))
    dfTumor <- as.data.frame(matrix(data=rgamma(n=150, shape=4, rate=3),
        nrow=15, ncol=10))
    rownames(dfTumor) <- paste0("patientT", 1:nrow(dfTumor))
    colnames(dfTumor) <- paste0("gene", 1:ncol(dfTumor))
    mmParams <- mixModelParams(dfNml, dfTumor)
    expect_true(is(mmParams, "data.frame")) #DataFrame
})
