#' # 2c-RR_RNA_test.R
#' 
#' - Scott Funkhouser and Ting Shen -
#' - STT855 LUAD Final Project -
#' 
#' ## Objective
#' 1. We want to test whether using the RNA seq data with a ridge regression estimator
#'	can produce a sensitive and specific estimator based on one trail of training and testing.

setwd("/mnt/home/funkhou9/Project-LUAD/scripts")
library(BGLR)
library(magrittr)
library(pROC)

#' Load data from 1-assemble_data.R
load("../data/processed/for_analysis/data_for_analysis.RData")

#' ### Ridge regression with rna data
#' Center and scale X.
X <- scale(mrna[, -1])

#' Outcome - binary tumor status.
Y <- mrna[, 1]

#' Randomly mask 200 outcomes for prediction.
mask <- sample(1:length(Y), size = 100)
Y[mask] <- NA

#' Check for missing values.
apply(X, 2, 
	  function(x) {
	  	any(is.na(x))
	  }) %>%
	sum()

#' Remove genes with any missing values.
idx <- apply(X, 2, 
	  		 function(x) {
	  			any(is.na(x))
	  		 })
X_full <- X[, !idx]

#' Set up linear predictor and fit using ridge regression.
ETA <- list(mrna = list(X = X_full, model = 'BL'))

#+ results = 'hide'
fm <- BGLR(y = Y,
	   	   ETA = ETA,
	   	   response_type = "ordinal",
	   	   nIter = 20000,
	   	   burnIn = 2000,
	   	   thin = 2,
	   	   saveAt = "../data/processed/for_analysis/posterior/mrna/BL-")

#' Inspect trace plots.
lambda <- scan("../data/processed/for_analysis/posterior/mrna/BL-ETA_mrna_lambda.dat")
varE <- scan("../data/processed/for_analysis/posterior/cnv/BL-varE.dat")
mu <- scan("../data/processed/for_analysis/posterior/cnv/BL-mu.dat")
thresholds <- scan("../data/processed/for_analysis/posterior/cnv/BL-thresholds.dat")

plot(lambda)

#' Check prediction accuracy and determine optimal classification threshold.
plot(fm$yHat[mask],
	 mrna[, 1][mask],
	 xlab = "Preditions",
	 ylab = "Tumor status")

#' Check ROC curve
roc <- roc(mrna[, 1][mask], fm$yHat[mask])
plot(roc, col = "green")

save.image(file = "../2c-RR_RNA_test.RData")


