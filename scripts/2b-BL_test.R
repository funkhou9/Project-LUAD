#' # 2a-RR_test,R
#' 
#' - Scott Funkhouser and Ting Shen -
#' - STT855 LUAD Final Project -
#' 
#' ## Objective
#' 1. Now that we have downloaded cnv and normalized rna-seq data,
#'	we will test if we can obtain uncorrelated samples from the posterior
#'	using ridge regression.

setwd("/mnt/home/funkhou9/Project-LUAD/scripts")
library(BGLR)
library(magrittr)
library(pROC)

#' Load data from 1-assemble_data.R
load("../data/processed/for_analysis/data_for_analysis.RData")

#' ### Ridge regression with cnv data
#' Center and scale X.
X <- scale(cnv[, -1])

#' Outcome - binary tumor status.
Y <- cnv[, 1]

#' Randomly mask 200 outcomes for prediction.
mask <- sample(1:length(Y), size = 200)
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

#' Set up linear predictor and fit using bayesian lasso.
ETA <- list(cnv = list(X = X_full, model = 'BL'))

#+ results = 'hide'
fm <- BGLR(y = Y,
	   	   ETA = ETA,
	   	   response_type = "ordinal",
	   	   nIter = 20000,
	   	   burnIn = 2000,
	   	   thin = 2,
	   	   saveAt = "../data/processed/for_analysis/posterior/cnv/BL-")

#' Inspect trace plots.
lambda <- scan("../data/processed/for_analysis/posterior/cnv/BL-ETA_cnv_lambda.dat")
varE <- scan("../data/processed/for_analysis/posterior/cnv/BL-varE.dat")
mu <- scan("../data/processed/for_analysis/posterior/cnv/BL-mu.dat")
thresholds <- scan("../data/processed/for_analysis/posterior/cnv/BL-thresholds.dat")

plot(lambda)
plot(thresholds)

#' ### Check prediction accuracy
roc <- roc(cnv[, 1][mask], fm$yHat[mask])

#' Determine optimal threshold for classification.
thresh <- which(roc$sensitivities == max(roc$sensitivities) & roc$specificities == max(roc$specificities))

#' Visualze prections with one classification threshold.
plot(fm$yHat[mask],
	 cnv[, 1][mask],
	 xlab = "Preditions",
	 ylab = "Tumor status",
	 col = c("red", "blue")[(cnv[, 1][mask] == 1) + 1])
abline(v = roc$thresholds[thresh][1])

#' Plot ROC curve using optimal thresholds(s).
roc <- roc(cnv[, 1][mask], fm$yHat[mask])
plot(roc, col = "blue")


save.image(file = "../2b-BL_test.RData")



