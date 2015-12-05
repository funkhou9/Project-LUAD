#' # 2a-RR_test.R
#' 
#' - Scott Funkhouser and Ting Shen
#' - STT855 LUAD Final Project
#' 
#' ## Objective
#' 1. We want to test whether prediction of tissue status - normal or tumor can be done using
#'	the 18 tumor CNVs identified in TCGA-LUAD paper 2014. And whether prediction accuracy can be
#'	increased by using all 20K CNV parameters.
#'
#' ### Table of contents
#' 1. [Prediction with candidate CNVs](#prediction-with-candidate-cnvs)
#' 2. [Prediction with all CNVs](#prediction-with-all-cnvs)

setwd("/mnt/home/funkhou9/Project-LUAD/scripts")
library(BGLR)
library(magrittr)
library(pROC)

#' Load data from 1-assemble_data.R
load("../data/processed/for_analysis/data_for_analysis.RData")
X <- scale(cnv[, -1])

#' Outcome - binary tumor status.
Y <- cnv[, 1]
Y_masked <- Y

#' Randomly mask 200 outcomes for prediction.
mask <- sample(1:length(Y), size = 200)
Y_masked[mask] <- NA

#' ### Prediction with candidate CNVs
#' Save names from candidate genes with copy number alterations. Identified with GISTIC.
sigcopy <- c("NKX2-1",
			 "TERT",
			 "MDM2",
			 "MCL1",
			 "KRAS",
			 "EGFR",
			 "MET",
			 "TERC",
			 "CCND1",
			 "CCND3",
			 "CDK4",
			 "ERBB2",
			 "CDK3",
			 "ZNF217",
			 "CDKN2A",
			 "PTPRD",
			 "KAT6A",
			 "JAK2")

X_sigcopy <- X[, colnames(X) %in% sigcopy]

#' Check for missing values.
idx <- apply(X_sigcopy, 2,
	  		 function(x) {
	  			any(is.na(x))
	  		 })

X_sigcopy_full <- X_sigcopy[, !idx]

#' Check if any filtering of X is necessary.
apply(X_sigcopy_full, 2, summary)

#' Eigenvector decomposition of G derived from X_sigcopy_full.
G <- X_sigcopy_full %*% t(X_sigcopy_full) / ncol(X_sigcopy_full)
evd <- eigen(G)
plot(evd$vector[, 1],
	 evd$vector[, 2],
	 col = c("blue", "red")[(Y == 1) + 1],
	 xlab = "PC1",
	 ylab = "PC2")
legend(0.1,
	   -0.15,
	   pch = c(1, 1),
	   c("Normal", "Tumor"),
	   col = c("blue", "red"))

#' Perform RR.
X_sigcopy_full_std <- X_sigcopy_full / sqrt(ncol(X_sigcopy_full))
ETA <- list(cnv = list(X = X_sigcopy_full_std, model = 'BRR'))

#+ results = 'hide'
fm <- BGLR(y = Y_masked,
	   	   ETA = ETA,
	   	   response_type = "ordinal",
	   	   nIter = 10000,
	   	   burnIn = 2000,
	   	   thin = 2,
	   	   saveAt = "../data/processed/for_analysis/posterior/cnv/RR-subset-")

#' Check trace plot of var(b).
cnv_varB <- scan("../data/processed/for_analysis/posterior/cnv/RR-subset-ETA_cnv_varB.dat")
plot(cnv_varB,
	 xlab = "Iteration",
	 ylab = "varB",
	 main = "Trace plot")

#' Check prediction accuracy.
plot(fm$yHat[mask],
	 Y[mask],
	 col = c("red", "blue")[(Y[mask] == 1) + 1],
	 xlab = "linear predictor (Xb)",
	 ylab = "tumor status (Y)",
	 main = "Predictions for TST")

roc <- roc(Y[mask], fm$yHat[mask])
plot(roc,
	 col = "darkgreen")
text(0.2, 0.2,
	 paste("AUC = ", round(roc$auc, 3), sep = ''),
	 cex = 2)


#' ### Prediction with all CNVs
#' Remove genes with any missing values.
idx <- apply(X, 2, 
	  		 function(x) {
	  			any(is.na(x))
	  		 })
X_full <- X[, !idx]

#' Eigenvector decomposition of G derived from X_full.
G_full <- X_full %*% t(X_full) / ncol(X_full)
evd_full <- eigen(G_full)
plot(evd_full$vector[, 1],
	 evd_full$vector[, 2],
	 col = c("blue", "red")[(Y == 1) + 1],
	 xlab = "PC1",
	 ylab = "PC2")
legend(0.1,
	   -0.15,
	   pch = c(1, 1),
	   c("Normal", "Tumor"),
	   col = c("blue", "red"))

#' Set up linear predictor and fit using ridge regression.
ETA <- list(cnv = list(X = X_full, model = 'BRR'))

#+ results = 'hide'
fm_full <- BGLR(y = Y_masked,
	   	   		ETA = ETA,
	   	   		response_type = "ordinal",
	   	   		nIter = 120000,
	   	   		burnIn = 2000,
	   	   		thin = 2,
	   	   		saveAt = "../data/processed/for_analysis/posterior/cnv/RR-")

#' Check trace plot of var(b).
cnv_varB_full <- scan("../data/processed/for_analysis/posterior/cnv/RR-ETA_cnv_varB.dat")
plot(cnv_varB_full,
	 xlab = "Iteration",
	 ylab = "varB",
	 main = "Trace plot")

#' Check prediction accuracy.
plot(fm_full$yHat[mask],
	 Y[mask],
	 col = c("red", "blue")[(Y[mask] == 1) + 1],
	 xlab = "linear predictor (Xb)",
	 ylab = "tumor status (Y)",
	 main = "Predictions for TST")

roc <- roc(Y[mask], fm_full$yHat[mask])
plot(roc,
	 col = "darkgreen")
text(0.2, 0.2,
	 paste("AUC = ", round(roc$auc, 3), sep = ''),
	 cex = 2)


save.image(file = "../2a-RR_test.RData")