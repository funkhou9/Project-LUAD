# 2c-RR_RNA_test.R

- Scott Funkhouser and Ting Shen
- STT855 LUAD Final Project

## Objective
1. We want to test whether using the RNA seq data with a ridge regression estimator
	can produce a sensitive and specific estimator based on one trail of training and testing.


```r
setwd("/mnt/home/funkhou9/Project-LUAD/scripts")
library(BGLR)
library(magrittr)
library(pROC)
```

```
## Type 'citation("pROC")' for a citation.
## 
## Attaching package: 'pROC'
## 
## The following objects are masked from 'package:stats':
## 
##     cov, smooth, var
```

Load data from 1-assemble_data.R


```r
load("../data/processed/for_analysis/data_for_analysis.RData")
```

### Ridge regression with rna data
Center and scale X.


```r
X <- scale(mrna[, -1])
```

Outcome - binary tumor status.


```r
Y <- mrna[, 1]
```

Randomly mask 200 outcomes for prediction.


```r
mask <- sample(1:length(Y), size = 100)
Y[mask] <- NA
```

Check for missing values.


```r
apply(X, 2, 
	  function(x) {
	  	any(is.na(x))
	  }) %>%
	sum()
```

```
## [1] 332
```

Remove genes with any missing values.


```r
idx <- apply(X, 2, 
	  		 function(x) {
	  			any(is.na(x))
	  		 })
X_full <- X[, !idx]
```

Set up linear predictor and fit using ridge regression.


```r
ETA <- list(mrna = list(X = X_full, model = 'BRR'))
```

```r
fm <- BGLR(y = Y,
	   	   ETA = ETA,
	   	   response_type = "ordinal",
	   	   nIter = 20000,
	   	   burnIn = 2000,
	   	   thin = 2,
	   	   saveAt = "../data/processed/for_analysis/posterior/mrna/RR-")
```

```
## Warning in sqrt(post_threshold2 - post_threshold^2): NaNs produced
```

Inspect trace plots.


```r
mrna_varB <- scan("../data/processed/for_analysis/posterior/mrna/RR-ETA_mrna_varB.dat")
varE <- scan("../data/processed/for_analysis/posterior/mrna/RR-varE.dat")
mu <- scan("../data/processed/for_analysis/posterior/mrna/RR-mu.dat")

plot(mrna_varB)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png) 

Check prediction accuracy and determine optimal classification threshold.


```r
plot(fm$yHat[mask],
	 mrna[, 1][mask],
	 xlab = "Preditions",
	 ylab = "Tumor status")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png) 

Check ROC curve


```r
roc <- roc(mrna[, 1][mask], fm$yHat[mask])
plot(roc, col = "green")
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png) 

```
## 
## Call:
## roc.default(response = mrna[, 1][mask], predictor = fm$yHat[mask])
## 
## Data: fm$yHat[mask] in 89 controls (mrna[, 1][mask] 0) < 11 cases (mrna[, 1][mask] 1).
## Area under the curve: 1
```

```r
save.image(file = "../2c-RR_RNA_test.RData")
```

