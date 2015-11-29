# 2a-RR_test,R

- Scott Funkhouser and Ting Shen
- STT855 LUAD Final Project

## Objective
1. Now that we have downloaded cnv and normalized rna-seq data,
	we will test if we can obtain uncorrelated samples from the posterior
	using ridge regression.


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

### Ridge regression with cnv data
Center and scale X.


```r
X <- scale(cnv[, -1])
```

Outcome - binary tumor status.


```r
Y <- cnv[, 1]
```

Randomly mask 200 outcomes for prediction.


```r
mask <- sample(1:length(Y), size = 200)
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
## [1] 4485
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
ETA <- list(cnv = list(X = X_full, model = 'BRR'))
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
cnv_varB <- scan("../data/processed/for_analysis/posterior/cnv/RR-ETA_cnv_varB.dat")
varE <- scan("../data/processed/for_analysis/posterior/cnv/RR-varE.dat")
mu <- scan("../data/processed/for_analysis/posterior/cnv/RR-mu.dat")

plot(cnv_varB)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png) 

Check prediction accuracy and determine optimal classification threshold.


```r
plot(fm$yHat[mask],
	 cnv[, 1][mask],
	 xlab = "Preditions",
	 ylab = "Tumor status")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png) 

Check ROC curve


```r
roc <- roc(cnv[, 1][mask], fm$yHat[mask])
plot(roc, col = "green")
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png) 

```
## 
## Call:
## roc.default(response = cnv[, 1][mask], predictor = fm$yHat[mask])
## 
## Data: fm$yHat[mask] in 99 controls (cnv[, 1][mask] 0) < 101 cases (cnv[, 1][mask] 1).
## Area under the curve: 0.9993
```

```r
save.image(file = "../2a-RR_test.RData")
```

