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

Set up linear predictor and fit using bayesian lasso.


```r
ETA <- list(cnv = list(X = X_full, model = 'BL'))
```

```r
fm <- BGLR(y = Y,
	   	   ETA = ETA,
	   	   response_type = "ordinal",
	   	   nIter = 20000,
	   	   burnIn = 2000,
	   	   thin = 2,
	   	   saveAt = "../data/processed/for_analysis/posterior/mrna/BL-")
```

```
## Warning in sqrt(post_threshold2 - post_threshold^2): NaNs produced
```

Inspect trace plots.


```r
lambda <- scan("../data/processed/for_analysis/posterior/cnv/BL-ETA_cnv_lambda.dat")
varE <- scan("../data/processed/for_analysis/posterior/cnv/BL-varE.dat")
mu <- scan("../data/processed/for_analysis/posterior/cnv/BL-mu.dat")
thresholds <- scan("../data/processed/for_analysis/posterior/cnv/BL-thresholds.dat")

plot(lambda)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png) 

```r
plot(thresholds)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-2.png) 

Check prediction accuracy and determine optimal classification threshold.


```r
plot(fm$yHat[mask],
	 cnv[, 1][mask],
	 col = c("red", "blue")[(fm$yHat[mask] > 0.5) + 1],
	 xlab = "Preditions",
	 ylab = "Tumor status")
abline(v = 0.5,
	   col = "darkgreen")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png) 

With an threshold 0.5, a perfect classifier is established - with a 100% TP (ture positive) rate and
	100% TN (true negative) rate.


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
## Area under the curve: 0.9847
```

```r
save.image(file = "../2b-BL_test.RData")
```

