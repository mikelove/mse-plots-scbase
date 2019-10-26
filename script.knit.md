---
title: "MSE plots"
author: "Michael Love"
---


```r
library(LoomExperiment)
library(ggplot2)
files <- list.files(".", ".loom")
dat <- list()
nms <- sub("deng.(.*).loom","\\1",files)
for (i in 1:4) {
  dat[[nms[i]]] <- import(files[i])
}
sapply(dat, dim)
```

```
##      full.em full.uq sim.em sim.uq
## [1,]   43346   43346  43346  43346
## [2,]     122     122    122    122
```

```r
lapply(dat, assayNames)
```

```
## $full.em
## [1] "matrix" "B"      "F"     
## 
## $full.uq
## [1] "matrix" "aln_B"  "aln_F"  "uniq_B" "uniq_F"
## 
## $sim.em
## [1] "matrix"   "B"        "F"        "lambda_k" "p_k"      "pi_bk"   
## [7] "pi_mk"    "pi_pk"   
## 
## $sim.uq
##  [1] "matrix"     "aln_B"      "aln_F"      "lambda_k"   "p_k"       
##  [6] "pi_bk"      "pi_mk"      "pi_pk"      "uniq_B"     "uniq_F"    
## [11] "uniq_total"
```

```r
mean(colSums(assay(dat[["sim.em"]], "matrix"))) # 147538
```

```
## [1] 147538.4
```


```r
# these will be the counts used on the x-axis
counts <- assay(dat[["sim.uq"]], "uniq_total")
uniq.total <- (assay(dat[["full.uq"]], "uniq_F") +
               assay(dat[["full.uq"]], "uniq_B"))
assay(dat[["full.uq"]], "uniq_total") <- uniq.total
```


```r
get.mse <- function(x, y, loc=NULL) {
  # loc is a binary matrix of which locations to use
  if (is.null(loc)) {
    loc <- apply(as.matrix(x), 2, function(z) !is.nan(z))
  }
  s <- rowSums(loc * (as.matrix(x)-as.matrix(y))^2, na.rm=TRUE)
  n <- rowSums(loc)
  s/n
}
calcDiffMSE <- function(nopp, pp, true) {
  mse.nopp <- get.mse(nopp, true)
  loc <- apply(as.matrix(nopp), 2, function(x) !is.nan(x))
  mse.pp <- get.mse(pp, true, loc)
  keep <- rowSums(!is.nan(nopp)) >= 15 & !is.nan(mse.nopp - mse.pp)
  print(sum(keep))
  d <- data.frame(log10ave=log10(cts),
                  nopp.minus.pp=(mse.nopp - mse.pp))
  d <- d[keep,]
  print(with(d, sum(nopp.minus.pp > .05)))
  print(with(d, sum(nopp.minus.pp < -.05)))
  print(with(d, mean(nopp.minus.pp)))
  d
}
```

EM version


```r
idx <- mcols(dat[["sim.em"]])$Rhat_ase < 1.1 &
       mcols(dat[["sim.em"]])$Selected > 0
table(idx)
```

```
## idx
## FALSE  TRUE 
## 36017  7329
```

```r
true <- (assay(dat[["full.em"]], "F") / assay(dat[["full.em"]], "matrix"))[idx,]
nopp <- (assay(dat[["sim.em"]], "F") / assay(dat[["sim.em"]], "matrix"))[idx,]
pp <- assay(dat[["sim.em"]], "p_k")[idx,]
cts <- rowMeans(counts[idx,])
d <- calcDiffMSE(nopp, pp, true)
```

```
## [1] 7195
## [1] 859
## [1] 121
## [1] 0.01686243
```

```r
ggplot(d, aes(x=log10ave,y=nopp.minus.pp)) + # EM
  geom_hex(bins=120) +
  geom_hline(yintercept=0, alpha=.6,color="grey") +
  ylim(-.4,.4) +
  scale_fill_gradient2(low="black",mid="red",high="white",midpoint=19)
```

```
## Warning: Removed 1 rows containing non-finite values (stat_binhex).
```

<img src="script_files/figure-html/em-1.png" width="672" />

uniq version


```r
idx <- mcols(dat[["sim.uq"]])$Rhat_ase < 1.1 &
       mcols(dat[["sim.uq"]])$Selected > 0
table(idx)
```

```
## idx
## FALSE  TRUE 
## 36478  6868
```

```r
true <- (assay(dat[["full.uq"]], "uniq_F") / assay(dat[["full.uq"]], "uniq_total"))[idx,]
nopp <- (assay(dat[["sim.uq"]], "uniq_F") / assay(dat[["sim.uq"]], "uniq_total"))[idx,]
pp <- assay(dat[["sim.uq"]], "p_k")[idx,]
cts <- rowMeans(counts[idx,])
d <- calcDiffMSE(nopp, pp, true)
```

```
## [1] 5759
## [1] 617
## [1] 3
## [1] 0.01759068
```

```r
ggplot(d, aes(x=log10ave,y=nopp.minus.pp)) + # uniq
  geom_hex(bins=80) +
  geom_hline(yintercept=0, alpha=.6,color="grey") +
  ylim(-.15,.15) +
  scale_fill_gradient2(low="black",mid="red",high="white",midpoint=11)
```

```
## Warning: Removed 9 rows containing non-finite values (stat_binhex).
```

<img src="script_files/figure-html/uniq-1.png" width="672" />
