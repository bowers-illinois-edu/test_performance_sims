## Load simulation data and parameters for simulation

library(here)
library(data.table)
library(manytestsr)

load(here::here("Data", "idat_not_equal_nb.rda"),verbose=TRUE)
load(here::here("Data", "bdat_not_equal_nb.rda"),verbose=TRUE)
load(here::here("Data", "idat_equal_nb.rda"),verbose=TRUE)
load(here::here("Data", "bdat_equal_nb.rda"),verbose=TRUE)
load(here::here("Data", "blocks_sampled.rda"),verbose=TRUE)

setkey(bdat_not_equal_nb, bC)
setkey(idat_not_equal_nb, bC)

icols <- c("y0","bC", "bF", "trtF", "trt", "covscluster", "covsplits", "vb1", "vb2", "vb3", "ub0")
bcols <- c("bC", "bF", "hwt", "covscluster", "covsplits")

## The unequal sized blocks
bdat_not_equal_B10 <- droplevels(bdat_not_equal_nb[.(blocks_not_equalB10), ..bcols])
bdat_not_equal_B25 <- droplevels(bdat_not_equal_nb[.(blocks_not_equalB25), ..bcols])
bdat_not_equal_B50 <- droplevels(bdat_not_equal_nb[.(blocks_not_equalB50), ..bcols])
bdat_not_equal_B75 <- droplevels(bdat_not_equal_nb[.(blocks_not_equalB75), ..bcols])
bdat_not_equal_B100 <- droplevels(bdat_not_equal_nb[.(blocks_not_equalB100), ..bcols])
bdat_not_equal_B1000 <- droplevels(bdat_not_equal_nb[.(blocks_not_equalB1000), ..bcols])
setkey(bdat_not_equal_B10, bC)
setkey(bdat_not_equal_B25, bC)
setkey(bdat_not_equal_B50, bC)
setkey(bdat_not_equal_B75, bC)
setkey(bdat_not_equal_B100, bC)
setkey(bdat_not_equal_B1000, bC)

idat_not_equal_B10 <- droplevels(idat_not_equal_nb[.(blocks_not_equalB10), ..icols, drop = TRUE])
idat_not_equal_B25 <- droplevels(idat_not_equal_nb[.(blocks_not_equalB25), ..icols, drop = TRUE])
idat_not_equal_B50 <- droplevels(idat_not_equal_nb[.(blocks_not_equalB50), ..icols, drop = TRUE])
idat_not_equal_B75 <- droplevels(idat_not_equal_nb[.(blocks_not_equalB75), ..icols, drop = TRUE])
idat_not_equal_B100 <- droplevels(idat_not_equal_nb[.(blocks_not_equalB100), ..icols, drop = TRUE])
idat_not_equal_B1000 <- droplevels(idat_not_equal_nb[.(blocks_not_equalB1000), ..icols, drop = TRUE])
setkey(idat_not_equal_B10, bC)
setkey(idat_not_equal_B25, bC)
setkey(idat_not_equal_B50, bC)
setkey(idat_not_equal_B75, bC)
setkey(idat_not_equal_B100, bC)
setkey(idat_not_equal_B1000, bC)

## The equal sized blocks
bdat_equal_B10 <- droplevels(bdat_equal_nb[.(blocks_equalB10), ..bcols])
bdat_equal_B25 <- droplevels(bdat_equal_nb[.(blocks_equalB25), ..bcols])
bdat_equal_B50 <- droplevels(bdat_equal_nb[.(blocks_equalB50), ..bcols])
bdat_equal_B75 <- droplevels(bdat_equal_nb[.(blocks_equalB75), ..bcols])
bdat_equal_B100 <- droplevels(bdat_equal_nb[.(blocks_equalB100), ..bcols])
bdat_equal_B1000 <- droplevels(bdat_equal_nb[.(blocks_equalB1000), ..bcols])
setkey(bdat_equal_B10, bC)
setkey(bdat_equal_B25, bC)
setkey(bdat_equal_B50, bC)
setkey(bdat_equal_B75, bC)
setkey(bdat_equal_B100, bC)
setkey(bdat_equal_B1000, bC)

idat_equal_B10 <- droplevels(idat_equal_nb[.(blocks_equalB10), ..icols, drop = TRUE])
idat_equal_B25 <- droplevels(idat_equal_nb[.(blocks_equalB25), ..icols, drop = TRUE])
idat_equal_B50 <- droplevels(idat_equal_nb[.(blocks_equalB50), ..icols, drop = TRUE])
idat_equal_B75 <- droplevels(idat_equal_nb[.(blocks_equalB75), ..icols, drop = TRUE])
idat_equal_B100 <- droplevels(idat_equal_nb[.(blocks_equalB100), ..icols, drop = TRUE])
idat_equal_B1000 <- droplevels(idat_equal_nb[.(blocks_equalB1000), ..icols, drop = TRUE])
setkey(idat_equal_B10, bC)
setkey(idat_equal_B25, bC)
setkey(idat_equal_B50, bC)
setkey(idat_equal_B75, bC)
setkey(idat_equal_B100, bC)
setkey(idat_equal_B1000, bC)

## Set up y0.
## It turns out that we lose control over y0 if we do not do it separately for each dataset:
## So this is redoing a lot of the work from make_test_data. Sigh.
set.seed(12345)
idat_equal_B10[, ub0new := rnorm(.N, mean = runif(1, min = 5, max = 20), sd = runif(1, min = 1, max = 3)), by = bC]
idat_equal_B10[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB10 <- idat_equal_B10[, .(y0ccR2 = summary(lm(y0new ~ covscluster))$r.squared), by = bC]
## We don't want to the outcome to be perfectly predicted by the covariate within block
## It can have a strong relationship but not perfect.
stopifnot(sum(checkB10$y0ccR2 > .99) < 1)

set.seed(12345)
idat_equal_B25[, ub0new := rnorm(.N, mean = runif(1, min = 1, max = 20), sd = runif(1, min = 1, max = 3)), by = bC]
idat_equal_B25[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB25 <- idat_equal_B25[, .(y0ccR2 = summary(lm(y0new ~ covscluster))$r.squared), by = bC]
stopifnot(sum(checkB25$y0ccR2 > .99) < 10)

set.seed(12345)
idat_equal_B50[, ub0new := rnorm(.N, mean = runif(1, min = 1, max = 20), sd = runif(1, min = 1, max = 3)), by = bC]
idat_equal_B50[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB50 <- idat_equal_B50[, .(y0ccR2 = summary(lm(y0new ~ covscluster))$r.squared), by = bC]
stopifnot(sum(checkB50$y0ccR2 > .99) < 10)

set.seed(12345)
idat_equal_B75[, ub0new := rnorm(.N, mean = runif(1, min = 1, max = 20), sd = runif(1, min = 1, max = 3)), by = bC]
idat_equal_B75[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB75 <- idat_equal_B75[, .(y0ccR2 = summary(lm(y0new ~ covscluster))$r.squared), by = bC]
stopifnot(sum(checkB75$y0ccR2 > .99) < 10)

set.seed(12345)
idat_equal_B100[, ub0new := rnorm(.N, mean = runif(1, min = 1, max = 20), sd = runif(1, min = 1, max = 3)), by = bC]
idat_equal_B100[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB100 <- idat_equal_B100[, .(y0ccR2 = summary(lm(y0new ~ covscluster))$r.squared), by = bC]
stopifnot(sum(checkB100$y0ccR2 > .99) < 10)

set.seed(12345)
idat_equal_B1000[, ub0new := rnorm(.N, mean = runif(1, min = 1, max = 20), sd = runif(1, min = 1, max = 3)), by = bC]
idat_equal_B1000[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB1000 <- idat_equal_B1000[, .(y0ccR2 = summary(lm(y0new ~ covscluster))$r.squared), by = bC]
stopifnot(sum(checkB1000$y0ccR2 > .99) < 10)

## Now for the unnot_equal sized blocks
set.seed(12345)
idat_not_equal_B10[, ub0new := rnorm(.N, mean = runif(1, min = 5, max = 20), sd = runif(1, min = 1, max = 3)), by = bC]
idat_not_equal_B10[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB10 <- idat_not_equal_B10[, .(y0ccR2 = summary(lm(y0new ~ covscluster))$r.squared), by = bC]
## We don't want to the outcome to be perfectly predicted by the covariate within block
## It can have a strong relationship but not perfect.
stopifnot(sum(checkB10$y0ccR2 > .99) < 1)

set.seed(12345)
idat_not_equal_B25[, ub0new := rnorm(.N, mean = runif(1, min = 1, max = 20), sd = runif(1, min = 1, max = 3)), by = bC]
idat_not_equal_B25[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB25 <- idat_not_equal_B25[, .(y0ccR2 = summary(lm(y0new ~ covscluster))$r.squared), by = bC]
stopifnot(sum(checkB25$y0ccR2 > .99) < 10)

set.seed(12345)
idat_not_equal_B50[, ub0new := rnorm(.N, mean = runif(1, min = 1, max = 20), sd = runif(1, min = 1, max = 3)), by = bC]
idat_not_equal_B50[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB50 <- idat_not_equal_B50[, .(y0ccR2 = summary(lm(y0new ~ covscluster))$r.squared), by = bC]
stopifnot(sum(checkB50$y0ccR2 > .99) < 10)

set.seed(12345)
idat_not_equal_B75[, ub0new := rnorm(.N, mean = runif(1, min = 1, max = 20), sd = runif(1, min = 1, max = 3)), by = bC]
idat_not_equal_B75[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB75 <- idat_not_equal_B75[, .(y0ccR2 = summary(lm(y0new ~ covscluster))$r.squared), by = bC]
stopifnot(sum(checkB75$y0ccR2 > .99) < 10)

set.seed(12345)
idat_not_equal_B100[, ub0new := rnorm(.N, mean = runif(1, min = 1, max = 20), sd = runif(1, min = 1, max = 3)), by = bC]
idat_not_equal_B100[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB100 <- idat_not_equal_B100[, .(y0ccR2 = summary(lm(y0new ~ covscluster))$r.squared), by = bC]
stopifnot(sum(checkB100$y0ccR2 > .99) < 10)

set.seed(12345)
idat_not_equal_B1000[, ub0new := rnorm(.N, mean = runif(1, min = 1, max = 20), sd = runif(1, min = 1, max = 3)), by = bC]
idat_not_equal_B1000[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB1000 <- idat_not_equal_B1000[, .(y0ccR2 = summary(lm(y0new ~ covscluster))$r.squared), by = bC]
stopifnot(sum(checkB1000$y0ccR2 > .99) < 10)

## Add a column that has the variance of y0new. Lower variance would be more
## precision and can be used for the splitLOO algorithm (i.e. lower variance
## should have more power). We cannot split based on hwt or nb because this is
## the same for all blocks. Using 1/var because I think splitLOO choose the
## maximum value of the splitby column

idat_equal_B10[,vary0new:=1/var(y0new),by=bC]
idat_equal_B100[,vary0new:=1/var(y0new),by=bC]
idat_equal_B1000[,vary0new:=1/var(y0new),by=bC]
idat_equal_B25[,vary0new:=1/var(y0new),by=bC]
idat_equal_B50[,vary0new:=1/var(y0new),by=bC]
idat_equal_B75[,vary0new:=1/var(y0new),by=bC]
idat_not_equal_B10[,vary0new:=1/var(y0new),by=bC]
idat_not_equal_B100[,vary0new:=1/var(y0new),by=bC]
idat_not_equal_B1000[,vary0new:=1/var(y0new),by=bC]
idat_not_equal_B25[,vary0new:=1/var(y0new),by=bC]
idat_not_equal_B50[,vary0new:=1/var(y0new),by=bC]
idat_not_equal_B75[,vary0new:=1/var(y0new),by=bC]

bdat_equal_B10       <- bdat_equal_B10      [idat_equal_B10[,.(vary0new=mean(vary0new)),by=bC]]
bdat_equal_B100      <- bdat_equal_B100     [ idat_equal_B100[,.(vary0new=mean(vary0new)),by=bC]]
bdat_equal_B1000     <- bdat_equal_B1000    [ idat_equal_B1000[,.(vary0new=mean(vary0new)),by=bC]]
bdat_equal_B25       <- bdat_equal_B25      [ idat_equal_B25[,.(vary0new=mean(vary0new)),by=bC]]
bdat_equal_B50       <- bdat_equal_B50      [ idat_equal_B50[,.(vary0new=mean(vary0new)),by=bC]]
bdat_equal_B75       <- bdat_equal_B75      [ idat_equal_B75[,.(vary0new=mean(vary0new)),by=bC]]
bdat_not_equal_B10   <- bdat_not_equal_B10  [ idat_not_equal_B10[,.(vary0new=mean(vary0new)),by=bC]]
bdat_not_equal_B100  <- bdat_not_equal_B100 [ idat_not_equal_B100[,.(vary0new=mean(vary0new)),by=bC]]
bdat_not_equal_B1000 <- bdat_not_equal_B1000[ idat_not_equal_B1000[,.(vary0new=mean(vary0new)),by=bC]]
bdat_not_equal_B25   <- bdat_not_equal_B25  [ idat_not_equal_B25[,.(vary0new=mean(vary0new)),by=bC]]
bdat_not_equal_B50   <- bdat_not_equal_B50  [ idat_not_equal_B50[,.(vary0new=mean(vary0new)),by=bC]]
bdat_not_equal_B75   <- bdat_not_equal_B75  [ idat_not_equal_B75[,.(vary0new=mean(vary0new)),by=bC]]

rm(list=ls(patt="objs"))
not_equal_objs_names <- ls(patt="not_equal")
equal_objs_names <- grep("not",ls(patt="_equal"),invert = TRUE, value=TRUE)
stopifnot(all.equal(length(not_equal_objs_names),length(equal_objs_names)))

save(list=not_equal_objs_names,file=here("Data","unequal_blocksize_data.rda"))
save(list=equal_objs_names,file=here("Data","equal_blocksize_data.rda"))

system("touch Data/block_datasets.done")

