## Load simulation data and parameters for simulation
library(here)
library(data.table)
library(manytestsr)

load(here::here("Data", "mtwrkdat.rda"))
load(here::here("Data", "blocks.rda"))
load(here::here("Analysis", "simparms.rda"))

setkey(mtbdat, blockC)
setkey(mtwrkdat, blockC)

icols <- c("YContNorm", "y0ContNorm", "block", "blockC", "blockF", "trtF", "trt", "covscluster", "covsplits", "vb1", "vb2", "vb3", "ub0")
bcols <- c("block", "blockC", "blockF", "hwt", "truemndiffb", "covscluster", "covsplits")

bdatB10 <- droplevels(mtbdat[.(blocksB10), ..bcols])
bdatB25 <- droplevels(mtbdat[.(blocksB25), ..bcols])
bdatB50 <- droplevels(mtbdat[.(blocksB50), ..bcols])
bdatB75 <- droplevels(mtbdat[.(blocksB75), ..bcols])
bdatB100 <- droplevels(mtbdat[.(blocksB100), ..bcols])
bdatB1000 <- droplevels(mtbdat[.(blocksB1000), ..bcols])
setkey(bdatB10, blockC)
setkey(bdatB25, blockC)
setkey(bdatB50, blockC)
setkey(bdatB75, blockC)
setkey(bdatB100, blockC)
setkey(bdatB1000, blockC)

idatB10 <- droplevels(mtwrkdat[.(blocksB10), ..icols, drop = TRUE])
idatB25 <- droplevels(mtwrkdat[.(blocksB25), ..icols, drop = TRUE])
idatB50 <- droplevels(mtwrkdat[.(blocksB50), ..icols, drop = TRUE])
idatB75 <- droplevels(mtwrkdat[.(blocksB75), ..icols, drop = TRUE])
idatB100 <- droplevels(mtwrkdat[.(blocksB100), ..icols, drop = TRUE])
idatB1000 <- droplevels(mtwrkdat[.(blocksB1000), ..icols, drop = TRUE])
setkey(idatB10, blockC)
setkey(idatB25, blockC)
setkey(idatB50, blockC)
setkey(idatB75, blockC)
setkey(idatB100, blockC)
setkey(idatB1000, blockC)

rm(mtwrkdat)

## Set up y0. It turns out that we lose control over y0 if we do not do it separately for each dataset:
set.seed(12345)
idatB10[, ub0new := rnorm(.N, mean = runif(1, min = 5, max = 20), sd = runif(1, min = 1, max = 3)), by = blockC] ## this works
idatB10[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB10 <- idatB10[, .(y0ccR2 = summary(lm(y0ContNorm ~ covscluster))$r.squared), by = blockC]
stopifnot(sum(checkB10$y0ccR2 > .99) < 10)

idatB25[, ub0new := rnorm(.N, mean = runif(1, min = 1, max = 20), sd = runif(1, min = 1, max = 3)), by = blockC] ## this works
idatB25[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB25 <- idatB25[, .(y0ccR2 = summary(lm(y0ContNorm ~ covscluster))$r.squared), by = blockC]
stopifnot(sum(checkB25$y0ccR2 > .99) < 10)

idatB50[, ub0new := rnorm(.N, mean = runif(1, min = 1, max = 20), sd = runif(1, min = 1, max = 3)), by = blockC] ## this works
idatB50[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB50 <- idatB50[, .(y0ccR2 = summary(lm(y0ContNorm ~ covscluster))$r.squared), by = blockC]
stopifnot(sum(checkB50$y0ccR2 > .99) < 10)

idatB75[, ub0new := rnorm(.N, mean = runif(1, min = 1, max = 20), sd = runif(1, min = 1, max = 3)), by = blockC] ## this works
idatB75[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB75 <- idatB75[, .(y0ccR2 = summary(lm(y0ContNorm ~ covscluster))$r.squared), by = blockC]
stopifnot(sum(checkB75$y0ccR2 > .99) < 10)

idatB100[, ub0new := rnorm(.N, mean = runif(1, min = 1, max = 20), sd = runif(1, min = 1, max = 3)), by = blockC] ## this works
idatB100[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB100 <- idatB100[, .(y0ccR2 = summary(lm(y0ContNorm ~ covscluster))$r.squared), by = blockC]
stopifnot(sum(checkB100$y0ccR2 > .99) < 10)

idatB1000[, ub0new := rnorm(.N, mean = runif(1, min = 1, max = 20), sd = runif(1, min = 1, max = 3)), by = blockC] ## this works
idatB1000[, y0new := {
  tmp <- .1 * sd(ub0new) * (vb1 == 1 & vb2 == 1) + .1 * sd(ub0new) * covscluster + ub0new
  ifelse(tmp < 0, abs(tmp), tmp)
}]
checkB1000 <- idatB1000[, .(y0ccR2 = summary(lm(y0ContNorm ~ covscluster))$r.squared), by = blockC]
stopifnot(sum(checkB1000$y0ccR2 > .99) < 10)
