## Make some testing data for use in the tests in this directory
## These data sets are the large ones and we will subsample them to assess performance with smaller numbers of blocks

library(randomizr)
library(data.table)
library(dplyr)
library(dtplyr)
library(fabricatr)
setDTthreads(1)

nblocks <- 1000
nunits <- 30000

set.seed(12345)
## For now, create a very simple canceling out arrangement and simple design setup
idat <- data.table(i = 1:nunits, b = rep(c(1:nblocks), length = nunits))
bdat <- data.table(b = 1:nblocks)
setkey(bdat, b)
setkey(idat, b)

## Make vb1, vb2, vb3 all nest within each other like states, counties, districts (with the individual blocks nested within all three)
## 2 levels at highest level
bdat[, vb1 := rep(1:2, length.out = .N)]
## 2 levels within each state
bdat[, vb2 := rep(1:2, length.out = .N), by = vb1]
bdat[, vb3 := rep(1:2, length.out = .N), by = interaction(vb1, vb2, lex.order = TRUE, drop = TRUE)]
ftable(bdat$vb1, bdat$vb2)
ftable(bdat$vb1, bdat$vb2, bdat$vb3)

## Now merge the block data onto the individual level data
## make potential outcome to control/status quo depend on a covariate at both the individual and block levels.
idat <- bdat[idat]

set.seed(12345)
idat[, vi1 := round(runif(.N, min = 0, max = 10)), by = b]
idat[, ub0 := draw_normal_icc(mean = 0, N = .N, clusters = b, ICC = .1)]
idat[, y0 := vi1 + vb1 + ub0 + rnorm(.N), by = b]
idat[, cc := .05 * sd(y0) * y0 + ub0, by = b]

## Check out the relationships between potential outcome to control and covariates
blah <- idat[, .(
  y0ccR2 = summary(lm(y0 ~ cc))$r.squared,
  y0R2 = summary(lm(y0 ~ vb1 + vi1 + cc))$r.squared,
  nb = .N,
  bary0 = mean(y0),
  barub0 = mean(ub0),
  sdy0 = sd(y0),
  sdub0 = sd(ub0),
  sdcc = sd(cc)
), by = b]
summary(blah)
stopifnot(sum(blah$y0ccR2 > .99) < 10)
blah[y0R2 > .80, ]
summary(blah[y0R2 > .80, ])
rm(blah)

idat[, tau := ifelse(b <= nblocks / 2, -5 * sd(y0), 5 * sd(y0))]
idat[, tauhomog := sd(y0) * 5]
idat[, taunormb := rnorm(.N, mean = (sd(y0) / b), sd = sd(y0) / 2), by = b]
idat[, y1 := y0 + tau]
idat[, y1homog := y0 + tauhomog]
idat[, y1normb := y0 + taunormb]
stopifnot(all.equal(mean(idat$y1 - idat$y0), 0))
## Randomly assign half to treatment
idat[, trt := complete_ra(N = nunits / nblocks), by = b]
testra <- idat[, sum(trt), by = b]
stopifnot(all(testra$V1 == (nunits / nblocks) / 2))
idat[, Y := trt * y1 + (1 - trt) * y0]
## Now ensure that the no effects case is really no effects
idat[, y0null := resid(lm(y0 ~ trt, data = .SD)), by = b]
idat[, y1null := y0null]
idat[, Ynull := trt * y1null + (1 - trt) * y0null]
stopifnot(all.equal(idat$Ynull, idat$y0null))
idat[, Yhomog := trt * y1homog + (1 - trt) * y0]
idat[, Ynormb := trt * y1normb + (1 - trt) * y0]
idat$trtF <- factor(idat$trt)
idat$bF <- factor(idat$b)
idat$gF <- interaction(idat$bF, idat$trtF)
## Created aligned versions of the  variables
idat[, Ymd := Y - mean(Y), by = b]
idat[, trtmd := trt - mean(trt), by = b]

## Observed Effects by Block
idat[, list(
  canceling = mean(Y[trt == 1] - Y[trt == 0]),
  null = mean(Ynull[trt == 1] - Ynull[trt == 0]),
  homog = mean(Yhomog[trt == 1] - Yhomog[trt == 0]),
  normb = mean(Ynormb[trt == 1] - Ynormb[trt == 0])
), by = b]

bdat_tmp <- idat[, .(
  nb = .N,
  nt = sum(trt),
  nc = sum(1 - trt),
  pb = mean(trt),
  covscluster = mean(cc)
), by = b]

bdat_tmp[, hwt := (2 * (nc * nt) / (nc + nt))]
setkey(bdat_tmp, b)
bdat <- bdat_tmp[bdat]
bdat$bF <- factor(bdat$b)
rm(bdat_tmp)

setkey(idat, bF)
setkey(bdat, bF)

#### Now make data with unequal sized blocks, but no block over nb=100,
## also unequal tau (causal effects)
bdat2 <- data.table(nb = rep(seq(4, 400, length = 100), length = nblocks), bF = factor(1:nblocks))
table(bdat2$nb)

## Make vb1, vb2, vb3 all nest within each other like states, counties, districts (with the individual blocks nested within all three)
## 2 levels at highest level
bdat2[, vb1 := rep(1:2, length.out = .N)]
## 2 levels within each state
bdat2[, vb2 := rep(1:2, length.out = .N), by = vb1]
bdat2[, vb3 := rep(1:2, length.out = .N), by = interaction(vb1, vb2, lex.order = TRUE, drop = TRUE)]
ftable(bdat2$vb1, bdat2$vb2)
ftable(bdat2$vb1, bdat2$vb2, bdat2$vb3)

setkey(bdat2, "bF")
setkey(idat, "bF")

## Make a dataset with unequal sized blocks
set.seed(12345)
idat2 <- data.table(bF = rep(as.character(bdat2$bF), bdat2$nb))
## Add the vb variables to the indiv level data
idat3 <- bdat2[idat2]
ftable(idat3$vb1, idat3$vb2, idat3$vb3)
idat3[, bF := as.factor(bF)]
idat3[, nb := .N, by = bF]
## Now giving it the same blocking names as idat for ease with testing
setkey(idat3, bF)
idat3 <- idat3[, b := as.numeric(as.character(bF))]
## Make more indiv level covariates
idat3[, c("v1", "v2", "v3", "v4") := lapply(1:4, function(i) {
  rnorm(.N, mean = .N / i, sd = 1 / i)
}), by = bF]

set.seed(12345)
idat3[, ub0 := draw_normal_icc(mean = 0, N = .N, clusters = b, ICC = .1)]
idat3[, y0 := v1 + vb1 + ub0 + rnorm(.N), by = b]
idat3[, cc := .05 * sd(y0) * y0 + ub0, by = b]

## Check out the relationships between potential outcome to control and covariates
blah <- idat3[, .(
  y0ccR2 = summary(lm(y0 ~ cc))$r.squared,
  y0R2 = summary(lm(y0 ~ vb1 + v1 + cc))$r.squared,
  nb = .N,
  bary0 = mean(y0),
  barub0 = mean(ub0),
  sdy0 = sd(y0),
  sdub0 = sd(ub0),
  sdcc = sd(cc)
), by = b]
summary(blah)
stopifnot(sum(blah$y0ccR2 > .99) < 10)
blah[y0R2 > .90, ]
blah[y0ccR2 > .90, ]
summary(blah[y0R2 > .90, ])
summary(blah[y0ccR2 > .90, ])

idat3[, tau := ifelse(b <= nblocks, -5 * sd(y0), 5 * sd(y0))]
idat3[, tauhomog := 5 * sd(y0)]
## Taunorm_inc increases the effect with the block number and also size of the block.
## We should be increasingly likely to detect effects for blocks with higher numbers.
idat3[, taunorm_inc := rnorm(.N, mean = b / nblocks * sd(y0), sd = .5), by = bF]
idat3[, taunorm_inc := ifelse(taunorm_inc < 0, 0, taunorm_inc)]
## taunorm_dec has effect size decreasing in block size
idat3[, taunorm_dec := rnorm(.N, mean = 4 / b * sd(y0), sd = .5), by = bF]
idat3[, taunorm_dec := ifelse(taunorm_dec < 0, 0, taunorm_dec)]
idat3[, sd(y0), by = b]
## tauv2 has no systematic relationship between block number or block size and effect size
idat3[, tauv2 := ifelse(v2 < quantile(v2, sample(c(.9, .1), size = 1)), sd(y0) * (taunorm_inc) + sd(y0) * runif(.N, min = .25, max = 2), 0), by = bF]
idat3[, lapply(.SD, mean), .SDcols = grep("tau|nb|y0", names(idat3), value = TRUE), by = b]
idat3[, y1 := y0 + tau]
idat3[, y1null := y0]
idat3[, y1homog := y0 + tauhomog]
set.seed(1235) ## set y1norm_inc=y0 for 1/4 of the blocks
nullbF <- sample(unique(levels(idat3$bF)), size = length(unique(idat3$bF)) / 4)
idat3[, y1norm_inc := ifelse(bF %in% nullbF, y0, y0 + taunorm_inc)]
idat3[, y1norm_dec := ifelse(bF %in% nullbF, y0, y0 + taunorm_dec)]
idat3[, y1tauv2 := y0 + tauv2]

## Randomly assign treatment but allow probability of assignment to vary by block
set.seed(1235)
tmp <- idat3[, .(nb = .N), by = bF]
mb <- round(tmp$nb * runif(nblocks, min = .3, max = .7))
idat3[, trt := block_ra(blocks = bF, block_m = mb)]
idat3[, Y := trt * y1 + (1 - trt) * y0]
idat3[, Ynull := trt * y1null + (1 - trt) * y0]
idat3[, Yhomog := trt * y1homog + (1 - trt) * y0]
idat3[, Ynorm_inc := trt * y1norm_inc + (1 - trt) * y0]
idat3[, Ynorm_dec := trt * y1norm_dec + (1 - trt) * y0]
idat3[, Ytauv2 := trt * y1tauv2 + (1 - trt) * y0]
idat3$trtF <- factor(idat3$trt)
idat3[, .(.N, mean(y1tauv2 - y0)), by = bF]
bdat3 <- idat3[, .(
  nb = .N, nt = sum(trt), nc = sum(1 - trt), pb = mean(trt),
  v1 = mean(v1), v2 = mean(v2), v3 = mean(v3), v4 = mean(v4),
  covscluster = mean(cc), barub0 = mean(ub0),
  ate_tau = mean(y1 - y0),
  ate_null = mean(y1null - y0),
  ate_homog = mean(y1homog - y0),
  ate_norm_inc = mean(y1norm_inc - y0),
  ate_norm_dec = mean(y1norm_dec - y0),
  ate_tauv2 = mean(y1tauv2 - y0),
  vb1 = unique(vb1),
  vb2 = unique(vb2),
  vb3 = unique(vb3)
), by = bF]
bdat3[, hwt := (2 * (nc * nt) / (nc + nt))]

###  A factor or categorial covariate for pre-determined splitting functions
bdat[, covsplits := interaction(vb1, vb2, vb3, lex.order = TRUE, drop = TRUE)]
bdat3[, covsplits := interaction(vb1, vb2, vb3, lex.order = TRUE, drop = TRUE)]
idat[, covsplits := interaction(vb1, vb2, vb3, lex.order = TRUE, drop = TRUE)]
idat3[, covsplits := interaction(vb1, vb2, vb3, lex.order = TRUE, drop = TRUE)]
## covscluster
idat[, covscluster := cc]
idat3[, covscluster := cc]

## Just comparing the two datasets: one with equal numbers of people per block and the other with heterogeneous block sizes
table(idat$bF)
table(idat3$bF)
sort(table(idat3$bF))
nrow(idat)
nrow(idat3)

## idat = "individual level data"
## bdat = "block level data" (mostly aggregated from the idat)
idat_equal_nb <- idat
idat_not_equal_nb <- idat3
bdat_equal_nb <- bdat
bdat_not_equal_nb <- bdat3

## Creating a character version of bF because I think keys in data.table are better as characters but I'm not sure.
idat_equal_nb[, bC := as.character(bF)]
bdat_equal_nb[, bC := as.character(bF)]
idat_not_equal_nb[, bC := as.character(bF)]
bdat_not_equal_nb[, bC := as.character(bF)]
setkey(idat_equal_nb, bC)
setkey(bdat_equal_nb, bC)
setkey(idat_not_equal_nb, bC)
setkey(bdat_not_equal_nb, bC)

save(idat_equal_nb, file = here::here("Data", "idat_equal_nb.rda"))
save(idat_not_equal_nb, file = here::here("Data", "idat_not_equal_nb.rda"))
save(bdat_equal_nb, file = here::here("Data", "bdat_equal_nb.rda"))
save(bdat_not_equal_nb, file = here::here("Data", "bdat_not_equal_nb.rda"))

system("touch Data/make_test_data.done")
