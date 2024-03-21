## Make a  set of parameters for simulations that can be shared across files

library(data.table)
library(manytestsr)

nblocks <- c(10, 25, 50, 75, 100, 1000)
taus <- c(0, .1, .25, .5, .75, 1, 2, 3)
prop_blocks_0 <- c(0, .1, .2, .5, .8, .9, 1)
pfn <- c("pOneway", "pWilcox", "pIndepDist")
splittingfns <- c("splitLOO", "splitEqualApprox", "splitCluster", "splitSpecifiedFactor", "splitSpecifiedFactorMulti")
afn <- c("alpha_investing", "alpha_saffron", "alpha_addis", "NULL")
padjmethods <- c("split", "fdr", "holm", "hommel")
stopsplitting <- c(TRUE, FALSE)

simparms <- expand.grid(
  nblocks = nblocks,
  tau_sizes = taus,
  prop_blocks_0 = prop_blocks_0,
  pfn = pfn,
  splitfn = splittingfns,
  afn = afn,
  p_adj_method = padjmethods,
  stopsplitting = stopsplitting,
  stringsAsFactors = FALSE
)
simparms$idatnm <- paste("idatB", simparms$nblocks, sep = "")
simparms$bdatnm <- paste("bdatB", simparms$nblocks, sep = "")
setDT(simparms)
dim(simparms)
simparms[, splitby := ifelse(splitfn %in% c("splitSpecifiedFactor", "splitSpecifiedFactorMulti"), "covsplits", "hwt")]
simparms[splitfn == "splitCluster", splitby := "covscluster"]
## Allow treatment effects to vary by a covariate called "covscluster" if we use one of the tau_fn_covariates... functions
simparms[, covars := "covscluster"] ## ifelse(splitfn=="splitCluster","covcluster","")]
simparms <- simparms[order(simparms$prop_blocks_0, decreasing = TRUE), ] ## assess fwer or tau=0 first##
## Keep only one row for each  prop_blocks_0=1 since all effects will be null
## simparms[, keeprows := (prop_blocks_0 < 1 | (prop_blocks_0 == 1 & tau_sizes == 0))]
simparms <- droplevels(simparms[(prop_blocks_0 < 1 | (prop_blocks_0 == 1 & tau_sizes == 0)), ])
head(simparms)
tail(simparms)

## Do not test the splitEqual or pWilcox approaches
simparms <- droplevels(simparms[splitfn != "splitEqual" & pfn != "pWilcox", ])
## simparms$reveal_and_test_fn=rep("reveal_po_and_test_siup",nrow(simparms))
## simparms <- rbind(simparms,c("NULL","NULL","NULL","reveal_po_and_test"))

## Do not iterate over splitters or  afns etc when p_adj_method is not "split"
simparms <- simparms[, splitfn := ifelse(p_adj_method != "split", "NULL", splitfn)]
simparms <- simparms[, afn := ifelse(p_adj_method != "split", "NULL", afn)]
simparms <- simparms[, splitby := ifelse(p_adj_method != "split", "NULL", splitby)]
## simparms <- simparms[,covars:=ifelse(p_adj_method!="split","NULL",covars)]

## Remove duplicated rows
## We should just have one row per combination of nblocks,tau_sizes, prop_blocks_0, stopsplitting for any given p_adjustment method
simparms[, dups := duplicated(simparms)]
table(simparms$dups)
sp1 <- simparms[!(dups), ]
tmp1 <- simparms[p_adj_method == "fdr" & nblocks == 10 & tau_sizes == .5 & prop_blocks_0 == 0 & pfn == "pOneway" & stopsplitting == TRUE, ]
tmp2 <- sp1[p_adj_method == "fdr" & nblocks == 10 & tau_sizes == .5 & prop_blocks_0 == 0 & pfn == "pOneway" & stopsplitting == TRUE, ]
tmp3 <- sp1[p_adj_method == "split" & nblocks == 10 & tau_sizes == .5 & prop_blocks_0 == 0 & pfn == "pOneway" & stopsplitting == TRUE, ]

stopifnot(nrow(sp1) < nrow(simparms))
stopifnot(nrow(tmp2) == 1)
stopifnot(tmp2 == tmp1[1, ])

thealpha <- .05

sp2 <- sp1[, dups := NULL]

simparms <- droplevels(sp2[(stopsplitting), ])

## Make the simulations faster --- finish a lot of stuff first so that we can prototype plots and tables and such.
simparms <- simparms[order(nblocks), ]
simparms[, idx := seq_len(.N)]
setkey(simparms, idx)
save(thealpha, simparms, file = here::here("Analysis", "simparms.rda"))
