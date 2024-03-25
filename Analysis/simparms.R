## Make a  set of parameters for simulations that can be shared across files

library(data.table)
library(manytestsr)

## these are all the possibilities. We will not use them all.
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
  stringsAsFactors = FALSE,
  block_sizes = c("equal","not_equal")
)
setDT(simparms)
dim(simparms)

simparms[, idatnm:=paste("idat_",block_sizes,"_B", simparms$nblocks, sep = "")]
simparms[, bdatnm:=paste("bdat_",block_sizes,"_B", simparms$nblocks, sep = "")]

## The specified factor splitters use a categorial covariate with a specific structure like state.district.school
simparms[, splitby := ifelse(splitfn %in% c("splitSpecifiedFactor", "splitSpecifiedFactorMulti"), "covsplits", "hwt")]
## There is no point of splitting off a single largest block (or block with largest statistical power)
simparms[, splitby := ifelse(splitfn=="splitLOO" & block_sizes=="equal", "vary0new", "hwt")]
simparms[splitfn == "splitCluster", splitby := "covscluster"]
## Allow treatment effects to vary by a covariate called "covscluster" if we use one of the tau_fn_covariates... functions
simparms[, covars := "covscluster"]

## Keep only one row for each  prop_blocks_0=1 since all effects will be null
simparms <- droplevels(simparms[(prop_blocks_0 < 1 | (prop_blocks_0 == 1 & tau_sizes == 0)), ])

## Do not test the splitEqual or pWilcox approaches
simparms <- droplevels(simparms[splitfn != "splitEqual" & pfn != "pWilcox", ])
## Do not iterate over splitters or  afns etc when p_adj_method is not "split"
simparms[, splitfn := ifelse(p_adj_method != "split", "NULL", splitfn)]
simparms[, afn := ifelse(p_adj_method != "split", "NULL", afn)]
simparms[, splitby := ifelse(p_adj_method != "split", "NULL", splitby)]

## Remove duplicated rows
## We should just have one row per combination of nblocks,tau_sizes, prop_blocks_0, stopsplitting for any given p_adjustment method and block_sizes argument
simparms[, dups := duplicated(simparms)]
table(simparms$dups)

sp1 <- droplevels(simparms[!(dups), ])
## Test that the duplication removal actually works
tmp1 <- simparms[p_adj_method == "fdr" & nblocks == 10 & tau_sizes == .5 & prop_blocks_0 == 0 & pfn == "pOneway" & stopsplitting == TRUE & block_sizes=="equal", ]

tmp2 <- sp1[p_adj_method == "fdr" & nblocks == 10 & tau_sizes == .5 & prop_blocks_0 == 0 & pfn == "pOneway" & stopsplitting == TRUE & block_sizes=="equal", ]

stopifnot(nrow(sp1) < nrow(simparms))
stopifnot(nrow(tmp2) == 1)
stopifnot(tmp2 == tmp1[1, ])

sp2 <- sp1[, dups := NULL]

## Also, we basically want to stop splitting when the splitby variable is
## constant. So, stop assessing what happens if you keep splitting at random.
simparms <- droplevels(sp2[(stopsplitting), ])

## Make the simulations faster --- finish a lot of small and fast stuff first
## so that we can prototype plots and tables and such.
setorder(simparms,-prop_blocks_0,nblocks)

simparms[, idx := seq_len(.N)]
setkey(simparms, idx)

head(simparms)
tail(simparms)

thealpha <- .05

save(thealpha, simparms, file = here::here("Analysis", "simparms.rda"))

