# specific-sims.R modified to use mclapply

suppressMessages(library(here))
library(data.table)
suppressMessages(library(manytestsr))
library(parallel)
setDTthreads(1)

## This loads the data files and simulation parameters
load(here("Data", "equal_blocksize_data.rda"))
load(here("Data", "unequal_blocksize_data.rda"))
load(here("Analysis", "simparms.rda"))
RNGkind("L'Ecuyer-CMRG")
set.seed(12345)
nsims <- 1000
simparms$nsims <- nsims

## For speed, break the jobs into those with fewer blocks versus more blocks.
## This should let us catch mistakes in code more quickly
## Also focus on smaller datasets for now and for now also just equal sized blocks

simparms <- droplevels(simparms[nblocks < 1000 & block_sizes == "equal" & splitfn != "splitLOO", ])

# testing just nblocks == 10
#simparms <- droplevels(simparms[nblocks == 10])
rm(list = ls(patt = "B1000"))
#rm(list = ls(patt = "B100"))
#rm(list = ls(patt = "B75"))
#rm(list = ls(patt = "B50"))
#rm(list = ls(patt = "B25"))

## Sometimes here() has problems with multi-node parallel processing
## So try another way outside of the loop
data_path <- file.path(here(), "Analysis/CSVS")

## Define the global objects that we want the whole cluster to know about
globals <- unique(as.vector(as.matrix(simparms[, .(idatnm, bdatnm)])))
globals <- unique(c(globals, "simparms", "nsims", "thealpha", "data_path"))
globals <- globals[globals != "NULL"]

## Define the function that will do the nsims simulations for each set of
## parameters defined in a row of simparms where i is the row of simparms

test_and_summarize_fn <- function(i) {
    starttime = Sys.time()
    # MVM: can this be hoisted?
    set.seed(12345)
    ## simparms is keyed by idx so this is the row with idx==i
    # simparms is a data.table
    # MVM - these are not actually equivalent syntax!
    # simparms[47,], simparms[47], simparms[idx == 47], simparms[.(47)]
    x <- simparms[i,]
    simidx <- x[["idx"]]
    ## nmsx0 <- names(x)
    ## nmsx <- grep('unequal',nmsx0,invert=TRUE,value=TRUE)
    ## xnm <- paste(x[,.SD,.SDcols=nmsx], collapse = "_")
    xnm <- paste(x, collapse = "_")
    message(xnm)
    ptm <- proc.time()
    p_sims_tab <- padj_test_fn(
        idat = get(x[["idatnm"]]),
        bdat = get(x[["bdatnm"]]),
        blockid = "bF",
        trtid = "trt",
        fmla = Y ~ trtF | bF,
        ybase = "y0new",
        prop_blocks_0 = x[["prop_blocks_0"]],
        tau_fn = tau_norm_covariate_outliers,
        tau_size = x[["tau_sizes"]],
        covariate = x[["covars"]],
        pfn = getFromNamespace(x[["pfn"]], ns = "manytestsr"),
        nsims = nsims,
        ncores = 1, ## parallelize over the  rows of simparms
        afn = ifelse(x[["afn"]] != "NULL", getFromNamespace(x[["afn"]], ns = "manytestsr"), "NULL"),
        splitfn = ifelse(x[["splitfn"]] != "NULL", getFromNamespace(x[["splitfn"]], ns = "manytestsr"), "NULL"),
        splitby = x[["splitby"]],
        p_adj_method = x[["p_adj_method"]],
        stop_splitby_constant = x[["stopsplitting"]]
    )
    etm <- proc.time() - ptm
    ## Since these simulations take a long time. Save them to disc as we go.
    res_rates <- p_sims_tab[, lapply(.SD, mean, na.rm = TRUE)]
    res <- cbind(x, res_rates, elapsed = etm["elapsed"])
    # fwrite(res, file = here::here("Analysis/CSVS", paste(xnm, ".csv", sep = "")))
    # MVM: commenting this out for profiling purposes.
    #fwrite(res, file = paste0(data_path, "/", xnm, ".csv", collapse = "", sep = ""))
    # MVM: this might be from start of the root process since mclapply forks, IIUC
    endtime = Sys.time()
    message('MVM,',i,',',simidx,',',starttime,',',endtime,',',etm["elapsed"])
    return(res)
}

#MVM: simparms is a data.table
#MVM: mclapply will coerce the first object into a list, so can't send
# a slice of rows... sigh.
p_sims_res <- mclapply(1:dim(simparms)[[1]], test_and_summarize_fn, mc.preschedule=FALSE, mc.cores=128)

# MVM: save or append to file.
save(p_sims_res, file = here::here("Analysis", "mclapply-test-small-equal-noLOO-nopresched-all.rda"))
