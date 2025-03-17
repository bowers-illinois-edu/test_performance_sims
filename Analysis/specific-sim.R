# siusims.R modified to take a command line argument for the
# specific line of simparms.[rda|csv] to use.
# R-cluster handling has been removed to make this a truly
# serial script.


suppressMessages(library(here))
library(data.table)
suppressMessages(library(manytestsr))
setDTthreads(0)
print(paste('MVM, getDTthreads:', getDTthreads()))

argv <- commandArgs(trailingOnly=TRUE)

rowNum = as.integer(argv[[1]])
print(paste('MVM, rowNum', rowNum))
## This loads the data files and simulation parameters
load(here("Data", "equal_blocksize_data.rda"))
load(here("Data", "unequal_blocksize_data.rda"))
load(here("Analysis", "simparms.rda"))
RNGkind("L'Ecuyer-CMRG")
set.seed(12345)
nsims <- 40 # originally 1000
simparms$nsims <- nsims

## For speed, break the jobs into those with fewer blocks versus more blocks.
## This should let us catch mistakes in code more quickly
## Also focus on smaller datasets for now and for now also just equal sized blocks
#simparms <- droplevels(simparms[nblocks < 1000 & block_sizes == "equal", ])
#rm(list = ls(patt = "B1000"))

## Sometimes here() has problems with multi-node parallel processing
## So try another way outside of the loop
data_path <- file.path(here(), "Analysis/CSVS")

## Define the global objects that we want the whole cluster to know about
globals <- unique(as.vector(as.matrix(simparms[, .(idatnm, bdatnm)])))
globals <- unique(c(globals, "simparms", "nsims", "thealpha", "data_path"))
globals <- globals[globals != "NULL"]

## Define the function that will do the nsims simulations for each set of
## parameters defined in a row of simparms where i is the row of simparms

test_and_summarize_fn <- function(x) {
    set.seed(12345)
    ## simparms is keyed by idx so this is the row with idx==i
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
    fwrite(res, file = paste0(data_path, "/", xnm, ".csv", collapse = "", sep = ""))
    return(res)
}

#myrow <- simparms[.(rowNum)]
Rprof(
      filename=here::here("Analysis",
                          paste0("profile-",rowNum,"-",nsims,".prof")),
      line.profiling=TRUE,
      gc.profiling=TRUE,
      filter.callframes=TRUE )
p_sims_res <- test_and_summarize_fn(simparms[simparms$idx == rowNum,])
Rprof(NULL)


# MVM: save or append to file.
save(p_sims_res, file = here::here("Analysis", paste0("mvm-test-",rowNum,".rda")))
