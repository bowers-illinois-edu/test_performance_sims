## Error rate and power simulations for the siu based approach
## Focusing on the equal block sized approach here

library(here)
library(data.table)
library(manytestsr)
setDTthreads(1)

## This loads the data files and simulation parameters
load(here("Data","equal_blocksize_data.rda"),verbose=TRUE)
load(here("Analysis","simparms.rda"),verbose=TRUE)
RNGkind("L'Ecuyer-CMRG")
set.seed(12345)
nsims <- 1000
simparms$nsims <- nsims

## For speed, break the jobs into those with fewer blocks versus more blocks.
## This should let us catch mistakes in code more quickly
## Also focus on smaller datasets for now and for now also just equal sized blocks
simparms <- droplevels(simparms[nblocks < 1000 & block_sizes=="equal", ])
rm(list = ls(patt = "B1000"))

## Now set up the cluster environment
## clusterswitch <- "keeling-socket"
## clusterswitch <- "keeling-future"
source(here::here("Analysis", "setup-clusters.R"), verbose = TRUE, echo = TRUE)

print(clusterswitch)
cl_exists <- exists("cl")
blah <- plan()
plan_exists <- !inherits(try(class(blah)),"try-error")

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
    set.seed(12345)
    x <- simparms[i, ]
    ## nmsx0 <- names(x)
    ##nmsx <- grep('unequal',nmsx0,invert=TRUE,value=TRUE)
    ##xnm <- paste(x[,.SD,.SDcols=nmsx], collapse = "_")
    xnm <- paste(x,collapse="_")
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

#if (clusterswitch == "keeling-socket") {
if (cl_exists & !plan_exists) {
    clusterEvalQ(cl, library(here))
    clusterEvalQ(cl, library(manytestsr))
    clusterEvalQ(cl, setDTthreads(1))
    clusterExport(cl, globals)
    clusterExport(cl, "test_and_summarize_fn")
}

## If the run has stopped at some point don't redo all of the past
## simulations. See the Makefile for code that collects the CSV names and
## figures out which rows of simparms have already been completed and which
## have not yet been done

## if there is no not_done_idx.rda file then use seq_len on simparms
if (file.exists(here("Analysis", "not_done_idx.rda"))) {
    load(here("Analysis", "not_done_idx.rda"), verbose = TRUE)
    theidx <- not_done_idx
} else {
    theidx <- seq_len(nrow(simparms))
}

## In general each row of simparms will require different amounts of time to
## run. So below we try different approaches to load balancing. Still not sure
## if any of them are really great at not leaving cores unused while the main
## process waits for things to finish.

#if (clusterswitch == "keeling=-socket") {
if (cl_exists & !plan_exists) {
    p_sims_res <- parLapplyLB(cl, theidx,
        fun = function(i) {
            test_and_summarize_fn(i)
        },chunk.size=1
    )
    stopCluster(cl)
}

# if (clusterswitch == "keeling=-future") {
if (plan_exists){
    ## here we are trying to load balance since each set of parameters can require a very different amount of time to run
    options(future.globals.maxSize = +Inf)
    options(future.wait.interval = 0L) # .001)
    options(future.globals.onReference = NULL)

    p_sims_res <- future_lapply(theidx,
        FUN = function(i) {
            setDTthreads(1)
            test_and_summarize_fn(i)
        }, future.seed = 0xBEEF,
        future.packages = c("here", "manytestsr"),
        future.globals = c("test_and_summarize_fn", globals),
        future.scheduling = Inf,
        future.chunk.size=NULL
    )

    plan(sequential)
}
##  This next for debugging
## p_sims_res <- mclapply(theidx,
##     FUN = function(i) {
##         message(i)
##         test_and_summarize_fn(i)
##     }, mc.cores = 12
## )

## p_sims_res <- lapply(theidx,
##     FUN = function(i) {
##         message(i)
##         test_and_summarize_fn(i)
##     }
## )

## Just save the expensive computation immediately using the CSVS.
## These next are not used since it takes days to run the simulation, at least for now.
## But we leave them here since maybe they will be used at some other time
save(p_sims_res, file = here::here("Analysis", "p_sims_res_equal_nb_2024.rda"))
#siusimsres_notB1000 <- rbindlist(p_sims_res, idcol = TRUE)
#save(siusimsres_notB1000, file = here::here("Analysis", "siusimsres-notB1000_2024.rda"))
