## Error rate and power simulations for the siu based approach

library(here)
library(data.table)
library(manytestsr)
setDTthreads(1)

## This loads the data files and simulation parameters
source(here::here("Analysis", "load_subsetdata.R"))
RNGkind("L'Ecuyer-CMRG")
# set.seed(12345) ## mac
set.seed(1234) ## keeling
nsims <- 1000
simparms$nsims <- nsims

## For speed, break the jobs into those with fewer blocks versus more blocks.
## This should let us catch mistakes in code more quickly
simparms <- droplevels(simparms[nblocks < 1000, ])
rm(list = ls(patt = "B1000"))

## clusterswitch <- "keeling-socket"
## clusterswitch <- "keeling-future"
source(here::here("Analysis", "setup-clusters.R"), verbose = TRUE, echo = TRUE)

print(clusterswitch)
## Sometimes here() has problems with multi-node parallel processing
## So try another way outside of the loop
data_path <- file.path(here(), "Analysis/CSVS")

## Define the global objects that we want the whole cluster to know about
globals <- unique(as.vector(as.matrix(simparms[, .(idatnm, bdatnm)])))
globals <- unique(c(globals, "simparms", "nsims", "thealpha", "data_path"))
globals <- globals[globals != "NULL"]

test_and_summarize_fn <- function(i) {
    set.seed(12345)
    x <- simparms[i, ]
    xnm <- paste(x, collapse = "_")
    message(xnm)
    ptm <- proc.time()
    p_sims_tab <- padj_test_fn(
        idat = get(x[["idatnm"]]),
        bdat = get(x[["bdatnm"]]),
        blockid = "blockF",
        trtid = "trt",
        fmla = Y ~ trtF | blockF,
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
    clusterEvalQ(cl, library(here))
    clusterEvalQ(cl, library(manytestsr))
    clusterEvalQ(cl, setDTthreads(1))
    clusterExport(cl, globals)
    clusterExport(cl, "test_and_summarize_fn")
#}

## wrap in if statement
## if there is no not_done_idx.rda file then use seq_len(...
if (file.exists(here("Analysis", "not_done_idx.rda"))) {
    load(here("Analysis", "not_done_idx.rda"), verbose = TRUE)
    theidx <- not_done_idx
} else {
    theidx <- seq_len(nrow(simparms))
}

#if (clusterswitch == "keeling=-socket") {
    p_sims_res <- parLapplyLB(cl, theidx,
        fun = function(i) {
            test_and_summarize_fn(i)
        },chunk.size=1
    )
    stopCluster(cl)
#}
if (clusterswitch == "keeling=-future") {
    # options(future.globals.maxSize = +Inf)
    options(future.wait.interval = 0L) # .001)
    options(future.globals.onReference = NULL)
    p_sims_res <- future_lapply(theidx,
        FUN = function(i) {
            setDTthreads(1)
            test_and_summarize_fn(i)
        }, future.seed = 0xBEEF,
        future.packages = c("here", "manytestsr"),
        future.globals = c("test_and_summarize_fn", globals),
        future.scheduling = Inf
    )
    plan(sequential)
}
## This next for debugging
## p_sims_res <- mclapply(theidx,
##    FUN = function(i) {
##        message(i)
##        test_and_summarize_fn(i)
##    }, mc.cores = 12
## )

## Just save the expensive computation immediately using the CSVS.
## These next are not used since it takes days to run the simulation, at least for now.
save(p_sims_res, file = here::here("Analysis", "p_sims_res_2024.rda"))
siusimsres_notB1000 <- rbindlist(p_sims_res, idcol = TRUE)
save(siusimsres_notB1000, file = here::here("Analysis", "siusimsres-notB1000_2024.rda"))
