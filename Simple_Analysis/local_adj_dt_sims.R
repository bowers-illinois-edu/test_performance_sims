## Assess performance of testing procedure with local adjustment using the Simes globaltest at each parent node

library(here)
library(data.table)
## Making threads =1 since We will be parallel over larger chunks below
setDTthreads(threads = 1)
library(dplyr)
library(dtplyr)
library(parallel)
library(TreeTestsSim)

machine_name <- Sys.getenv("MACHINE")
message(machine_name)
print(machine_name)

load(here("Simple_Analysis", "sim_parms.rda"), verbose = TRUE)
nrow(sim_parms)

data_path <- file.path(here(), "Simple_Analysis/CSVS_latest")
nsims <- 10000

## Get number of cores perhaps from an environment variable
if (Sys.getenv("CORES") == "" & !exists("numcores")) {
  # numcores <- detectCores() ## as.numeric(system("sysctl hw.ncpu | awk '{print $2}'",intern=TRUE))
  ncores <- future::availableCores()
  print(ncores)
} else {
  if (!exists("ncores")) {
    ncores <- as.numeric(Sys.getenv("CORES")[[1]])
  }
}
print(ncores)

## This next comes from Simple_Analysis/what_has_been_done.R
## if there is no not_done_idx.rda file then use seq_len on simparms
if (file.exists(here("Simple_Analysis", "not_done_idx.rda"))) {
  load(here("Simple_Analysis", "not_done_idx.rda"), verbose = TRUE)
  ## since we took away some rows of simparms above, we might by accident refer
  ## to an empty row if we were not careful here.
  theidx <- not_done_idx[not_done_idx %in% sim_parms$idx]
} else {
  theidx <- sim_parms$idx # seq_len(nrow(simparms))
}
str(theidx)

## Not duplicate work: send most to campus cluster and some to keeling because the CampusCluster times out after 3 days.
## Trying to keep this simple and one machine/node and use mclapply rather than
## future_apply or some job array via slurm which might be a good idea but basically annoying because of shipping objects around.

if (length(theidx) == nrow(sim_parms)) {
 ## machine_name <- Sys.getenv("MACHINE")
  if (machine_name == "CampusCluster") {
    # theidx <- theidx[1:floor(length(theidx) * .75)]
    theidx <- theidx[floor(length(theidx) * .5):length(theidx)]
  }
  if (machine_name == "Keeling") {
    theidx <- theidx[1:floor(length(theidx) * .5)]
    # theidx <- theidx[floor(length(theidx) * .75):length(theidx)]
  }
}
str(theidx)
#theidx <- theidx[1:12]
## Parallelizing the outer loop not the inner loop for most of the sims, but not all
res <- mclapply(theidx, function(i) {
##res <- lapply(theidx, function(i) {
  set.seed(12345) ## same seed for each set of parms
  ## Using the key
  parms <- sim_parms[.(i)]
  message(paste(c(i, parms[1, ]), collapse = " "))
  ptm <- proc.time()
  res <- simulate_many_runs_DT(
    n_sim = nsims,
    k = parms$k,
    max_level = parms$l,
    t = parms$prop_tau_nonzero,
    alpha = .05,
    N_total = parms$total_nodes * 100,
    beta_base = parms$beta_base,
    adj_effN = parms$adj_effN,
    local_adj_p_fn = getFromNamespace(parms[["local_adj_fn"]], ns = "TreeTestsSim"),
    global_adj = "hommel",
    alpha_method = parms$alpha_fn,
    return_details = FALSE,
    final_global_adj = parms$final_adj_method,
    ##multicore = TRUE
    multicore = FALSE
  )
  etm <- proc.time() - ptm
  parms[, names(res) := as.list(res)]
  parms$time <- etm["elapsed"]
  parms$sims <- nsims
  filename <- paste(data_path, "/sim_", paste(parms[1, 1:8], collapse = "_"), ".csv", collapse = "", sep = "")
  ## Since these simulations take a long time. Save them to disc as we go.
  data.table::fwrite(parms, file = filename)
  message(paste(c(parms[1, 1:8], parms$time), collapse = " "))
  return(parms)
 ## })
}, mc.cores = ncores, mc.preschedule = FALSE, mc.set.seed = TRUE)

## This all takes so long, that not worth actually saving.
## save(res, file = here("Simple_Analysis", paste("simple_latest_results_dt", "_", machine_name, ".rda", collapse = "")))
