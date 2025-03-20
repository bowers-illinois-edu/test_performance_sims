## Assess performance of testing procedure with local adjustment using the Simes globaltest at each parent node

library(here)
library(data.table)
## Making threads =1 since We will be parallel over larger chunks below
setDTthreads(threads = 1)
library(dplyr)
library(dtplyr)
library(parallel)
library(TreeTestsSim)

alpha_methods <- c("fixed", "fixed_k_adj", "adaptive_k_adj", "spending", "investing")
final_adj_methods <- c("none", "fdr", "fwer")
local_adj_methods <- c("local_hommel_all_ps", "local_unadj_all_ps")
adj_effN <- c(TRUE, FALSE)

sim_parms <- expand.grid(
  k = sort(unique(c(seq(2, 20, 2), 50, 100))),
  l = sort(unique(c(seq(2, 20, 2), 50, 100))),
  prop_tau_nonzero = c(0, .1, .5, .9, 1),
  adj_effN = adj_effN,
  local_adj_fn = local_adj_methods,
  alpha_fn = alpha_methods,
  final_adj_method = final_adj_methods
)
sim_parms$total_nodes <- with(sim_parms, ((k^(l + 1)) - 1) / (k - 1))
sim_parms$local_adj_fn <- as.character(sim_parms$local_adj_fn)

## Restrict the simulation to cases where the number of nodes is not gigantic
sim_parms <- sim_parms %>%
  filter(total_nodes < 5e+6) %>%
  arrange(total_nodes) %>%
  as.data.table()
nrow(sim_parms)
sim_parms[, idx := seq_len(.N)]
setkey(sim_parms, idx)
## Don't know why these are factor and not character from the expand.grid call
sim_parms[, alpha_fn := as.character(alpha_fn)]
sim_parms[, final_adj_method := as.character(final_adj_method)]

save(sim_parms, file = here("Simple_Analysis", "sim_parms.rda"))

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


## We are using beta_params = c(.1,1) (very high power at each
## node, all non-null nodes same high power) This last is very unrealistic
## since we are going to actually split the data at each node. So that N at
## each parent is always larger than N at any child.

## Setting N_total very large just means that the power of the rbeta continues
## to go down as the tree grows rather than stabilizing at something like 1.

# mean(rbeta(10000, .1, 1) <= .05)
# ## Notice that this is runif
# mean(rbeta(10000, 1, 1) <= .05)
# ## This is slightly more powerful than the dist under the null
# mean(rbeta(10000, .9, 1) <= .05)
# mean(rbeta(10000, .5, 1) <= .05)
# ## This is closer to the cap
# mean(rbeta(10000, 1.64, 1) <= .05)

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

## Not duplicate work: send most to campus cluster and some to keeling because the CampusCluster times out after 3 days.
## Trying to keep this simple and one machine/node and use mclapply rather than
## future_apply or some job array via slurm which might be a good idea but basically annoying because of shipping objects around.

if (length(theidx) == nrow(sim_parms)) {
  machine_name <- Sys.getenv("MACHINE")
  if (machine_name == "CampusCluster") {
    theidx <- theidx[1:floor(length(theidx) * .75)]
  }
  if (machine_name == "Keeling") {
    theidx <- theidx[floor(length(theidx) * .75):length(theidx)]
  }
}

## Parallelizing the outer loop not the inner loop
res <- mclapply(theidx, function(i) {
  ## res <- lapply(seq_len(nrow(sim_parms)), function(i) {
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
    beta_base = .1,
    adj_effN = parms$adj_effN,
    local_adj_p_fn = getFromNamespace(parms[["local_adj_fn"]], ns = "TreeTestsSim"),
    global_adj = "hommel",
    alpha_method = parms$alpha_fn,
    return_details = FALSE,
    final_global_adj = x$final_adj_method,
    alpha_method = x$alpha_method,
    multicore = FALSE
  )
  etm <- proc.time() - ptm
  parms[, names(res) := as.list(res)]
  parms$time <- etm["elapsed"]
  parms$sims <- nsims
  filename <- paste(data_path, "/sim_", paste(parms[1, 1:7], collapse = "_"), ".csv", collapse = "", sep = "")
  ## Since these simulations take a long time. Save them to disc as we go.
  data.table::fwrite(parms, file = filename)
  message(paste(c(parms[1, 1:7], parms$time), collapse = " "))
  return(parms)
  # })
}, mc.cores = ncores, mc.preschedule = FALSE, mc.set.seed = TRUE)

save(res, file = here("Simple_Analysis", paste("simple_latest_results_dt", "_", MACHINE, ".rda", collapse = "")))
