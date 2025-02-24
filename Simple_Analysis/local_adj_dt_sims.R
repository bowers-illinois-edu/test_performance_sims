## Assess performance of testing procedure with local adjustment using the Simes globaltest at each parent node

library(here)
library(data.table)
## Making threads =1 since We will be parallel over larger chunks below
setDTthreads(threads = 1)
library(dplyr)
library(parallel)
library(TreeTestsSim)

sim_parms <- expand.grid(
  k = sort(unique(c(seq(2, 20, 2), 50, 100))),
  l = sort(unique(c(seq(2, 20, 2), 50, 100))),
  prop_tau_nonzero = c(0, .1, .5, .9, 1),
  adj_effN = c(TRUE, FALSE),
  local_adj_fn = c("local_simes", "local_hommel_all_ps", "local_unadj_all_ps"),
  alpha_fn = c("fixed", "spending", "investing")
)
sim_parms$total_nodes <- with(sim_parms, (k^l - 1) / (k <- 1))
sim_parms$local_adj_fn <- as.character(sim_parms$local_adj_fn)

## Restrict the simulation to cases where the number of nodes is not gigantic
sim_parms <- sim_parms %>%
  filter(total_nodes < 5e+6) %>%
  as.data.table()
nrow(sim_parms)

data_path <- file.path(here(), "Simple_Analysis/CSVS_latest")
nsims <- 10000
ncores <- future::availableCores()

## We are using alpha = .05 and beta_params = c(.1,1) (very high power at each
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

## nsims <- 10
res <- mclapply(seq_len(nrow(sim_parms)), function(i) {
  ## res <- lapply(seq_len(nrow(sim_parms)), function(i) {
  set.seed(12345) ## same seed for each set of parms
  parms <- sim_parms[i, ]
  message(paste(i, parms[1, ], collapse = " "))
  ptm <- proc.time()
  res <- simulate_many_runs_DT(
    n_sim = nsims,
    k = parms$k,
    max_level = parms$l,
    t = parms$prop_tau_nonzero,
    alpha = .05, N_total = parms$total_nodes * 100,
    beta_base = .1,
    adj_effN = parms$adj_effN,
    local_adj_p_fn = getFromNamespace(parms[["local_adj_fn"]], ns = "TreeTestsSim"),
    global_adj = "hommel",
    return_details = FALSE
  )
  etm <- proc.time() - ptm
  parms$fwer <- mean(res)
  parms$time <- etm["elapsed"]
  parms$sims <- nsims
  parms[, names(res) := as.list(res)]
  filename <- paste(data_path, "/sim_", paste(parms[1, 1:6], collapse = "_"), ".csv", collapse = "", sep = "")
  ## Since these simulations take a long time. Save them to disc as we go.
  data.table::fwrite(parms, file = filename)
  message(paste(c(parms[1, 1:6], parms$time), collapse = " "))
  return(parms)
  # })
}, mc.cores = ncores - 1, mc.set.seed = TRUE)

save(res, file = here("Simple_Analysis", "simple_latest_results_dt.rda"))
