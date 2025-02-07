## Assess performance of testing procedure with local adjustment using the Simes globaltest at each parent node

library(here)
library(dplyr)
library(parallel)
library(TreeTestsSim)

sim_parms <- expand.grid(
  k = sort(unique(c(seq(2, 20, 1), 100))),
  l = sort(unique(c(seq(2, 20, 2), 100))),
  prop_tau_nonzero = c(0, .1, .5, .9, 1)
)
sim_parms$total_nodes <- with(sim_parms, (k^l - 1) / (k <- 1))
## Restrict the simulation to cases where the number of nodes is not gigantic
sim_parms <- sim_parms %>% filter(total_nodes < 5e+6)

data_path <- file.path(here(), "Simple_Analysis/CSVS_adj")
nsims <- 1000
ncores <- future::availableCores()

## We are using alpha = .05 and beta_params = c(.1,1) (very high power at each
## node, all non-null nodes same high power) This last is very unrealistic
## since we are going to actually split the data at each node. So that N at
## each parent is always larger than N at any child.

res <- mclapply(1:nrow(sim_parms), function(i) {
  set.seed(12345) ## same seed for each set of parms
  parms <- sim_parms[i, ]
  message(paste(parms[1, ], collapse = " "))
  ptm <- proc.time()
  res <- simulate_many_runs(n_sim = nsims, k = parms$k, L = parms$l, t = parms$prop_tau_nonzero)
  etm <- proc.time() - ptm
  parms$fwer <- mean(res)
  parms$time <- etm["elapsed"]
  parms$sims <- nsims
  filename <- paste(data_path, "/sim_", paste(parms[1, 1:3], collapse = "_"), ".csv", collapse = "", sep = "")
  ## Since these simulations take a long time. Save them to disc as we go.
  data.table::fwrite(parms, file = filename)
  message(paste(c(parms[1, 1:3], parms$time), collapse = " "))
  return(parms)
}, mc.cores = ncores - 1, mc.set.seed = TRUE)

save(res, file = here("Simple_Analysis", "simple_adj_results.rda"))
