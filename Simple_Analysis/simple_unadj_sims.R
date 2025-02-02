## Assess performance of testing procedure without adjustment other
## than following the monotonicity rule

library(here)
source(here("Simple_Analysis", "simple_p_draws_fns.R"))

dfs_test <- function(k, L, prop_tau_nonzero, beta_a, beta_b) {
  res <- generate_tree_data_dfs(k = k, L = L, prop_tau_nonzero = prop_tau_nonzero, beta_a = beta_a, beta_b = beta_b, p0_min = 0)
  any(res$p[!res$has_false_hyp] <= .05)
}


sim_parms <- expand.grid(
  k = sort(unique(c(seq(2, 20, 1), 100))),
  l = sort(unique(c(seq(2, 20, 2), 100))),
  prop_tau_nonzero = c(0, .1, .5, .9, 1)
)
sim_parms$total_nodes <- with(sim_parms, (k^l - 1) / (k <- 1))
## Restrict the simulation to cases where the number of nodes is not gigantic
sim_parms <- sim_parms %>% filter(total_nodes < 5e+6)

data_path <- file.path(here(), "Simple_Analysis/CSVS_unadj")
nsims <- 1000
library(parallel)
ncores <- future::availableCores()

res <- mclapply(1:nrow(sim_parms), function(i) {
  parms <- sim_parms[i, ]
  message(paste(parms[1, ], collapse = " "))
  ## The commented code here is to compare the bfs versus dfs approaches. The DFS approaches are much faster in general
  ## ptm_bfs <- proc.time()
  ## res_bfs <- replicate(nsims, bfs_test(k = parms$k, L = parms$l, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1))
  ## etm_bfs <- proc.time() - ptm_bfs
  ptm_dfs <- proc.time()
  res_dfs <- replicate(nsims, dfs_test(k = parms$k, L = parms$l, prop_tau_nonzero = parms$prop_tau_nonzero, beta_a = .1, beta_b = 1))
  etm_dfs <- proc.time() - ptm_dfs
  ## fpr_bfs <- mean(res_bfs)
  fpr_dfs <- mean(res_dfs)
  ## parms$fpr_bfs <- fpr_bfs
  parms$fpr_dfs <- fpr_dfs
  ## parms$time_bfs <- etm_bfs["elapsed"]
  parms$time_dfs <- etm_dfs["elapsed"]
  parms$sims <- nsims
  filename <- paste(data_path, "/sim_", paste(parms[1, 1:3], collapse = "_"), ".csv", collapse = "", sep = "")
  ## Since these simulations take a long time. Save them to disc as we go.
  data.table::fwrite(parms, file = filename)
  return(parms)
}, mc.cores = ncores - 1, mc.set.seed = TRUE)

save(res, file = here("Simple_Analysis", "simple_raw_unadj_results.rda"))
