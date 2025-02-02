## Assess performance of splitting procedures with an adjustment at each node like the Simes global procedure.

library(here)
source(here("Analysis", "simple_p_draws_fns.R"))

## Example:
## If node is TRUE, this means draw from beta_b (i.e. true treatment effect)
set.seed(101)
df_dfs_half_true <- generate_tree_data_dfs(k = 4, L = 3, prop_tau_nonzero = 0.5, beta_a = 0.2, beta_b = 1, p0_min = 0)
df_dfs_half_true %>% arrange(node_id)
set.seed(101)
df_dfs_all_null <- generate_tree_data_dfs(k = 4, L = 3, prop_tau_nonzero = 0, beta_a = 0.2, beta_b = 1, p0_min = 0)
df_dfs_all_null %>% arrange(node_id)

### Compare to the BFS version:
set.seed(101)
df_bfs_half_true <- generate_tree_data_bfs(k = 4, L = 3, prop_tau_nonzero = 0.5, beta_a = 0.2, beta_b = 1, p0_min = 0)
df_bfs_half_true %>% arrange(node_id)
set.seed(101)
df_bfs_all_null <- generate_tree_data_bfs(k = 4, L = 3, prop_tau_nonzero = 0, beta_a = 0.2, beta_b = 1, p0_min = 0)
df_bfs_all_null %>% arrange(node_id)

## The two approaches to setting up the tree are the same.
stopifnot(nrow(df_bfs_half_true) == nrow(df_dfs_half_true))
stopifnot(nrow(df_bfs_all_null) == nrow(df_dfs_all_null))
stopifnot(df_bfs_half_true$label == df_dfs_half_true$label)
stopifnot(df_bfs_all_null$label == df_dfs_all_null$label)
stopifnot(df_bfs_half_true$has_false_hyp == df_dfs_half_true$has_false_hyp)
stopifnot(df_bfs_all_null$has_false_hyp == df_dfs_all_null$has_false_hyp)

## The two approaches yield different results in regards the p-values because they draw in different orders and even given the same seed the numbers come out different.
## the p_0 is the same and p_1 is the same
df_all_null <- full_join(df_bfs_all_null, df_dfs_all_null, by = c("label"), suffix = c("_bfs", "_dfs"))
df_all_null %>% dplyr::select(node_id_bfs, label, has_false_hyp_bfs, p_bfs, p_dfs)

df_half_true <- full_join(df_bfs_half_true, df_dfs_half_true, by = c("label"), suffix = c("_bfs", "_dfs"))
df_half_true %>% dplyr::select(node_id_bfs, label, has_false_hyp_bfs, p_bfs, p_dfs)

## for example, p=.03141 for lv_2 (2nd node of level 2 --- the first level
## under the root using bfs) and for lv_1_1 (the first node of the second level
## under the root, level 3). I don't have reason to believe that either is more
## or less valid at this point. Try them out in simulations.

bfs_test <- function(k, L, prop_tau_nonzero, beta_a, beta_b) {
  res <- generate_tree_data_bfs(k = k, L = L, prop_tau_nonzero = prop_tau_nonzero, beta_a = beta_a, beta_b = beta_b, p0_min = 0)
  any(res$p[!res$has_false_hyp] <= .05)
}

dfs_test <- function(k, L, prop_tau_nonzero, beta_a, beta_b) {
  res <- generate_tree_data_dfs(k = k, L = L, prop_tau_nonzero = prop_tau_nonzero, beta_a = beta_a, beta_b = beta_b, p0_min = 0)
  any(res$p[!res$has_false_hyp] <= .05)
}

## Using 1000 sims for speed
sim_se <- 2 * sqrt(.05 * (1 - .05) / 1000)


## So, BFS versus DFS doesn't matter when all are null.
set.seed(12345)
fpr_bfs <- replicate(1000, bfs_test(k = 3, L = 4, prop_tau_nonzero = 0, beta_a = .2, beta_b = 1))
mean(fpr_bfs)
fpr_dfs <- replicate(1000, dfs_test(k = 3, L = 4, prop_tau_nonzero = 0, beta_a = .2, beta_b = 1))
mean(fpr_dfs)
stopifnot(abs(mean(fpr_bfs) - mean(fpr_dfs)) < sim_se)

## Also looks correct when we have a mix of true and false nulls
set.seed(12345)
fpr_bfs_some_null <- replicate(1000, bfs_test(k = 3, L = 4, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1))
mean(fpr_bfs_some_null)
fpr_dfs_some_null <- replicate(1000, dfs_test(k = 3, L = 4, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1))
mean(fpr_dfs_some_null)
stopifnot(abs(mean(fpr_bfs_some_null) - mean(fpr_dfs_some_null)) < sim_se)

## Also the return the same kind of results when we have many more nodes per
## level. (i.e. we have a false postive rate control problem)

set.seed(12345)
## We know that the total number of nodes in an k-ary tree is given by the following:
# (k^L - 1) / (k - 1)
(10^2 - 1) / (10 - 1)
fpr_bfs_some_null_k10_l2 <- replicate(1000, bfs_test(k = 10, L = 2, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1))
mean(fpr_bfs_some_null_k10_l2)
fpr_dfs_some_null_k10_l2 <- replicate(1000, dfs_test(k = 10, L = 2, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1))
mean(fpr_dfs_some_null_k10_l2)
stopifnot(abs(mean(fpr_bfs_some_null_k10_l2) - mean(fpr_dfs_some_null_k10_l2)) < sim_se)

# (k^L - 1) / (k - 1)
(100^2 - 1) / (100 - 1)
fpr_bfs_some_null_k100_l2 <- replicate(1000, bfs_test(k = 100, L = 2, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1))
mean(fpr_bfs_some_null_k100_l2)
fpr_dfs_some_null_k100_l2 <- replicate(1000, dfs_test(k = 100, L = 2, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1))
mean(fpr_dfs_some_null_k100_l2)
stopifnot(abs(mean(fpr_bfs_some_null_k100_l2) - mean(fpr_dfs_some_null_k100_l2)) < sim_se)

## This next is slow but is here to see about timing of the dfs versus bfs approach
## Unclear exactly what is  going on here with the "system" verus "user" time here.
(2^10 - 1) / (2 - 1)
system.time(
  fpr_bfs_some_null_k2_l10 <- replicate(1000, bfs_test(k = 2, L = 10, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1))
)
mean(fpr_bfs_some_null_k2_l10)
system.time(
  fpr_dfs_some_null_k2_l10 <- replicate(1000, dfs_test(k = 2, L = 10, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1))
)
mean(fpr_dfs_some_null_k2_l10)
stopifnot(abs(mean(fpr_bfs_some_null_k2_l10) - mean(fpr_dfs_some_null_k2_l10)) < sim_se)

## We want a fairly powerful test --- that rejects the null of no effects
## frequently when it is false.

mean(rbeta(1000, .1, 1) <= .05)
mean(rbeta(1000, .2, 1) <= .05)


sim_parms <- expand.grid(
  k = sort(unique(c(seq(2, 20, 1), 100))),
  l = sort(unique(c(seq(2, 20, 2), 100))),
  prop_tau_nonzero = c(0, .1, .5, .9, 1)
)
sim_parms$total_nodes <- with(sim_parms, (k^l - 1) / (k <- 1))
sim_parms <- sim_parms %>% filter(total_nodes < 5e+6)

data_path <- file.path(here(), "CSVS")
nsims <- 1000
library(parallel)
ncores <- future::availableCores()

res <- mclapply(1:nrow(sim_parms), function(i) {
  parms <- sim_parms[i, ]
  message(paste(parms[1, ], collapse = " "))
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

save(res, file = "res_weak.rda")

## Assuming that we are using the CSVs
load(here("Analysis", "simple_sims_results.rda"), verbose = TRUE)

head(simp_simsres)
summary(simp_simsres$sims)

## with(simp_simsres,table(time_bfs>time_dfs))

simp_simsres$k_vs_l <- with(simp_simsres, k / l)
## simp_simsres$bfs_vs_dfs <- with(simp_simsres,time_bfs - time_dfs)

## simp_simsres %>% select(-file) %>% arrange(k_vs_l)
## simp_simsres %>% filter(abs(bfs_vs_dfs) > 1) %>% select(-file) %>% arrange(k_vs_l)
## simp_simsres %>% select(-file) %>% arrange(total_nodes,l,k)
## ## It seems like dfs is faster when the number of nodes is large.

## Roughly how long might it take for 1000 sims using dfs?
sum(simp_simsres$time_dfs * 100) / 60

res_dts <- lapply(res, function(lst) {
  tmp <- as.data.frame(lst)
  names(tmp)[4] <- "FWER"
  return(tmp)
})
res_dt <- bind_rows(res_dts, .id = "sim_id")



## Start here: I think that BFS is better until l > 10  for small k.
library(microbenchmark)

(10^4 - 1) / (10 - 1)
(4^10 - 1) / (4 - 1)

benches <- microbenchmark(
  bfs_k10_l2 = bfs_test(k = 10, L = 2, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1),
  dfs_k10_l2 = dfs_test(k = 10, L = 2, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1),
  bfs_k2_l10 = bfs_test(k = 2, L = 10, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1),
  dfs_k2_l10 = dfs_test(k = 2, L = 10, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1),
  bfs_k4_l4 = bfs_test(k = 4, L = 4, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1),
  dfs_k4_l4 = dfs_test(k = 4, L = 4, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1),
  bfs_k4_l6 = bfs_test(k = 4, L = 6, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1),
  dfs_k4_l6 = dfs_test(k = 4, L = 6, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1),
  bfs_k100_l2 = bfs_test(k = 100, L = 2, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1),
  dfs_k100_l2 = dfs_test(k = 100, L = 2, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1),
  bfs_k20_l2 = bfs_test(k = 20, L = 2, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1),
  dfs_k20_l2 = dfs_test(k = 20, L = 2, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1),
  bfs_k2_l15 = bfs_test(k = 2, L = 15, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1),
  dfs_k2_l15 = dfs_test(k = 2, L = 15, prop_tau_nonzero = .5, beta_a = .2, beta_b = 1),
  times = 10L
)
benches
