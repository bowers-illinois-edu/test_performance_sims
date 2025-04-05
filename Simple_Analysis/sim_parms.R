## Make the set of simulation parameters

library(here)
library(dplyr)
library(data.table)
library(dtplyr)

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

