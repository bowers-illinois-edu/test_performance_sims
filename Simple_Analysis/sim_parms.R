## Make the set of simulation parameters

library(here)
library(dplyr)
library(data.table)
library(dtplyr)

alpha_methods <- c("fixed", "fixed_k_adj", "adaptive_k_adj", "spending", "investing")
final_adj_methods <- c("none", "fdr", "fwer")
local_adj_methods <- c("local_hommel_all_ps", "local_bh_all_ps","local_unadj_all_ps")
adj_effN <- c(TRUE, FALSE)
beta_bases <- c(.001,.01,.05,.1,.2)
theks <- sort(unique(c(seq(2, 10, 2), 50, 100)))
thels <- sort(unique(c(seq(2, 10, 2), 50, 100)))
theeffects <- c(0, .1, .5, .9, 1)

sim_parms <- expand.grid(
  k = theks,
  l = thels,
  prop_tau_nonzero = theeffects,
  adj_effN = adj_effN,
  local_adj_fn = local_adj_methods,
  alpha_fn = alpha_methods,
  final_adj_method = final_adj_methods,
  beta_base = beta_bases
)
sim_parms$total_nodes <- with(sim_parms, ((k^(l + 1)) - 1) / (k - 1))
sim_parms$local_adj_fn <- as.character(sim_parms$local_adj_fn)

## Restrict the simulation to cases where the number of nodes is not gigantic
sim_parms <- sim_parms %>%
  filter(total_nodes < 1e+6) %>%
  arrange(total_nodes) %>%
  as.data.table()
nrow(sim_parms)
sim_parms[, idx := seq_len(.N)]
setkey(sim_parms, idx)
## Don't know why these are factor and not character from the expand.grid call
sim_parms[, alpha_fn := as.character(alpha_fn)]
sim_parms[, final_adj_method := as.character(final_adj_method)]

with(sim_parms,table(k,l))

#### Looking at the times required for the last runs
#### Not super clear but not recording whether run on campus cluster (100+ cores) versus keeling (24 cores)
##simp_simsres_latest$mins  <- simp_simsres_latest$time/60
##simp_simsres_latest %>%
##  group_by(k,l) %>%
##  summarize(min(mins),median(mins),mean(mins),max(mins)) %>%
##  print(n=100)
##
##
##simp_simsres_latest %>% filter(mins > 120) %>%
##  select(k,l,mins,prop_tau_nonzero,adj_effN,local_adj_fn,alpha_fn,final_adj_method,total_nodes) %>%
##  arrange(mins)
##
#### Let's try to keep the simulations less than an hour per:
##slow_sims <- simp_simsres_latest %>% filter(mins > 60 & k %in% theks & l %in% thels) %>%
##  select(k,l,mins,prop_tau_nonzero,adj_effN,local_adj_fn,alpha_fn,final_adj_method,total_nodes)
##
##with(slow_sims,table(k,l))
##with(sim_parms,table(k,l))
##
##slow_sims %>%  arrange(mins)

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

