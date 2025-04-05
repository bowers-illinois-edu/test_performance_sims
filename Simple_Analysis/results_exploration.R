## Explore and interpret the results of the simulations. Here focusing on the
## version without adjustment at each node.

library(here)
library(dplyr)
library(dtplyr)
library(data.table)
library(conflicted)
library(ggplot2)
library(viridis)
library(xtable)
conflicts_prefer(dplyr::filter)

# load(here("Simple_Analysis", "simple_sims_unadj_results.rda"), verbose = TRUE)
# load(here("Simple_Analysis", "simple_sims_adj_results.rda"), verbose = TRUE)
load(here("Simple_Analysis", "simple_sims_latest_results.rda"), verbose = TRUE)
load(here("Simple_Analysis", "sim_parms.rda"))

## Look at progress of the simulations
with(sim_parms, table(k, l))
with(simp_simsres_latest, table(k, l))

## we use alpha=.05 below
## and maximum simulation error is this
sim_se <- 2 * sqrt(.5 * (1 - .5) / unique(simp_simsres_latest$sims))
##  Since alpha=.05, we want to be less than this number:
.05 + sim_se

## using the dtplyr approach for readablilty for those not used to data.table syntax

################### WEAK CONTROL OF THE FWER ############
## Check on weak control: looks good. No serious departures from nominal
## control when all H0 is true

key_char <- c("false_error", "power", "num_leaves_tested", "leaf_power", "num_nodes_tested")

summary_fn <- function(x) {
  if (all(is.na(x))) {
    return(data.table(summary_stat = c("min", "median", "max"), value = c(NA, NA, NA)))
  } else {
    return(data.table(
      summary_stat = c("min", "median", "max"),
      value = quantile(x, c(0, 0.5, 1), na.rm = TRUE)
    ))
  }
}

## For those simulations where the null of no effect is true, summarize the key
## characteristics across all of the tree sizes by type of adjustment and rule.

res_all_true <- simp_simsres_latest[prop_tau_nonzero == 0,
  {
    res_list <- lapply(key_char, function(nm) {
      res0 <- summary_fn(get(nm))
      res0[, variable := nm]
      return(res0)
    })
    rbindlist(res_list)
  },
  by = c("alpha_fn", "adj_effN", "local_adj_fn", "final_adj_method")
]

## We have weak control regardless of the local or global adjustment etc.
max_fwer_weak_control <- res_all_true %>%
  filter(variable %in% c("false_error") & summary_stat == "max")
stopifnot(max_fwer_weak_control$value <= .05 + sim_se)

## There is not appreciable variation in the control of the FWER here:
res_all_true %>%
  filter(variable %in% c("false_error") & summary_stat == "max") %>%
  select(value) %>%
  summary()

## The simplest rules involve fixed alpha and no local adjustment and no global final adjustment
## Notice that we do not even need data splitting

basic_weak_control_res <- simp_simsres_latest %>%
  filter(prop_tau_nonzero == 0 & alpha_fn == "fixed" &
    local_adj_fn == "local_unadj_all_ps" & final_adj_method == "none") %>%
  select(-file & one_of(c("k", "l", "adj_effN", "total_nodes", "num_leaves", key_char))) %>%
  arrange(k, l)

basic_weak_control_res

weak_control_sim_tab0 <- basic_weak_control_res %>%
  group_by(k) %>%
  summarize(
    min_l = min(l),
    max_l = max(l),
    max_fwer = max(false_error),
    max_nodes_tested = max(num_nodes_tested),
    min_total_nodes = min(total_nodes),
    max_total_nodes = max(total_nodes),
    min_leaves = as.integer(min(num_leaves)),
    max_leaves = as.integer(max(num_leaves))
  )

names(weak_control_sim_tab0) <- c("$k$", "Min", "Max", "FWER", "Tests", "Min", "Max", "Min", "Max")
# weak_control_sim_tab0[["Max Leaves"]] <- as.integer(weak_control_sim_tab0[["Max Leaves"]])

weak_control_sim_tab <- xtable(weak_control_sim_tab0, digits = 3)

caption(weak_control_sim_tab) <- "This table show weak-control of the
  family-wise error rate across 10,000 simulations for hypotheses using
  $\\alpha=.05$ on $k$-ary trees with $k$ nodes per level where the hypothesis
  of no effects is true for all nodes. Each row summarizes the results of
  simulations for the trees with a given $k$ and between min and max levels.
  Maximum average false positive rates shown in the 'Max FWER' column. The
  maximum nodes tested are shown in 'Max Tests'. And information about the
  range of trees is also shown: the total number of nodes in the trees and the
  total number of terminal nodes or leaves."

label(weak_control_sim_tab) <- "tab:weak_control_sim"

second_row <- paste(paste(names(weak_control_sim_tab0), collapse = " & "), "\\\\", collapse = "")
first_row <- "&\\multicolumn{2}{c}{Levels} &Max&Max& \\multicolumn{2}{c}{Nodes} & \\multicolumn{2}{c}{Leaves}\\\\"

align(weak_control_sim_tab) <- xalign(weak_control_sim_tab0)
# align(weak_control_sim_tab)["Max Tests"] <- "b{1cm}"
# align(weak_control_sim_tab)["Max FWER"] <- "b{1cm}"


print(weak_control_sim_tab,
  add.to.row = list(
    pos = list(0, 0),
    command = c(first_row, second_row)
  ),
  file = here("Paper", "weak_control_sim_tab.tex"),
  include.rownames = FALSE,
  include.colnames = FALSE,
  booktabs = TRUE
)

## When all hyps are false then of course no false errors.

vars2 <- c("power", "leaf_power", "num_leaves_tested", "num_leaves", "bottom_up_power", "false_error")

res_all_false <- simp_simsres_latest[prop_tau_nonzero == 1,
  {
    res_list <- lapply(vars2, function(nm) {
      res0 <- summary_fn(get(nm))
      res0[, variable := nm]
      return(res0)
    })
    rbindlist(res_list)
  },
  by = c("alpha_fn", "adj_effN", "local_adj_fn")
]

## No false errors possible given the simulation
res_all_false %>%
  filter(variable == "false_error") %>%
  summary(mean(value))


## Few big power differences across the approaches because the power of the
## individual tests is more or less the same across all of them (interesting
## that the adj_effN doesn't distinguish between them, though.)

res_all_false %>%
  filter(variable == "power" & summary_stat == "max")

res_all_false %>%
  filter(variable == "leaf_power" & summary_stat == "max")

## The top down approach tests fewer leaves but tends to nearly always
## reject the false leaves compared to the bottom up approach which fails to
## reject a lot (comparing leaf_power to bottom_up_power here).

## For example, when k=16 and l=4, there are a total of 65536 leaves (all of
## which tested by the bottom up approach) but we only test (on average) 4692
## of them. We reject false nulls at the leaf level 91% of the time.

simp_simsres_latest %>%
  filter(prop_tau_nonzero == 1) %>%
  select(one_of(c("k", "l", vars2)))

## Now what about strong control?
## here is where we need some local and/or global adjustments.
## excluding the k=2, l=2 example because when prop_tau_nonzero==.1 in that case there are no leaves set as a non-null effect
some_tau_res <- simp_simsres_latest %>%
  filter(prop_tau_nonzero > 0 & prop_tau_nonzero < 1 & (k != 2 & l != 2))

## Some definite failures to control FWER here
some_tau_res %>%
  group_by(alpha_fn, adj_effN, prop_tau_nonzero, local_adj_fn, final_adj_method) %>%
  summarize(
    mean_fwer = mean(false_error),
    max_fwer = max(false_error),
    mean_bottom_up_fwer = mean(bottom_up_false_error),
    max_bottom_up_fwer = max(bottom_up_false_error),
    mean_power = mean(power),
    min_power = min(power),
    mean_leaf_power = mean(leaf_power),
    min_leaf_power = mean(leaf_power)
  ) %>%
  print(n = 200)

## So there are multiple ways to control the FWER here and they don't look so different in power
some_tau_res %>%
  filter(false_error <= .05 + sim_se) %>%
  droplevels() %>%
  group_by(alpha_fn, adj_effN, prop_tau_nonzero, local_adj_fn, final_adj_method) %>%
  summarize(
    mean_fwer = mean(false_error),
    max_fwer = max(false_error),
    mean_bottom_up_fwer = mean(bottom_up_false_error),
    max_bottom_up_fwer = max(bottom_up_false_error),
    mean_power = mean(power),
    min_power = min(power),
    mean_leaf_power = mean(leaf_power),
    min_leaf_power = mean(leaf_power)
  ) %>%
  print(n = 200)

## We probably don't want to have both local_hommel plus conservative alpha
## adjustment at each step let alone plus final overall adjustment. Also use
## the data splitting approach.

some_tau_res %>%
  filter(local_adj_fn == "local_hommel_all_ps" &
    !(alpha_fn %in% c("adaptive_k_adj", "fixed_k_adj")) & adj_effN & final_adj_method == "none") %>%
  group_by(alpha_fn, prop_tau_nonzero, local_adj_fn) %>%
  summarize(
    mean_fwer = mean(false_error),
    max_fwer = max(false_error),
    mean_bottom_up_fwer = mean(bottom_up_false_error),
    max_bottom_up_fwer = max(bottom_up_false_error),
    mean_power = mean(power),
    min_power = min(power),
    mean_leaf_power = mean(leaf_power),
    min_leaf_power = mean(leaf_power)
  ) %>%
  print(n = 100)

## Make a data set of this (with simulated data splitting)
## Allow uncontrolled FWER to learn where we see it

strong_ctrl_res <- simp_simsres_latest %>%
  filter(adj_effN == TRUE & local_adj_fn == "local_hommel_all_ps" &
    (prop_tau_nonzero > 0 & prop_tau_nonzero < 1)) %>%
  as_tibble()

## Here we didn't get to test any leaves. The bottom up approach does 1024 tests with very low power

## But the top down approach stops early --- it tells us **correctly** that at
## least one of the leaves of a rejected parent is nonnul
## but it also says that we cannot detect the specific leaf
##

## This is just a summary of what has been run
with(strong_ctrl_res, table(k, l))

## When it does test a nonnull leaf, it always rejects
summary(strong_ctrl_res$num_leaves_tested)
## nearly always
summary(strong_ctrl_res$leaf_power)

strong_ctrl_res %>%
  filter(leaf_power < .5) %>%
  select(k, l, prop_tau_nonzero, false_error, power, num_nodes_tested, leaf_power, num_leaves, num_leaves_tested, bottom_up_false_error, bottom_up_power)

## Often stops after just a few nodes: (these are means, TODO record more than mean, like max and mid and median)
## rarely (on average) tests leaves
summary(strong_ctrl_res$num_nodes_tested)

strong_ctrl_res %>%
  filter(k == 6, l == 8) %>%
  select(prop_tau_nonzero, false_error, power, num_nodes_tested, leaf_power, num_leaves, num_leaves_tested, bottom_up_false_error, bottom_up_power)

strong_ctrl_res %>%
  filter(k == 8, l == 6) %>%
  select(prop_tau_nonzero, false_error, power, num_nodes_tested, leaf_power, num_leaves, num_leaves_tested, bottom_up_false_error, bottom_up_power)

strong_ctrl_res %>%
  filter(num_leaves_tested > 0) %>%
  select(k, l, prop_tau_nonzero, false_error, power, num_nodes_tested, leaf_power, num_leaves, num_leaves_tested, bottom_up_false_error, bottom_up_power) %>%
  arrange(k, l, prop_tau_nonzero)


plot_dat <- rbind(strong_ctrl_res, strong_ctrl_res)
plot_dat <- plot_dat %>% mutate(top_down = rep(c(1, 0), each = nrow(strong_ctrl_res)))
plot_dat <- plot_dat %>% mutate(
  fwer = ifelse(top_down == 1, false_error, bottom_up_false_error),
  num_disc = ifelse(top_down == 1, power, bottom_up_power)
)

g_tmp <- ggplot(data = plot_dat, aes(x = k, y = num_disc, groups = factor(l), color = factor(l))) +
  geom_line() +
  facet_wrap(~ prop_tau_nonzero + top_down)

print(g_tmp)

## Filling the contours at the alpha levels but plus sim error (sim_se)
g_fill_unadj <- ggplot(data = simp_simsres_unadj, aes(x = k, y = l, z = fpr_dfs)) +
  # geom_contour_filled(breaks = c(0, .05, .1, .2, .3, .4, .5, 1)) +
  geom_contour_filled(breaks = c(0, .05 + sim_se, .1 + sim_se, .2, .3, .4, .5, 1)) +
  labs(
    x = "Number of nodes per level",
    fill = "FWER"
  ) +
  facet_wrap(~prop_tau_nonzero) +
  #  scale_fill_viridis(discrete = TRUE) +
  scale_y_continuous(
    name = "Number of Levels"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16)
  ) +
  coord_cartesian(clip = "off")

print(g_fill_unadj)


## Filling the contours at the alpha levels but plus sim error (sim_se)
g_fill_adj <- ggplot(data = simp_simsres_adj, aes(x = k, y = l, z = fwer)) +
  # geom_contour_filled(breaks = c(0, .05, .1, .2, .3, .4, .5, 1)) +
  geom_contour_filled(breaks = c(0, .05 + sim_se, .1 + sim_se, .2, .3, .4, .5, 1)) +
  labs(
    x = "Number of nodes per level",
    fill = "FWER"
  ) +
  facet_wrap(~prop_tau_nonzero) +
  #  scale_fill_viridis(discrete = TRUE) +
  scale_y_continuous(
    name = "Number of Levels"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16)
  ) +
  coord_cartesian(clip = "off")

print(g_fill_adj)
## When there are no treatment effects (i.e. all hypotheses of no effects are
## true), then we control at all kinds of graphs (prop_tau_nonzero=0)

## When all nodes have a treatment effect (i.e. all hyps of no effects are
## false) then we controll all kind of graphcs (prop_tau_nonzero=1)

## So this is weak control of the FWER

## We also see some nominal FWERS --- at low k (number of nodes) and low l
## (number of levels). But as nodes and levels increase, we start making too many
## errors.

## Below we see something special: even in mixed treatment effects scenarios,
## FWER can remain below some number depending on a function of number of nodes
## and number of levels. So, for example, with fewer than 5 levels, and 2 nodes
## per level, the FWER is fairly well controlled.

## This is also evident in the contour plots --- there is a kind of smooth
## function that we are seeing here with a kind of limit or boundary when k is very large.

g_lines_nodes_unadj <- ggplot(
  data = simp_simsres_unadj,
  aes(groups = l, color = factor(l), x = k, y = fpr_dfs)
) +
  facet_wrap(~prop_tau_nonzero) +
  geom_line() +
  geom_hline(yintercept = .05)

g_lines_nodes_unadj

g_lines_nodes_adj <- ggplot(
  data = simp_simsres_adj,
  aes(groups = l, color = factor(l), x = k, y = fwer)
) +
  facet_wrap(~prop_tau_nonzero) +
  geom_line() +
  geom_hline(yintercept = .05)

g_lines_nodes_adj


g_lines_levels_unadj <- ggplot(
  data = simp_simsres_unadj,
  aes(groups = k, color = factor(k), x = l, y = fpr_dfs)
) +
  geom_line() +
  facet_wrap(~prop_tau_nonzero) +
  geom_hline(yintercept = .05)

g_lines_levels_unadj


g_lines_levels_adj <- ggplot(
  data = simp_simsres_adj,
  aes(groups = k, color = factor(k), x = l, y = fwer)
) +
  geom_line() +
  facet_wrap(~prop_tau_nonzero) +
  geom_hline(yintercept = .05)

g_lines_levels_adj


## Look at the latest

### First, weak control

simp_simsres_latest %>%
  filter(prop_tau_nonzero == 0) %>%
  group_by(adj_effN, local_adj_fn) %>%
  reframe(summary(fwer))

## with no possible false positives (prop_tau_nonzero==1) or all true null
## (prop_tau_nonzero=0) we have control at all levels
simp_simsres_unadj %>%
  filter(prop_tau_nonzero %in% c(0, 1)) %>%
  select(-file) %>%
  reframe(summary(fpr_dfs))
#  summary(fpr_dfs)
# 1       0.00000000
# 2       0.00000000
# 3       0.01750000
# 4       0.02488636
# 5       0.04800000
# 6       0.06800000


simp_simsres_unadj %>%
  filter(prop_tau_nonzero > 0 & prop_tau_nonzero < 1) %>%
  select(-file) %>%
  arrange(fpr_dfs)

simp_simsres_adj %>%
  filter(prop_tau_nonzero > 0 & prop_tau_nonzero < 1) %>%
  select(-file) %>%
  arrange(fwer)

## What is the function that is f(k,l,prop_tau_nonzero) -> fwer (i.e. fdr_dfs) here alpha=.05
simp_simsres_unadj %>%
  # filter(prop_tau_nonzero > 0 & prop_tau_nonzero < 1 & k < 4 & l < 15) %>%
  filter(prop_tau_nonzero == .5 & k < 4 & l < 15) %>%
  select(-file) %>%
  select(k, l, prop_tau_nonzero, fpr_dfs) %>%
  arrange(prop_tau_nonzero, k, l, fpr_dfs)



system("touch Simple_Analysis/simple_results_exploration.done")
