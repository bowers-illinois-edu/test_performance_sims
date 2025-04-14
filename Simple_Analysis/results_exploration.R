## Explore and interpret the results of the simulations. Here focusing on the
## version without adjustment at each node.

library(here)
library(tidyverse)
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

## Make the table with the non-data splitting versions
weak_control_sim_tab0 <- basic_weak_control_res %>%
  filter(adj_effN == FALSE) %>%
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
first_row <- "&\\multicolumn{2}{c}{Levels} &Max&Avg. Max& \\multicolumn{2}{c}{Nodes} & \\multicolumn{2}{c}{Leaves}\\\\"

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


####################################################
## Now what about strong control?
## here is where we need some local and/or global adjustments.
## excluding the k=2, l=2 example because when prop_tau_nonzero==.1 in that case there are no leaves set as a non-null effect
some_tau_res <- simp_simsres_latest %>%
  filter(prop_tau_nonzero > 0 & prop_tau_nonzero < 1)

## The simple approach does not control FWER in a strong sense
some_tau_naive <- simp_simsres_latest %>%
  filter(prop_tau_nonzero == .5 &
    alpha_fn == "fixed" & local_adj_fn == "local_unadj_all_ps" & final_adj_method == "none") %>%
  select(-file & one_of(c("k", "l", "adj_effN", "total_nodes", "num_leaves", "prop_tau_nonzero", key_char)))

some_tau_naive

strong_control_naive_tab0 <- some_tau_naive %>%
  group_by(k, adj_effN) %>%
  summarize(
    min_l = min(l),
    max_l = max(l),
    min_fwer = min(false_error),
    max_fwer = max(false_error),
    max_nodes_tested = max(num_nodes_tested),
    min_total_nodes = min(total_nodes),
    max_total_nodes = max(total_nodes),
    min_leaves = as.integer(min(num_leaves)),
    max_leaves = as.integer(max(num_leaves))
  ) %>%
  arrange(adj_effN, k)

print(strong_control_naive_tab0, n = 100)

names(strong_control_naive_tab0) <- c("$k$", "Splitting", "Min", "Max", "Min", "Max", "Tests", "Min", "Max", "Min", "Max")
# strong_control_naive_tab0[["Max Leaves"]] <- as.integer(strong_control_naive_tab0[["Max Leaves"]])

strong_control_naive_tab <- xtable(strong_control_naive_tab0, digits = 2)

caption(strong_control_naive_tab) <- "This table shows  lack of strong-control
  of the family-wise error rate across 10,000 simulations for hypotheses using
  $\\alpha=.05$ on $k$-ary trees with $k$ nodes per level where the hypothesis
  of no effects is false for half of the leaves. Any ancestor of a leaf where
  the null hypothesis of no effects is false is also false by construction.
  Each row summarizes the results of simulations for the trees with a given $k$
  and between min and max levels. Minimum and maximum average false positive
  rates shown in the 'Min' and 'Max' FWER columns. The average maximum nodes
  tested across the 10,000 simulatiosn are shown in 'Avg. Max Tests'. And
  information about the range of trees is also shown: the total number of nodes
  in the trees and the total number of terminal nodes or leaves."

label(strong_control_naive_tab) <- "tab:strong_control_naive"

second_row <- paste(paste(names(strong_control_naive_tab0), collapse = " & "), "\\\\", collapse = "")
first_row <- "&Data&\\multicolumn{2}{c}{Levels} &\\multicolumn{2}{c}{FWER} & Avg. Max &
  \\multicolumn{2}{c}{Nodes} & \\multicolumn{2}{c}{Leaves}\\\\"

align(strong_control_naive_tab) <- xalign(strong_control_naive_tab0)

print(strong_control_naive_tab,
  add.to.row = list(
    pos = list(0, 0),
    command = c(first_row, second_row)
  ),
  file = here("Paper", "strong_control_naive_tab.tex"),
  include.rownames = FALSE,
  include.colnames = FALSE,
  booktabs = TRUE
)


### Start Here TODO

##### Focus on what does control the FWER in a strong sense. For now only
### focusing on prop_tau_nonzero==.5 just to keep it simple.


## Approaches that do control the FWER in the strong sense
some_tau_control0 <- simp_simsres_latest %>%
  filter(prop_tau_nonzero == .5 & false_error <= (.05 + sim_se)) %>%
  select(-file & one_of(c("k", "l", "adj_effN", "alpha_fn", "local_adj_fn", "final_adj_method", "total_nodes", "num_leaves", "prop_tau_nonzero", "bottom_up_power", key_char)))

some_tau_control0

## Focus first on single approaches
with(some_tau_control0, table(alpha_fn, local_adj_fn, final_adj_method))

## Try to focus on the application of single approaches
some_tau_control1 <- some_tau_control0 %>%
  filter((alpha_fn != "fixed" & final_adj_method == "none" & local_adj_fn == "local_unadj_all_ps") |
    (alpha_fn == "fixed" & final_adj_method == "none" & local_adj_fn == "local_hommel_all_ps") |
    (alpha_fn == "fixed" & final_adj_method != "none" & local_adj_fn == "local_unadj_all_ps"))

some_tau_control1
with(some_tau_control1, table(alpha_fn, local_adj_fn, final_adj_method))

strong_control_tab0 <- some_tau_control1 %>%
  group_by(k, adj_effN, alpha_fn, local_adj_fn, final_adj_method) %>%
  summarize(
    min_l = min(l),
    max_l = max(l),
    min_fwer = min(false_error),
    max_fwer = max(false_error),
    min_power = min(power),
    max_power = max(power),
    min_leaf_power = min(leaf_power, na.rm = TRUE),
    max_leaf_power = max(leaf_power, na.rm = TRUE),
    max_nodes_tested = max(num_nodes_tested),
    min_total_nodes = min(total_nodes),
    max_total_nodes = max(total_nodes),
    min_leaves_tested = min(num_leaves_tested),
    max_leaves_tested = max(num_leaves_tested),
    min_leaves = as.integer(min(num_leaves)),
    max_leaves = as.integer(max(num_leaves)),
    min_bottom_up_power = min(bottom_up_power),
    max_bottom_up_power = max(bottom_up_power)
  ) %>%
  ungroup() %>%
  arrange(adj_effN, k)

## this next from stackoverflow
## df_transpose <- function(df) {
##   first_name <- colnames(df)[1]
##   temp <-
##     df %>%
##     tidyr::pivot_longer(-1) %>%
##     tidyr::pivot_wider(names_from = 1, values_from = value)
##   colnames(temp)[1] <- first_name
##   return(temp)
## }

## Look at the parts
### One approach is to adjust all of the p-values at a given level (this is the local hommel approach)
## tmp <- strong_control_tab0 %>%
##   filter(adj_effN == TRUE & final_adj_method == "none" & alpha_fn == "fixed") %>%
##   select(-adj_effN, -final_adj_method, -alpha_fn, -local_adj_fn, -max_power, -max_leaf_power) %>%
##   df_transpose() %>%
##   column_to_rownames(var = "k") %>%
##   zapsmall() %>%
##   select(`2`, `4`, `10`, `20`, `100`)
## tmp

local_hommel <- strong_control_tab0 %>%
  filter(local_adj_fn == "local_hommel_all_ps" & adj_effN == TRUE) %>%
  select(
    k, min_l, max_l, max_fwer, min_power, max_nodes_tested, min_leaf_power,
    max_leaves_tested, min_bottom_up_power, max_leaves # , alpha_fn
  ) %>%
  mutate(bottom_up_detected_leaves = min_bottom_up_power * (max_leaves * .5))
local_hommel
# # A tibble: 12 × 11
#        k min_l max_l max_fwer min_power max_nodes_tested min_leaf_power max_leaves_tested min_bottom_up_power max_leaves detected_leaves
#    <int> <int> <int>    <dbl>     <dbl>            <dbl>          <dbl>             <dbl>               <dbl>      <int>           <dbl>
#  1     2     2    18   0.0327     0.738             3.29              1             0.478            0.000437     262144           57.3
#  2     4     2     8   0.0317     0.730             2.91              1             0.451            0.000877      65536           28.7
#  3     6     2     6   0.0346     0.739             2.73              1             0.412            0.00104       46656           24.2
#  4     8     2     6   0.0317     0.732             2.70              1             0.455            0.000437     262144           57.3
#  5    10     2     6   0.0296     0.741             2.64              1             0.454            0.000224    1000000          112.
#  6    12     2     4   0.027      0.748             2.56              1             0.456            0.00156       20736           16.1
#  7    14     2     4   0.0269     0.737             2.49              1             0.467            0.00114       38416           21.9
#  8    16     2     4   0.0251     0.737             2.37              1             0.450            0.000876      65536           28.7
#  9    18     2     4   0.0189     0.731             2.28              1             0.438            0.000693     104976           36.4
# 10    20     2     4   0.0177     0.732             2.14              1             0.389            0.000558     160000           44.7
# 11    50     2     2   0.0175     0.74              2.39              1             0.630            0.00450        2500            5.62
# 12   100     2     2   0.0269     0.739             3.25              1             1.21             0.00225       10000           11.3


global_approach <- strong_control_tab0 %>%
  filter(local_adj_fn == "local_unadj_all_ps" & adj_effN == TRUE & final_adj_method != "none" & alpha_fn == "fixed") %>%
  select(
    final_adj_method, k, min_l, max_l, max_fwer, min_power, max_nodes_tested, min_leaf_power,
    max_leaves_tested, min_bottom_up_power, max_leaves
  ) %>%
  mutate(bottom_up_detected_leaves = min_bottom_up_power * (max_leaves * .5))
global_approach


global_approach_raw <- simp_simsres_latest %>%
  filter(final_adj_method != "none" & alpha_fn == "fixed" & prop_tau_nonzero > 0 & local_adj_fn == "local_unadj_all_ps" & adj_effN == TRUE) %>%
  select(-file & one_of(c(
    "k", "l", "adj_effN", "alpha_fn", "local_adj_fn",
    "final_adj_method", "total_nodes", "num_leaves", "prop_tau_nonzero", "bottom_up_power", key_char
  ))) %>%
  arrange(false_error)

global_approach_raw %>% filter(false_error <= (.05 + sim_se))
global_approach_raw %>% filter(false_error > (.05 + sim_se))

##   df_transpose() %>%
##   column_to_rownames(var = "k") %>%
##   zapsmall() %>%
##   select(`2`, `4`, `10`, `20`, `100`)
names(strong_control_tab0) <- c("$k$", "Splitting", "Min", "Max", "Min", "Max", "Tests", "Min", "Max", "Min", "Max")
# strong_control_tab0[["Max Leaves"]] <- as.integer(strong_control_tab0[["Max Leaves"]])

strong_control_tab <- xtable(strong_control_tab0, digits = 2)

caption(strong_control_tab) <- "This table shows  lack of strong-control
  of the family-wise error rate across 10,000 simulations for hypotheses using
  $\\alpha=.05$ on $k$-ary trees with $k$ nodes per level where the hypothesis
  of no effects is false for half of the leaves. Any ancestor of a leaf where
  the null hypothesis of no effects is false is also false by construction.
  Each row summarizes the results of simulations for the trees with a given $k$
  and between min and max levels. Minimum and maximum average false positive
  rates shown in the 'Min' and 'Max' FWER columns. The average maximum nodes
  tested across the 10,000 simulatiosn are shown in 'Avg. Max Tests'. And
  information about the range of trees is also shown: the total number of nodes
  in the trees and the total number of terminal nodes or leaves."

label(strong_control_tab) <- "tab:strong_control"

second_row <- paste(paste(names(strong_control_tab0), collapse = " & "), "\\\\", collapse = "")
first_row <- "&Data&\\multicolumn{2}{c}{Levels} &\\multicolumn{2}{c}{FWER} & Avg. Max &
  \\multicolumn{2}{c}{Nodes} & \\multicolumn{2}{c}{Leaves}\\\\"

align(strong_control_tab) <- xalign(strong_control_tab0)

print(strong_control_tab,
  add.to.row = list(
    pos = list(0, 0),
    command = c(first_row, second_row)
  ),
  file = here("Paper", "strong_control_tab.tex"),
  include.rownames = FALSE,
  include.colnames = FALSE,
  booktabs = TRUE
)



## The simple approach does have some cases with control and mostly requires data splitting.

some_tau %>%
  filter(false_error < .05 + sim_se) %>%
  arrange(k, l, prop_tau_nonzero)


## there are multiple ways to control the FWER here and they don't look so different in power
tmp <- some_tau_res %>%
  filter(false_error <= .05 + sim_se & adj_effN == FALSE) %>%
  droplevels() %>%
  group_by(alpha_fn, prop_tau_nonzero, local_adj_fn, final_adj_method) %>%
  summarize(
    mean_fwer = mean(false_error),
    max_fwer = max(false_error),
    #    mean_bottom_up_fwer = mean(bottom_up_false_error),
    #    max_bottom_up_fwer = max(bottom_up_false_error),
    mean_power = mean(power),
    min_power = min(power),
    mean_leaf_power = mean(leaf_power),
    min_leaf_power = mean(leaf_power)
  )

print(tmp, n = 200)

tmp %>% filter(final_adj_method == "fdr" & local_adj_fn == "local_unadj_all_ps" & alpha_fn == "fixed")
tmp %>% filter(final_adj_method == "fwer" & local_adj_fn == "local_unadj_all_ps" & alpha_fn == "fixed")
tmp %>% filter(final_adj_method == "none" & local_adj_fn == "local_unadj_all_ps" & alpha_fn == "fixed")

## We probably don't want to have both local_hommel plus conservative alpha
## adjustment at each step let alone plus final overall adjustment. Also use
## the data splitting approach.

some_tau_res %>%
  filter((local_adj_fn == "local_hommel_all_ps" & alpha_fn == "fixed") | (local_adj_fn != "local_hommel_all_ps" & alpha_fn != "fixed") &
    adj_effN == TRUE) %>%
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
  filter(adj_effN == TRUE & (prop_tau_nonzero > 0 & prop_tau_nonzero < 1)) %>%
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
  filter(k == 6, l == 6) %>%
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
