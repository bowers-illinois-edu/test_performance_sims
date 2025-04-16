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

caption(weak_control_sim_tab) <- "This table shows weak-control of the
  family-wise error rate across 10,000 simulations using $\\alpha=.05$ on
  $k$-ary trees with $k$ nodes per level where the hypothesis of no effects is
  true for all nodes. Each row summarizes the results of simulations for the
  trees with a given $k$ and between the minimum and maximum number of levels
  for that number of nodes per level. Maximum average false positive rates are
  shown in the 'Max FWER' column. The maximum number of nodes tested are shown
  in 'Max Tests'. And information about the range of trees is also shown: the
  total number of nodes in the trees and the total number of terminal nodes or
  leaves."

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

caption(strong_control_naive_tab) <- "This table shows lack of strong-control
  of the family-wise error rate across 10,000 simulations for hypotheses using
  $\\alpha=.05$ on $k$-ary trees with $k$ nodes per level where the hypothesis
  of no effects is false for half of the leaves. Any ancestor of a leaf where
  the null hypothesis of no effects is false is also false by construction. No
  local or global adjustments to the $p$-values or rejection thresholds were
  applied for the simulations in this table. Each row summarizes the results of
  simulations for the trees with a given $k$ and between min and max levels.
  Minimum and maximum average false positive rates across all of the levels for
  a given $k$ are shown in the 'Min' and 'Max' FWER columns. The average
  maximum nodes tested across the 10,000 simulatiosn are shown in 'Avg. Max
  Tests'. And information about the range of trees is also shown: the total
  number of nodes in the trees and the total number of terminal nodes or
  leaves."

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
  booktabs = TRUE,
  table.placement = "H"
)

some_tau_naive %>%
  filter(k == 2) %>%
  arrange(adj_effN, k, l)
# Index: <prop_tau_nonzero>
#         k     l prop_tau_nonzero adj_effN total_nodes num_nodes_tested false_error  power num_leaves_tested leaf_power num_leaves
#     <int> <int>            <num>   <lgcl>       <int>            <num>       <num>  <num>             <num>      <num>      <num>
#  1:     2     2              0.5    FALSE           7           2.7354      0.0490 0.7410            0.6429          1          4
#  2:     2     4              0.5    FALSE          31           5.7143      0.0751 0.7455            0.7358          1         16
#  3:     2     6              0.5    FALSE         127          13.2254      0.1404 0.7386            1.1086          1         64
#  4:     2     8              0.5    FALSE         511          26.7850      0.2081 0.7367            2.0409          1        256
#  5:     2    10              0.5    FALSE        2047          55.5525      0.3244 0.7422            3.9855          1       1024
#  6:     2    12              0.5    FALSE        8191         103.5610      0.4173 0.7401            7.0598          1       4096
#  7:     2    14              0.5    FALSE       32767         191.1142      0.4811 0.7403           12.6456          1      16384
#  8:     2    16              0.5    FALSE      131071         343.4361      0.5274 0.7498           21.7593          1      65536
#  9:     2    18              0.5    FALSE      524287         607.7383      0.5519 0.7389           37.3772          1     262144
# 10:     2     2              0.5     TRUE           7           2.5191      0.0452 0.7373            0.5693          1          4
# 11:     2     4              0.5     TRUE          31           3.5107      0.0346 0.7364            0.2367          1         16
# 12:     2     6              0.5     TRUE         127           4.1352      0.0097 0.7477            0.0266          1         64
# 13:     2     8              0.5     TRUE         511           4.2031      0.0009 0.7456            0.0018          1        256
# 14:     2    10              0.5     TRUE        2047           4.1435      0.0000 0.7380            0.0001          1       1024
# 15:     2    12              0.5     TRUE        8191           4.1718      0.0000 0.7418            0.0000         NA       4096
# 16:     2    14              0.5     TRUE       32767           4.1442      0.0000 0.7468            0.0000         NA      16384
# 17:     2    16              0.5     TRUE      131071           4.1641      0.0000 0.7391            0.0000         NA      65536
# 18:     2    18              0.5     TRUE      524287           4.1380      0.0000 0.7319            0.0000         NA     262144


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
  filter(final_adj_method != "none" & alpha_fn == "fixed" & (prop_tau_nonzero > 0 & prop_tau_nonzero < 1) & local_adj_fn == "local_unadj_all_ps" & adj_effN == TRUE) %>%
  select(-file & one_of(c(
    "k", "l", "adj_effN", "alpha_fn", "local_adj_fn",
    "final_adj_method", "total_nodes", "num_leaves", "prop_tau_nonzero", "bottom_up_power", key_char
  ))) %>%
  arrange(false_error)

global_approach_raw %>% filter(false_error <= (.05 + sim_se))
global_approach_raw %>% filter(false_error > (.05 + sim_se))

global_control_tab0 <- global_approach_raw %>%
  group_by(k, final_adj_method) %>%
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
  select(-max_leaf_power) %>%
  arrange(k)

print(global_control_tab0, n = 100)
# # A tibble: 24 × 18
#        k final_adj_method min_l max_l min_fwer max_fwer min_power max_power min_leaf_power max_nodes_tested min_total_nodes max_total_nodes min_leaves_tested max_leaves_tested min_leaves max_leaves min_bottom_up_power max_bottom_up_power
#    <int> <chr>            <int> <int>    <dbl>    <dbl>     <dbl>     <dbl>          <dbl>            <dbl>           <int>           <int>             <dbl>             <dbl>      <int>      <int>               <dbl>               <dbl>
#  1     2 fdr                  2    18   0         0.069    NA        NA                  1             4.20               7          524287            0                  0.810          4     262144           NA                   NA
#  2     2 fwer                 2    18   0         0.069    NA        NA                  1             4.20               7          524287            0                  0.810          4     262144           NA                   NA
#  3     4 fdr                  2     8   0.0016    0.162     0.730     0.742              1             6.59              21           87381            0.0086             1.70          16      65536            0.000873             0.106
#  4     4 fwer                 2     8   0.0016    0.162     0.730     0.742              1             6.59              21           87381            0.0086             1.70          16      65536            0.000873             0.106
#  5     6 fdr                  2     6   0.0364    0.289     0.736     0.747              1            12.3               43           55987            0.346              2.70          36      46656            0.00102              0.0385
#  6     6 fwer                 2     6   0.0364    0.289     0.736     0.747              1            12.3               43           55987            0.346              2.70          36      46656            0.00102              0.0385
#  7     8 fdr                  2     6   0.0974    0.481     0.732     0.741              1            29.0               73          299593            0.954              7.94          64     262144            0.000437             0.0287
#  8     8 fwer                 2     6   0.0974    0.481     0.732     0.741              1            29.0               73          299593            0.954              7.94          64     262144            0.000437             0.0287
#  9    10 fdr                  2     6   0.124     0.622     0.741     0.746              1            66.2              111         1111111            1.48              26.3          100    1000000            0.000223             0.0226
# 10    10 fwer                 2     4   0.124     0.622     0.741     0.746              1            32.9              111           11111            1.48              16.2          100      10000            0.00224              0.0226
# 11    12 fdr                  2     4   0.166     0.680     0.74      0.748              1            52.3              157           22621            2.05              29.7          144      20736            0.00156              0.0188
# 12    12 fwer                 2     6   0.166     0.680     0.737     0.748              1           140.               157         3257437            2.05              69.8          144    2985984            0.000129             0.0188
# 13    14 fdr                  2     4   0.190     0.706     0.737     0.747              1            79.4              211           41371            2.61              49.8          196      38416            0.00114              0.0162
# 14    14 fwer                 2     4   0.190     0.706     0.737     0.747              1            79.4              211           41371            2.61              49.8          196      38416            0.00114              0.0162
# 15    16 fdr                  2     4   0.226     0.725     0.737     0.749              1           115.               273           69905            3.25              77.4          256      65536            0.000874             0.0142
# 16    16 fwer                 2     4   0.226     0.725     0.737     0.749              1           115.               273           69905            3.25              77.4          256      65536            0.000874             0.0142
# 17    18 fdr                  2     4   0.260     0.724     0.733     0.743              1           159.               343          111151            3.81             113.           324     104976            0.000689             0.0126
# 18    18 fwer                 2     4   0.260     0.724     0.733     0.743              1           159.               343          111151            3.81             113.           324     104976            0.000689             0.0126
# 19    20 fdr                  2     4   0.284     0.729     0.732     0.747              1           211.               421          168421            4.19             156.           400     160000            0.000558             0.0113
# 20    20 fwer                 2     4   0.284     0.729     0.732     0.747              1           211.               421          168421            4.19             156.           400     160000            0.000558             0.0113
# 21    50 fdr                  2     2   0.669     0.738     0.74      0.74               1            89.5             2551            2551           26.3               80.6         2500       2500            0.00442              0.00450
# 22    50 fwer                 2     2   0.669     0.738     0.74      0.74               1            89.5             2551            2551           26.3               80.6         2500       2500            0.00442              0.00450
# 23   100 fdr                  2     2   0.737     0.739     0.739     0.739              1           339.             10101           10101          105.               322.         10000      10000            0.00223              0.00225
# 24   100 fwer                 2     2   0.737     0.739     0.739     0.739              1           339.             10101           10101          105.               322.         10000      10000            0.00223              0.00225
# # A tibble: 24 × 18
#        k final_adj_method min_l max_l min_fwer max_fwer min_power max_power min_leaf_power max_nodes_tested min_total_nodes max_total_nodes min_leaves_tested max_leaves_tested min_leaves max_leaves min_bottom_up_power max_bottom_up_power
#    <int> <chr>            <int> <int>    <dbl>    <dbl>     <dbl>     <dbl>          <dbl>            <dbl>           <int>           <int>             <dbl>             <dbl>      <int>      <int>               <dbl>               <dbl>
#  1     2 fdr                  2    18   0         0.069    NA        NA                  1             4.20               7          524287            0                  0.810          4     262144           NA                    NA
#  2     2 fwer                 2    18   0         0.069    NA        NA                  1             4.20               7          524287            0                  0.810          4     262144           NA                    NA
#  3     4 fdr                  2     8   0.0016    0.162     0.730     0.742              1             6.59              21           87381            0.0086             1.70          16      65536            0.000873              0.106
#  4     4 fwer                 2     8   0.0016    0.162     0.730     0.742              1             6.59              21           87381            0.0086             1.70          16      65536            0.000873              0.106
#  5     6 fdr                  2     6   0.0364    0.289     0.736     0.747              1            12.3               43           55987            0.346              2.70          36      46656            0.00102               0.0385
#  6     6 fwer                 2     6   0.0364    0.289     0.736     0.747              1            12.3               43           55987            0.346              2.70          36      46656            0.00102               0.0385
#  7     8 fdr                  2     6   0.0974    0.481     0.732     0.741              1            29.0               73          299593            0.954              7.94          64     262144            0.000437              0.0287
#  8     8 fwer                 2     6   0.0974    0.481     0.732     0.741              1            29.0               73          299593            0.954              7.94          64     262144            0.000437              0.0287
#  9    10 fdr                  2     6   0.124     0.622     0.741     0.746              1            66.2              111         1111111            1.48              26.3          100    1000000            0.000223              0.0226
# 10    10 fwer                 2     4   0.124     0.622     0.741     0.746              1            32.9              111           11111            1.48              16.2          100      10000            0.00224               0.0226
# # ℹ 14 more rows
# # ℹ Use `print(n = ...)` to see more rows


hommel_raw <- simp_simsres_latest %>%
  filter(final_adj_method == "none" & alpha_fn == "fixed" & (prop_tau_nonzero > 0 & prop_tau_nonzero < 1) & local_adj_fn == "local_hommel_all_ps" & adj_effN == TRUE) %>%
  select(-file & one_of(c(
    "k", "l", "adj_effN", "alpha_fn", "local_adj_fn",
    "final_adj_method", "total_nodes", "num_leaves", "prop_tau_nonzero", "bottom_up_power", key_char
  ))) %>%
  arrange(false_error)

hommel_control_tab0 <- hommel_raw %>%
  group_by(k) %>%
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
  arrange(k)

print(hommel_control_tab0)
# # A tibble: 12 × 18
#        k min_l max_l min_fwer max_fwer min_power max_power min_leaf_power max_leaf_power max_nodes_tested min_total_nodes max_total_nodes min_leaves_tested max_leaves_tested min_leaves max_leaves min_bottom_up_power max_bottom_up_power
#    <int> <int> <int>    <dbl>    <dbl>     <dbl>     <dbl>          <dbl>          <dbl>            <dbl>           <int>           <int>             <dbl>             <dbl>      <int>      <int>               <dbl>               <dbl>
#  1     2     2    18   0        0.0538    NA        NA                  1              1             3.32               7          524287            0                  0.686          4     262144           NA                   NA
#  2     4     2     8   0        0.0342     0.730     0.749              1              1             3.08              21           87381            0                  0.834         16      65536            0.000877             0.106
#  3     6     2     6   0        0.0424     0.739     0.751              1              1             3.02              43           55987            0                  0.696         36      46656            0.00102              0.0385
#  4     8     2     6   0        0.0391     0.732     0.741              1              1             3.03              73          299593            0                  0.778         64     262144            0.000437             0.0287
#  5    10     2     6   0        0.0459     0.741     0.747              1              1             2.99             111         1111111            0                  0.800        100    1000000            0.000224             0.0226
#  6    12     2     4   0.0002   0.0418     0.748     0.749              1              1             2.91             157           22621            0.0017             0.804        144      20736            0.00155              0.0188
#  7    14     2     4   0.0001   0.0397     0.737     0.739              1              1             2.83             211           41371            0.0016             0.810        196      38416            0.00113              0.0162
#  8    16     2     4   0        0.0377     0.737     0.742              1              1             2.71             273           69905            0.0007             0.789        256      65536            0.000876             0.0142
#  9    18     2     4   0        0.0325     0.731     0.733              1              1             2.61             343          111151            0.0005             0.769        324     104976            0.000693             0.0126
# 10    20     2     4   0.0001   0.0295     0.732     0.748              1              1             2.44             421          168421            0.0008             0.694        400     160000            0.000556             0.0113
# 11    50     2     2   0.0046   0.0301     0.74      0.74               1              1             2.88            2551            2551            0.150              1.12        2500       2500            0.00442              0.00450
# 12   100     2     2   0.0069   0.0451     0.739     0.739              1              1             4.17           10101           10101            0.282              2.13       10000      10000            0.00223              0.00225

alpha_adaptive_k_raw <- simp_simsres_latest %>%
  filter(final_adj_method == "none" & alpha_fn == "adaptive_k_adj" & (prop_tau_nonzero > 0 & prop_tau_nonzero < 1) & local_adj_fn == "local_unadj_all_ps" & adj_effN == TRUE) %>%
  select(-file & one_of(c(
    "k", "l", "adj_effN", "alpha_fn", "local_adj_fn",
    "final_adj_method", "total_nodes", "num_leaves", "prop_tau_nonzero", "bottom_up_power", key_char
  ))) %>%
  arrange(false_error)

adaptive_k_control_tab0 <- alpha_adaptive_k_raw %>%
  group_by(k) %>%
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
  select(-max_leaf_power, -min_leaves_tested, -max_bottom_up_power) %>%
  arrange(k)

print(adaptive_k_control_tab0, n = 100)
# # A tibble: 12 × 15
#        k min_l max_l min_fwer max_fwer min_power max_power min_leaf_power max_nodes_tested min_total_nodes max_total_nodes max_leaves_tested min_leaves max_leaves min_bottom_up_power
#    <int> <int> <int>    <dbl>    <dbl>     <dbl>     <dbl>          <dbl>            <dbl>           <int>           <int>             <dbl>      <int>      <int>               <dbl>
#  1     2     2    18   0        0.0502    NA        NA              0.776             3.77               7          524287             0.678          4     262144           NA
#  2     4     2     8   0.0008   0.0422     0.727     0.750          0.945             3.78              21           87381             0.866         16      65536            0.000876
#  3     6     2     6   0.0026   0.0628     0.739     0.750          0.993             4.03              43           55987             0.867         36      46656            0.00103
#  4     8     2     6   0.0029   0.0734     0.732     0.738          0.991             4.61              73          299593             1.08          64     262144            0.000436
#  5    10     2     6   0.0069   0.0891     0.735     0.741          0.997             5.31             111         1111111             1.23         100    1000000            0.000223
#  6    12     2     4   0.0115   0.0907     0.748     0.750          0.998             5.19             157           22621             1.35         144      20736            0.00155
#  7    14     2     4   0.0141   0.09       0.737     0.743          0.999             5.49             211           41371             1.38         196      38416            0.00113
#  8    16     2     4   0.0159   0.098      0.737     0.744          1.00              5.46             273           69905             1.45         256      65536            0.000869
#  9    18     2     4   0.0156   0.0975     0.733     0.742          0.999             5.64             343          111151             1.50         324     104976            0.000692
# 10    20     2     4   0.0163   0.104      0.732     0.739          0.999             5.83             421          168421             1.70         400     160000            0.000560
# 11    50     2     2   0.0322   0.154      0.74      0.74           1.00              5.43            2551            2551             3.27        2500       2500            0.00442
# 12   100     2     2   0.0627   0.222      0.739     0.739          1                10.6            10101           10101             7.97       10000      10000            0.00223

table(simp_simsres_latest$alpha_fn)
#
# adaptive_k_adj          fixed    fixed_k_adj      investing       spending
#           1996           2022           2003           1982           1993


fixed_k_raw <- simp_simsres_latest %>%
  filter(final_adj_method == "none" & alpha_fn == "fixed_k_adj" & (prop_tau_nonzero > 0 & prop_tau_nonzero < 1) & local_adj_fn == "local_unadj_all_ps" & adj_effN == TRUE) %>%
  select(-file & one_of(c(
    "k", "l", "adj_effN", "alpha_fn", "local_adj_fn",
    "final_adj_method", "total_nodes", "num_leaves", "prop_tau_nonzero", "bottom_up_power", key_char
  ))) %>%
  arrange(false_error)

fixed_k_control_tab0 <- fixed_k_raw %>%
  group_by(k) %>%
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
  select(-max_leaf_power, -min_leaves_tested, -max_bottom_up_power) %>%
  arrange(k)

print(fixed_k_control_tab0, n = 100)

investing_raw <- simp_simsres_latest %>%
  filter(final_adj_method == "none" & alpha_fn == "investing" & (prop_tau_nonzero > 0 & prop_tau_nonzero < 1) & local_adj_fn == "local_unadj_all_ps" & adj_effN == TRUE) %>%
  select(-file & one_of(c(
    "k", "l", "adj_effN", "alpha_fn", "local_adj_fn",
    "final_adj_method", "total_nodes", "num_leaves", "prop_tau_nonzero", "bottom_up_power", key_char
  ))) %>%
  arrange(false_error)

investing_control_tab0 <- investing_raw %>%
  group_by(k) %>%
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
  select(-max_leaf_power, -min_leaves_tested, -max_bottom_up_power) %>%
  arrange(k)

print(investing_control_tab0, n = 100)



spending_raw <- simp_simsres_latest %>%
  filter(final_adj_method == "none" & alpha_fn == "spending" & (prop_tau_nonzero > 0 & prop_tau_nonzero < 1) & local_adj_fn == "local_unadj_all_ps" & adj_effN == TRUE) %>%
  select(-file & one_of(c(
    "k", "l", "adj_effN", "alpha_fn", "local_adj_fn",
    "final_adj_method", "total_nodes", "num_leaves", "prop_tau_nonzero", "bottom_up_power", key_char
  ))) %>%
  arrange(false_error)

spending_control_tab0 <- spending_raw %>%
  group_by(k) %>%
  summarize(
    min_l = min(l), max_l = max(l),
    min_fwer = min(false_error), max_fwer = max(false_error),
    min_power = min(power), max_power = max(power),
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
  select(-max_leaf_power, -min_leaves_tested, -max_bottom_up_power) %>%
  arrange(k)

print(spending_control_tab0, n = 100)

## Looking at combinations with hommel

hommel_plus_raw <- simp_simsres_latest %>%
  filter(final_adj_method == "none" & (prop_tau_nonzero > 0 & prop_tau_nonzero < 1) & local_adj_fn == "local_hommel_all_ps" & adj_effN == TRUE) %>%
  select(-file & one_of(c(
    "k", "l", "adj_effN", "alpha_fn", "local_adj_fn",
    "final_adj_method", "total_nodes", "num_leaves", "prop_tau_nonzero", "bottom_up_power", key_char
  ))) %>%
  arrange(false_error)

hommel_plus_control_tab0 <- hommel_plus_raw %>%
  group_by(k, alpha_fn) %>%
  summarize(
    min_l = min(l), max_l = max(l),
    min_fwer = min(false_error), max_fwer = max(false_error),
    min_power = min(power), max_power = max(power),
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
  select(-max_leaf_power, -min_leaves_tested, -max_bottom_up_power, -max_total_nodes, -min_total_nodes, -max_leaves, -min_bottom_up_power) %>%
  arrange(k, alpha_fn)

print(hommel_plus_control_tab0, n = 100)

## So, only use alpha adjustments along with hommell. Something weird about k=2, l=18 (I suspect) with the investing adjustment

## What about with global?

global_plus_raw <- simp_simsres_latest %>%
  filter(final_adj_method != "none" & (prop_tau_nonzero > 0 & prop_tau_nonzero < 1) &
    local_adj_fn == "local_hommel_all_ps" & adj_effN == TRUE & alpha_fn == "fixed") %>%
  select(-file & one_of(c(
    "k", "l", "adj_effN", "alpha_fn", "local_adj_fn",
    "final_adj_method", "total_nodes", "num_leaves", "prop_tau_nonzero", "bottom_up_power", key_char
  ))) %>%
  arrange(false_error)

global_plus_control_tab0 <- global_plus_raw %>%
  group_by(k, final_adj_method, alpha_fn) %>%
  summarize(
    min_l = min(l), max_l = max(l),
    min_fwer = min(false_error), max_fwer = max(false_error),
    min_power = min(power), max_power = max(power),
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
  select(-max_leaf_power, -min_leaves_tested, -max_bottom_up_power, -max_total_nodes, -min_total_nodes, -max_leaves, -min_bottom_up_power) %>%
  arrange(k, final_adj_method, alpha_fn)

## Looks great but not clear better power than hommel alone.
print(global_plus_control_tab0, n = 100)

## TODO
## Choose combinations that can be simulated more quickly to assess different power.
## Is there a combination that has best power while controlling FWER? (say, hommel+investing+global?)
