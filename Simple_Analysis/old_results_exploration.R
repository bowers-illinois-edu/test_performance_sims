
## TODO: Question is whether hommel + spending or investing etc.. has more power.

##### OLD STUFF BELOW
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
