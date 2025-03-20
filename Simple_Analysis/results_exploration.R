## Explore and interpret the results of the simulations. Here focusing on the
## version without adjustment at each node.

library(here)
library(dplyr)
library(dtplyr)
library(data.table)
library(conflicted)
library(ggplot2)
library(viridis)
conflicts_prefer(dplyr::filter)

# load(here("Simple_Analysis", "simple_sims_unadj_results.rda"), verbose = TRUE)
# load(here("Simple_Analysis", "simple_sims_adj_results.rda"), verbose = TRUE)
load(here("Simple_Analysis", "simple_sims_latest_results.rda"), verbose = TRUE)

## we use alpha=.05 below
## and simulation error is this (assuming 1000 sims)
sim_se <- 2 * sqrt(.05 * (1 - .05) / 1000)
## We want to be less that this number:
.05 + sim_se

## simp_simsres_unadj <- simp_simsres_unadj %>%
##   filter(!is.na(prop_tau_nonzero)) %>%
##   droplevels() %>% as.data.table()

## using the dtplyr approach for readablilty for those not used to data.table syntax

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

res_all_true <- simp_simsres_latest[prop_tau_nonzero == 0,
  {
    res_list <- lapply(key_char, function(nm) {
      res0 <- summary_fn(get(nm))
      res0[, variable := nm]
      return(res0)
    })
    rbindlist(res_list)
  },
  by = c("adj_effN", "local_adj_fn")
]

## We have weak control
res_all_true %>% filter(variable %in% c("false_error"))

## When all hyps are false then of course no false errors.

vars2 <- c("power", "leaf_power", "num_leaves_tested", "num_leaves", "bottom_up_power")
res_all_false <- simp_simsres_latest[prop_tau_nonzero == 1,
  {
    res_list <- lapply(vars2, function(nm) {
      res0 <- summary_fn(get(nm))
      res0[, variable := nm]
      return(res0)
    })
    rbindlist(res_list)
  },
  by = c("adj_effN", "local_adj_fn")
]

res_all_false

## So, the top down approach tests fewer leaves but tends to nearly always
## reject the false leaves compared to the bottom up approach which fails to
## reject a lot (comparing leaf_power to bottom_up_power here)

## Now what about strong control? Without both simulating the idea of splitting
## (adj_effN=TRUE) and local adjustment, we do not have strong control. And, in
## fact, we need direct p-value adjustment (local_hommel_all_ps) rather than
## the Simes local/global test.

simp_simsres_latest %>%
  filter(prop_tau_nonzero > 0 & prop_tau_nonzero < 1) %>%
  group_by(adj_effN, prop_tau_nonzero, local_adj_fn) %>%
  summarize(
    mean_fwer = mean(false_error),
    max_fwer = max(false_error),
    mean_bottom_up_fwer = mean(bottom_up_false_error),
    max_bottom_up_fwer = max(bottom_up_false_error)
  )

## So, our winner is adj_effN=TRUE and local_hommel_all_ps  and we should
## compare this against the bottom up hommel in terms of power.

simp_simsres_latest %>%
  filter(adj_effN == TRUE & local_adj_fn == "local_hommel_all_ps") %>%
  group_by(adj_effN, prop_tau_nonzero, local_adj_fn) %>%
  summarize(
    mean_fwer = mean(false_error),
    max_fwer = max(false_error),
    mean_bottom_up_fwer = mean(bottom_up_false_error),
    max_bottom_up_fwer = max(bottom_up_false_error)
  )

simp_simsres_latest %>%
  filter(prop_tau_nonzero > 0 & adj_effN == TRUE & local_adj_fn == "local_hommel_all_ps") %>%
  group_by(adj_effN, prop_tau_nonzero, local_adj_fn) %>%
  summarize(
    min_power = min(power, na.rm = TRUE),
    med_power = median(power, na.rm = TRUE),
    max_power = max(power, na.rm = TRUE),
    min_leaf_power = min(leaf_power, na.rm = TRUE),
    med_leaf_power = median(leaf_power, na.rm = TRUE),
    max_leaf_power = max(leaf_power, na.rm = TRUE),
    min_bot_up_power = min(bottom_up_power, na.rm = TRUE),
    med_bot_up_power = median(bottom_up_power, na.rm = TRUE),
    max_bot_up_power = max(bottom_up_power, na.rm = TRUE)
  )


## Make a data set of this:

strong_ctrl_res <- simp_simsres_latest %>%
  filter(adj_effN == TRUE & local_adj_fn == "local_hommel_all_ps" &
    (prop_tau_nonzero > 0 & prop_tau_nonzero < 1)) %>%
  as_tibble()

## Here we didn't get to test any leaves. The bottom up approach does 1024 tests with very low power
## But the top down approach stops early --- it tells us **correctly** that at least one of the leaves of a rejected parent is nonnull
## but it also says that we cannot detect the specific leaf
##
with(strong_ctrl_res, table(k, l))
## When it does test a nonnull leaf, it always rejects
summary(strong_ctrl_res$num_leaves_tested)
summary(strong_ctrl_res$leaf_power)
table(strong_ctrl_res$leaf_power, exclude = c())

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
