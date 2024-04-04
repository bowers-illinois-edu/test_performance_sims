## Show power differences between the approaches

library(tidyverse)
library(here)
library(data.table)
library(dtplyr)
library(dplyr)
library(ggthemes)
library(RColorBrewer)
library(manytestsr)
library(xtable)
library(ggh4x) ## for facet_grid2

load(here::here("figures", "simsres_plotting.rda"))

names(simsres)
head(simsres)

## Standard error of the simulations
nsims <- 1000
std_err_of_sim <- 2 * sqrt(.05 * (1 - .05) / nsims)

table(simsres$prop_blocks_0, exclude = c())
table(simsres$prop_blocks_0_labels, exclude = c())
table(simsres$tau_labels, exclude = c())
table(simsres$splitfn)

## First, compare top-down and bottom-up when alpha is fixed (so that the
## bottom up approaches are targeting either FDR or FWER (holm, hommel).

## We do not show the prespecified splitters here because those are not chosen
## for power: in those cases a hypothetical statistical consultant to a
## decision maker already knows the nesting and wants to test as much as
## possible (but also to protect the FWER). So increasing power would have to
## be done via some other means such as covariance adjustment.

## Just as an example, choose middle effect sizes to look at the power of the pre-specified appraoches --- which is pretty low here
simsres %>%
  filter(prop_blocks_0 == .5 & afn == "fixed (a=.05)" & tau_sizes == 1 & nblocks %in% c(10, 100) & pfn == "pIndepDist") %>%
  select(splitfnF, nblocks, true_pos_prop, false_pos_prop) %>%
  arrange(nblocks, splitfnF)


plotdat_alpha_fixed <- simsres %>%
  filter(prop_blocks_0 < 1 &
    nblocks < 1000 &
    pfn == "pIndepDist" &
    afn == "fixed (a=.05)" &
    tau_sizes > 0 &
    prop_blocks_0 %in% c(0, 0.1, 0.5, 0.9)) %>%
  droplevels() %>%
  as.data.frame()

g_pwr_id_alpha_fixed <- ggplot(
  data = plotdat_alpha_fixed,
  aes(x = nblocks, y = true_pos_prop, color = splitfnF, shape = split_vs_not)
) +
  geom_line() +
  geom_point() +
  facet_grid2(
    rows = vars(prop_blocks_0_labels), cols = vars(tau_labels),
    labeller = labeller(tau_labels = label_parsed, prop_blocks_0_labels = label_parsed),
    scales = "free_y",
    axes = "y",
    remove_labels = "x",
    independent = "y"
  ) +
  labs(
    color = "",
    shape = "",
    y = "True Positive Proportion (Power)",
    x = "Number of Hypotheses (Blocks)"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9)) +
  theme(legend.position = "bottom", legend.background = element_blank()) +
  scale_color_colorblind() +
  scale_shape_manual(values = c(3, 16, 17))
print(g_pwr_id_alpha_fixed)

ggsave(g_pwr_id_alpha_fixed, file = here::here("figures", "power_alpha_fixed.pdf"), width = 15, height = 5)


## Some facts

max_power_alpha_fixed <- plotdat_alpha_fixed %>%
  group_by(nblocks, prop_blocks_0) %>%
  summarize(
    top_down = max(true_pos_prop[split_vs_not == "Top-Down"]),
    holm_hommell = max(true_pos_prop[p_adj_method %in% c("holm", "hommel")]),
    fdr = max(true_pos_prop[p_adj_method == "fdr"]),
    top_down_vs_fwer = top_down / holm_hommell,
    top_down_vs_fdr = top_down / fdr
  )
max_power_alpha_fixed
#  > max_power_alpha_fixed
#  # A tibble: 20 × 7
#  # Groups:   nblocks [5]
#     nblocks prop_blocks_0 top_down holm_hommell    fdr top_down_vs_fwer top_down_vs_fdr
#       <int>         <dbl>    <dbl>        <dbl>  <dbl>            <dbl>           <dbl>
#   1      10           0     0.834        0.724  0.809              1.15            1.03
#   2      10           0.1   0.745        0.635  0.719              1.17            1.04
#   3      10           0.5   0.366        0.314  0.341              1.17            1.07
#   4      10           0.9   0.1          0.1    0.1                1               1
#   5      25           0     0.805        0.649  0.788              1.24            1.02
#   6      25           0.1   0.732        0.597  0.712              1.23            1.03
#   7      25           0.5   0.418        0.333  0.400              1.26            1.05
#   8      25           0.9   0.0788       0.0195 0.0271             4.04            2.91
#   9      50           0     0.876        0.728  0.868              1.20            1.01
#  10      50           0.1   0.775        0.595  0.762              1.30            1.02
#  11      50           0.5   0.415        0.320  0.397              1.30            1.05
#  12      50           0.9   0.0830       0.0727 0.0789             1.14            1.05
#  13      75           0     0.829        0.644  0.818              1.29            1.01
#  14      75           0.1   0.751        0.564  0.735              1.33            1.02
#  15      75           0.5   0.456        0.370  0.442              1.23            1.03
#  16      75           0.9   0.0839       0.0805 0.0818             1.04            1.03
#  17     100           0     0.801        0.601  0.784              1.33            1.02
#  18     100           0.1   0.723        0.544  0.704              1.33            1.03
#  19     100           0.5   0.385        0.286  0.353              1.35            1.09
#  20     100           0.9   0.0819       0.0767 0.0803             1.07            1.02

#################
## Does allowing alpha to vary increase power?
table(simsres$afn)
table(simsres$split_vs_not)
table(simsres$afn, simsres$splitfn)
table(simsres$afn, simsres$p_adj_method, exclude = c())

simsres <- simsres %>% mutate(split_adj_fn = ifelse(p_adj_method == "split", splitfn, p_adj_method))
with(simsres, table(split_adj_fn, p_adj_method, exclude = c()))

alpha_power_compared <- simsres %>%
  filter(prop_blocks_0 < 1 &
    nblocks < 1000 &
    pfn == "pIndepDist" &
    tau_sizes > 0) %>%
  droplevels() %>%
  group_by(split_adj_fn, afn, tau_sizes, prop_blocks_0, nblocks) %>%
  summarize(
    # fixed_vs_addis = true_pos_prop[afn == "fixed (a=.05)" & p_adj_method=="split"] - true_pos_prop[afn == "alpha_addis" & p_adj_method=="split"],
    # fixed_vs_saffron = true_pos_prop[afn == "fixed (a=.05)" & p_adj_method=="split"] - true_pos_prop[afn == "alpha_saffron"& p_adj_method=="split" ],
    # fixed_vs_investing = true_pos_prop[afn == "fixed (a=.05)" & p_adj_method=="split" ] - true_pos_prop[afn == "alpha_investing" & p_adj_method=="split"],
    max_true_prop = max(true_pos_prop),
    max_false_pos_prop = max(false_pos_prop),
    max_false_disc_prop = max(false_disc_prop)
  ) %>%
  arrange(split_adj_fn, nblocks, tau_sizes, prop_blocks_0) %>%
  ungroup() %>%
  as.data.frame()

## alpha_power_compared %>%
##  filter(tau_sizes == .5 & prop_blocks_0 == .5) %>%
##  group_by(split_adj_fn, nblocks) %>%
##  reframe(vs_saffron = median(fixed_vs_saffron),
##    vs_investing = median(fixed_vs_investing),
##    vs_addis = median(fixed_vs_addis),
##    maxpwd = max(max_true_prop),
##    maxfwer = max(max_false_pos_prop),
##    maxfdr = max(max_false_disc_prop)) %>%
##  print(n = 100)

library(ggrepel)
library(lemon) ## for reposition_legend

plotdat_alpha_varying <- simsres %>%
  filter(prop_blocks_0 < 1 &
    nblocks < 1000 &
    pfn == "pIndepDist" &
    tau_sizes == .5 &
    prop_blocks_0 == .5 & split_vs_not == "Top-Down") %>%
  droplevels() %>%
  group_by(nblocks, splitfnF, afn) %>%
  mutate(
    max_fwer = max(false_pos_prop),
    max_fdr = max(false_disc_prop)
  ) %>%
  ungroup() %>%
  as.data.frame()

loo_label_dat <- plotdat_alpha_varying %>%
  filter(splitfn == "splitLOO") %>%
  mutate(max_fwer = ifelse(max_fwer > .05, max_fwer, NA)) %>%
  droplevels()


bottom_up_comparison <- simsres %>%
  filter(prop_blocks_0 < 1 & nblocks < 1000 & pfn == "pIndepDist" & tau_sizes == .5 & prop_blocks_0 == .5 &
    splitfn != "splitSpecifiedFactorMulti" & split_vs_not == "Bottom-Up") %>%
  droplevels() %>%
  group_by(nblocks, splitfnF) %>%
  select(nblocks, splitfnF, p_adj_method, true_pos_prop, false_pos_prop, false_disc_prop) %>%
  ungroup() %>%
  arrange(nblocks, splitfnF) %>%
  as.data.frame()
## bottom_up_comparison$splitfnF <- rep("Equal Sample Size Split", nrow(bottom_up_comparison))

bottom_up_comparison

bottom_up_comparison %>%
  group_by(splitfnF) %>%
  summarize(
    max_fwer = max(false_pos_prop),
    max_fdr = max(false_disc_prop), max_power = max(true_pos_prop)
  )

g_pwr_id_alpha_varying <- ggplot(
  data = plotdat_alpha_varying,
  aes(x = nblocks, y = true_pos_prop, shape = afn, color = afn, label = round(max_fwer, 3))
) +
  geom_line() +
  geom_point() +
  # geom_line(data = bottom_up_comparison, aes(x = nblocks, y = true_pos_prop, group = p_adj_method, color = p_adj_method)) +
  facet_wrap(~splitfnF,
    scales = "free"
  ) +
  geom_text_repel(data = loo_label_dat, color = "black") +
  labs(
    shape = "",
    color = "",
    y = "True Positive Proportion (Power)",
    x = "Number of Hypotheses (Blocks)"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9)) +
  theme(legend.position = "bottom", legend.background = element_blank()) +
  scale_shape_manual(values = c(1, 2, 3, 4, 16, 17, 18))
g_pwr_id_alpha_varying





# https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
shift_legend3 <- function(p) {
  pnls <- cowplot::plot_to_gtable(p) %>%
    gtable::gtable_filter("panel") %>%
    with(setNames(grobs, layout$name)) %>%
    purrr::keep(~ identical(.x, zeroGrob()))

  if (length(pnls) == 0) stop("No empty facets in the plot")

  lemon::reposition_legend(p, "center", panel = names(pnls))
}

# g2 <- shift_legend3(g_pwr_id_alpha_varying)

ggsave(g_pwr_id_alpha_varying, file = here("figures", "power_alpha_varying.pdf"), width = 8, height = 8)
