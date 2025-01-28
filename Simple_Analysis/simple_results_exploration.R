library(here)
library(tidyverse)
library(viridis)

load(here("Analysis", "simple_sims_results.rda"), verbose = TRUE)

## we use alpha=.05 below
## and simulation error is
sim_se <- 2 * sqrt(.05 * (1 - .05) / 1000)

simp_simsres <- simp_simsres %>%
  filter(!is.na(prop_tau_nonzero)) %>%
  droplevels()

g_fill <- ggplot(data = simp_simsres, aes(x = k, y = l, z = fpr_dfs)) +
  # geom_contour_filled(breaks = c(0, .05, .1, .2, .3, .4, .5, 1)) +
  geom_contour_filled(breaks = c(0, .05 + .014, .1 + .014, .2, .3, .4, .5, 1)) +
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

print(g_fill)


g_lines_nodes <- ggplot(
  data = simp_simsres,
  aes(groups = l, color = factor(l), x = k, y = fpr_dfs)
) +
  facet_wrap(~prop_tau_nonzero) +
  geom_line() +
  geom_hline(yintercept = .05)

g_lines_nodes


g_lines_levels <- ggplot(
  data = simp_simsres,
  aes(groups = k, color = factor(k), x = l, y = fpr_dfs)
) +
  geom_line() +
  facet_wrap(~prop_tau_nonzero) +
  geom_hline(yintercept = .05)

g_lines_levels


## with no possible false positives (prop_tau_nonzero==1) or all true null
## (prop_tau_nonzero=0) we have control at all levels
simp_simsres %>%
  filter(prop_tau_nonzero %in% c(0, 1)) %>%
  select(-file) %>%
  reframe(summary(fpr_dfs))


simp_simsres %>%
  filter(prop_tau_nonzero > 0 & prop_tau_nonzero < 1) %>%
  select(-file) %>%
  arrange(fpr_dfs)
