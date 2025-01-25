library(here)
library(tidyverse)
library(viridis)

load(here("Analysis", "simple_sims_results.rda"), verbose = TRUE)

## we use alpha=.05 below

g_fill <- ggplot(data = simp_simsres, aes(x = k, y = l, z = fpr_dfs)) +
  geom_contour_filled(breaks = c(0, .05, .1, .2, .3, .4, .5, 1)) +
  labs(
    x = "Number of nodes per level",
    fill = "FWER"
  ) +
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
  geom_line() +
  geom_hline(yintercept = .05)

g_lines_nodes


g_lines_levels <- ggplot(
  data = simp_simsres,
  aes(groups = k, color = factor(k), x = l, y = fpr_dfs)
) +
  geom_line() +
  geom_hline(yintercept = .05)

g_lines_levels

## Seems like as long as we have only 2 splits we are always ok.

simp_simsres %>% filter(k == 2)
#   V1 k  l total_nodes fpr_dfs time_dfs sims               file
# 1  0 2  2           3   0.026    0.738 1000      sim_2_2_3.csv
# 2 10 2  4          15   0.012    1.283 1000     sim_2_4_15.csv
# 3 13 2  6          63   0.020    3.791 1000     sim_2_6_63.csv
# 4 15 2  8         255   0.012   11.811 1000    sim_2_8_255.csv
# 5 20 2 10        1023   0.005   45.877 1000  sim_2_10_1023.csv
# 6 25 2 12        4095   0.003  192.909 1000  sim_2_12_4095.csv
# 7 28 2 14       16383   0.002  740.287 1000 sim_2_14_16383.csv
# 8 31 2 16       65535   0.002 2913.331 1000 sim_2_16_65535.csv

## But the same is not for levels
simp_simsres %>% filter(l == 2)
#    V1   k l total_nodes fpr_dfs time_dfs sims               file
# 1   0   2 2           3   0.026    0.738 1000      sim_2_2_3.csv
# 2   1   4 2          15   0.038    0.752 1000     sim_4_2_15.csv
# 3   2   8 2          63   0.072    0.811 1000     sim_8_2_63.csv
# 4   3   6 2          35   0.060    0.848 1000     sim_6_2_35.csv
# 5   4  10 2          99   0.098    0.850 1000    sim_10_2_99.csv
# 6   5  16 2         255   0.159    0.842 1000   sim_16_2_255.csv
# 7   6  12 2         143   0.125    0.866 1000   sim_12_2_143.csv
# 8   7  14 2         195   0.126    0.899 1000   sim_14_2_195.csv
# 9   8  18 2         323   0.165    0.958 1000   sim_18_2_323.csv
# 10  9  20 2         399   0.200    0.984 1000   sim_20_2_399.csv
# 11 11 100 2        9999   0.447    2.202 1000 sim_100_2_9999.csv

So, basicaly 
