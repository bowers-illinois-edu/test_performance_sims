library(here)
library(abind)
library(tidyverse)

load(here("test_l10_n2_bigsims.rda"), verbose = TRUE)


get_max_fpr <- function(parent_lv, child_lv, obj = test_l6_n2) {
  res <- obj %>%
    group_by(.data[[paste0("sim_id_level_", parent_lv)]]) %>%
    summarize(mean_fpr = mean(.data[[paste0("node_fpr_", child_lv)]])) %>%
    ungroup()
  ## This is the maximum false positive rate at a given level in the tree
  return(max(res$mean_fpr))
}

options(scipen = 10)
set.seed(12345)
num_levs <- 10
str(test_l10_n2)
head(test_l10_n2)
test_l10_n2 %>% select(starts_with("p"), sim_id_level_0, sim_id_level_1)
## Only the root level tests <= alpha in p_0
summary(test_l10_n2$p_0)
## We only keep the p_0 if <= alpha.
length(unique(test_l10_n2$p_0))
unique(test_l10_n2$node_fpr_0)
stopifnot(all.equal(length(unique(test_l10_n2$p_0)) / 100000, unique(test_l10_n2$node_fpr_0)))
lvs <- seq(0, num_levs - 1)
res_l10_n2_tmp <- sapply(lvs, function(parent_lv) {
  get_max_fpr(parent_lv = parent_lv, child_lv = parent_lv + 1, obj = test_l10_n2)
})
res_l10_n2 <- c(unique(test_l10_n2$node_fpr_0), res_l10_n2_tmp)
names(res_l10_n2) <- seq(0, num_levs)
res_l10_n2
