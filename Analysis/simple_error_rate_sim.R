## Basic simulations of error rates but using directly generated p-values

library(here)
library(abind)
library(tidyverse)

#' Simulate P-values on a tree
#'
#' Assuming no treatement effects at all. Generate sims (say, 1000) p-values at
#' the root level. For any p<=.05 generate p-values for child nodes --- where
#' each node has a true null of no effects and where the minimum p-value is the
#' minimum p-value <= alpha and the maximum is 1. Do this for each level. Here,
#' assuming the same number of nodes per level throughout the tree below the
#' root level.
#' value

sim_ps <- function(num_levels, nodes_per_level, alpha, sims) {
  ## For now, we assume that there are no false nulls. This is a completely no effects world.
  ## should make .05 of values less than .05
  get_child_ps <- function(parent_p, depth) {
    node_ps_lst <- lapply(1:nodes_per_level, function(node_num) {
      thep <- runif(sims, min = parent_p, max = 1)
      thep_le_alpha <- max(0, thep[thep <= alpha])
      res <- vector(length = 3L)
      names(res) <- c(
        paste("node_fpr_", depth, sep = ""),
        paste("p_", depth, sep = ""),
        paste("p_", depth - 1, sep = "")
      )
      res[1:3] <- c(
        mean(thep <= alpha),
        min(thep),
        parent_p
      )
      return(res)
    })
    return(node_ps_lst)
  }

  ## Start with the root node.
  p_0 <- runif(sims, min = 0, max = 1)
  p_0_le_alpha <- p_0[p_0 <= alpha]
  node_fpr_0 <- mean(p_0 <= alpha)
  sim_id_level_0 <- as.character(1:length(p_0_le_alpha))

  if (node_fpr_0 == 0) {
    message("No Root FPR", "nodes per level", nodes_per_level, "num levels", num_levels)
    return(tibble(node_fpr_0 = 0, p_0 = min(p_0), p_p = min(p_0), sim_id_level_0 = 0))
  }

  ## for now, require at least a root node and children nodes
  stopifnot(num_levels > 1)

  ## Then do the first level after the root by hand
  p_1s <- lapply(p_0_le_alpha, function(parent_p) {
    get_child_ps(parent_p = parent_p, depth = 1)
  })
  stopifnot(all.equal(length(p_1s), length(p_0_le_alpha)))
  p_1_dt <- dplyr::bind_rows(p_1s)
  p_1_dt$sim_id_level_0 <- rep(sim_id_level_0, each = nodes_per_level)
  p_1_dt$sim_id_level_1 <- paste(p_1_dt$sim_id_level_0, rep(1:nodes_per_level), sep = ".")
  stopifnot(all.equal(nrow(p_1_dt), nodes_per_level * length(p_0_le_alpha)))

  ## Setup a list to receive the data from the other levels
  res_list <- vector("list", length = num_levels)
  res_list[[1]] <- p_1_dt

  # Notice that the minimum p-value for the children at level 1 are highly
  # correlated with the parent p-value<=0  because the minimum p-value for the
  # children is the parent p-value.

  # with(p_1_dt,plot(p_0,p_1,xlim=c(0,.07),ylim=c(0,.07))); abline(0,1)

  library(parallel)
  library(future)
  library(future.apply)
  ncores <- future::availableCores() - 1
  for (i in 2:num_levels) {
    res <- mclapply(res_list[[i - 1]][[paste0("p_", i - 1)]], function(parent_p) {
      get_child_ps(parent_p = parent_p, depth = i)
    }, mc.cores = ncores)
    res_dt <- dplyr::bind_rows(res)
    res_dt[[paste0("sim_id_level_", i - 1)]] <- rep(res_list[[i - 1]][[paste0("sim_id_level_", i - 1)]], each = nodes_per_level)
    res_dt[[paste0("sim_id_level_", i)]] <- paste(res_dt[[paste0("sim_id_level_", i - 1)]], rep(1:nodes_per_level), sep = ".")
    res_list[[i]] <- res_dt
    rm(res, res_dt)
  }

  # p_dt <- left_join(p_2_dt, p_1_dt, by = join_by(p1 == min_p_le_alpha))
  my_join <- function(xdat, ydat) {
    xdat_id <- grep("^sim_id_level_[0-9]", names(xdat), value = TRUE)
    ydat_id <- grep("^sim_id_level_[0-9]", names(ydat), value = TRUE)
    stopifnot(nrow(xdat) >= nrow(ydat))
    left_join(xdat, ydat, by = intersect(xdat_id, ydat_id), suffix = c("", ".y"))
  }

  p_dt <- purrr::reduce(res_list[length(res_list):1], my_join)
  ## Just a few tests
  ## stopifnot(all(p_dt$p_3 >= p_dt$p_2))
  ## stopifnot(all(p_dt$p_2 >= p_dt$p_1))
  stopifnot(all(p_dt$p_1 >= p_dt$p_0))
  stopifnot(all.equal(p_dt$p_1.y, p_dt$p_1))
  # stopifnot(all.equal(p_dt$p_9.y,p_dt$p_9))
  ## Drop all columns that are duplicated
  p_dt <- p_dt %>% select(-ends_with(".y"))

  p_dt$node_fpr_0 <- node_fpr_0

  return(p_dt)
  ## End
}

get_max_fpr <- function(parent_lv, child_lv, obj = test_l6_n2) {
  res <- obj %>%
    group_by(.data[[paste0("sim_id_level_", parent_lv)]]) %>%
    summarize(mean_fpr = mean(.data[[paste0("node_fpr_", child_lv)]])) %>%
    ungroup()
  ## This is the maximum false positive rate at a given level in the tree
  return(max(res$mean_fpr))
}

#### Now use the simulation function

set.seed(12345)
num_levs <- 20
test_l10_n2 <- sim_ps(num_levels = num_levs, nodes_per_level = 2, alpha = .05, sims = 1000)
save(test_l10_n2, file = "test_l10_n2.rda")
str(test_l10_n2)
head(test_l10_n2)
test_l10_n2 %>% select(starts_with("p"), sim_id_level_0, sim_id_level_1)
## Only the root level tests <= alpha in p_0
summary(test_l10_n2$p_0)
## We only keep the p_0 if <= alpha.
length(unique(test_l10_n2$p_0))
unique(test_l10_n2$node_fpr_0)
stopifnot(all.equal(length(unique(test_l10_n2$p_0)) / 10000, unique(test_l10_n2$node_fpr_0)))
lvs <- seq(0, num_levs - 1)
res_l10_n2_tmp <- sapply(lvs, function(parent_lv) {
  get_max_fpr(parent_lv = parent_lv, child_lv = parent_lv + 1, obj = test_l10_n2)
})
res_l10_n2 <- c(unique(test_l10_n2$node_fpr_0), res_l10_n2_tmp)
names(res_l10_n2) <- seq(0, num_levs)
res_l10_n2

## Now with many more nodes per level
set.seed(12345)
num_levs <- 2
test_l3_n20 <- sim_ps(num_levels = num_levs, nodes_per_level = 20, alpha = .05, sims = 10000)
save(test_l3_n20, file = "test_l3_n20.rda")
str(test_l3_n20)
head(test_l3_n20)
lvs_l3 <- seq(0, num_levs - 1)
res_l3_n20_tmp <- sapply(lvs_l3, function(parent_lv) {
  get_max_fpr(parent_lv = parent_lv, child_lv = parent_lv + 1, obj = test_l3_n20)
})
res_l3_n20 <- c(unique(test_l3_n20$node_fpr_0), res_l3_n20_tmp)
names(res_l3_n20) <- seq(0, num_levs)
res_l3_n20


set.seed(12345)
test_l3_n100 <- sim_ps(num_levels = 2, nodes_per_level = 100, alpha = .05, sims = 10000)
save(test_l3_n100, file = "test_l3_n100.rda")
str(test_l3_n100)
head(test_l3_n100)
lvs_l3 <- 0:1
res_l3_n100_tmp <- sapply(lvs_l3, function(parent_lv) {
  get_max_fpr(parent_lv = parent_lv, child_lv = parent_lv + 1, obj = test_l3_n100)
})
res_l3_n100 <- c(unique(test_l3_n100$node_fpr_0), res_l3_n100_tmp)
names(res_l3_n100) <- 0:2
res_l3_n100

## Now try to see if there are patterns when we increase number of nodes per level
node_sim_parms <- sort(unique(c(seq(2, 10, 1), seq(10, 100, 10))))

set.seed(12345)
res_l3 <- lapply(node_sim_parms, function(n_nodes) {
  message(n_nodes)
  sim_ps(num_levels = 2, nodes_per_level = n_nodes, alpha = .05, sims = 10000)
})
save(res_l3, file = "res_l3.rda")

tests_l3 <- sapply(res_l3, function(obj) {
  res <- sapply(0:1, function(parent_lv) {
    tmp <- get_max_fpr(parent_lv = parent_lv, child_lv = parent_lv + 1, obj = obj)
    stopifnot(length(unique(obj$node_fpr_0)) == 1)
    return(tmp)
  })
  c(unique(obj$node_fpr_0), res)
})

colnames(tests_l3) <- node_sim_parms
rownames(tests_l3) <- 0:2

tests_l3

res_l5 <- lapply(node_sim_parms, function(n_nodes) {
  message(n_nodes)
  sim_ps(num_levels = 5, nodes_per_level = n_nodes, alpha = .05, sims = 10000)
})
save(res_l5, file = "res_l5.rda")

## Next, try strong control.
###

m <- 20 # Total number of hypotheses
pi0 <- 0.9 # Proportion of true null hypotheses
alpha <- 0.05 # Significance level
n_sim <- 10000 # Number of simulations
# Number of true and false nulls
m0 <- round(m * pi0)
m1 <- m - m0
p_values_null <- runif(m0, 0, 1)
p_values_alternative <- rbeta(m1, 0.1, 1)
p_values <- c(p_values_null, p_values_alternative)
