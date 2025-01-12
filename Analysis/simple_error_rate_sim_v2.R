## Basic simulations of error rates but using directly generated p-values

library(here)
library(abind)
library(tidyverse)
library(parallel)

#' Simulate P-values on a tree
#'
#' Assuming no treatement effects at all. Generate one p-value at the root
#' level. Then generate J p-values, one for each of j=1...J child nodes ---
#' where each node has a true null of no effects and where the minimum reported
#' p-value for each node under the root node is the minimum p-value <= alpha
#' and the maximum is 1. Do this for each level. Here, assuming the same number
#' of nodes per level throughout the tree below the root level.
#' For now, we assume that there are no false nulls. This is a completely no effects world.
#' should make .05 of values less than .05
sim_ps <- function(num_levels, nodes_per_level, alpha, return_only_errors = FALSE, test_limit = 1e+6) {
  stopifnot(nodes_per_level >= 2)
  stopifnot(num_levels >= 2)

  ## This generates p-values from a uniform distribution (i.e. all hypotheses are true)
  ## for each parent node (ex. given p_0 and anumber of nodes at level 1=2,
  ## this generates two p_values labeled "p_1" with both of their minimum of
  ## p_0)
  get_child_ps <- function(parent_p, num_nodes = nodes_per_level, depth) {
    node_ps <- sapply(1:num_nodes, function(i) {
      thep <- runif(1, min = parent_p, max = 1)
      res <- c(thep, parent_p)
      names(res) <- c(
        paste("p_", depth, sep = ""),
        paste("p_", depth - 1, sep = "")
      )
      return(res)
    })
    return(t(node_ps))
  }

  ## Start with the root node.
  p_0 <- runif(1, min = 0, max = 1)
  sim_id_level_0 <- "0"

  ## for now, require at least a root node and 1 set of child nodes
  stopifnot(num_levels > 0)

  ## Setup a list to receive the data from the other levels
  res_list <- vector("list", length = num_levels)

  ## Then do the first level after the root by hand
  p_1_dt <- as_tibble(get_child_ps(parent_p = p_0, depth = 1))
  p_1_dt$sim_id_level_0 <- rep(sim_id_level_0, each = nodes_per_level)
  p_1_dt$sim_id_level_1 <- paste(p_1_dt$sim_id_level_0, rep(1:nodes_per_level), sep = ".")

  if (!return_only_errors) {
    ## Setup a list to receive the data from the other levels
    res_list <- vector("list", length = num_levels)
    res_list[[1]] <- p_1_dt

    for (i in 2:num_levels) {
      res <- lapply(res_list[[i - 1]][[paste0("p_", i - 1)]], function(parent_p) {
        as_tibble(get_child_ps(parent_p = parent_p, depth = i))
      })
      res_dt <- dplyr::bind_rows(res)
      res_dt[[paste0("sim_id_level_", i - 1)]] <- rep(res_list[[i - 1]][[paste0("sim_id_level_", i - 1)]],
        each = nodes_per_level
      )
      res_dt[[paste0("sim_id_level_", i)]] <- paste(res_dt[[paste0("sim_id_level_", i - 1)]],
        rep(1:nodes_per_level),
        sep = "."
      )
      res_list[[i]] <- res_dt
      rm(res, res_dt)
    }

    ## A function to allow the use of the reduce strategy to merge multiple datasets stored in a list together in order.
    my_join <- function(xdat, ydat) {
      xdat_id <- grep("^sim_id_level_[0-9]", names(xdat), value = TRUE)
      ydat_id <- grep("^sim_id_level_[0-9]", names(ydat), value = TRUE)
      stopifnot(nrow(xdat) >= nrow(ydat))
      left_join(xdat, ydat, by = intersect(xdat_id, ydat_id), suffix = c("", ".y"))
    }

    p_dt <- purrr::reduce(res_list[length(res_list):1], my_join)
    ## Just a few tests
    stopifnot(all(p_dt$p_1 >= p_dt$p_0))
    stopifnot(all.equal(p_dt$p_1.y, p_dt$p_1))
    ## Drop all columns that are duplicated
    p_dt <- p_dt %>% select(-ends_with(".y"))

    return(p_dt)
  } else {
    tot_tests <- sum(nodes_per_level^seq(0, num_levels))
    stopifnot("You have asked for too many tests" = tot_tests < test_limit)

    fp_vec <- rep(NA, length = num_levels + 1)
    names(fp_vec) <- paste0("lv_", seq(0, num_levels))
    fp_vec["lv_0"] <- as.numeric(p_0 <= alpha)
    fp_vec["lv_1"] <- as.numeric(any(p_1_dt$p_1 <= alpha))

    for (i in 2:num_levels) {
      res <- lapply(res_list[[i - 1]][[paste0("p_", i - 1)]], function(parent_p) {
        as_tibble(get_child_ps(parent_p = parent_p, depth = i))
      })
      res_dt <- dplyr::bind_rows(res)
      fp_vec[names(fp_vec)[i + 1]] <- as.numeric(any(res_dt[[paste0("p_", i)]] <= alpha))
      rm(res, res_dt)
    }
    return(fp_vec)
  }
  ## End
}

get_max_fpr <- function(parent_lv, child_lv, obj, alpha = .05) {
  ## Return the proportion of the p-values that are less than .05 under a given node
  res <- obj %>%
    group_by(.data[[paste0("sim_id_level_", parent_lv)]]) %>%
    summarize(
      fp = mean(.data[[paste0("p_", child_lv)]] <= alpha),
      n_fp = sum(.data[[paste0("p_", child_lv)]] <= alpha),
      n_tests = length(unique(.data[[paste0("sim_id_level_", child_lv)]])),
      max_p = max(.data[[paste0("p_", child_lv)]]),
    ) %>%
    ungroup()

  ## res0 <- obj %>%
  ##  group_by(.data[[paste0("sim_id_level_", parent_lv)]]) %>%
  ##  reframe(mean_fpr = unique(.data[[paste0("p_", child_lv)]])) %>%
  ##  ungroup()

  ## This is the maximum false positive rate across nodes for a given level in the tree
  ## return(max(res$mean_fpr))
  return(res)
}

#### Now use the simulation function

## set.seed(12345)
## num_levs <- 10
## nodes_per_lev <- 2
## alpha <- .1
## nsims <- 1000000
## sims_l10_n2_lst <- mclapply(1:nsims, function(i) {
##   res <- sim_ps(num_levels = num_levs, nodes_per_level = nodes_per_lev, alpha = alpha, return_only_errors = TRUE)
##   return(res)
## },
## mc.cores = 10
## )
## sims_l10_n2_dt <- bind_rows(sims_l10_n2_lst, .id = "rep_id")
## summary(sims_l10_n2_dt)
## ## So, FWER for the whole procedure is proportion of p<=alpha across all reps:
## fwer_l10_n2 <- mean(unlist(sims_l10_n2_lst))
## fwer_l10_n2

## So this seems to have weak control when nodes_per_lev=2 and levels=10, now
## lets look across nodes per level and levels my prior is that as nodes per
## level increase the FWER will increase. But recall that since the fpr at node
## level 0 is bounded and controlled, I'm not sure how fast it will increase
## (i.e. some results suggest it is a function of \alpha^depth not sure about these though).

nodes_per_level <- sort(unique(c(seq(2, 10, 1), seq(10, 100, 10))))
number_levels <- seq(2, 20)
test_limit <- 1e+9

sim_parms <- expand.grid(num_levels = number_levels, num_nodes = nodes_per_level)
sim_parms <- sim_parms %>%
  rowwise() %>%
  mutate(tot_tests = sum(num_nodes^seq(0, num_levels))) %>%
  ungroup()
dim(sim_parms)
sim_parms_limited <- sim_parms %>% filter(tot_tests <= test_limit)
dim(sim_parms_limited)

fwer_sim_ps <- function(num_levels, num_nodes, nsims, alpha, test_limit = 1e+8) {
  res <- replicate(nsims, sim_ps(
    num_levels = num_levels,
    nodes_per_level = num_nodes, alpha = alpha, return_only_errors = TRUE, test_limit = test_limit
  ))
  ## The proportion of p_values <= alpha across the simulations for a given number of levels and number of nodes
  fwer <- mean(res)
  return(fwer)
}

blah <- fwer_sim_ps(num_levels = sim_parms$num_levels[1], num_nodes = sim_parms$num_nodes[1], nsims = 1000, alpha = .1)

set.seed(12345)
nsims <- 10000
alpha <- .1
res_lst <- mclapply(1:nrow(sim_parms_limited), function(i) {
  nl <- sim_parms_limited$num_levels[i]
  nn <- sim_parms_limited$num_nodes[i]
  message(nl, "-", nn, "tot tests =", sum(nn^seq(0, nl)))
  fwer_sim_ps(num_levels = nl, num_nodes = nn, nsims = nsims, alpha = .1, test_limit = test_limit)
}, mc.cores = 10)
sim_parms_limited$res_vec <- simplify2array(res_lst)

## Notice control of FWER for all combinations
head(sim_parms_limited)
summary(sim_parms_limited)
quantile(sim_parms_limited$res_vec, seq(0, 1, .1))
mean(sim_parms_limited$res_vec <= .1)

## Notice only one test with 19 or 20 levels.
sim_parms_limited %>%
  filter(num_levels > 15) %>%
  arrange(num_levels, num_nodes)

library(viridis)

## We good control here. Yes we have increasing FWER in both number of nodes and number of levels, but
## no sign of approaching .1.
## In fact, with a binary tree, we see the basic FWER is very low.-- the more levels the lower the FWER at low levels of nodes.
## And the lines do not cross (or almost not at all)

g_lines <- ggplot(
  data = sim_parms_limited,
  aes(groups = num_levels, color = factor(num_levels), x = num_nodes, y = res_vec)
) +
  geom_line()

g_lines

g_res <- ggplot(data = sim_parms_limited, aes(x = num_levels, y = num_nodes, z = res_vec)) +
  geom_contour_filled(breaks = c(0, .1, 1)) +
  labs(
    x = "Number of levels",
    y = "Number of nodes per level",
    fill = "FWER"
  ) +
  geom_text(label = round(res_vec, 3), colour = "white", fontface = "bold") +
  scale_fill_viridis(discrete = TRUE) +
  # theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16)
  )

g_res

#### older below when we were using the raw p-values

summary(sims_l10_n2_dt$p_0)
mean(sims_l10_n2_dt$p_0 <= alpha)

## Focusing on the cases where p_0 <= .05
sims_l10_n2_dt_sig <- sims_l10_n2_dt %>%
  filter(p_0 <= alpha) %>%
  droplevels()

get_max_fpr(parent_lv = 1, child_lv = 2, obj = sims_l10_n2_dt, alpha = .1)

## So the overall question is proportion of p_0 <= alpha across simulations (this is the simple FPR for p_0) but also:
## whether any of the 2^1 p_1 <= alpha (number of tests/nodes/child nodes at a given level is num_nodes_per_level^level_number)
## whether:
## any of the 2^2 p_2 <=alpha (the proportion across simulations should be less than alpha): ((# p_2 <= alpha)/4)/# sims is the FWER by level 2.
##

tmp <- obj %>%
  group_by(rep_id) %>%
  reframe(p_1 = unique(p_1)) %>%
  group_by(rep_id) %>%
  summarize(
    fp_p_1 = mean(p_1 <= alpha),
    n_p_1_le_alpha = sum(p_1 <= alpha),
    n_tests = n()
  )
stopifnot(length(unique(tmp$n_tests)) == 1)
table(tmp$fp_p_1)
table(tmp$n_p_1_le_alpha)
## Out of 1000 sims, only 8/1000 have any false positives at level 1.
## If this tree stopped at level 1, then we would have an FWER of number of fp at level 0+number of fp at level 1.
mean(tmp$n_p_1_le_alpha)

obj %>%
  filter(rep_id == 191) %>%
  reframe(p_1 = unique(p_1))

obj %>%
  group_by(rep_id) %>%
  reframe(p_1 = unique(p_1))

res_l10_n2_tmp <- lapply(seq(1, num_levs - 1), function(parent_lv) {
  get_max_fpr(
    parent_lv = parent_lv, child_lv = parent_lv + 1,
    obj = sims_l10_n2_dt_sig
  )
})

res_l10_n2 <- c(mean(sims_l10_n2_dt$p_0 <= alpha), res_l10_n2_tmp)
names(res_l10_n2) <- seq(0, num_levs)
res_l10_n2


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
