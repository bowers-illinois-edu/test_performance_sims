## Basic simulations of error rates but using directly generated p-values

library(abind)

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

  p_0 <- runif(sims, min = 0, max = 1)
  p_0_le_alpha <- p_0[p_0 <= alpha]
  node_fpr_0 <- mean(p_0 <= alpha)

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

  if (node_fpr_0 == 0) {
    return(c(node_fpr_0 = 0, p_0 = min(p_0), p_0 = min(p_0)))
  }

  ## for now, require at least a root node and children nodes
  stopifnot(num_levels > 1)

  p_1s <- lapply(p_0_le_alpha, function(parent_p) {
    get_child_ps(parent_p = parent_p, depth = 1)
  })
  stopifnot(all.equal(length(p_1s), length(p_0_le_alpha)))
  p_1_dt <- dplyr::bind_rows(p_1s)
  stopifnot(all.equal(nrow(p_1_dt), 2 * length(p_0_le_alpha)))

  res_list <- vector("list", length = num_levels)
  res_list[[1]] <- p_1_dt

  for (i in 2:num_levels) {
    res <- lapply(res_list[[i - 1]][[paste0("p_", i - 1)]], function(parent_p) {
      get_child_ps(parent_p = parent_p, depth = i)
    })
    res_list[[i]] <- dplyr::bind_rows(res)
  }

  # p_dt <- left_join(p_2_dt, p_1_dt, by = join_by(p1 == min_p_le_alpha))
  my_join <- function(xdat, ydat) {
    xdat_p_nms <- grep("^p_[0-9]", names(xdat), value = TRUE)
    ydat_p_nms <- grep("^p_[0-9]", names(ydat), value = TRUE)
    stopifnot(nrow(xdat) >= nrow(ydat))
    left_join(xdat, ydat, by = intersect(xdat_p_nms, ydat_p_nms))
  }

  p_dt <- purrr::reduce(res_list[length(res_list):1], my_join)
  # p_dt1 <- purrr::accumulate(res_list[length(res_list):1], my_join)
  stopifnot(all(p_dt$p_3 >= p_dt$p_2))
  stopifnot(all(p_dt$p_2 >= p_dt$p_1))
  stopifnot(all(p_dt$p_1 >= p_dt$p_0))
  p_dt$node_fpr_0 <- node_fpr_0

  return(p_dt)
  ## End
}

set.seed(123)
tmp <- sim_ps(num_levels = 3, nodes_per_level = 2, alpha = .05, sims = 10000)
str(tmp)
head(tmp)
## Only the root level tests <= alpha in p_0
summary(tmp$p_0)

fwer_level_1 <- tmp %>%
  group_by(p_0) %>%
  summarize(fpr = mean(node_fpr_1))

summary(fwer_level_1)

fwer_level_2 <- tmp %>%
  group_by(p_1) %>%
  summarize(fpr = mean(node_fpr_2))

summary(fwer_level_2)


fwer_level_3 <- tmp %>%
  group_by(p_2) %>%
  summarize(fpr = mean(node_fpr_3))

summary(fwer_level_3)
