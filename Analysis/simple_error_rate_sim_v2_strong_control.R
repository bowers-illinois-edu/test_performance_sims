## The idea here is to directly generate p_values from either a uniform
## distribution (for hypothese that are true ) and p_values from a beta
## distribution for hypotheses that are false.

## We know that the total number of tests in a tree is
## sum(num_nodes^seq(0,num_levels)). For example, in a tree with 2 nodes per
## level and 3 levels we have 2^0+2^1+2^2+2^3=1+2+4+8=15 tests. Now, which
## tests should be false? We could randomly choose some proportion prop_false
## out of the total. Some of them will be at the higher levels and some the
## lower levels. We could also do this by flipping a coin for each tests: is
## this a true test versus a false test.

## Alternatively we could set prop_false for the proportion false among the
## leave node hypotheses (which are num_nodes^num_levels in number). And then
## any parent_p between that child/leaf node and the root is also false (i.e.
## draws from rbeta rather than runif).


library(here)
library(abind)
library(tidyverse)
library(parallel)


# Function to generate Beta random variables with minimum 'a' and maximum 'c'
beta_with_min <- function(n, alpha, beta, a, c) {
  # Generate standard Beta random variables
  x <- rbeta(n, shape1 = alpha, shape2 = beta)
  # Transform to desired range
  y <- a + (c - a) * x
  return(y)
}

# Example usage
set.seed(123)
samples <- beta_with_min(
  n = 1000, # Number of samples
  alpha = .1, # Shape parameter 1
  beta = 1, # Shape parameter 2
  a = 0, # Minimum value
  c = 1 # Maximum value
)
summary(samples)
# plot(density(samples))

#' Assign false hypotheses in a k-ary tree, propagate upward, and label the nodes.
#'
#' The children of node i are the indices \{\, (i-1)\cdot k + 2,\; (i-1)\cdot k + 3,\; \dots,\; i\cdot k + 1 \},
#' as long as they do not exceed the total number of nodes.
#' Leaves occupy the last k^{L-1} indices of \{1,\dots,N\}, where
#' N = \frac{k^L - 1}{k - 1}.
#'
#' @param nodes_per_level k: branching factor (each non-leaf node has k children)
#' @param num_levels L: total number of levels (root is level 1, leaves are level L)
#' @param prop_false numeric in [0, 1], the proportion of leaf nodes to mark as TRUE
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item node_id: integer ID (1-based) in array representation
#'     \item label: hierarchical character label (e.g. "lv_0", "lv_1", "lv_1_1", ...)
#'     \item has_false_hyp: logical indicating if node is TRUE or FALSE
#'   }
#'
#' @examples
#' set.seed(123)
#' df <- assign_false_H_with_labels(2, 3, 0.4)
#' print(df)
assign_false_H_with_labels <- function(nodes_per_level, num_levels, prop_false) {
  k <- as.integer(nodes_per_level)
  L <- as.integer(num_levels)

  if (k < 2) {
    stop("'nodes_per_level' (k) must be at least 2.")
  }
  if (L < 1) {
    stop("'num_levels' must be at least 1.")
  }
  if (prop_false < 0 || prop_false > 1) {
    stop("'prop_false' must be between 0 and 1.")
  }

  # Total number of nodes in a complete k-ary tree of depth L:
  # N = (k^L - 1) / (k - 1), root is level 1, leaves are level L
  total_nodes <- (k^L - 1) / (k - 1)
  total_nodes <- as.integer(total_nodes)

  # Number of leaves
  num_leaves <- k^(L - 1)

  # Number of leaves to mark as TRUE
  num_selected <- floor(prop_false * num_leaves)

  # Randomly pick which leaves are TRUE
  selected_leaves <- sample.int(num_leaves, num_selected)

  # Initialize the logical vector
  has_false_hyp <- logical(total_nodes)

  # Leaves occupy the last 'num_leaves' indices of [1..total_nodes]
  leaf_start <- total_nodes - num_leaves + 1
  has_false_hyp[leaf_start + selected_leaves - 1] <- TRUE

  # -------------------------------
  # 2) Propagate TRUE upward.
  #
  # Using standard array logic for a k-ary tree (1-based):
  #   If has_false_hyp[i] == TRUE and i > 1, then set parent to TRUE.
  #   The parent of node i (i>1) is floor((i - 2)/k) + 1.
  # -------------------------------
  for (i in seq.int(total_nodes, 2)) {
    if (has_false_hyp[i]) {
      parent <- floor((i - 2) / k) + 1L
      has_false_hyp[parent] <- TRUE
    }
  }

  # -------------------------------
  # 3) Label the nodes in BFS order:
  #    - Node 1 is the root => "lv_0".
  #    - The children of "lv_0" => "lv_1", "lv_2", ..., "lv_k".
  #    - Then the children of "lv_1" => "lv_1_1", "lv_1_2", ...
  #    - etc.
  # -------------------------------
  labels <- character(total_nodes)
  labels[1] <- "lv_0"

  # We'll do a queue-based BFS. For each node i:
  #   children are from k*(i-1)+2  to  k*i+1
  #   (clipped to total_nodes)
  queue <- 1L

  while (length(queue) > 0) {
    node <- queue[1]
    queue <- queue[-1]

    child_start <- (node - 1) * k + 2
    child_end <- node * k + 1

    if (child_start <= total_nodes) {
      child_end <- min(child_end, total_nodes)
      # label each child
      for (child in seq.int(child_start, child_end)) {
        # offset among siblings => 1..k
        offset <- child - child_start + 1
        # If parent is "lv_0", children => "lv_1", "lv_2", ...
        if (labels[node] == "lv_0") {
          labels[child] <- paste0("lv_", offset)
        } else {
          # children => "lv_1_1", etc.
          labels[child] <- paste0(labels[node], "_", offset)
        }
      }
      # Enqueue children
      queue <- c(queue, child_start:child_end)
    }
  }

  # Make a data frame
  out <- data.frame(
    node_id = seq_len(total_nodes),
    label = labels,
    has_false_hyp = has_false_hyp,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  return(out)
}

set.seed(123)
df <- assign_false_H_with_labels(nodes_per_level = 2, num_levels = 3, prop_false = 0.4)
print(df)

df <- assign_false_H_with_labels(nodes_per_level = 3, num_levels = 4, prop_false = 0.3)
print(df)

library(tidygraph)
library(ggraph)


#' Plot a k-ary tree (from assign_false_H_with_labels) using ggraph.
#'
#' @param df A data.frame with columns:
#'   \itemize{
#'     \item \code{node_id} (integer) unique 1-based IDs.
#'     \item \code{label} (character) human-readable labels, e.g. "lv_0", "lv_1", "lv_1_1".
#'     \item \code{has_false_hyp} (logical) whether node is TRUE or FALSE.
#'   }
#' @param nodes_per_level The branching factor k (same \code{k} you used when creating \code{df}).
#' @param layout Character string specifying the layout in \code{ggraph}.
#'   Defaults to \code{"tree"}. Other options include \code{"dendrogram"},
#'   \code{"circlepack"}, etc.
#' @param flip_tree Logical. If \code{TRUE}, flip the tree vertically (root at top).
#'   Defaults to \code{FALSE}.
#'
#' @return A \code{ggplot} object which, when printed, displays the tree.
#'
#' @examples
#' # Suppose df is from assign_false_H_with_labels(k=2, L=3, prop_false=0.4)
#' # plot_kary_tree(df, nodes_per_level=2)
plot_kary_tree <- function(df, nodes_per_level, layout = "tree", flip_tree = FALSE) {
  if (!requireNamespace("tidygraph", quietly = TRUE)) {
    stop("Package 'tidygraph' is required but is not installed.")
  }
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("Package 'ggraph' is required but is not installed.")
  }

  k <- as.integer(nodes_per_level)
  total_nodes <- nrow(df)

  # Build a parent->child edge list from the standard array logic for a k-ary tree:
  # For node i > 1, parent(i) = floor((i - 2)/k) + 1
  edges <- data.frame(from = numeric(0), to = numeric(0))
  for (i in 2:total_nodes) {
    parent_i <- floor((i - 2) / k) + 1L
    edges <- rbind(edges, data.frame(from = parent_i, to = i))
  }

  # Construct a tbl_graph
  # We'll rename df columns for clarity in the graph
  graph <- tidygraph::tbl_graph(
    nodes = df,
    edges = edges,
    directed = TRUE
  )

  # Plot with ggraph
  p <- ggraph(graph, layout = layout) +
    geom_edge_link(
      arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
      end_cap = circle(3, "mm"),
      start_cap = circle(3, "mm"),
      colour = "gray60"
    ) +
    geom_node_point(aes(color = has_false_hyp), size = 5) +
    geom_node_text(aes(label = label), vjust = -0.8, size = 3.5) +
    # Map TRUE -> red, FALSE -> blue:
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
    theme_minimal()

  # Optionally flip the tree (put root at top)
  if (flip_tree) {
    p <- p + coord_flip()
  }

  return(p)
}

# 1) Build your k-ary tree data frame:
set.seed(123)
df <- assign_false_H_with_labels(nodes_per_level = 3, num_levels = 4, prop_false = 0.4)

# 2) Plot it:
p <- plot_kary_tree(df, nodes_per_level = 3)
print(p)


# 1) Build your k-ary tree data frame:
set.seed(123)
df <- assign_false_H_with_labels(nodes_per_level = 5, num_levels = 3, prop_false = 0.1)

# 2) Plot it:
p <- plot_kary_tree(df, nodes_per_level = 5)
print(p)



#' Simulate P-values on a tree where some nodes have treatment effects
#' (hypotheses are false) and others have true nulls (no treatment effects)
#'
#' Specify Generate one p-value at the root
#' level. Then generate J p-values, one for each of j=1...J child nodes ---
#' where each node has a true null of no effects and where the minimum reported
#' p-value for each node under the root node is the minimum p-value <= alpha
#' and the maximum is 1. Do this for each level. Here, assuming the same number
#' of nodes per level throughout the tree below the root level. For now, we
#' assume that there are no false nulls. This is a completely no effects world.
#' should make .05 of values less than .05
sim_ps <- function(num_levels, nodes_per_level, alpha, return_only_errors =
                     FALSE, test_limit = 1e+6) {
  stopifnot(nodes_per_level >= 2)
  stopifnot(num_levels >= 2)

  ## Number of nodes in a k-ary tree: canonical formula has k nodes per level and n levels
  ## tot_tests_1 <- sum(nodes_per_level^seq(0, (num_levels-1)))
  tot_tests_2 <- ((nodes_per_level^num_levels) - 1) / (nodes_per_level - 1)
  stopifnot(tot_tests_1 == tot_tests_2)

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
    tot_tests_1 <- sum(nodes_per_level^seq(0, num_levels))
    tot_tests_2 <- ((nodes_per_level^num_levels) - 1) / (nodes_per_level - 1)
    stopifnot(tot_tests_1 == tot_tests_2)

    stopifnot("You have asked for too many tests" = tot_tests_2 < test_limit)

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
