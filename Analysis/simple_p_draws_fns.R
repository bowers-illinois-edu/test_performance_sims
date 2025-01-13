library(tidyverse)
library(tidygraph)
library(ggraph)

#' Draw Beta random variables in [themin, themax].
#'
#' @param n number of draws
#' @param alpha,beta shape parameters for Beta
#' @param themin,themax numeric boundaries of range
#'
#' @return A numeric vector of length n
beta_with_min <- function(n, alpha = 0.1, beta = 1, themin = 0, themax = 1) {
  x <- stats::rbeta(n, shape1 = alpha, shape2 = beta)
  y <- themin + (themax - themin) * x
  return(y)
}

#' Build a k-ary tree of depth L, and randomly assign some leaves as TRUE.
#' Then propagate TRUE up to the root. Also label the nodes.
#'
#' @param k integer, the branching factor (number of children per node).
#' @param L integer, the number of levels (root is level 1, leaves are level L).
#' @param prop_false numeric in [0,1], fraction of leaves to mark as TRUE.
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item node_id (1-based ID in array representation)
#'     \item label (character, e.g. "lv_0", "lv_1", "lv_1_1", etc.)
#'     \item has_false_hyp (logical, TRUE or FALSE)
#'   }
assign_false_H_with_labels <- function(k, L, prop_false = 0.3) {
  if (k < 1) stop("k must be >= 1")
  if (L < 1) stop("L must be >= 1")
  if (prop_false < 0 || prop_false > 1) stop("prop_false must be in [0,1]")

  # Total number of nodes
  total_nodes <- (k^L - 1) / (k - 1)
  total_nodes <- as.integer(total_nodes)

  # Number of leaves
  num_leaves <- k^(L - 1)

  # Mark some leaves as TRUE
  num_selected <- floor(prop_false * num_leaves)
  selected_leaves <- sample.int(num_leaves, size = num_selected, replace = FALSE)

  # Initialize
  has_false_hyp <- logical(total_nodes) # all FALSE
  leaf_start <- total_nodes - num_leaves + 1
  has_false_hyp[leaf_start + selected_leaves - 1] <- TRUE

  # Propagate TRUE up to the root
  for (i in seq.int(total_nodes, 2)) {
    if (has_false_hyp[i]) {
      parent_i <- floor((i - 2) / k) + 1L
      has_false_hyp[parent_i] <- TRUE
    }
  }

  # Labels
  labels <- character(total_nodes)
  labels[1] <- "lv_0"

  # BFS labeling
  queue <- 1L
  while (length(queue) > 0) {
    node <- queue[1]
    queue <- queue[-1]

    child_start <- (node - 1) * k + 2
    child_end <- node * k + 1
    if (child_start <= total_nodes) {
      child_end <- min(child_end, total_nodes)
      for (child in seq(child_start, child_end)) {
        offset <- child - child_start + 1
        if (labels[node] == "lv_0") {
          labels[child] <- paste0("lv_", offset)
        } else {
          labels[child] <- paste0(labels[node], "_", offset)
        }
      }
      queue <- c(queue, child_start:child_end)
    }
  }

  # Return data frame
  df <- data.frame(
    node_id = seq_len(total_nodes),
    label = labels,
    has_false_hyp = has_false_hyp,
    stringsAsFactors = FALSE
  )
  rownames(df) <- NULL
  return(df)
}

#' BFS Approach: Given a k-ary tree data frame (from assign_false_H_with_labels),
#' draw a random value for each node, with the min determined by parent's value.
#'
#' - If a node is TRUE, we draw from Beta_with_min().
#' - If a node is FALSE, we draw from Uniform().
#'
#' The root uses themin=0 (or user-provided p0_min) for the initial draw range.
#'
#' @param tree_df data.frame with columns (node_id, label, has_false_hyp).
#' @param k branching factor.
#' @param alpha,beta shape parameters for Beta.
#' @param p0_min numeric, the minimum for the root node's random draw (default 0).
#' @return The same data frame, but with an extra column "rand_value".
draw_values_along_tree <- function(tree_df,
                                   k,
                                   alpha = 0.1,
                                   beta = 1,
                                   p0_min = 0) {
  n <- nrow(tree_df)

  rand_value <- numeric(n)

  # BFS queue
  queue <- 1L # root is node_id=1
  # Root draw
  if (tree_df$has_false_hyp[1]) {
    rand_value[1] <- beta_with_min(
      n = 1, alpha = alpha, beta = beta,
      themin = p0_min, themax = 1
    )
  } else {
    rand_value[1] <- runif(1, min = p0_min, max = 1)
  }

  while (length(queue) > 0) {
    node <- queue[1]
    queue <- queue[-1]

    val_parent <- rand_value[node]
    child_start <- (node - 1) * k + 2
    child_end <- node * k + 1

    if (child_start <= n) {
      child_end <- min(child_end, n)

      for (child in seq(child_start, child_end)) {
        if (tree_df$has_false_hyp[child]) {
          # Beta
          rand_value[child] <- beta_with_min(
            n = 1, alpha = alpha, beta = beta,
            themin = val_parent, themax = 1
          )
        } else {
          # Uniform
          rand_value[child] <- runif(1, min = val_parent, max = 1)
        }
      }
      # Enqueue children
      queue <- c(queue, child_start:child_end)
    }
  }

  # Add column to data frame
  out <- tree_df
  out$rand_value <- rand_value
  return(out)
}

#' BFS Approach: Generate a k-ary tree with T/F assignment and random draws
#' according to the rules described in the question.
#'
#' @param k branching factor
#' @param L number of levels
#' @param prop_false fraction of leaves to assign as TRUE
#' @param alpha,beta shape parameters for Beta
#' @param p0_min numeric, the min for the root node's range (default=0).
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item node_id
#'     \item label
#'     \item has_false_hyp
#'     \item rand_value
#'   }
#'
#' @examples
#' set.seed(123)
#' df <- generate_tree_data(k = 2, L = 3, prop_false = 0.4)
#' df
generate_tree_data_bfs <- function(k, L, prop_false = 0.3,
                                   alpha = 0.1, beta = 1,
                                   p0_min = 0) {
  # Step 1 & 2: build tree, assign T/F
  tree_df <- assign_false_H_with_labels(k, L, prop_false)

  # Step 3 & 4: draw random values
  out <- draw_values_along_tree(tree_df, k, alpha, beta, p0_min)
  return(out)
}

## Example:
## If node is TRUE, this means draw from Beta (i.e. true treatment effect)
set.seed(101)
df_bfs <- generate_tree_data_bfs(k = 4, L = 3, prop_false = 0.6, alpha = 0.5, beta = 2, p0_min = 0)
df_bfs %>% arrange(node_id)

# test-generate_tree_data.R
library(testthat)

test_that("generate_tree_data produces correct shape of output", {
  set.seed(123)
  k <- 4
  L <- 3
  prop_false <- 0.4

  df_test <- generate_tree_data_bfs(
    k = k, L = L, prop_false = prop_false,
    alpha = 0.1, beta = 1, p0_min = 0
  )

  # total nodes
  total_nodes <- (k^L - 1) / (k - 1)
  expect_equal(nrow(df_test), total_nodes)
  expect_setequal(names(df_test), c("node_id", "label", "has_false_hyp", "rand_value"))

  # root is node_id=1
  # if any leaf is TRUE, the root is forced to TRUE
  # but if prop_false=0, we could have a root=FALSE
  # For this test, with prop_false=0.4,  let's see if root ended up being TRUE:
  # (We won't rely on it, but let's just check that it's indeed a logical)
  expect_type(df$has_false_hyp, "logical")

  # check rand_value is numeric, all in [0,1]
  expect_type(df$rand_value, "double")
  expect_true(all(df$rand_value >= 0 & df$rand_value <= 1))

  # check that if node is T, child is drawn from Beta with min = parent's val
  # if node is F, child is from Uniform with min = parent's val
  # We won't do a super rigorous test, but we can do quick checks here to ensure that the child p is >= parent p:

  for (i in 2:nrow(df)) {
    parent_i <- floor((i - 2) / k) + 1
    # child's value >= parent's value
    expect_true(df$rand_value[i] >= df$rand_value[parent_i],
      info = paste0("Child rand_value < parent at i=", i)
    )
  }
})


#' Assign T/F to leaves in a k-ary tree (prop_false fraction = TRUE),
#' propagate TRUE to the root, and label nodes with DFS.
#'
#' @param k integer, branching factor (# children per non-leaf node).
#' @param L integer, number of levels (root is level 1, leaves are level L).
#' @param prop_false numeric in [0,1], fraction of leaves to set as TRUE.
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item \code{node_id} (1-based index in array form)
#'     \item \code{label} (character, e.g. "lv_0", "lv_0_1", "lv_0_1_1", etc.)
#'     \item \code{has_false_hyp} (logical, TRUE/FALSE)
#'   }
assign_false_H_with_labels_DFS <- function(k, L, prop_false) {
  if (k < 1) stop("k must be >= 1.")
  if (L < 1) stop("L must be >= 1.")
  if (prop_false < 0 || prop_false > 1) {
    stop("prop_false must be in [0,1].")
  }

  # Total number of nodes in a complete k-ary tree
  total_nodes <- (k^L - 1) / (k - 1)
  total_nodes <- as.integer(total_nodes)

  # Number of leaves
  num_leaves <- k^(L - 1)

  # Randomly assign some leaves as TRUE
  num_selected <- floor(prop_false * num_leaves)
  selected_leaves <- sample.int(num_leaves, size = num_selected, replace = FALSE)

  # Logical vector: all start as FALSE
  has_false_hyp <- logical(total_nodes)

  # Leaves occupy the last num_leaves indices: [leaf_start..total_nodes]
  leaf_start <- total_nodes - num_leaves + 1
  # Mark the selected leaves
  has_false_hyp[leaf_start + selected_leaves - 1] <- TRUE

  # Propagate TRUE upward to the root (bottom-up)
  # parent(i) = floor((i - 2)/k) + 1
  for (i in seq.int(total_nodes, 2)) {
    if (has_false_hyp[i]) {
      parent_i <- floor((i - 2) / k) + 1L
      has_false_hyp[parent_i] <- TRUE
    }
  }

  # Now we label the nodes with a DFS approach.
  labels <- character(total_nodes)
  # Root is node 1
  # We'll push (node, label_for_node) onto a stack
  stack <- list(list(node = 1L, label = "lv_0"))

  while (length(stack) > 0) {
    # Pop the top
    top <- stack[[length(stack)]]
    stack <- stack[-length(stack)]

    node <- top$node
    lbl <- top$label

    # Assign label
    labels[node] <- lbl

    # Identify children in array form:
    child_start <- (node - 1) * k + 2
    child_end <- node * k + 1
    if (child_start <= total_nodes) {
      child_end <- min(child_end, total_nodes)
      # We'll push children in reverse order so that child_start is visited last
      # if we want them in ascending order, or vice versa.
      # The order typically doesn't matter for labeling, as long as it's consistent.
      for (child in seq(child_end, child_start, by = -1)) {
        offset <- child - child_start + 1
        if (lbl == "lv_0") {
          child_label <- paste0("lv_", offset)
        } else {
          child_label <- paste0(lbl, "_", offset)
        }
        stack[[length(stack) + 1]] <- list(node = child, label = child_label)
      }
    }
  }

  df <- data.frame(
    node_id = seq_len(total_nodes),
    label = labels,
    has_false_hyp = has_false_hyp,
    stringsAsFactors = FALSE
  )
  return(df)
}


#' Draw random values along the k-ary tree using a DFS approach.
#'
#' - If node is TRUE => Beta distribution in [parent_val, 1].
#' - If node is FALSE => Uniform distribution in [parent_val, 1].
#'
#' @param tree_df data.frame with columns (node_id, label, has_false_hyp).
#' @param k branching factor
#' @param alpha,beta shape parameters for Beta
#' @param p0_min numeric, the minimum for the root node's range
#'
#' @return The same data frame + a new numeric column "rand_value".
draw_values_along_tree_DFS <- function(tree_df,
                                       k,
                                       alpha = 0.1,
                                       beta = 1,
                                       p0_min = 0) {
  n <- nrow(tree_df)
  rand_value <- numeric(n)

  # We'll use a stack of (node, parent_val)
  # For the root (node=1), parent_val = p0_min
  stack <- list(list(node = 1L, parent_val = p0_min))

  while (length(stack) > 0) {
    top <- stack[[length(stack)]]
    stack <- stack[-length(stack)]

    node <- top$node
    parent_val <- top$parent_val

    # Draw the random value for this node
    if (tree_df$has_false_hyp[node]) {
      # Beta
      rand_value[node] <- beta_with_min(
        n      = 1,
        alpha  = alpha,
        beta   = beta,
        themin = parent_val,
        themax = 1
      )
    } else {
      # Uniform
      rand_value[node] <- runif(1, min = parent_val, max = 1)
    }

    # Now push children, if any
    child_start <- (node - 1) * k + 2
    child_end <- node * k + 1

    if (child_start <= n) {
      child_end <- min(child_end, n)
      for (child in seq(child_end, child_start, by = -1)) {
        stack[[length(stack) + 1]] <- list(
          node = child,
          parent_val = rand_value[node]
        )
      }
    }
  }

  out <- tree_df
  out$rand_value <- rand_value
  return(out)
}


#' Generate data from a k-ary tree using DFS:
#' 1) Assign T/F at leaves, propagate up, label by DFS.
#' 2) Draw random values by DFS from root to leaves.
#'
#' @param k branching factor
#' @param L number of levels
#' @param prop_false fraction of leaves to be TRUE
#' @param alpha,beta shape parameters for Beta
#' @param p0_min numeric, min for the root's distribution (default=0)
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item node_id
#'     \item label
#'     \item has_false_hyp
#'     \item rand_value
#'   }
#'
#' @examples
#' set.seed(42)
#' df <- generate_tree_data_DFS(k = 2, L = 3, prop_false = 0.5, alpha = 0.5, beta = 2, p0_min = 0)
#' df
generate_tree_data_DFS <- function(k, L, prop_false,
                                   alpha = 0.1,
                                   beta = 1,
                                   p0_min = 0) {
  # Step 1: Build tree & assign T/F, label nodes by DFS
  tree_df <- assign_false_H_with_labels_DFS(k, L, prop_false)

  # Step 2: Draw random values by DFS
  out <- draw_values_along_tree_DFS(tree_df, k, alpha, beta, p0_min)
  return(out)
}

# test-generate_tree_data_DFS.R
# library(testthat)

test_that("DFS-based tree generation works for a small example", {
  set.seed(101)
  k <- 2
  L <- 3
  prop_false <- 0.4

  df <- generate_tree_data_DFS(
    k          = k,
    L          = L,
    prop_false = prop_false,
    alpha      = 0.5,
    beta       = 2,
    p0_min     = 0
  )

  # total nodes
  total_nodes <- (k^L - 1) / (k - 1)
  expect_equal(nrow(df), total_nodes)
  expect_true(all(c("node_id", "label", "has_false_hyp", "rand_value") %in% names(df)))

  # rand_value in [0,1]
  expect_true(all(df$rand_value >= 0 & df$rand_value <= 1))

  # If child is TRUE => Beta in [parent_val,1]. If child is FALSE => Uniform in [parent_val,1].
  # We'll do a quick check that child's value >= parent's value.
  for (i in 2:total_nodes) {
    parent_i <- floor((i - 2) / k) + 1
    expect_true(df$rand_value[i] >= df$rand_value[parent_i],
      info = paste("Child rand_value < parent's at node", i)
    )
  }
})


## Example:
## If node is TRUE, this means draw from Beta (i.e. true treatment effect)
set.seed(101)
df_bfs <- generate_tree_data_bfs(k = 4, L = 3, prop_false = 0.4, alpha = 0.1, beta = 1, p0_min = 0)
df_bfs %>% arrange(node_id)

set.seed(101)
df_dfs <- generate_tree_data_DFS(k = 4, L = 3, prop_false = 0.4, alpha = 0.1, beta = 1, p0_min = 0)
df_dfs %>% arrange(node_id)

stopifnot(nrow(df_bfs) == nrow(df_dfs))
stopifnot(df_bfs$label == df_dfs$label)
stopifnot(df_bfs$has_false_hyp == df_dfs$has_false_hyp)


## The two approaches yield different results
df_big <- left_join(df_bfs, df_dfs, by = c("label"), suffix = c("_bfs", "_dfs"))
df_big

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
    geom_node_text(aes(label = round(rand_value, 2)), vjust = 0.8, size = 3.5) +
    # Map TRUE -> red, FALSE -> blue:
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )

  # Optionally flip the tree (put root at top)
  if (flip_tree) {
    p <- p + coord_flip()
  }

  return(p)
}


p_bfs <- plot_kary_tree(df_bfs, nodes_per_level = 4)
print(p_bfs)
p_dfs <- plot_kary_tree(df_dfs, nodes_per_level = 4)
print(p_dfs)
