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

source(here("Analysis", "simple_p_draws_fns.R"))

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
