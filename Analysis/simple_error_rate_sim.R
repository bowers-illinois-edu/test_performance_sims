## Basic simulations of error rates but using directly generated p-values

library(abind)

sim_ps <- function(num_levels, nodes_per_level, alpha, sims) {
  ## For now, we assume that there are no false nulls. This is a completely no effects world.
  set.seed(125)
  p_0 <- runif(sims) ## should make 5% of values less than .05
  p_0_le_alpha <- p_0[p_0 <= alpha]
  fpr_0 <- mean(p_0 <= alpha)
  if (fpr_0 == 0) {
    return(p_0)
  }

  p_1s <- lapply(p_0_le_alpha, function(p0) {
    node_ps_lst <- lapply(1:nodes_per_level, function(node_num) {
      thep <- runif(sims, min = p0, max = 1)
      thep_le_alpha <- max(0, thep[thep <= alpha])
      ## return(matrix(cbind(rep(p0, length(thep)), thep), ncol = 2))
      return(c(
        node_fpr = mean(thep <= alpha),
        min_p_le_alpha = min(thep_le_alpha),
        pre_p = p0
      ))
    })
    res <- simplify2array(node_ps_lst)
  })



  ## End
}
