## Set up simsres data for use by the multiple figure making files.
library(here)
library(data.table)

load(here::here("Analysis", "simsres.rda"))

## What simulations do we have available?
table(simsres$nblocks)

## Missing alpha functions are fixed alpha approaches
simsres[afn == "", afn := "fixed (a=.05)"]

simsres$splitfnF <- factor(interaction(simsres$splitfn, simsres$p_adj_method, drop = TRUE),
  labels = c(
    "All Blocks (FDR adjust)",
    "All Blocks (Holm adjust)",
    "All Blocks (Hommel adjust)",
    "K-Means Split on Prognostic Covariate",
    "Equal Sample Size Split",
    "Split at Largest Block Precision Weight",
    "Pre-Specified Covariate Splits Binary",
    "Pre-Specified Covariate Splits Multi"
  )
)

simsres[, split_vs_not := factor(fifelse(splitfn != "", "Top-Down", "Bottom-Up"))]

## Check the recoding
table(simsres$splitfnF, interaction(simsres$splitfn, simsres$p_adj_method, drop = TRUE), exclude = c())
table(simsres$split_vs_not, interaction(simsres$splitfn, simsres$p_adj_method, drop = TRUE), exclude = c())

## # relabel certain variables for figures
simsres$tau_labels <- paste(expression(bar(tau)), "==", simsres$tau_sizes)
simsres$prop_blocks_0_labels <- paste(expression(prop ~ tau[bi] == 0), ":", simsres$prop_blocks_0)

save(simsres, file = here::here("figures", "simsres_plotting.rda"))
