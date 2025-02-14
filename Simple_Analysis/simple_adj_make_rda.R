## Convert the CSV file of simulation results to an R data frame
library(data.table)
library(here)

simp_simsres_adj <- fread(here::here("Simple_Analysis", "simple_sims_adj_merged_results.csv"))

save(simp_simsres_adj, file = here::here("Simple_Analysis", "simple_sims_adj_results.rda"))
