## Convert the CSV file of simulation results to an R data frame
library(data.table)
library(here)

simp_simsres <- fread(here::here("Analysis", "simple_sims_merged_results.csv"))

save(simp_simsres, file = here::here("Analysis", "simple_sims_results.rda"))
