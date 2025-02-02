## Convert the CSV file of simulation results to an R data frame
library(data.table)
library(here)

simp_simsres_unadj <- fread(here::here("Simple_Analysis", "simple_sims_unadj_merged_results.csv"))

save(simp_simsres_unadj, file = here::here("Simple_Analysis", "simple_sims_unadj_results.rda"))
