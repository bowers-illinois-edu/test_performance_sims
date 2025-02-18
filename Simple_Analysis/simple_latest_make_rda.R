## Convert the CSV file of simulation results to an R data frame
library(data.table)
library(here)

simp_simsres_latest <- fread(here::here("Simple_Analysis", "simple_sims_latest_merged_results.csv"))

save(simp_simsres_latest, file = here::here("Simple_Analysis", "simple_sims_latest_results.rda"))
