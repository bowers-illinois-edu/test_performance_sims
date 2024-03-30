## Convert the CSV file of simulation results to an R data frame
library(data.table)
library(here)

simsres <- fread(here::here("Analysis", "merged_results.csv"))

save(simsres, file = here::here("Analysis", "simsres.rda"))
