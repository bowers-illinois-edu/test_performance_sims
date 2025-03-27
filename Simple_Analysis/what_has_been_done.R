# Remember to change 1000 to nsims below
library(data.table)
library(stringi)
library(here)
load(here("Simple_Analysis", "sim_parms.rda"))

## Relies on simparms being keyed to idx
stopifnot(key(sim_parms) == "idx")

sim_parms[, nms := do.call(paste, c("sim", .SD[, 1:7], sep = "_"))]

done <- stri_replace_all_fixed(list.files(path = here("Simple_Analysis/CSVS_latest"), patt = "[0-9]*csv"),
  replacement = "", pattern = ".csv", vectorize_all = TRUE
)

idx <- sim_parms$idx # seq_len(nrow(simparms))
done_idx <- sim_parms$idx[sim_parms$nms %in% done]
not_done_idx <- setdiff(idx, done_idx)

stopifnot(all(not_done_idx %in% sim_parms$idx))

save(not_done_idx, file = here("Simple_Analysis", "not_done_idx.rda"))
