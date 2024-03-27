# Remember to change 1000 to nsims below
library(data.table)
library(stringi)
library(here)
load(here("Analysis", "simparms.rda"))

## Relies on simparms being keyed to idx
stopifnot(key(simparms) == "idx")

## We have been doing 1000 sims
simparms[, nms := do.call(paste, c(.SD, "1000", sep = "_"))]

done <- stri_replace_all_fixed(list.files(path = here("Analysis/CSVS"), patt = "[0-9]*csv"), replacement = "", pattern = ".csv", vectorize_all = TRUE)

idx <- simparms$idx # seq_len(nrow(simparms))
done_idx <- simparms$idx[simparms$nms %in% done]
not_done_idx <- setdiff(idx, done_idx)
save(not_done_idx, file = here("Analysis", "not_done_idx.rda"))
