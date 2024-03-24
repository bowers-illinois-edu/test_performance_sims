## Make sets of blocks fixed for simulations.

## Because we have some nested
## structure, we build up from the bottom --- sampling first 10 blocks, and
## then 25-10, etc.. adding the original blocks to the newly sampled blocks
## each time

library(here)
library(data.table)
load(here::here("Data", "bdat_equal_nb.rda"), verbose = TRUE)
load(here::here("Data", "bdat_not_equal_nb.rda"), verbose = TRUE)

## I'm not sure about uses factor versus character variables as keys in data.table
bdat_equal_nb[, bC := as.character(bF)]
bdat_not_equal_nb[, bC := as.character(bF)]
setkey(bdat_equal_nb, bC)
setkey(bdat_not_equal_nb, bC)

## Numbers of blocks: 10,25,50,75,100

## The nested structure is block within vb3 within vb2 within vb1. So, first we
## focus only on two levels within each so that we can have some structured
## testing possible with only 10 blocks We then add levels as the number of
## blocks that we sample increases

set.seed(12345)
blocks_not_equalB10 <- bdat_not_equal_nb[vb1 %in% c(1, 2) & vb2 %in% c(1, 2) & vb3 %in% c(1, 2), sample(bC, size = 10)]
stopifnot(length(unique(blocks_not_equalB10)) == 10)
## Now make a 25 block dataset by adding to the original 10 block dataset
## Point is to minimize simulation error between datasets
blocks_not_equalB25samp <- bdat_not_equal_nb[!blocks_not_equalB10, ][vb1 %in% c(1, 2) & vb2 %in% c(1, 2) & vb3 %in% c(1, 2), sample(bC, 25 - 10)]
blocks_not_equalB25 <- c(blocks_not_equalB10, blocks_not_equalB25samp)
stopifnot(length(unique(blocks_not_equalB25)) == 25)
blocks_not_equalB50 <- c(blocks_not_equalB25, bdat_not_equal_nb[!blocks_not_equalB25, ][vb1 %in% c(1, 2, 3) & vb2 %in% c(1, 2, 3) & vb3 %in% c(1, 2, 3), sample(bC, 50 - 25)])
stopifnot(length(unique(blocks_not_equalB50)) == 50)
blocks_not_equalB75 <- c(blocks_not_equalB50, bdat_not_equal_nb[!blocks_not_equalB50, ][vb1 %in% c(1, 2, 3) & vb2 %in% c(1, 2, 3) & vb3 %in% c(1, 2, 3), sample(bC, 75 - 50)])
stopifnot(length(unique(blocks_not_equalB75)) == 75)
blocks_not_equalB100 <- c(blocks_not_equalB75, bdat_not_equal_nb[!blocks_not_equalB75, ][vb1 %in% c(1, 2, 3) & vb2 %in% c(1, 2, 3) & vb3 %in% c(1, 2, 3), sample(bC, 100 - 75)])
stopifnot(length(unique(blocks_not_equalB100)) == 100)
blocks_not_equalB1000 <- c(blocks_not_equalB100, bdat_not_equal_nb[!blocks_not_equalB100, ][vb1 %in% c(1, 2, 3) & vb2 %in% c(1, 2, 3, 4) & vb3 %in% c(1, 2, 3, 4), sample(bC, 1000 - 100)])
stopifnot(length(unique(blocks_not_equalB1000)) == 1000)


set.seed(12345)
blocks_equalB10 <- bdat_equal_nb[vb1 %in% c(1, 2) & vb2 %in% c(1, 2) & vb3 %in% c(1, 2), sample(bC, size = 10)]
stopifnot(length(unique(blocks_equalB10)) == 10)
## Now make a 25 block dataset by adding to the original 10 block dataset
## Point is to minimize simulation error between datasets
blocks_equalB25samp <- bdat_equal_nb[!blocks_equalB10, ][vb1 %in% c(1, 2) & vb2 %in% c(1, 2) & vb3 %in% c(1, 2), sample(bC, 25 - 10)]
blocks_equalB25 <- c(blocks_equalB10, blocks_equalB25samp)
stopifnot(length(unique(blocks_equalB25)) == 25)
blocks_equalB50 <- c(blocks_equalB25, bdat_equal_nb[!blocks_equalB25, ][vb1 %in% c(1, 2, 3) & vb2 %in% c(1, 2, 3) & vb3 %in% c(1, 2, 3), sample(bC, 50 - 25)])
stopifnot(length(unique(blocks_equalB50)) == 50)
blocks_equalB75 <- c(blocks_equalB50, bdat_equal_nb[!blocks_equalB50, ][vb1 %in% c(1, 2, 3) & vb2 %in% c(1, 2, 3) & vb3 %in% c(1, 2, 3), sample(bC, 75 - 50)])
stopifnot(length(unique(blocks_equalB75)) == 75)
blocks_equalB100 <- c(blocks_equalB75, bdat_equal_nb[!blocks_equalB75, ][vb1 %in% c(1, 2, 3) & vb2 %in% c(1, 2, 3) & vb3 %in% c(1, 2, 3), sample(bC, 100 - 75)])
stopifnot(length(unique(blocks_equalB100)) == 100)
blocks_equalB1000 <- c(blocks_equalB100, bdat_equal_nb[!blocks_equalB100, ][vb1 %in% c(1, 2, 3) & vb2 %in% c(1, 2, 3, 4) & vb3 %in% c(1, 2, 3, 4), sample(bC, 1000 - 100)])
stopifnot(length(unique(blocks_equalB1000)) == 1000)

save(blocks_not_equalB1000, blocks_not_equalB100, blocks_not_equalB10,
  blocks_not_equalB25, blocks_not_equalB50, blocks_not_equalB75,
  blocks_equalB1000, blocks_equalB100, blocks_equalB10, blocks_equalB25,
  blocks_equalB50, blocks_equalB75, file = here::here("Data",
    "blocks_sampled.rda"))
