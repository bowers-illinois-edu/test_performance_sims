## Show that all of the different adjustment methods control the FWER
## Compare control in the weak sense (when the sharp null of no effects is true for all units)
## and in the strong sense (when the sharp null of no effects is true for some but not all units)

library(tidyverse)
library(here)
library(data.table)
library(dtplyr)
library(dplyr)
library(ggthemes)
library(RColorBrewer)
library(manytestsr)
library(xtable)

load(here::here("figures", "simsres_plotting.rda"))

names(simsres)
head(simsres)
nrow(simsres)

## Standard error of the simulations
nsims <- 1000
std_err_of_sim <- 2 * sqrt(.05 * (1 - .05) / nsims)

table(simsres$prop_blocks_0, exclude = c())
table(simsres$nblocks, exclude = c())

## Look at amount of time spent per sim on average

sims_time_summary <- simsres %>% group_by(nblocks,prop_blocks_0,splitfn,afn) %>%
  reframe(mintime=min(elapsed),
          pct25time=quantile(elapsed,.25),
          meantime=mean(elapsed),
          mediantime=median(elapsed),
          pct75time=quantile(elapsed,.75),
          maxtime=max(elapsed)) %>% 
as.data.table()

sims_time_summary0 <- simsres %>% group_by(nblocks) %>%
  reframe(mintime=min(elapsed),
          pct25time=quantile(elapsed,.25),
          meantime=mean(elapsed),
          mediantime=median(elapsed),
          pct75time=quantile(elapsed,.75),
          maxtime=max(elapsed)) %>% 
as.data.table()

sims_time_summary %>% arrange(mediantime)

time_plot <- ggplot(sims_time_summary,aes(x=nblocks,y=maxtime,color=prop_blocks_0))+
  geom_point() +
  facet_wrap(~splitfn+afn,scales="free")
time_plot



## Look at the simulation results in two types Weak and Strong:
## FDR control is also shown

### The sharp null of no effects is true for all units (and all blocks)
simsres_alltrue <- droplevels(simsres[prop_blocks_0 == 1 & nblocks < 1000, ])

## Here we are looking for strong control of the FWER: when the null is true, reject rarely.
## True nulls will be mixed with false nulls (i.e. in these sims, effects=0 versus effects not 0)
### At least one block as a true effect and at least has no effect (i.e. a "true" null)
simsres_sometrue <- droplevels(simsres[(prop_blocks_0 < 1 & prop_blocks_0 > 0) & nblocks < 1000, ])


## See: https://copyprogramming.com/howto/why-is-controlling-fdr-less-stringent-than-controlling-fwer
## Notice that all approaches control the FWER in a WEAK form (when the
## sharp null is true)

## "Controlling the familywise error rate in the weak sense then implies the
## control of the probability of declaring an effect for at least one outcome
## variable, when there is no effct on any variable." (page 14, Bretz et al
## 2011)

## I think that controlling the FWER in a strong sense means that no more than
## .05  of the true hypotheses are rejected. (This is not the same as . 05 of
## the rejected hypotheses are falsely rejected --- this is the fdr).

## "If, for a given multiple comparison procedure, the Type I error rate is
## controlled under any partial configuration of true and false null
## hypotheses, the error control is denoted as strong." (page 14, Bretz et al
## 2011): maximum p for all tests of true null hypotheses "In the previous
## example, controlling the familywise error rate in the strong sense implies
## the control of the probability of declaring an effect for at least one
## outcome variable, regardless of the effect sizes for any of the outcome
## variables." (page 14--15, bretz et al 2011)

## Page 25 of Bretz cites Marcus  et al 1976 as saying that the closure
## principle can be proven to control the FWER in the strong sense."For the
## final inference, an elementary null hypothesis Hi is rejected if and only if
## all intersection hypotheses implying Hi are rejected by their individual
## tests at local level ↵, too. It can be shown that the above procedure
## controls the familywise error rate strongly at level ↵ (Marcus et al. 1976).
## "

## I am not sure whether we can guarantee strong FWER control always. So this is more exploratory here.

############## WEAK CONTROL of the FWER ########

### Reminder: the false positive proportion is the proportion of discoveries
### (whether of a single block or a group of blocks) occurring when all of the
### effects are 0. The false discovery proportion is the proportion of
### rejections of the true null out of all of the rejections.

stopifnot(all(simsres_alltrue$false_pos_prop < (.05 + std_err_of_sim)))
summary(simsres_alltrue$false_pos_prop)

## make a table for the appendix (maybe. This is a lot actually.)
simsres_alltrue %>%
  arrange(splitfnF) %>%
  select(splitfnF, nblocks, pfn, afn, false_pos_prop, false_disc_prop)

## So maybe collapse across nblocks
weak_control_tab_dt <- simsres_alltrue %>%
  group_by(splitfnF, afn) %>%
  droplevels() %>%
  summarize(
    max_fwer = max(false_pos_prop),
    max_fdr = max(false_disc_prop)
  ) %>%
  ungroup() %>%
  arrange(splitfnF, afn)

# %>% print(n = 100)

weak_control_tab_dt$afn <- gsub("_", " ", weak_control_tab_dt$afn)

print(weak_control_tab_dt, n = 100)

load(here("Data", "equal_blocksize_data.rda"), verbose = TRUE)
load(here("Analysis", "simparms.rda"), verbose = TRUE)
table(table(idat_equal_B100$bF))
table(table(idat_equal_B10$bF))

weak_control_tab <- xtable(weak_control_tab_dt,
  caption = "All approaches control FWER and FDR when the sharp null of no
  effects is true. Entries are the maximum FWER or FDR over simulations
  involving between 10 (N=200) and 100 (N=2000) blocks with all blocks
having 20 units and half of the units assigned
to treatment at random. In this simulation the null
  hypothesis of no effects is true for all units and all blocks ($\\tau_{i,b} =
  0$ for all $i$ and $b$). Twice the simulation error for 1000 simulations is
  roughly .02.",
  label = "tab:weak_fwer_control",
  auto = TRUE,
  # align = "llp{6ex}cccp{4ex}p{4ex}p{4ex}p{4ex}p{4ex}p{4ex}cccc",
  digits = 3
)

colnames(weak_control_tab) <- c("", "", "Max FWER", "Max FDR")

print(weak_control_tab,
  type = "latex",
  include.rownames = FALSE,
  include.colnames = TRUE,
  # add.to.row = add_to_row,
  #  sanitize.text.function = force,
  booktabs = TRUE,
  floating = TRUE,
  table.placement = "H",
  file = here("figures", "sim_fwer_weak.tex")
)
# End

############# STRONG CONTROL OF THE FWER ############

## From
## https://en.wikipedia.org/wiki/Family-wise_error_rate#Controlling_procedures
## A procedure controls the FWER in the weak sense if the FWER control at level � {\displaystyle \alpha \,\!} is guaranteed only when all null hypotheses are true (i.e. when � 0 = � {\displaystyle m_{0}=m}, meaning the "global null hypothesis" is true)

## A procedure controls the FWER in the strong sense if the FWER control at
## level  {\displaystyle \alpha \,\!} is guaranteed for any configuration of true and
## non-true null hypotheses (whether the global null hypothesis is true or not)

## So, we can calculate strong control of the null when we have some blocks
## with a true sharp null of no effects and some blocks with true sharp null of
## some effects.

## Mostly doesn't reject true null but sometimes does
summary(simsres_sometrue$false_pos_prop)

not_strong_ctrl <- simsres_sometrue %>%
  filter(false_pos_prop > (.05 + 2 * std_err_of_sim)) %>%
  droplevels() %>%
  select(
    splitfn, nblocks, afn, false_pos_prop, false_disc_prop, true_disc_prop, tau_sizes, prop_blocks_0, nreject, naccept,
    prop_reject, prop_accept, tot_true0, tot_truenot0, tot_reject_true0, tot_not_reject_true0, tot_reject_truenot0,
    correct_pos_effect_prop, correct_pos_nulls_prop, prop_not_reject_true0
  ) %>%
  arrange(nblocks, tau_sizes)


## The issue is with the leave one out splitter.
with(not_strong_ctrl, table(splitfn))
## And this is for the approach to FDR not FWER anyway. So, splitLOO has strong control as far as we can see.
with(not_strong_ctrl, table(afn, exclude = c()))
## Problem gets worse as we get more blocks (all of the issues with nblocks=10 are within 2 * SE of Simulation)
with(not_strong_ctrl, table(nblocks))
with(not_strong_ctrl, table(tau_sizes))

## Interpreting the results:
not_strong_ctrl


## This from calc_errs()
##  deterrs <- detobj[, .(
##    nreject = sum(hitb),
##    naccept = sum(1 - hitb),
##    prop_reject = mean(hitb),
##    prop_accept = mean(1 - hitb),
##    tot_true0 = sum(true0),
##    tot_truenot0 = sum(truenot0),
##    tot_reject_true0 = sum(hitb * true0),
##    tot_not_reject_true0 = sum((1 - hitb) * true0),
##    tot_reject_truenot0 = sum(hitb * truenot0),
##    tot_not_reject_truenot0 = sum((1 - hitb) * truenot0),
##    # or 1-prop_reject
##    # Proportion of the total blocks that have an effect where we detect that effect:
##    correct_pos_effect_prop = sum(hitb * truenot0) / max(1, sum(truenot0)),
##    correct_pos_nulls_prop = sum(hitb * true0) / max(1, sum(true0)),
##    true_pos_prop = mean(hitb * truenot0),
##    ## Proportion of rejections of a true null out of the total number of tests
##    false_pos_prop = mean(hitb * true0),
##    # If we do not reject and all blocks are truly null, then we have no error.
##    prop_not_reject_true0 = mean((1 - hitb) * true0),
##    # If we do not reject/detect and at least one of the blocks actually has an effect, we have
##    # a false negative error --- a failure to detect the truth: proportion accept when true is not 0
##    false_neg_prop = mean((1 - hitb) * truenot0),
##    # Now look at false and true discoveries: false rejections as a proportion of rejections
##    false_disc_prop = sum(hitb * true0) / max(1, sum(hitb)),
##    true_disc_prop = sum(hitb * truenot0) / max(1, sum(hitb)),
##    # Failure to reject when the null is not 0 out of the total failures to reject
##    false_nondisc_prop = sum((1 - hitb) * truenot0) / max(1, sum(1 - hitb)),
##    ## Failure to reject when the null is true (we want this) out of total accepts/failures to reject
##    true_nondisc_prop = sum((1 - hitb) * true0) / max(1, sum(1 - hitb))
##    # meangrpsize = mean(grpsize),
##    # medgrpsize = median(grpsize)
##  )]

strong_control_tab_dt <- simsres_sometrue %>%
  group_by(splitfnF, afn) %>%
  droplevels() %>%
  summarize(
    max_fwer = max(false_pos_prop),
    max_fdr = max(false_disc_prop)
  ) %>%
  ungroup() %>%
  arrange(splitfnF, afn)

print(strong_control_tab_dt, n = 100)

strong_control_tab_dt$afn <- gsub("_", " ", strong_control_tab_dt$afn)

## We see strong control almost always --- we expect it by proof for the fixed
## alpha and bottom up methods
strong_control_tab_dt %>% filter(afn == "fixed (a=.05)")
## We include FDR control here but we are not actually aiming to control it.
## Interesting that the structured alpha investing style approaches to FDR
## control for the closed testing procedures do a better job of FWER control
## than FDR control. Do they show a loss in power?

## Also interesting the the split one out has poor FWER and FDR control given
## the alpha investing procedures. (Again, we don't have big expectations about
## this. We present it to spur future work/investigation.)

source(here("Analysis", "load_subsetdata.R"))

nb_B10 <- table(idatB10$blockF)
nb_B100 <- table(idatB100$blockF)
nb_range <- range(c(nb_B10, nb_B100))

pb_B10 <- idatB10[, .(pb = mean(trt)), by = blockF]
pb_B100 <- idatB100[, .(pb = mean(trt)), by = blockF]
summary(pb_B10)
summary(pb_B100)

strong_control_tab <- xtable(strong_control_tab_dt,
  caption = "Control of the FWER and FDR when at least some blocks have
  $\\bar{\\tau}_b \\ne 0 $ and some blocks have $\\bar{\\tau}_b = 0$. Fixed
  $\\alpha$ splitting rules control the FWER in this strong sense. Maximum FWER
  and FDR taken over simulations with between 10 blocks (N=568) and 100 blocks
  (N=13065) with block sizes that vary from 5 to 3147 units and between .4 and
  .6 of the units within the blocks assigned to treatment at random.  Twice the
  simulation error for 1000 simulations is roughly .02.",
  label = "tab:strong_fwer_control",
  auto = TRUE,
  # align = "llp{6ex}cccp{4ex}p{4ex}p{4ex}p{4ex}p{4ex}p{4ex}cccc",
  digits = 3
)

colnames(strong_control_tab) <- c("", "", "Max FWER", "Max FDR")

print.xtable(strong_control_tab,
  type = "latex",
  include.rownames = FALSE,
  include.colnames = TRUE,
  # add.to.row = add_to_row,
  sanitize.text.function = function(x) {
    x
  },
  booktabs = TRUE,
  floating = TRUE,
  table.placement = "H",
  file = here("figures", "sim_fwer_strong.tex")
)
