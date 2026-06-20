# run_all.R — reproduce the full COL14 analysis pipeline.
#
# Usage (with the COL14 repository, or any folder inside it, as the working
# directory):
#   source("run_all.R")     # from an interactive R session
#   Rscript run_all.R       # from a shell
#
# WARNING: this runs every model from scratch — 4 chains x 100,000 iterations
# for the double-exponential growth model and for each proxy model — and can
# take a long time. If you only need to inspect the results, load the
# precomputed objects in results/ instead (see results/README.md).
#
# Requires the `here` package and the packages listed in R/README.md.

library(here)

orig_wd <- getwd()
on.exit(setwd(orig_wd), add = TRUE)   # restore wd even if a script setwd()s

message("COL14 repository root: ", here::here())

# 1. Core demographic model — must run first: it estimates the change point
#    whose 95% HPDI upper bound (3808 cal BP) defines the window the proxy
#    models use.
source(here::here("R", "col14_mcmc_v12.R"))

# 2-3. Bin-level proxy models.
source(here::here("R", "eljunco.R"))
source(here::here("R", "xrfpallcacocha.R"))

message("Done. Fresh outputs are under outputs/.")
