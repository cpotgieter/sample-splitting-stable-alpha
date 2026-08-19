# load_all.R

# Sources the methodology functions in dependency order.
#
# Usage, from the repository root:
#
#   source("load_all.R")
#
# Makes available:
#
#   wls.ests.tvec()     two-sample quantile-based WLS scale estimator
#   run.split()         main repeated-splitting estimator
#   run.split.mult()    cumulative estimates for B = 1, ..., n.splits
#   reduced.split.cov() (B1,B2) covariance reconstruction
#   reduced.split.se()  full bootstrap SE pipeline
#   compute_optimal()   inverse-covariance optimal combination
#   rstable_new()       stable RNG (needed only for simulated examples)

.sp_dir <- if (exists(".sp_dir")) .sp_dir else "."

source(file.path(.sp_dir, "R", "dependencies.R"))
source(file.path(.sp_dir, "R", "stable_rng.R"))
source(file.path(.sp_dir, "R", "estimation_functions.R"))
source(file.path(.sp_dir, "R", "bootstrap_functions.R"))
source(file.path(.sp_dir, "R", "reduced_split_bootstrap.R"))
source(file.path(.sp_dir, "R", "combination_functions.R"))

rm(.sp_dir)

message("Sample-splitting functions loaded.")
