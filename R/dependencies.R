# Package dependencies for the sample-splitting estimator.
#
# Originally taken from the header of split_sample_functions.R. Matrix and
# kernelboot have been dropped: nearPD(), vech()/invvech() and ruvk() were
# used only by split.boot(), which is no longer part of this repository.
#
#   ks           kde() in wls.ests.tvec()
#   stabledist   stableMode() / .om(), used by rstable_new() when pm != 0
#   StableEstim  comparison estimators in scripts/02 (RUN_COMPARISONS)

library(stabledist)
library(StableEstim)
library(ks)
