# Inverse-covariance optimal linear combination of estimators.
# Source: stable_wls_cov_sim.R (lines 107-116), unmodified except
# for removal of two leading spaces of indentation per line.

compute_optimal <- function(covmat, est) {
  ones       <- rep(1, ncol(covmat))
  inv.cov    <- solve(covmat)
  w.opt      <- inv.cov %*% ones
  w.opt      <- as.vector(w.opt / sum(w.opt))
  est.opt    <- sum(w.opt * est)
  var.opt    <- 1 / as.numeric(t(ones) %*% inv.cov %*% ones)
  list(weight = w.opt, estimate = est.opt, var = var.opt)
}



# Simple data-driven combinations of the WLS estimators (the ones reported in
# the paper's XRP application). `covmat` is the (bootstrap) covariance matrix
# of `est`; both should be restricted to the WLS estimators.
#
#   Avg     equal weights
#   InvVar  weights proportional to 1 / variance (diagonal of covmat only)
#   MinVar  the single estimator with the smallest estimated variance
combine_simple <- function(covmat, est) {
  v <- diag(covmat)
  w.inv <- (1 / v) / sum(1 / v)
  list(Avg    = list(estimate = mean(est),
                     var      = sum(covmat) / length(est)^2),
       InvVar = list(estimate = sum(w.inv * est),
                     var      = as.numeric(t(w.inv) %*% covmat %*% w.inv)),
       MinVar = list(estimate = est[which.min(v)],
                     var      = min(v)))
}
