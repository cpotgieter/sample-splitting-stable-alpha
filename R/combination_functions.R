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

