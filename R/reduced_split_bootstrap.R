# ---------------------------------------------------------------------------
# Reduced-split bootstrap inference.
#
# NOTE ON PROVENANCE: unlike the other files under R/, this file is NEW code.
# It was not extracted from the original scripts. It packages the covariance
# reconstruction that appears INLINE in XRE_analysis.R (lines 103-140) and in
# stable_wls_cov_sim.R (lines 51-100) so that both scripts in scripts/ can
# share one implementation.
#
# The arithmetic is unchanged from those inline blocks:
#
#   G = (B2*C2 - B1*C1) / (B2 - B1)
#   S = B1*C1 - (B1-1)*G
#
# averaged over all pairs B1 < B2, then combined as
#
#   cov = (1/B.full)*Sigma_avg + ((B.full-1)/B.full)*Gamma_avg
#
# This is the procedure that produced the standard errors reported in the
# paper.
# ---------------------------------------------------------------------------


reduced.split.cov <- function(boot.est.splits, B.full, B.min = 2) {
  # boot.est.splits : array of dim c(n.boot, n.est, max.split), containing the
  #                   cumulative estimates after 1,...,max.split splits for each
  #                   bootstrap replicate, already on the scale of interest.
  # B.full          : the number of splits used for the point estimate.
  # B.min           : smallest B1 to use. Defaults to 2; B = 1 is excluded
  #                   because the reconstruction is unstable there.

  max.split <- dim(boot.est.splits)[3]
  n.est     <- dim(boot.est.splits)[2]

  stopifnot(max.split >= B.min + 1, B.full >= 1)

  Gamma_sum <- Sigma_sum <- matrix(0, n.est, n.est)
  n.pairs   <- 0

  for (B1 in B.min:(max.split - 1)) {
    for (B2 in (B1 + 1):max.split) {
      C1 <- cov(boot.est.splits[, , B1])
      C2 <- cov(boot.est.splits[, , B2])
      G  <- (B2 * C2 - B1 * C1) / (B2 - B1)
      S  <- B1 * C1 - (B1 - 1) * G
      Gamma_sum <- Gamma_sum + G
      Sigma_sum <- Sigma_sum + S
      n.pairs   <- n.pairs + 1
    }
  }

  Gamma_avg <- Gamma_sum / n.pairs
  Sigma_avg <- Sigma_sum / n.pairs

  cov.reconstructed <- (1 / B.full) * Sigma_avg +
                       (B.full - 1) / B.full * Gamma_avg

  list(cov     = cov.reconstructed,
       Sigma   = Sigma_avg,
       Gamma   = Gamma_avg,
       n.pairs = n.pairs,
       B.full  = B.full)
}


reduced.split.se <- function(Z, n.boot, max.split, B.full, p,
                             t_ratio, t_wls,
                             alpha.scale = TRUE,
                             B.min = 2,
                             verbose = TRUE,
                             pd.tol = 1e-8) {
  # Full inference pipeline:
  #   resample Z -> run.split.mult() -> (optional) alpha transform
  #   -> reduced.split.cov() -> standard errors.
  #
  # alpha.scale = TRUE applies log(2)/log(.) to the cumulative sigma estimates
  # BEFORE computing covariances, so the returned covariance is on the alpha
  # scale. This matches XRE_analysis.R.
  
  n     <- length(Z)
  n.est <- length(t_ratio) + length(t_wls)
  
  boot.est.splits <- array(NA_real_, dim = c(n.boot, n.est, max.split))
  
  for (b in seq_len(n.boot)) {
    Z.boot <- sample(Z, size = n, replace = TRUE)
    
    bsplits <- run.split.mult(Z.boot,
                              n.splits = max.split,
                              p        = p,
                              t_ratio  = t_ratio,
                              t_wls    = t_wls)
    
    for (k in seq_len(max.split)) {
      boot.est.splits[b, , k] <- bsplits[[k]][seq_len(n.est)]
    }
    
    if (verbose && (b %% 25 == 0 || b == n.boot)) {
      cat(sprintf("\r  bootstrap replicate %d / %d", b, n.boot))
      utils::flush.console()
    }
  }
  if (verbose) cat("\n")
  
  if (alpha.scale) {
    boot.est.splits <- log(2) / log(boot.est.splits)
  }
  
  n.bad <- sum(!is.finite(boot.est.splits))
  if (n.bad > 0) {
    warning(sprintf(
      "%d of %d bootstrap split estimates were non-finite (%.2f%%). %s",
      n.bad, length(boot.est.splits),
      100 * n.bad / length(boot.est.splits),
      "These propagate as NA into cov(); inspect before trusting the SEs."))
  }
  
  out <- reduced.split.cov(boot.est.splits, B.full = B.full, B.min = B.min)
  
  # Check whether the reconstructed covariance matrix is positive definite.
  # If needed, repair it by flooring eigenvalues at a scale-aware threshold.
  eig <- eigen(out$cov, symmetric = TRUE)
  
  scale <- mean(diag(out$cov))
  eig.floor <- pd.tol * scale
  
  if (min(eig$values) <= eig.floor) {
    eig.values.fixed <- pmax(eig$values, eig.floor)
    
    cov.fixed <- eig$vectors %*%
      diag(eig.values.fixed, nrow = length(eig.values.fixed)) %*%
      t(eig$vectors)
    
    dimnames(cov.fixed) <- dimnames(out$cov)
    
    if (verbose) {
      cat(sprintf(
        "  covariance repair applied: minimum eigenvalue %.3e -> %.3e\n",
        min(eig$values), eig.floor
      ))
    }
    
    out$cov.raw <- out$cov
    out$cov <- cov.fixed
    out$cov.repaired <- TRUE
  } else {
    out$cov.repaired <- FALSE
  }
  
  out$se              <- sqrt(diag(out$cov))
  out$boot.est.splits <- boot.est.splits
  out$n               <- n
  out$n.boot          <- n.boot
  out$alpha.scale     <- alpha.scale
  
  out
}
