# Core estimation: quantile-based WLS scale estimator and the
# repeated-splitting estimator.

# select_bandwidth(x)
#
# Plug-in bandwidth ks::hpi() for the kernel density estimate of Y. Falls back
# to Silverman's rule of thumb (stats::bw.nrd0) if hpi() errors or returns a
# non-finite / non-positive value (can happen for very heavy tails, e.g.
# alpha = 0.2), and finally to a small fraction of the data range so that
# kde() never receives h = 0.
select_bandwidth <- function(x) {
  h <- tryCatch(ks::hpi(x), error = function(e) NA_real_)
  if (!is.finite(h) || h <= 0) h <- stats::bw.nrd0(x)
  if (!is.finite(h) || h <= 0) {
    spread <- diff(range(x))
    h <- if (spread > 0) spread / 1000 else 1e-6
  }
  h
}

# wls.ests.tvec(X, Y, t0, weights = c("diag", "gls"))
#
# Slope of the weighted regression of the Y-quantiles on the X-quantiles
# (with intercept) at probabilities t0. It estimates sigma = 2^(1/alpha).
#
# weights:
#   "diag" (default) density weighting, W = diag(fY^2). This is the estimator
#                    used in the paper.
#   "gls"            full GLS weights W = D C^-1 D, with D = diag(fY) and
#                    C_ij = min(t_i, t_j) - t_i t_j (Brownian-bridge covariance).
#                    NOTE: in our simulations "gls" generally performs
#                    subpar relative to "diag" (higher RMSE in roughly 90% of
#                    settings, and increasingly so for larger k), so it is
#                    kept only for comparison.
wls.ests.tvec <- function(X, Y, t0, weights = c("diag", "gls")) {

  weights <- match.arg(weights)
  k0 <- length(t0)

  # Sample quantiles of X and Y at the probabilities t0
  qX <- as.numeric(quantile(X, probs = t0))
  qY <- as.numeric(quantile(Y, probs = t0))

  # Density of Y at its quantiles: kernel density estimate on the asinh scale
  # (plug-in bandwidth), back-transformed by the Jacobian 1 / sqrt(1 + y^2).
  aY <- asinh(Y)
  fY <- kde(aY,
            eval.points = asinh(qY),
            h = select_bandwidth(aY))$estimate / sqrt(1 + qY^2)

  if (weights == "diag") {
    # Closed-form weighted slope with weights fY^2. The data are rescaled
    # first (the slope is scale-invariant) so that extremely heavy-tailed
    # samples (alpha = 0.2, |Y| up to ~1e100) do not overflow.
    s <- max(abs(qY))
    if (!is.finite(s) || s <= 0) s <- 1
    w  <- (fY / max(fY))^2
    x  <- qX / s
    y  <- qY / s
    sw <- sum(w)
    xb <- sum(w * x) / sw
    yb <- sum(w * y) / sw
    den <- sum(w * (x - xb)^2)
    if (!is.finite(den) || den <= 0) return(NA_real_)
    return(sum(w * (x - xb) * (y - yb)) / den)
  }

  # "gls": W = D C^-1 D, normal equations solved directly
  C <- outer(t0, t0, FUN = pmin) - outer(t0, t0, FUN = `*`)
  W <- outer(fY, fY, FUN = "*") * solve(C)
  Q <- cbind(1, qX)
  tryCatch(solve(t(Q) %*% W %*% Q, t(Q) %*% W %*% qY)[2],
           error = function(e) NA_real_)
}

run.split <- function(Z, n.splits, p, t_ratio, t_wls, weights = c("diag", "gls")) {
  weights <- match.arg(weights)   # see wls.ests.tvec(); "gls" generally performs subpar
  
  # Number of ratio and wls values to compute
  n_ratio <- length(t_ratio)
  n_wls <- length(t_wls)
  
  # Initialize a matrix to store estimates from each bootstrap iteration.
  # First n_ratio columns for the ratio estimates and next n_wls columns for the wls estimates.
  par.store <- matrix(NA, nrow = n.splits, ncol = n_ratio + n_wls)
  
  n <- length(Z)
  # Pre-generate bootstrap masks for each split in a matrix form.
  B.vec <- matrix(rbinom(n.splits * n, 1, p), nrow = n.splits, ncol = n)
  
  for (b in 1:n.splits) {
    # Extract the b-th bootstrap binary mask.
    B <- B.vec[b, ]
    
    # Subset Z into X (where B==1) and Z0 (where B==0)
    X <- Z[B == 1]
    Z0 <- Z[B == 0]
    
    # Compute Y as the sum of all unique pairwise combinations from Z0 using outer.
    # Only compute the upper triangular part of the outer sum matrix.
    if (length(Z0) >= 2) {
      Y.full <- outer(Z0, Z0, "+")
      Y <- Y.full[upper.tri(Y.full)]
    } else {
      # If there are fewer than 2 elements in Z0, Y is set to NA (or handle as needed).
      Y <- NA
    }
    
    # Compute ratio estimates for each element in t_ratio.
    ratio_estimates <- sapply(1:n_ratio, function(k) {
      diff(quantile(Y, probs = t_ratio[[k]])) / diff(quantile(X, probs = t_ratio[[k]]))
    })
    
    # Compute the weighted least squares estimates for each element in t_wls.
    wls_estimates <- sapply(1:n_wls, function(k) {
      wls.ests.tvec(X, Y, t_wls[[k]], weights = weights)
    })
    
    # Store the combined estimates: first ratios, then wls estimates.
    par.store[b, ] <- c(ratio_estimates, wls_estimates)
  }
  
  # Compute the mean estimates over all bootstrap iterations.
  par.out <- colMeans(par.store, na.rm = TRUE)
  
  # Assign names: ratio estimates as "ratio1", "ratio2", ... and wls estimates as "wls1", "wls2", ...
  names(par.out) <- c(paste0("ratio", 1:n_ratio), paste0("wls", 1:n_wls))
  
  # Return a list with the ratio estimates and wls estimates.
  return(list(ratio = par.out[1:n_ratio],
              wls = par.out[(n_ratio + 1):(n_ratio + n_wls)]))
}
