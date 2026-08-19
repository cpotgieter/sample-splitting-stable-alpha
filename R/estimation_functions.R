# Core estimation: quantile-based WLS scale estimator and the
# repeated-splitting estimator.
# Source: split_sample_functions.R (lines 51-150), unmodified.

wls.ests.tvec <- function(X, Y, t0) {
  
  k0 <- length(t0)
  
  # Construct a k0-by-k0 matrix where each row is the vector t
  # t.mat <- matrix(rep(t0, k0), nrow = k0)
  
  # Compute the covariance structure (Sigma matrix)
  # Sigma <- pmin(t.mat, t(t.mat)) - t.mat * t(t.mat)
  Sigma <- outer(t0, t0, FUN = pmin) - outer(t0, t0, FUN = `*`)
  
  # Compute sample quantiles for X and Y at the probabilities in t
  qX <- as.numeric(quantile(X, probs = t0))
  qY <- as.numeric(quantile(Y, probs = t0))
  
  # Estimate density fY at the transformed quantiles using kernel density estimation.
  # The transformation is applied as log(y + sqrt(y^2 + 1)) and then back-transformed.
  fY <- kde(asinh(Y), 
            eval.points = asinh(qY),
            h = bw.nrd(asinh(Y)))$estimate / sqrt(1 + qY^2)
  
  # Create a k0-by-k0 matrix where each column is the density fY vector
  #fY.mat <- matrix(rep(fY, k0), nrow = k0)
  
  # Compute the weight matrix: product of densities divided by the corresponding Sigma element.
  # w <- fY.mat * t(fY.mat) / Sigma
  w <- outer(fY, fY, FUN = "*") / Sigma
  
  # Create a design matrix for regression with an intercept and qX as predictor.
  qX <- matrix(c(rep(1, k0), qX), nrow = k0, byrow = FALSE)
  
  # Compute the weighted least squares estimator using the closed form solution:
  #wlm <- solve(t(qX) %*% w %*% qX) %*% (t(qX) %*% w %*% qY)
  
  A <- t(qX) %*% w %*% qX
  b <- t(qX) %*% w %*% qY
  
  # Solve the linear system A * wlm = b directly.
  wlm <- solve(A, b)
  
  # Return the slope (second coefficient) as the estimator
  return(wlm[2])
}

run.split <- function(Z, n.splits, p, t_ratio, t_wls) {
  
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
      wls.ests.tvec(X, Y, t_wls[[k]])
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
