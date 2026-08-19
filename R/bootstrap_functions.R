# Cumulative split estimates: same repeated-splitting procedure as
# run.split(), but returning the running mean after 1, ..., n.splits
# splits. This is the input to the reduced-split bootstrap.
# Source: split_sample_functions.R (lines 152-202), unmodified.

run.split.mult <- function(Z, n.splits, p, t_ratio, t_wls) {
  # Number of ratio and wls values to compute
  n_ratio <- length(t_ratio)
  n_wls <- length(t_wls)
  
  # Initialize a matrix to store estimates from each bootstrap iteration.
  # First n_ratio columns for the ratio estimates and next n_wls columns for the wls estimates.
  par.store <- matrix(NA, nrow = n.splits, ncol = n_ratio + n_wls)
  
  n <- length(Z)
  # Pre-generate bootstrap masks for all iterations.
  B.vec <- matrix(rbinom(n.splits * n, 1, p), nrow = n.splits, ncol = n)
  
  for (b in 1:n.splits) {
    # Extract the b-th bootstrap mask.
    B <- B.vec[b, ]
    
    # Subset Z into X (where B == 1) and Z0 (where B == 0)
    X <- Z[B == 1]
    Z0 <- Z[B == 0]
    
    # Compute Y as the sum of all unique pairwise combinations from Z0 using outer.
    if (length(Z0) >= 2) {
      Y.full <- outer(Z0, Z0, "+")
      Y <- Y.full[upper.tri(Y.full)]
    } else {
      # If there are fewer than 2 elements in Z0, set Y to NA (or handle as needed).
      Y <- NA
    }
    
    # Compute ratio estimates for each element in t_ratio.
    ratio_estimates <- sapply(1:n_ratio, function(k) {
      diff(quantile(Y, probs = t_ratio[[k]])) / diff(quantile(X, probs = t_ratio[[k]]))
    })
    
    # Compute the weighted least squares estimator for each element in t_wls.
    wls_estimates <- sapply(1:n_wls, function(k) {
      wls.ests.tvec(X, Y, t_wls[[k]])
    })
    
    # Store the combined estimates: first ratios, then wls estimates.
    par.store[b, ] <- c(ratio_estimates, wls_estimates)
  }
  
  # Compute the cumulative mean estimates over the first 1 to n.splits rounds.
  par.out <- lapply(seq_len(n.splits), function(k) {
    colMeans(par.store[1:k, , drop = FALSE], na.rm = TRUE)
  })
  
  return(par.out)
}
