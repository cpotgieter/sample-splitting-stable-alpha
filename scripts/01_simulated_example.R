# 01_simulated_example.R
#
# A complete, self-contained walkthrough of the method on simulated data,
# where the true tail index is known:
#
#   1. generate stable data with known alpha / sigma
#   2. set the tuning vectors
#   3. estimate sigma by repeated splitting, transform to alpha if desired
#   4. bootstrap standard errors by the reduced-split procedure
#   5. combine the estimators using the inverse-covariance optimal weights
#
# Run from the repository root:
#
#   source("scripts/01_simulated_example.R")
#

source("load_all.R")

set.seed(20260818)


#  1. Generate data 

alpha_true <- 0.40
beta_true  <- 0.75
n          <- 500

sigma_true <- 2^(1 / alpha_true)   # the population value the estimator targets

Z <- rstable_new(n, alpha = alpha_true, beta = beta_true)

#  2. Tuning vectors 

# One ratio estimator at the quartiles, and three WLS estimators using
# k = 5, 10, 16 interior quantile levels. 

t_ratio <- list(c(1/4, 3/4)) # IQR-based estimator (very robust)

t_wls <- lapply(c(5, 10, 16), function(m)
  round(seq(1/(m + 1), m/(m + 1), length.out = m), 3)
)

est.names <- c("iqr_ratio", "wls05", "wls10", "wls16")


#  3. Point estimation 

B.full <- 250   # number of splits averaged over

sigma.hat <- run.split(Z,
                       n.splits = B.full,
                       p        = 0.5,
                       t_ratio  = t_ratio,
                       t_wls    = t_wls)

sigma.hat <- unlist(sigma.hat)
alpha.hat        <- log(2) / log(sigma.hat)
names(sigma.hat) <- est.names
names(alpha.hat) <- est.names

print(round(cbind(
  true = c(sigma_true, alpha_true),
  rbind(
    sigma.hat = sigma.hat,
    alpha.hat = alpha.hat
  )
), 4))


#  4. Bootstrap standard errors 

# The reduced-split idea: running B = 250 splits inside every bootstrap
# replicate is prohibitive. Instead run only max.split = 25, record the
# cumulative estimate after 1,...,25 splits, and extrapolate the covariance
# out to B = 250 using the Sigma/Gamma decomposition.

n.boot    <- 200
max.split <- 25

## Reconstructed covariance of sigma.hat:
boot.out.sigma <- reduced.split.se(Z,
                                   n.boot    = n.boot,
                                   max.split = max.split,
                                   B.full    = B.full,
                                   p         = 0.5,
                                   t_ratio   = t_ratio,
                                   t_wls     = t_wls,
                                   alpha.scale = FALSE) # If FALSE, returns sigma scale

se.hat.sigma        <- boot.out.sigma$se
names(se.hat.sigma) <- est.names

# Estimated standard errors
print(se.hat.sigma)

## Reconstructed covariance of alpha.hat:
boot.out.alpha <- reduced.split.se(Z,
                                   n.boot    = n.boot,
                                   max.split = max.split,
                                   B.full    = B.full,
                                   p         = 0.5,
                                   t_ratio   = t_ratio,
                                   t_wls     = t_wls,
                                   alpha.scale = TRUE) # If TRUE, returns alpha scale

se.hat.alpha        <- boot.out.alpha$se
names(se.hat.alpha) <- est.names

# Estimated standard errors
print(se.hat.alpha)

#  5. Optimal combination (illustrated for alpha scale) 

opt.all <- compute_optimal(boot.out.alpha$cov, alpha.hat)
opt.wls <- compute_optimal(boot.out.alpha$cov[2:4, 2:4], alpha.hat[2:4])


#  6. Results 

results <- data.frame(
  Estimator = c(est.names, "optimal (all)", "optimal (WLS only)"),
  Estimate  = c(alpha.hat, opt.all$estimate, opt.wls$estimate),
  SE        = c(se.hat.alpha, sqrt(opt.all$var), sqrt(opt.wls$var)),
  stringsAsFactors = FALSE
)

print(format(results, digits = 4), row.names = FALSE)


#  7. How the variance falls off with B -

# Var(alpha.hat_B) = Sigma/B + (B-1)/B * Gamma. The Gamma term is the floor
# that cannot be removed by taking more splits.

Bgrid   <- 1:B.full
var.wls10 <- boot.out.alpha$Sigma[3, 3] / Bgrid +
             (Bgrid - 1) / Bgrid * boot.out.alpha$Gamma[3, 3]

plot(Bgrid, sqrt(var.wls10), type = "l", log = "x",
     xlab = "Number of splits (B)", ylab = "SE of alpha.hat",
     main = "Reduced-split variance extrapolation (WLS_10)")
abline(h = sqrt(boot.out.alpha$Gamma[3, 3]), lty = 3, col = "grey40")
abline(v = max.split, lty = 2, col = "grey40")

cat("\nSE floor as B -> Inf (WLS_10):",
    round(sqrt(boot.out.alpha$Gamma[3, 3]), 4), "\n")
