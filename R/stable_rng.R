# rstable_new(): stable random variate generation with corrected
# handling of the alpha = 1 case.
# Source: split_sample_functions.R (lines 7-48), unmodified.

rstable_new <- function(n, alpha, beta, gamma = 1, delta = 0, pm = 0) {
  # Base function comes from library(stabledist)
  # but did not adequately adjust for alpha = 1
  stopifnot(0 < alpha, alpha <= 2, length(alpha) == 1,
            -1 <= beta, beta <= 1, length(beta) == 1,
            0 <= gamma, length(pm) == 1, pm %in% 0:2)
  
  if (pm == 1) {
    delta <- delta + beta * gamma * .om(gamma, alpha)
  } else if (pm == 2) {
    delta <- delta - alpha^(-1/alpha) * gamma * stableMode(alpha, beta)
    gamma <- alpha^(-1/alpha) * gamma
  }
  
  # Generate common random variables
  theta <- pi * (runif(n) - 1/2)  # U ~ Uniform(-pi/2, pi/2)
  w <- -log(runif(n))             # W ~ Exponential(1)
  
  result <- if (alpha == 1) {
    if (beta == 0) {
      # Special case: symmetric Cauchy
      rcauchy(n)
    } else {
      # Use the standard limit formula for alpha == 1 with beta != 0:
      # X = (2/pi) * [ (pi/2 + beta*theta)*tan(theta) - beta * log((pi/2*w*cos(theta))/(pi/2+beta*theta)) ]
      (2/pi) * ((pi/2 + beta * theta) * tan(theta) -
                  beta * log((pi/2 * w * cos(theta))/(pi/2 + beta * theta)))
    }
  } else {
    # Standard CMS (or Weron-type) algorithm for alpha != 1
    pi2 <- pi / 2
    b.tan.pa <- beta * tan(pi2 * alpha)
    theta0 <- min(max(-pi2, atan(b.tan.pa) / alpha), pi2)
    c <- (1 + b.tan.pa^2)^(1/(2 * alpha))
    a.tht <- alpha * (theta + theta0)
    r <- (c * sin(a.tht) / (cos(theta))^(1/alpha)) *
      (cos(theta - a.tht) / w)^((1 - alpha)/alpha)
    r - b.tan.pa
  }
  
  result * gamma + delta
}
