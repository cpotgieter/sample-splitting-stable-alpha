# ===========================================================================
# 02_xrp_application.R
#
# Applies the sample-splitting estimator to XRP (Ripple) daily log-returns,
# Aug 2013 - May 2025, split into an early and a late subperiod.
#
# Run from the repository root:
#
#   source("scripts/02_xrp_application.R")
#
# RUNTIME WARNING. The cost of one split is driven by the pairwise-sum
# construction, which is O(n^2): each subperiod has n ~ 1437, so every split
# builds roughly 258,000 pairwise sums and runs a KDE over them. At the paper
# settings (n.boot = 500, max.split = 25) this script takes on the order of an
# hour per subperiod. Set QUICK <- TRUE for a ~5 minute sanity run.
# ===========================================================================

source("load_all.R")

QUICK <- TRUE     # < lower computational cost version

if (QUICK) {
  B.full <- 100; n.boot <- 50;  max.split <- 15
} else {
  B.full <- 500; n.boot <- 500; max.split <- 25
}

set.seed(20250818)


#  1. Import and prepare the data -

csv.path <- file.path("data", "ripple_2013-08-05_2025-05-27.csv")
stopifnot(file.exists(csv.path))

# fileEncoding handles the UTF-8 byte-order mark on the header line.
data_raw <- read.csv(csv.path, header = TRUE, fileEncoding = "UTF-8-BOM")

# The file is ordered most-recent-first, so reverse before differencing.
n_raw <- nrow(data_raw)
date  <- as.Date(data_raw$End[n_raw:1])
price <- data_raw$Close[n_raw:1]

data <- diff(log(price))   # daily log-returns
date <- date[-1]           # drop first date to match
n    <- length(data)


cat("\n XRP log-returns \n")
cat(sprintf("%d returns, %s to %s\n",
            n, format(min(date)), format(max(date))))

# Outer thirds: an early and a late regime, symmetric by construction.
one_third   <- floor(n / 3)
left_index  <- one_third
right_index <- n - one_third + 1

data_old <- data[1:left_index]
data_new <- data[right_index:n]

periods <- list(
  old = list(x = data_old,
             label = sprintf("%s to %s",
                             format(date[1], "%b %Y"),
                             format(date[left_index], "%b %Y"))),
  new = list(x = data_new,
             label = sprintf("%s to %s",
                             format(date[right_index], "%b %Y"),
                             format(date[n], "%b %Y")))
)

cat(sprintf("early period: n = %d  (%s)\n", length(data_old), periods$old$label))
cat(sprintf("late  period: n = %d  (%s)\n", length(data_new), periods$new$label))


#  2. Series plot with the two subperiods shaded 

op <- par(mfrow = c(1, 1), mar = c(5, 4.5, 2, 1))
plot(date, data, type = "n", xlab = "", ylab = "XRP log-return",
     xaxt = "n", xaxs = "i")
usr <- par("usr")
rect(min(date), usr[3], date[left_index], usr[4], col = "lightgray", border = NA)
rect(date[right_index], usr[3], max(date), usr[4], col = "lightgray", border = NA)
lines(date, data)
abline(v = date[c(left_index, right_index)], lty = 3, lwd = 1.2)
box()
ticks  <- date[c(1, left_index, right_index, n)]
axis(1, at = ticks, labels = FALSE)
text(x = ticks, y = usr[3] - 0.05 * diff(usr[3:4]),
     labels = format(ticks, "%b %e\n%Y"),
     srt = 45, adj = 1, xpd = TRUE, cex = 0.8)
par(op)


#  3. Stability diagnostic: Q-Q of returns vs self-convolution 

op <- par(mfrow = c(1, 2), mar = c(3, 3, 2, 1), oma = c(2.5, 2.5, 1.5, 1))
for (nm in names(periods)) {
  x1 <- periods[[nm]]$x
  x2 <- apply(combn(x1, 2), 2, sum)
  qqplot(x1, x2, pch = 19, cex = 0.5, cex.main = 0.9,
         main = periods[[nm]]$label, xlab = "", ylab = "")
  pp <- seq(0.01, 0.99, length.out = 50)
  abline(lm(quantile(x2, pp) ~ quantile(x1, pp)), lty = 2)
}
mtext("XRP log-returns", side = 1, outer = TRUE, line = 1.2)
mtext("Empirical Self-Convolution", side = 2, outer = TRUE, line = 1.2)
par(op)


#  4. Tuning vectors 

t_ratio <- list(c(1/4, 3/4))

t_wls <- lapply(c(5, 10, 16), function(m)
  round(seq(1/(m + 1), m/(m + 1), length.out = m), 3)
)

est.names <- c("iqr_ratio", "wls05", "wls10", "wls16")


#  5. Estimation and inference for each subperiod 

analyse <- function(x, label) {

  cat("\n===============================================================\n")
  cat("  ", label, "   (n = ", length(x), ")\n", sep = "")
  cat("===============================================================\n")

  # Point estimate
  cat("Estimating with B =", B.full, "splits ...\n")
  sigma.hat        <- run.split(x, n.splits = B.full, p = 0.5,
                                t_ratio = t_ratio, t_wls = t_wls)
  alpha.hat        <- log(2) / log(unlist(sigma.hat))
  names(alpha.hat) <- est.names

  # Reduced-split bootstrap standard errors
  cat("Bootstrapping:", n.boot, "replicates x", max.split, "splits ...\n")
  boot.out <- reduced.split.se(x,
                               n.boot    = n.boot,
                               max.split = max.split,
                               B.full    = B.full,
                               p         = 0.5,
                               t_ratio   = t_ratio,
                               t_wls     = t_wls)

  se.hat        <- boot.out$se
  names(se.hat) <- est.names

  # Optimally combined estimators
  opt.all <- compute_optimal(boot.out$cov, alpha.hat)
  opt.wls <- compute_optimal(boot.out$cov[2:4, 2:4], alpha.hat[2:4])

  tab <- data.frame(
    Estimator = c(est.names, "optimal (all)", "optimal (WLS only)"),
    Estimate  = c(alpha.hat, opt.all$estimate, opt.wls$estimate),
    SE        = c(se.hat, sqrt(opt.all$var), sqrt(opt.wls$var)),
    stringsAsFactors = FALSE
  )
  tab$Lower <- tab$Estimate - 1.96 * tab$SE
  tab$Upper <- tab$Estimate + 1.96 * tab$SE

  cat("\n")
  print(format(tab, digits = 4), row.names = FALSE)

  list(label     = label,
       n         = length(x),
       alpha.hat = alpha.hat,
       cov       = boot.out$cov,
       Sigma     = boot.out$Sigma,
       Gamma     = boot.out$Gamma,
       opt.all   = opt.all,
       opt.wls   = opt.wls,
       table     = tab)
}

res <- lapply(names(periods),
              function(nm) analyse(periods[[nm]]$x, periods[[nm]]$label))
names(res) <- names(periods)


#  6. Side-by-side summary 

summary.tab <- data.frame(
  Estimator = res$old$table$Estimator,
  Early     = sprintf("%.4f (%.4f)", res$old$table$Estimate, res$old$table$SE),
  Late      = sprintf("%.4f (%.4f)", res$new$table$Estimate, res$new$table$SE),
  stringsAsFactors = FALSE
)

cat("\n===============================================================\n")
cat("  alpha estimates, standard errors in parentheses\n")
cat("===============================================================\n\n")
print(summary.tab, row.names = FALSE)


# Informal comparison of the two periods using the optimally combined
# estimator. The subperiods are disjoint, so treating them as independent is
# reasonable, but note this ignores any estimation of the change point and the
# fact that both were chosen after seeing the series.

d.est <- res$new$opt.all$estimate - res$old$opt.all$estimate
d.se  <- sqrt(res$new$opt.all$var + res$old$opt.all$var)

cat(sprintf("\nLate minus early (optimal combination): %.4f  (SE %.4f, z = %.2f)\n",
            d.est, d.se, d.est / d.se))
cat("Interpretation is left to the reader; see the caveat above.\n")


#  7. Optional comparison with standard estimators 

# Set to TRUE to add the StableEstim comparators. These are slow, especially
# the bootstrapped SEs and the ML fit.

RUN_COMPARISONS <- FALSE

if (RUN_COMPARISONS) {
  for (nm in names(periods)) {
    x <- periods[[nm]]$x

    Kout0 <- KoutParametersEstim(x)$Estim$par
    Kout  <- Kout0[1]
    McC   <- IGParametersEstim(x)[1]

    trad.ests <- array(NA, dim = c(n.boot, 2))
    for (b in seq_len(n.boot)) {
      Z.boot <- sample(x, size = length(x), replace = TRUE)
      trad.ests[b, ] <- c(KoutParametersEstim(Z.boot)$Estim$par[1],
                          IGParametersEstim(Z.boot)[1])
    }

    alpha.mle <- Estim(EstimMethod = "ML", data = x, ComputeCov = TRUE)

    cat("\n", periods[[nm]]$label, "comparators \n")
    print(round(data.frame(
      Estimator = c("Koutrouvelis", "McCulloch", "MLE"),
      Estimate  = c(Kout, McC, alpha.mle@par[1]),
      SE        = c(sqrt(diag(cov(trad.ests))), sqrt(alpha.mle@vcov[1, 1]))
    )[-1], 4))
  }
}


#  8. Save 

if (!dir.exists("output")) dir.create("output")
save(res, summary.tab, B.full, n.boot, max.split, QUICK,
     file = file.path("output", "xrp_results.RData"))

cat("\nSaved to output/xrp_results.RData\n")
if (QUICK) cat("NOTE: QUICK = TRUE. These are not the paper's settings.\n\n")
