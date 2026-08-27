# Sample-Splitting Estimation of the Stable Tail Index

R code accompanying the paper on estimating the (stable) tail index $\alpha$ or
transformed scale $\sigma = 2^{1/\alpha}$ using repeated sample splitting, and a
two-sample location-scale formulation. Also implements a reduced-split
bootstrap for inference. For the current manuscript, see arXiv:1705.09840.

This repo contains the methodology functions plus two runnable scripts: a
worked example that illustrates the methods on simulated data, and the
application to XRP log-returns in the paper.

---

## Structure

```
sample-splitting/
├── README.md
├── load_all.R                         # sources everything under R/
├── R/
│   ├── dependencies.R
│   ├── stable_rng.R                   # rstable_new()
│   ├── estimation_functions.R         # wls.ests.tvec(), run.split()
│   ├── bootstrap_functions.R          # run.split.mult()
│   ├── reduced_split_bootstrap.R      # reduced.split.cov/se()
│   └── combination_functions.R        # compute_optimal()
├── scripts/
│   ├── 01_simulated_example.R         # simulate -> estimate -> SE
│   └── 02_xrp_application.R           # real data -> estimate -> SE
├── data/
│   └── ripple_2013-08-05_2025-05-27.csv
└── simulation_results/
    ├── sim_combined.rds                # estimator-level simulation results for each config
    └── rmse_sigma.csv                  # RMSE summaries on sigma scale
    └── rmse_alpha.csv                  # RMSE summaries on alpha scale
```

---

## Installation

```r
install.packages(c("stabledist", "StableEstim", "ks"))
```

`kde` comes from **ks**; `stableMode` from **stabledist**, used by
`rstable_new()` only when `pm != 0`; and the comparison estimators
(`KoutParametersEstim`, `IGParametersEstim`, `Estim`) from **StableEstim**,
used only in the optional comparison block of `scripts/02`.

---

## Quick start

```r
source("load_all.R")

t_ratio <- list(c(1/4, 3/4))
t_wls   <- lapply(c(5, 10, 16), function(m)
  round(seq(1/(m + 1), m/(m + 1), length.out = m), 3))

Z <- rstable_new(400, alpha = 1.2, beta = 0.25)

# Point estimate
sigma.hat <- run.split(Z, n.splits = 250, p = 0.5, t_ratio, t_wls)
sigma.hat <- unlist(sigma.hat)
alpha.hat <- log(2) / log(sigma.hat)

# Standard errors for alpha.hat
boot.out <- reduced.split.se(Z, n.boot = 200, max.split = 25, B.full = 250,
                             p = 0.5, t_ratio = t_ratio, t_wls = t_wls)
boot.out$se
```

The estimator works on the $\sigma$ scale; the transformation
$\hat\alpha = \log 2 / \log \hat\sigma$ is applied explicitly by the caller.
Pass `alpha.scale = FALSE` to `reduced.split.se()` to get the covariance on the
$\sigma$ scale instead.

---

## The two scripts

### `scripts/01_simulated_example.R`

End-to-end walkthrough with a known truth. Generates stable data with
`rstable_new()`, estimates $\sigma$ and $\alpha$ with four tuning choices,
computes reduced-split bootstrap standard errors on both scales, and forms the
inverse-covariance optimal combination on the $\alpha$ scale. It also shows a
plot of how SE extrapolation works in the reduced-split bootstrap, including
the floor that additional splits cannot remove.

Runtime roughly 2–5 minutes at $n = 500$, $B = 250$, 200 bootstrap replicates.
Note that the bootstrap is run twice, once per scale.

### `scripts/02_xrp_application.R`

Reads the CSV from `data/`, converts to daily log-returns, and analyses the
outer thirds of the sample (1,437 returns each) as an early and a late regime.
Produces the shaded series plot and the two Q-Q panels, then estimates and
bootstraps each subperiod, reports both subperiods side by side with 95%
intervals, and saves to `output/xrp_results.RData`.

The script defaults to `QUICK <- TRUE` (a ~5 minute sanity run). **Set
`QUICK <- FALSE` for the paper's settings** ($B = 500$, 500 bootstrap
replicates, `max.split = 25`).

---

## Method

**`wls.ests.tvec(X, Y, t0)`** — two-sample WLS scale estimator. Builds the
Brownian-bridge covariance $\Sigma_{ij} = \min(t_i,t_j) - t_i t_j$, estimates
the density of `Y` by KDE on the `asinh` scale with bandwidth `bw.nrd`, forms
the weight matrix, and returns the WLS slope.

**`run.split(Z, n.splits, p, t_ratio, t_wls)`** — main estimator. For each of
`n.splits` Bernoulli($p$) splits, forms the baseline sample `X` and all unique
pairwise sums `Y` from the complement, computes the ratio estimator at each
element of `t_ratio` and the WLS estimator at each element of `t_wls`, and
averages across splits.

**`run.split.mult(...)`** — same, but returns the *cumulative* mean after
$1,\dots,B$ splits. This is what makes the reduced-split bootstrap possible.

**`reduced.split.se(...)`** — the inference procedure. Running $B = 500$ splits
inside every bootstrap replicate is prohibitive, so it runs only
`max.split` (25) and extrapolates using

$$\mathrm{Var}(\hat\alpha_B) = \frac{\Sigma}{B} + \frac{B-1}{B}\Gamma,$$

solving for $(\Sigma, \Gamma)$ from each pair $B_1 < B_2$ and averaging over
all pairs. With `alpha.scale = TRUE` (the default) the $\alpha$ transform is
applied *before* computing covariances, so the returned covariance is on the
$\alpha$ scale. The reconstruction is not guaranteed positive definite — the
subtraction $B_1 C_{B_1} - (B_1-1)\Gamma$ can produce small negative
eigenvalues — so eigenvalues are floored at `pd.tol` times the mean diagonal
and `$cov.repaired` records whether that was necessary.

**`compute_optimal(covmat, est)`** — inverse-covariance optimal linear
combination. Weights are unconstrained and may be negative when estimators are
strongly correlated.

```
Z ──► run.split() ──► sigma.hat ──► alpha.hat = log(2)/log(sigma.hat)
             │
             └── calls wls.ests.tvec() on each split

Z ──► reduced.split.se() ──► cov, SE ──► compute_optimal()
             │
             └── calls run.split.mult() on each bootstrap resample
```

---

## Simulation results

The `simulation_results/` directory contains the saved output from the
simulation studies reported in the paper, including additional RMSE summaries
not shown in the main text. The `.rds` file contains the individual estimator
values for each simulation configuration and replication, while the `.csv`
files contain the corresponding RMSE summaries used for reporting and
comparison.
