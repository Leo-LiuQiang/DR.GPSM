
# Doubly-Robust Generalized Propensity Score Matching

## Overview

**DR.GPSM** (Doubly-Robust Generalized Propensity Score Matching) is an
R package implementing the *Balance-based Generalized Propensity Score
Matching* (BGPSM) framework proposed by Dr. Mingrui Zhang and Dr. Yuan
Liu.

This method addresses the challenge of estimating **average treatment
effects (ATE)** in **observational studies with multi-valued
treatments**, combining:

- **Generalized Propensity Score (GPS)** modeling, to account for
  confounding;
- **Matching** strategies across multiple treatment groups, to construct
  balanced cohorts;
- **Doubly-robust adjustment** using optional outcome models to improve
  estimation efficiency and reduce bias.

The package is designed for researchers in **causal inference**,
**epidemiology**, **health services research**, and **social sciences**
who need flexible and statistically rigorous tools for ATE estimation
with **2+ treatment groups**.

------------------------------------------------------------------------

## Key Features

- **Generalized Propensity Score (GPS) modeling** using:

  - Logistic multinomial regression (`logit`)
  - Random Forest (`rf`)
  - Gradient Boosted Models (`gbm`)
  - Generalized Additive Models (`gam`)

- **Balance-based matching algorithm** that constructs matched sets
  across all treatment levels simultaneously.

- **Doubly-robust estimation** with optional outcome regression (`none`,
  `lm`, `rf`) for bias correction.

- **Visualization utilities** for GPS distribution (overlay histograms).

- Designed for both **binary** and **multi-valued** treatments.

------------------------------------------------------------------------

## Background

In observational studies, treatments are not randomly assigned, leading
to potential **confounding bias** when estimating the causal effect of
treatment on outcomes.

**Imbens (2000)** generalized the propensity score framework to
multi-valued treatments. Since then, multiple approaches (matching,
weighting, regression, clustering) have been developed to estimate ATE
in these settings.

**Matching-based strategies** (Lechner, 2001; Yang et al., 2016) are
attractive because they create *balanced pseudo-cohorts* that emulate
randomized trials. However, traditional pairwise or trio-based matching
approaches suffer from:

- **Inconsistent target populations** across pairwise contrasts;
- **Computational burden** when the number of treatment groups grows.

The **Balance-based GPS Matching** method by Zhang and Liu (2022)
improves upon these limitations by:

- Matching across **all treatment levels simultaneously**;
- Introducing a **balance-based stopping rule** to prune poor matches;
- Allowing **doubly-robust bias correction** via outcome regression.

------------------------------------------------------------------------

## Functions

| Function            | Description                                                                    |
|---------------------|--------------------------------------------------------------------------------|
| `dr_gpsm()`         | **Main API**: Doubly-robust GPS matching estimator for ATE.                    |
| `gps_pre_process()` | Fits the GPS model and augments the dataset with `gps_` and `loggps_` columns. |
| `gps_matching()`    | Core matching engine that computes pairwise contrasts and bootstrapped CIs.    |
| `build_contrast()`  | Generates all pairwise contrast matrices for multiple treatments.              |
| `gps_histogram()`   | Plots overlaid histograms of GPS distributions for all treatment levels.       |

------------------------------------------------------------------------

## Methodology

**Causal Framework**

- **Treatments:** $T \in \{ t_1, t_2, ..., t_Z \}$  
- **Outcome:** Continuous response $Y$  
- **Covariates:** $X$

Assumptions: - **Strong unconfoundedness:**
$T \perp (Y(t_1), ..., Y(t_Z)) \mid X$  
- **Positivity:** $P(T = t \mid X = x) > \eta$

Under these assumptions, the **average treatment effect** is
identifiable: $$
\tau(t,t') = E[Y(t) - Y(t')]
$$

**BGPSM Algorithm Highlights** 1. Estimate GPS:
$r(t,x) = P(T = t \mid X = x)$ 2. Construct **matching sets** using
nearest neighbor search on the GPS vector. 3. Apply a **balance-based
stopping rule** to exclude poor matches. 4. Optionally fit **outcome
regression models** (linear, RF) for **doubly-robust** bias correction.
5. Compute ATEs with bootstrap confidence intervals.

------------------------------------------------------------------------

## Installation

``` r
# From GitHub (recommended)
devtools::install_github("Leo-LiuQiang/DR.GPSM")
```

## Usage

``` r
library(DR.GPSM)

# Simulated example
set.seed(123)
n <- 100
dat <- data.frame(
  trt = factor(sample(c("A","B","C"), n, replace=TRUE)),
  x1 = rnorm(n),
  x2 = runif(n),
  y  = rnorm(n)
)

# Run doubly-robust GPS matching
res <- dr_gpsm(
  data = dat,
  treatment = 1,
  treatment_ref = "A",
  covariate = 2:3,
  outcome = 4,
  gps_model = "logit",
  outcome_model = "lm",
  folds = 2,
  nboot = 50
)

res$estimate

# Visualization of generalized propensity score overlap
gps_df <- gps_pre_process(dat, treatment=1, covariate=2:3, gps_model="logit")
p <- gps_histogram(gps_df)
print(p)
```
