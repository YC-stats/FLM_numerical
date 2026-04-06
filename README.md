# FLM_numerical

**Functional Linear Model (FLM) for Electricity Demand Forecasting**  
*R implementation by Yan Cui (YC-stats)*

---

## Overview

This repository implements a **Functional Linear Model (FLM)** for forecasting electricity demand. Rather than treating a day's electricity load profile (or temperature curve) as a vector of isolated hourly measurements, FLM treats it as a **continuous function** — capturing the full shape and dynamics of intraday patterns. This approach is grounded in Functional Data Analysis (FDA), a branch of statistics designed for data that are naturally expressed as curves or surfaces.

The core idea: given a functional predictor X(t) (e.g., yesterday's 24-hour load curve), predict a scalar or functional response Y (e.g., tomorrow's total demand or peak load) via the integral relationship:

```
Y = α + ∫ β(t) X(t) dt + ε
```

where β(t) is an unknown **coefficient function** estimated from data.

---

## Repository Structure

```
FLM_numerical/
├── app/          # Real-data application: electricity demand forecasting
├── simu/         # Simulation studies: benchmarking and method evaluation
└── README.md
```

### `app/` — Real Data Application

Contains R scripts that apply the FLM to real electricity demand data. Typical tasks performed here include:

- **Data loading and preprocessing** — reading in hourly/daily electricity consumption time series and covariates (temperature, calendar effects, etc.)
- **Functional data representation** — converting discrete time-series observations into functional objects using basis expansion (e.g., B-splines or Fourier bases)
- **FPCA (Functional Principal Component Analysis)** — decomposing the functional predictors into orthogonal principal components for dimension reduction
- **FLM fitting** — estimating the coefficient function β(t) by regressing the response onto the leading functional principal components
- **Forecasting and evaluation** — generating out-of-sample predictions and computing error metrics (e.g., RMSE, MAE)

### `simu/` — Simulation Studies

Contains R scripts for controlled numerical experiments to evaluate the FLM estimator's performance. Typical tasks include:

- **Data generation** — simulating functional covariates from a Gaussian process with a specified covariance structure and generating scalar responses according to a known β(t)
- **Method comparison** — benchmarking the FLM against baseline models under varying conditions
- **Sensitivity analysis** — testing the effect of sample size, measurement noise, and number of principal components on estimation accuracy
- **Performance metrics** — computing Average Mean Squared Prediction Error (AMSPE) across Monte Carlo replications

---

## Statistical Background

### Functional Linear Model

The FLM with a scalar response is defined as:

```
Yᵢ = α + ∫ β(t) Xᵢ(t) dt + εᵢ,   i = 1, ..., n
```

- **Xᵢ(t)**: functional predictor (e.g., intraday temperature or lagged load curve) observed on a domain T
- **β(t)**: coefficient function to be estimated — its shape reveals which time periods of X most influence Y
- **α**: scalar intercept
- **εᵢ**: i.i.d. zero-mean error

### FPCA-Based Estimation

The coefficient function β(t) is estimated using a two-step FPCA approach:

1. **FPCA on X**: Decompose each Xᵢ(t) into functional principal components (FPCs) φ_k(t) with scores ξᵢₖ = ∫ Xᵢ(t) φ_k(t) dt
2. **Scalar regression**: Regress Y on the leading FPC scores {ξᵢ₁, ..., ξᵢₖ} using ordinary least squares
3. **Reconstruct β(t)**: Back-project regression coefficients into the functional domain: β̂(t) = Σₖ b̂ₖ φ_k(t)

This approach is computationally efficient and statistically principled, relying on results from Yao et al. (2005) and Hall & Horowitz (2007).

---

## Requirements

**R version**: 4.0 or higher recommended

**Key R packages:**

| Package | Purpose |
|---|---|
| `fdapace` | FPCA and FLM estimation |
| `fda` | Functional data objects and basis representations |
| `refund` | Penalised regression for functional data |
| `ggplot2` | Visualisation |
| `dplyr` / `tidyr` | Data wrangling |
| `MASS` | Simulation utilities |

Install all dependencies at once:

```r
install.packages(c("fdapace", "fda", "refund", "ggplot2", "dplyr", "tidyr", "MASS"))
```

---

## Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/YC-stats/FLM_numerical.git
cd FLM_numerical
```

### 2. Run the simulation study

```r
# From the R console or RStudio, with working directory set to repo root:
source("simu/main_simu.R")
```

This will generate simulated functional predictors, fit the FLM, and report prediction errors across replications.

### 3. Run the real-data application

```r
source("app/main_app.R")
```

This will load the electricity data, fit the FLM, and output forecasts and evaluation metrics.

---

## Key Concepts for New Users

**Why use a Functional Linear Model instead of a standard regression?**

Electricity demand curves are inherently functional — they are smooth, continuous processes sampled at hourly intervals. Treating each hour as a separate predictor in a standard regression ignores the temporal smoothness and wastes degrees of freedom. The FLM instead:

- Preserves the **continuous-time structure** of the data
- Produces an interpretable **coefficient function** β(t) showing the influence of each time point
- Is **statistically efficient**, especially with many measurement points and limited observations

**What does β(t) tell us?**

In an electricity demand context, if X(t) represents yesterday's temperature curve and Y is today's total demand, then β(t) reveals at which hours of the day temperature most strongly drives demand — e.g., afternoon peaks in summer due to air conditioning.

---

## References

- Ramsay, J.O. and Silverman, B.W. (2005). *Functional Data Analysis* (2nd ed.). Springer.
- Yao, F., Müller, H.G., Wang, J.L. (2005). Functional linear regression analysis for longitudinal data. *Annals of Statistics*, 33, 2873–2903.
- Hall, P. and Horowitz, J.L. (2007). Methodology and convergence rates for functional linear regression. *The Annals of Statistics*, 35(1), 70–91.
- Hsing, T. and Eubank, R. (2015). *Theoretical Foundations of Functional Data Analysis, with an Introduction to Linear Operators*. Wiley.

---

## Author

**Yan Cui** — Statistician  
GitHub: [@YC-stats](https://github.com/YC-stats)

---

## License

This repository is publicly available for research and educational use. Please cite appropriately if this code contributes to published work.

