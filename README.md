# Semiparametric Models with Errors-in-Variables under Right-Censoring

This repository contains simulation code for comparing estimation methods in semiparametric models when covariates are measured with error and responses are right-censored.

## Overview

The code implements and compares different estimation approaches for semiparametric models of the form:

```
Y = X'β + m(T) + ε
```

where:
- `Y` is the response (subject to right-censoring)
- `X` are error-free covariates 
- `T` is a covariate measured with error (we observe `W = T + U`)
- `m(T)` is an unknown smooth function
- `β` are regression coefficients

## Key Features

- Measurement Error Correction: Implements deconvolution methods to handle errors-in-variables (Core contribution)
- Censoring Handling: Uses Buckley-James method for right-censored responses
- Multiple Estimators: Compares Full Model (FM), Submodel (SM), Pretest (PT), Shrinkage (S), and Positive Shrinkage (PS) estimators
- Flexible Smoothing: Supports both kernel smoothing and local polynomial methods
- High-Dimensional Support: Handles both low and high-dimensional settings

## Methods Implemented

### Nonparametric Estimation
- **Naive Kernel Smoothing**: Standard kernel methods ignoring measurement error
- **Deconvolution Kernel Smoothing**: Corrects for measurement error using deconvolution kernels
- **Local Polynomial Smoothing**: Alternative to kernel methods with both naive and deconvolution versions

### Parametric Estimation
- **Full Model (FM)**: Uses all covariates (Ridge regression in high-dim)
- **Submodel (SM)**: Variable selection via AIC (Lasso in high-dim)
- **Pretest (PT)**: Choose FM or SM based on hypothesis test
- **Shrinkage (S)**: James-Stein type shrinkage between FM and SM
- **Positive Shrinkage (PS)**: Positive-part shrinkage estimator

## Usage

### Basic Setup

```r
# Required packages
library(survival)
library(glmnet)

# Source the main functions
source("simulation_functions.R")
```

### Running Simulations

```r
# Low-dimensional case (n=100, p=10)
results_low <- run_simulation(
  n_sim = 100,
  n = 100, 
  p = 10,
  sigma_u = 0.1,          # measurement error variance
  censoring_rate = 0.2,   # proportion censored
  high_dim = FALSE
)

# High-dimensional case (n=100, p=300)  
results_high <- run_simulation(
  n_sim = 100,
  n = 100,
  p = 300,
  sigma_u = 0.1,
  censoring_rate = 0.2,
  high_dim = TRUE
)

# Display results
format_results(results_low, "Low-dimensional Results")
format_results(results_high, "High-dimensional Results")
```

### Key Parameters

- `n_sim`: Number of simulation replications
- `n`: Sample size
- `p`: Number of covariates
- `sigma_u`: Standard deviation of measurement error
- `sigma_eps`: Error term standard deviation  
- `censoring_rate`: Proportion of censored observations
- `h_range`: Range for bandwidth selection
- `high_dim`: Whether to use high-dimensional methods (Ridge/Lasso)
- `smooth_method`: Either "kernel" or "local_poly"

## Data Generation

The simulation generates data where:
- True covariate `T ~ Uniform(0,1)`
- Measurement error `U ~ N(0, σ²_u)`
- Observed covariate `W = T + U`
- Nonparametric function `m(t) = sin(2πt)`
- First 5 regression coefficients equal 1, rest are 0
- Censoring is generated to achieve specified censoring rate

## Performance Metrics

- **MSE**: Mean squared error for parametric component β
- **MISE**: Mean integrated squared error for nonparametric component m(·)
- **Relative Efficiency**: Comparison between estimators

## Key Functions

- `generate_data()`: Generates simulation data with measurement error and censoring
- `buckley_james()`: Handles right-censored responses
- `deconvolution_kernel()`: Implements deconvolution for measurement error correction
- `kernel_smoothing()`: Kernel-based nonparametric estimation
- `local_poly_smooth()`: Local polynomial alternative
- `estimate_parametric()`: Parametric component estimation
- `compute_shrinkage_estimators()`: Implements PT, S, and PS estimators
- `run_simulation()`: Main simulation runner
- `format_results()`: Pretty-prints results

## Expected Results

The deconvolution methods should generally outperform naive approaches, especially when measurement error is substantial. The positive shrinkage (PS) estimator typically shows the best performance among the parametric estimators.

## Dependencies

- R ≥ 3.6.0
- `survival` package for Kaplan-Meier estimation
- `glmnet` package for Ridge/Lasso regression

## Citation

If you use this code, please cite the associated papers:
- Aydın, D., Yılmaz, E., Chamidah, N., & Lestari, B. (2024). Right-censored partially linear regression model with error in variables: application with carotid endarterectomy dataset. The International Journal of Biostatistics, 20(1), 245-278.
- Aydın, D., Yılmaz, E., Chamidah, N., Lestari, B., & Budiantara, I. N. (2025). Right-censored nonparametric regression with measurement error. Metrika, 88(2), 183-214. 



## Contact
ersin.yilmaz@mu.edu.tr
