# Replication for CAR Calibration

This repository contains an R/Rcpp implementation for simulation studies and empirical analysis related to the paper entitled **“Integrating Heterogeneous Information in Randomized Experiments: A Unified Calibration Framework”** by Ma, Wu, and Zhang (2026).

## Attribution

**Important:** Parts of the implementation in this repository are adapted from code associated with **Tu, Ma, and Liu (2024)**. This repo reorganizes and extends that code for running the simulations/figures in this project.

Full reference:

Tu, F., Ma, W. & Liu, H. (2024), “A unified framework for covariate adjustment under stratified randomisation”, *Stat* 13(4), e70016.

If you use this repository, please make sure to **cite Tu, Ma, and Liu (2024)** in addition to citing this repository.

## Repository structure

- `run_model1.R`, `run_model2.R`, `run_model3.R`, `run_model4.R`: main simulation entry points for Models 1–4.
- `Model1.R`–`Model4.R`: data-generating processes (DGPs) for each model.
- `calibration_pipeline.R`: core calibration/estimation pipeline (e.g., `cal_main()`, `sdim()`, baseline estimators).
- `cal_estimator.R`: computes the calibration estimator by solving for calibration weights (uses `generate_objective_func.cpp`).
- `generate_objective_func.cpp`: RcppArmadillo code to build objective matrices used by calibration weight solvers.
- `vestimator.cpp`: Rcpp code used by the CAR assignment / helper routines.
- `createFolds_by_strat.R`: creates sample-splitting folds by splitting within each strata × treatment cell, then merging the per-cell folds into final folds.
- `estimators.R`: estimators adapted from Tu, Ma, and Liu (2024).
- `tables.R`: true value calculations for the models considered in simulations.
- `real_data.R`: empirical application using data in `Data/`.
- `Data/`: input datasets for the empirical analysis.
- `simulation_results/`: simulation outputs (CSV/PDF).
- `real_data_results/`: empirical outputs (CSV/PDF).

## Requirements

- R (tested on Linux; any recent R should work)
- R packages used across scripts include (not exhaustive):
  - `Rcpp`, `RcppArmadillo`, `MASS`, `dplyr`
  - `randomForest`, `caret`, `glmnet`, `kernlab`, `nnet`, `gbm`, `rpart`
  - `caretEnsemble`, `elasticnet`, `np`, `earth`, `neuralnet`
  - `quantreg`, `pls`, `mda`, `splines`
  - `ggplot2`, `haven`

If you see a “package not found” error, install the missing package via `install.packages("<pkg>")`.

## Quick start

Note: Run the commands below inside an **R session** started in the repository root (so that relative paths like `Data/` and `simulation_results/` resolve correctly). Do not run `source("...")` in a bash shell.

### 1) Run simulations (Models 1–4)

From an R session started in the repository root:

```r
source("run_model1.R")
source("run_model2.R")
source("run_model3.R")
source("run_model4.R")
```

Each script:
- Ensures the output folder `simulation_results/` exists.
- Produces boxplot PDFs like:
  - `simulation_results/YYYYMMDD_Model1_n500_SRS.pdf`
- Produces a summary CSV like:
  - `simulation_results/YYYYMMDD_Model1.csv`

You can control the simulation size by editing `S` (number of repetitions) and `sample_sizes` inside the corresponding `run_model*.R`.

### 2) Run the empirical analysis

The empirical script expects the `.dta` files under `Data/` to be present.

```r
source("real_data.R")
```

Outputs:
- `real_data_results/empirical.csv`
- `real_data_results/ATE_estimates_in_Uganda.pdf`
- `real_data_results/ATE_estimates_in_Malawi.pdf`

## Notes on compilation (Rcpp)

Several scripts call `Rcpp::sourceCpp()` / `sourceCpp()` to compile C++ code at runtime:
- `vestimator.cpp`
- `generate_objective_func.cpp`

On Linux this typically works out-of-the-box if you have a working compiler toolchain. If compilation fails, install system build tools (e.g., `g++`, `make`) and ensure your R is configured for package compilation.

## Reproducibility

Most scripts set:

```r
RNGkind("L'Ecuyer-CMRG")
set.seed(0)
```

so that parallel-capable RNG streams are reproducible.

## License

No explicit license file is included. Please treat this code as research code and ensure that any reused components adapted from Tu, Ma, and Liu (2024) comply with the original code’s licensing/terms.

