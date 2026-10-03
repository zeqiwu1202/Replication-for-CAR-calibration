# Integrating Heterogeneous Information in Randomized Experiments: A Unified Calibration Framework

Replication code for **Integrating Heterogeneous
Information in Randomized Experiments: A Unified Calibration Framework**,
Ma, Wu, and Zhang (2026). Preprint: <https://arxiv.org/abs/2603.07055>.

## Repository layout

```text
src/core.R                       Shared methods and simulation functions
simulation/run_simulation.R      Run simulations and summarize results
simulation/README.md             Simulation settings and instructions
simulation/results/              Simulation outputs (created when running)
empirical/run_empirical.R        Eight-estimator Uganda/Malawi empirical analysis
empirical/Data/                  Empirical input datasets
empirical/results/reproduced/    Empirical outputs (created when running)
```

## Requirements

- R with `MASS`, `randomForest`, and `neuralnet` for simulations.
- Empirical analysis uses `haven` and `randomForest`.

Install the required packages in R:

```r
install.packages(c("MASS", "randomForest", "neuralnet", "haven"))
```

Verified with R 4.4.3, MASS 7.3-61, randomForest 4.7-1.2, neuralnet 1.44.2,
and haven 2.5.5.

## Simulation study

- Models 1--4, 30 covariates, sample sizes 500 and 1500, 1000 replications.
- SRS, stratified permuted blocks of size 6 (SBR), and Pocock--Simon minimization
  with biased-coin probability 0.75 (PS).
- Two stratum-by-treatment-balanced outer folds, and five balanced inner folds
  for Super Learner/stacking.
- 24 model/design/sample-size combinations, with 22 methods per combination:
  quadratic/quartic calibration with RF, NN, RF+NN, grouped RF,
  grouped NN, and RF+linear libraries; AIPW RF/NN; AIPW Super Learner/stacking
  with RF+NN, grouped RF, and grouped NN libraries; full-sample `sdim`; and
  full-sample stratum-by-arm OLS (`lin`) without sample splitting.

Random forests use 50 trees. Neural networks have two hidden units and scale
outcomes using the training sample. Calibration uses the fold-specific
`v1+v2-v3` variance decomposition with correction `n_fk/(n_fk-rank_fk-1)`.
Cross-fitted AIPW uses the CAR variance within each fold and the fold's overall
realized treatment share. `sdim` and `lin` use the full sample.

Run the simulations from the repository root:

```bash
Rscript simulation/run_simulation.R
```

Or run it in an R session with the working directory set to the repository root:

```r
source("simulation/run_simulation.R")
```

The settings at the top of `simulation/run_simulation.R` specify the models, designs,
sample sizes, replication count, fold counts, worker count, and run name.
Defaults run all 24 combinations with 1000 replications and 100 workers.
Change `workers` to control parallelism or restrict `models`, `designs`, and
`sample_sizes` to run a subset. Changing `replications` allows a shorter run.
The script uses R's `parallel` package on Linux, macOS, and Windows; setting
`workers <- 1L` runs sequentially.

Results go to `simulation/results/<run_name>/`. Every run recomputes the selected
replications and overwrites their result files. Use a different `run_name` when
changing the replication count, fold counts, or implementation, to keep runs separate.
The script writes per-combination results plus the combined `summary.csv`.
The default run name is `formal_simulation`.
Progress messages appear in the R console.
The repository includes the [simulation summary](simulation/results/formal_simulation/summary.csv).

## Empirical application

The empirical analysis compares eight methods for Uganda and Malawi:
`sdim`, quadratic and quartic calibration with `X`, `info`, or `X_info`, plus
`aipw_info`. External outcome predictions come from 40-tree random forests
trained in the other country, using seed 0 with `L'Ecuyer-CMRG`. Estimation uses
the full sample without cross-fitting. All estimators use the shared methods
in `src/core.R`.

```bash
Rscript empirical/run_empirical.R
# Optional paths are resolved relative to the caller's working directory:
Rscript empirical/run_empirical.R --data-dir empirical/Data \
  --output-dir empirical/results/reproduced
```

Or run `source("empirical/run_empirical.R")` from the repository root in R.

The script writes `estimates.csv`, with one row per country and estimator.
The default output directory is `empirical/results/reproduced/`.
The repository includes the [empirical estimates](empirical/results/reproduced/estimates.csv).

[empirical/README.md](empirical/README.md) documents the data, preprocessing,
estimators, and output columns.

## Attribution

Parts of the implementation were adapted from Tu, Ma, and Liu (2024),
*A unified framework for covariate adjustment under stratified randomisation*,
**Stat**, 13(4), e70016. Cite both that paper and Ma, Wu, and Zhang (2026) when
using this replication code.

The empirical datasets are from the replication package for Dupas, Karlan,
Robinson, and Ubfal (2018), *Banking the Unbanked? Evidence from Three Countries*,
**American Economic Journal: Applied Economics**, 10(2), 257--297.
Source: <https://www.openicpsr.org/openicpsr/project/116346/version/V1/view>.
