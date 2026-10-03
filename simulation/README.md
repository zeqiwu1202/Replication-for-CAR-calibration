# Simulation

Run from the repository root:

```bash
Rscript simulation/run_simulation.R
```

Or from this directory:

```bash
Rscript run_simulation.R
```

In an R session, use `source("simulation/run_simulation.R")` from the repository
root. Edit the settings at the top of the script to choose models, designs,
sample sizes, replication count, folds, workers, and run name. The default run
name is `formal_simulation`, with two outer folds, five inner folds, and 100
workers. Setting `workers <- 1L` runs sequentially.

The script sources the shared implementation in `../src/core.R`. Each replication
is computed on every run and saved, overwriting its existing result file. Results
include per-combination summaries and a combined `summary.csv`.

## Study and outputs

Models 1--4, designs SRS/SBR/PS, sample sizes 500 and 1500, and 1000 replications
per combination give 24 combinations. Each combination includes 22 methods.

Runs create `simulation/results/<run_name>/` from the repository root, or
`results/<run_name>/` within this directory. With the default settings, outputs
are written to `results/formal_simulation/`:

- `<design>/n<n>/<model>/replications/rep_<id>.csv`: individual replications.
- `<design>/n<n>/<model>/replications_combined.csv`: combined replications.
- `<design>/n<n>/<model>/summary_<model>.csv`: per-combination summary.
- `summary.csv`: summaries for all combinations.

Summaries contain bias, standard deviation, mean standard error, RMSE, coverage,
the Monte Carlo standard error of coverage, median absolute error, and replication
counts. A completed full run produces 528 summary rows (24 combinations
times 22 methods).
