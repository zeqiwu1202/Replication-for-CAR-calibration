# Uganda and Malawi empirical application

The analysis estimates the effect of access to a bank account on savings in
Uganda and Malawi. Each country supplies random-forest outcome predictions for
the other country's calibration estimators.

Run from the repository root:

```bash
Rscript empirical/run_empirical.R
```

`--data-dir` and `--output-dir` accept custom paths relative to the caller's
working directory.

## Data and preprocessing

The two input files are from the replication package for Dupas, Karlan, Robinson,
and Ubfal (2018), *Banking the Unbanked? Evidence from Three Countries*,
**American Economic Journal: Applied Economics**, 10(2), 257--297.
[Replication package](https://www.openicpsr.org/openicpsr/project/116346/version/V1/view).

- `Data/Data_Uganda/glaseu_four_rounds.dta`: Uganda survey data.
- `Data/Data_Malawi/glasem_four_rounds.dta`: Malawi survey data.

The outcome is the sum of nine savings components at the first follow-up
(`mv1_amountsaved_resp_*`): home, bank, SACCO, ROSCA, friend, mobile, shop, leader,
and farm group. Each component is winsorized at the 99th percentile before
summation. The total is converted using factors 0.000368169 for Uganda and
0.005928101 for Malawi. The baseline covariate is `b_amountsaved_resp_tot2`.
Rows with missing analysis variables and strata with at most six observations
are excluded from the estimation sample.
The resulting samples contain 2115 observations in 37 strata for Uganda and
1987 observations in 67 strata for Malawi.

## Estimators

Estimation uses the full sample without cross-fitting. Both the empirical and
simulation scripts use the shared methods in `../src/core.R`.

- `sdim`: stratified difference in means, without covariate adjustment.
- `cal2_X`, `cal2_info`, `cal2_X_info`: quadratic calibration.
- `cal4_X`, `cal4_info`, `cal4_X_info`: quartic calibration.
- `aipw_info`: AIPW with the other country's RF outcome predictions.

For calibration, `X` is the untransformed baseline savings covariate, `info`
contains the external control and treatment predictions, and `X_info` combines
all three columns. External RF models use 40 trees and seed 0 with
`L'Ecuyer-CMRG`. Models trained on Malawi first supply predictions for Uganda;
models trained on Uganda then supply predictions for Malawi.

The six calibration methods call `calibration_fold()` and share the same variance
formula for a given variable library, with within-stratum effective rank and
correction `n_k/(n_k-rank_k-1)`. `sdim` and `aipw_info` call `estimate_car()`.
The required packages are `haven` and `randomForest`.

## Results

The script writes `results/reproduced/estimates.csv` by default. The table has
16 rows, one per country and estimator, with country, method, ATE, SE, normal
95% confidence limits, sample size, and number of strata.
The CSV retains full precision; console output shows three decimal places.

Custom paths can be supplied with:

```bash
Rscript empirical/run_empirical.R --data-dir empirical/Data \
  --output-dir empirical/results/reproduced
```

The script can also be run with `source("empirical/run_empirical.R")` from the
repository root in R.
