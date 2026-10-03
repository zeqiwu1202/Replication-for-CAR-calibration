#!/usr/bin/env Rscript
# Empirical application for Uganda and Malawi.

seed <- 0L
rf_trees <- 40L

script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
args <- if (is.null(script_file)) commandArgs(trailingOnly = TRUE) else character()
if (is.null(script_file)) {
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else
    if (file.exists("run_empirical.R")) "run_empirical.R" else
      "empirical/run_empirical.R"
}
script_dir <- dirname(normalizePath(script_file, mustWork = TRUE))
data_dir <- file.path(script_dir, "Data")
output_dir <- file.path(script_dir, "results", "reproduced")
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--help") {
    cat("Usage: Rscript empirical/run_empirical.R [--data-dir PATH] [--output-dir PATH]\n")
    quit(status = 0L)
  }
  if (!(args[i] %in% c("--data-dir", "--output-dir")) || i == length(args)) {
    stop("Expected --data-dir PATH or --output-dir PATH.")
  }
  if (args[i] == "--data-dir") data_dir <- args[i + 1L]
  if (args[i] == "--output-dir") output_dir <- args[i + 1L]
  i <- i + 2L
}
# Relative paths use the caller's working directory.
data_dir <- normalizePath(data_dir, mustWork = TRUE)

packages <- c("haven", "randomForest")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
source(file.path(script_dir, "..", "src", "core.R"))
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)

read_country <- function(country, conversion, stratum_count) {
  path <- file.path(data_dir, paste0("Data_", country),
                    if (country == "Uganda") "glaseu_four_rounds.dta" else
                      "glasem_four_rounds.dta")
  survey <- haven::read_dta(path)
  savings <- grep("^mv1_amountsaved_resp_(home|bank|sacco|rosca|friend|mobile|shop|leader|farmgroup)$",
                  names(survey), value = TRUE)
  for (name in savings) {
    cutoff <- quantile(survey[[name]], 0.99, na.rm = TRUE)
    survey[[name]] <- ifelse(survey[[name]] > cutoff & !is.na(survey[[name]]),
                            cutoff, survey[[name]])
  }
  other_cash <- survey$mv1_amountsaved_resp_shop +
    survey$mv1_amountsaved_resp_leader + survey$mv1_amountsaved_resp_farmgroup
  total <- survey$mv1_amountsaved_resp_home + survey$mv1_amountsaved_resp_bank +
    survey$mv1_amountsaved_resp_sacco + survey$mv1_amountsaved_resp_rosca +
    survey$mv1_amountsaved_resp_friend + survey$mv1_amountsaved_resp_mobile +
    other_cash
  strata <- numeric(nrow(survey))
  for (k in seq_len(stratum_count)) {
    strata[which(survey[[paste0("strata", k)]] == 1)] <- k
  }
  na.omit(data.frame(Y = total * conversion, X = survey$b_amountsaved_resp_tot2,
                     A = survey$treated, S = strata))
}

external_predictions <- function(target, external) {
  fit1 <- randomForest::randomForest(Y ~ X, data = external[external$A == 1, ],
                                     ntree = rf_trees)
  fit0 <- randomForest::randomForest(Y ~ X, data = external[external$A == 0, ],
                                     ntree = rf_trees)
  cbind(info0 = as.numeric(predict(fit0, newdata = target)),
        info1 = as.numeric(predict(fit1, newdata = target)))
}

estimate_country <- function(data, info, country) {
  counts <- table(data$S)
  keep <- data$S %in% names(counts[counts > 6])
  data <- data[keep, ]
  info <- info[keep, , drop = FALSE]
  n <- nrow(data)
  strata <- sort(unique(data$S))
  eval <- list(Y = data$Y, A = data$A, S = data$S)
  proxies <- list(X = cbind(X = data$X), info = info,
                  X_info = cbind(X = data$X, info))
  sdim <- estimate_car(eval)
  estimates <- c(sdim = sdim$estimate)
  variances <- c(sdim = sdim$variance)

  for (degree in c(2L, 4L)) {
    discrepancy <- if (degree == 2L) "quadratic" else "quartic"
    for (library in names(proxies)) {
      method <- paste0("cal", degree, "_", library)
      fit <- calibration_fold(eval, proxies[[library]], discrepancy)
      estimates[method] <- fit$estimate
      variances[method] <- fit$variance$variance
    }
  }

  aipw_info <- estimate_car(eval, info[, "info0"], info[, "info1"])
  estimates["aipw_info"] <- aipw_info$estimate
  variances["aipw_info"] <- aipw_info$variance

  se <- sqrt(unname(variances))
  data.frame(country = country, method = names(estimates),
             estimate = unname(estimates), se = se,
             lower = unname(estimates) - 1.96 * se,
             upper = unname(estimates) + 1.96 * se,
             n = n, strata = length(strata))
}

countries <- list(Uganda = read_country("Uganda", 0.000368169, 41L),
                  Malawi = read_country("Malawi", 0.005928101, 78L))
# Train on Malawi to predict Uganda, then train on Uganda to predict Malawi.
information <- list(
  Uganda = external_predictions(countries$Uganda, countries$Malawi),
  Malawi = external_predictions(countries$Malawi, countries$Uganda)
)
results <- do.call(rbind, lapply(names(countries), function(country) {
  estimate_country(countries[[country]], information[[country]], country)
}))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(output_dir, "estimates.csv"), row.names = FALSE)

display <- results
display[c("estimate", "se", "lower", "upper")] <-
  lapply(display[c("estimate", "se", "lower", "upper")], round, digits = 3L)
print(display, row.names = FALSE)
cat("Results: ", file.path(output_dir, "estimates.csv"), "\n", sep = "")
