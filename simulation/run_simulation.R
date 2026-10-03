# Simulation settings
models <- paste0("Model", 1:4)
designs <- c("SRS", "SBR", "PS")
sample_sizes <- c(500L, 1500L)
replications <- 1000L
outer_folds <- 2L
inner_folds <- 5L
workers <- 100L
run_name <- "formal_simulation"

script_file <- if (sys.nframe() > 0L) sys.frame(1)$ofile else NULL
if (is.null(script_file)) {
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else
    if (file.exists("run_simulation.R")) "run_simulation.R" else "simulation/run_simulation.R"
}
simulation_dir <- dirname(normalizePath(script_file))
root <- normalizePath(file.path(simulation_dir, ".."))
source(file.path(root, "src", "core.R"))

run_one <- function(job, output_dir, outer_folds, inner_folds, run_name) {
  directory <- file.path(output_dir, job$design, paste0("n", job$n), job$model)
  path <- file.path(directory, "replications", sprintf("rep_%03d.csv", job$replication))
  rows <- tryCatch(
    run_replication(job$model, job$replication, design = job$design, n = job$n,
                    outer_folds = outer_folds, inner_folds = inner_folds),
    error = function(e) empty_replication_rows(job$model, job$replication, job$design, job$n)
  )
  rows$run_id <- run_name
  atomic_write_csv(rows, path)
  invisible(NULL)
}

run_simulation <- function(models, designs, sample_sizes, replications,
                           outer_folds, inner_folds, workers, run_name, output_dir) {
  require_simulation_packages()
  positive_integer <- function(x) {
    is.numeric(x) && length(x) > 0L && all(is.finite(x)) && all(x >= 1 & x == floor(x))
  }
  if (!length(models) || anyNA(models) || anyDuplicated(models) ||
      !all(models %in% paste0("Model", 1:4))) stop("Choose models from Model1--Model4")
  if (!length(designs) || anyNA(designs) || anyDuplicated(designs) ||
      !all(designs %in% SIMULATION_DESIGNS)) stop("Choose designs from SRS, SBR, PS")
  if (!positive_integer(sample_sizes) || anyDuplicated(sample_sizes)) stop("Invalid sample sizes")
  if (length(replications) != 1L || !positive_integer(replications)) stop("Invalid replication count")
  if (length(workers) != 1L || !positive_integer(workers)) stop("Invalid worker count")
  if (length(outer_folds) != 1L || !positive_integer(outer_folds) || outer_folds < 2L ||
      length(inner_folds) != 1L || !positive_integer(inner_folds) || inner_folds < 2L) {
    stop("Fold counts must be integers of at least two")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1,
             VECLIB_MAXIMUM_THREADS = 1, NUMEXPR_NUM_THREADS = 1)

  grid <- expand.grid(design = designs, n = sample_sizes, model = models,
                      stringsAsFactors = FALSE)
  jobs <- lapply(seq_len(nrow(grid)), function(i) {
    cbind(grid[rep(i, replications), , drop = FALSE], replication = seq_len(replications))
  })
  jobs <- do.call(rbind, jobs)
  jobs <- split(jobs, seq_len(nrow(jobs)))
  workers <- min(as.integer(workers), length(jobs))
  message(sprintf("Running %d combinations, %d replications each, %d methods, %d workers",
                  nrow(grid), replications, nrow(method_catalog()), workers))

  if (workers > 1L) {
    cluster <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function(path, libraries) {
      .libPaths(libraries)
      source(path)
      RNGkind("Mersenne-Twister", "Inversion", "Rejection")
      NULL
    }, file.path(root, "src", "core.R"), .libPaths())
    invisible(parallel::parLapplyLB(cluster, jobs, run_one, output_dir = output_dir,
                                   outer_folds = outer_folds, inner_folds = inner_folds,
                                   run_name = run_name))
  } else {
    invisible(lapply(jobs, run_one, output_dir = output_dir, outer_folds = outer_folds,
                      inner_folds = inner_folds, run_name = run_name))
  }

  summaries <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    cell <- grid[i, ]
    directory <- file.path(output_dir, cell$design, paste0("n", cell$n), cell$model)
    result <- aggregate_replications(file.path(directory, "replications"), cell$model,
                                     cell$design, cell$n, replications)
    atomic_write_csv(result$raw, file.path(directory, "replications_combined.csv"))
    atomic_write_csv(result$summary, file.path(directory, paste0("summary_", cell$model, ".csv")))
    summaries[[i]] <- result$summary
    message(sprintf("Summarized %s / n%d / %s", cell$design, cell$n, cell$model))
  }
  atomic_write_csv(do.call(rbind, summaries), file.path(output_dir, "summary.csv"))
  message("Simulation results: ", output_dir)
  invisible(output_dir)
}

run_simulation(models, designs, sample_sizes, replications, outer_folds, inner_folds,
                workers, run_name, file.path(simulation_dir, "results", run_name))
