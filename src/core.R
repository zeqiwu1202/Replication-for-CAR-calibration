SVD_TOL <- 1e-10
BALANCE_TOL <- 1e-8
SIMULATION_DESIGNS <- c("SRS", "SBR", "PS")

require_simulation_packages <- function() {
  pkgs <- c("MASS", "randomForest", "neuralnet")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  invisible(pkgs)
}

seed_int <- function(...) {
  x <- as.double(unlist(list(...)))
  weights <- seq_along(x) * 104729
  as.integer((sum((x + 37) * weights) %% 2147483000) + 1)
}

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.csv(x, tmp, row.names = FALSE, na = "NA")
  if (!file.rename(tmp, path)) stop("Atomic rename failed for ", path)
  invisible(path)
}

pinv_svd <- function(x, tol = SVD_TOL) {
  x <- as.matrix(x)
  if (!length(x)) return(list(inverse = x, rank = 0L))
  sx <- svd(x, nu = min(dim(x)), nv = min(dim(x)))
  scale <- if (length(sx$d)) max(sx$d) else 0
  cutoff <- tol * max(1, scale)
  keep <- sx$d > cutoff
  inv <- matrix(0, nrow = ncol(x), ncol = nrow(x))
  if (any(keep)) {
    inv <- sx$v[, keep, drop = FALSE] %*%
      (t(sx$u[, keep, drop = FALSE]) / sx$d[keep])
  }
  list(inverse = inv, rank = sum(keep))
}

equi_cov <- function(q, rho) {
  out <- matrix(rho, q, q)
  diag(out) <- 1
  out
}

toeplitz_cov <- function(q, rho) toeplitz(rho^(0:(q - 1)))

model_parameters <- function(model) {
  switch(model,
    Model1 = list(alpha0 = 1, alpha1 = 4, beta0 = c(75, 35, 125, 80),
                  beta1 = c(100, 80, 60, 40), sigma0 = 1, sigma1 = 3),
    Model2 = list(alpha0 = -3, alpha1 = 0, beta0 = c(10, 24, 15, 20),
                  beta1 = c(20, 17, 10), sigma0 = 1, sigma1 = 3),
    Model3 = list(alpha0 = 5, alpha1 = 2, beta0 = c(42, 83),
                  beta1 = c(30, 75), sigma0 = 1, sigma1 = 3),
    Model4 = list(alpha0 = 5, alpha1 = 5, beta0 = c(20, 30, 50),
                  beta1 = c(20, 30, 65), sigma0 = 1, sigma1 = 3),
    stop("Unknown DGP: ", model)
  )
}

normal_square_reciprocal_moment <- function(scale) {
  a <- abs(as.numeric(scale))
  a2 <- a^2
  small <- a <= 0.05
  out <- 1 / 3 - a2 / 9 + a2^2 / 9 - 5 * a2^3 / 27 + 35 * a2^4 / 81
  if (any(!small)) {
    q <- sqrt(3 / 2) / a[!small]
    log_erfcx <- q^2 + log(2) + stats::pnorm(-q * sqrt(2), log.p = TRUE)
    out[!small] <- sqrt(base::pi / 2) / (sqrt(3) * a[!small]) * exp(log_erfcx)
  }
  out
}

model2_fixed_moments <- local({
  cache <- new.env(parent = emptyenv())
  function() {
    if (exists("moments", envir = cache, inherits = FALSE)) {
      return(get("moments", envir = cache, inherits = FALSE))
    }
    joint_expectation <- function(fun) {
      stats::integrate(function(x1) {
        vapply(x1, function(x1_value) {
          stats::integrate(function(x2) fun(x1_value, x2) / 4, -2, 2,
                           rel.tol = 1e-11)$value
        }, numeric(1)) * stats::dbeta(x1, 3, 4)
      }, 0, 1, rel.tol = 1e-11)$value
    }
    ans <- c(
      exp_moment = joint_expectation(
        function(x1, x2) exp((x1 + x2)^2 / 18)
      ),
      reciprocal_square_moment = joint_expectation(
        function(x1, x2) normal_square_reciprocal_moment(abs(x1 - x2) / 3)
      )
    )
    assign("moments", ans, envir = cache)
    ans
  }
})

population_tau <- function(model) {
  z <- model_parameters(model)
  if (model == "Model1") {
    return(z$alpha1 - z$alpha0 + (z$beta1[1] - z$beta0[1]) * 3 / 7 +
             (z$beta1[4] - z$beta0[4]) * 3.8)
  }
  if (model == "Model2") {
    moments <- model2_fixed_moments()
    elog <- stats::integrate(function(x) log1p(x) * stats::dbeta(x, 3, 4),
                             0, 1, rel.tol = 1e-12)$value
    eexp <- stats::integrate(function(x) exp(x + 2) * stats::dbeta(x, 3, 4),
                             0, 1, rel.tol = 1e-12)$value
    erecip <- stats::integrate(function(x) stats::dbeta(x, 3, 4) / (x + 1),
                               0, 1, rel.tol = 1e-12)$value
    ex2sq <- 4 / 3
    return(z$alpha1 - z$alpha0 + z$beta1[1] * eexp + z$beta1[2] * erecip +
             z$beta1[3] * ex2sq - z$beta0[1] * elog - z$beta0[2] * ex2sq -
             z$beta0[3] * moments[["exp_moment"]] -
             z$beta0[4] * moments[["reciprocal_square_moment"]])
  }
  if (model == "Model3") {
    term0 <- stats::integrate(function(x) {
      cval <- x + 2
      inner <- x * (1 - cval / 4 * log((cval + 2) / pmax(cval - 2, .Machine$double.xmin)))
      inner * stats::dbeta(x, 3, 4)
    }, 0, 1, subdivisions = 2000L, rel.tol = 1e-11)$value
    exp_beta <- stats::integrate(function(x) exp(-x) * stats::dbeta(x, 3, 4),
                                 0, 1, rel.tol = 1e-12)$value
    term1 <- (4 / 3) * exp(-2) * exp_beta
    return(z$alpha1 - z$alpha0 + z$beta1[1] + z$beta1[2] * term1 - z$beta0[1] * term0)
  }
  if (model == "Model4") {
    elog <- stats::integrate(function(x) log1p(x) * stats::dbeta(x, 3, 4),
                             0, 1, rel.tol = 1e-12)$value
    eexp <- (exp(2) - exp(-2)) / 4
    return(z$alpha1 - z$alpha0 + 0.5 * z$beta1[3] * eexp -
             0.5 * z$beta0[3] * elog)
  }
  stop("Unknown DGP: ", model)
}

stratification_frame <- function(S) {
  if (is.data.frame(S)) out <- S
  else if (is.matrix(S)) out <- as.data.frame(S)
  else out <- data.frame(S = S)
  if (!nrow(out) || anyNA(out)) stop("Stratification variables must be nonempty and complete")
  names(out) <- paste0("S", seq_len(ncol(out)))
  out
}

stratum_factor <- function(S) {
  frame <- stratification_frame(S)
  do.call(interaction, c(lapply(frame, as.factor),
                         list(drop = TRUE, lex.order = TRUE)))
}

assignment_counts <- function(S, A) {
  table(stratum_factor(S), factor(as.integer(A), levels = 0:1),
        dnn = c("stratum", "arm"))
}

assignment_is_valid <- function(S, A, minimum_cell = 10L) {
  counts <- assignment_counts(S, A)
  ncol(counts) == 2L && all(counts >= minimum_cell)
}

sbr_assignment <- function(S, block_size = 6L) {
  if (block_size != 6L) stop("The SBR design uses block size 6")
  key <- stratum_factor(S)
  A <- integer(length(key))
  template <- c(rep(0L, block_size / 2L), rep(1L, block_size / 2L))
  for (level in levels(key)) {
    idx <- which(key == level)
    for (start in seq.int(1L, length(idx), by = block_size)) {
      block_idx <- idx[start:min(start + block_size - 1L, length(idx))]
      A[block_idx] <- sample(template, block_size, replace = FALSE)[seq_along(block_idx)]
    }
  }
  A
}

pocock_simon_scores <- function(S, A_previous, subject_index) {
  frame <- stratification_frame(S)
  if (subject_index < 1L || subject_index > nrow(frame)) stop("Invalid subject index")
  if (length(A_previous) != subject_index - 1L) stop("A_previous has the wrong length")
  scores <- numeric(2L)
  for (candidate in 0:1) {
    total <- 0
    for (j in seq_len(ncol(frame))) {
      previous_same <- if (subject_index == 1L) logical(0) else
        frame[seq_len(subject_index - 1L), j] == frame[subject_index, j]
      imbalance <- if (subject_index == 1L || !any(previous_same)) 0L else
        sum(ifelse(A_previous[previous_same] == 1L, 1L, -1L))
      total <- total + abs(imbalance + if (candidate == 1L) 1L else -1L)
    }
    scores[candidate + 1L] <- total
  }
  names(scores) <- c("control", "treatment")
  scores
}

pocock_simon_assignment <- function(S, biased_coin = 0.75) {
  if (!isTRUE(all.equal(as.numeric(biased_coin), 0.75))) {
    stop("The Pocock--Simon design uses biased-coin probability 0.75")
  }
  n <- nrow(stratification_frame(S))
  A <- integer(n)
  for (i in seq_len(n)) {
    scores <- pocock_simon_scores(S, A[seq_len(i - 1L)], i)
    treatment_probability <- if (scores[["treatment"]] < scores[["control"]]) {
      biased_coin
    } else if (scores[["treatment"]] > scores[["control"]]) {
      1 - biased_coin
    } else {
      0.5
    }
    A[i] <- stats::rbinom(1L, 1L, treatment_probability)
  }
  A
}

draw_assignment <- function(S, design = "SRS", seed = NULL, pi = 0.5,
                            block_size = 6L, biased_coin = 0.75,
                            minimum_cell = 10L) {
  design <- match.arg(design, SIMULATION_DESIGNS)
  if (!is.null(seed)) set.seed(as.integer(seed))
  for (attempt in seq_len(1000L)) {
    A <- switch(design,
      SRS = stats::rbinom(nrow(stratification_frame(S)), 1L, pi),
      SBR = sbr_assignment(S, block_size = block_size),
      PS = pocock_simon_assignment(S, biased_coin = biased_coin)
    )
    if (assignment_is_valid(S, A, minimum_cell = minimum_cell)) return(as.integer(A))
  }
  stop("Could not draw a valid ", design, " assignment with non-small stratum-arm cells")
}

replication_seeds <- function(model, n, rep_id, design) {
  model_index <- match(model, paste0("Model", 1:4))
  design_index <- match(design, SIMULATION_DESIGNS)
  if (is.na(model_index) || is.na(design_index)) stop("Unknown model or design")
  data_seed <- seed_int(20260812, model_index, as.integer(n), as.integer(rep_id), 11L)
  assignment_seed <- seed_int(data_seed, design_index, 17L)
  list(data = data_seed, assignment = assignment_seed)
}

generate_dgp <- function(model, n = 1000L, p = 30L, design = "SRS",
                         data_seed = NULL, assignment_seed = NULL) {
  if (p != 30L) stop("The simulation uses p=30")
  design <- match.arg(design, SIMULATION_DESIGNS)
  if (!is.null(data_seed)) set.seed(as.integer(data_seed))
  z <- model_parameters(model)
  X1 <- stats::rbeta(n, 3, 4)
  X2 <- stats::runif(n, -2, 2)

  if (model == "Model1") {
    X3 <- sample(c(1, -1), n, replace = TRUE)
    X4 <- sample(c(3, 5), n, replace = TRUE, prob = c(0.6, 0.4))
    Xadd <- MASS::mvrnorm(n, rep(0, p - 4L), equi_cov(p - 4L, 0.2))
    X <- cbind(X1, X2, X3, X4, Xadd)
    S <- sample(1:4, n, replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2))
    g0 <- z$alpha0 + X[, 1:4, drop = FALSE] %*% z$beta0
    g1 <- z$alpha1 + X[, 1:4, drop = FALSE] %*% z$beta1
    eps0 <- stats::rnorm(n, sd = z$sigma0)
    eps1 <- stats::rnorm(n, sd = z$sigma1)
  } else if (model == "Model2") {
    Z3 <- stats::rnorm(n)
    Z4 <- stats::rnorm(n)
    X3 <- (X1 + X2) * Z3 / 3
    X4 <- (X1 - X2) * Z4 / 3
    Xrest <- MASS::mvrnorm(n, rep(0, p - 4L), equi_cov(p - 4L, 0.2))
    Xrest[, 1:5] <- Xrest[, 1:5, drop = FALSE] * X1
    Xrest[, 6:10] <- Xrest[, 6:10, drop = FALSE] * X2
    X <- cbind(X1, X2, X3, X4, Xrest)
    S <- sample(1:4, n, replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2))
    g0 <- z$alpha0 + z$beta0[1] * log1p(X[, 1]) + z$beta0[2] * X[, 2]^2 +
      z$beta0[3] * exp(X[, 3]) + z$beta0[4] / (X[, 4]^2 + 3)
    g1 <- z$alpha1 + z$beta1[1] * exp(X[, 1] + 2) +
      z$beta1[2] / (X[, 1] + 1) + z$beta1[3] * X[, 2]^2
    eps0 <- stats::rnorm(n, sd = z$sigma0)
    eps1 <- stats::rnorm(n, sd = z$sigma1)
  } else if (model == "Model3") {
    X3 <- stats::rnorm(n)
    X4 <- stats::runif(n, 0, 2)
    Xadd <- MASS::mvrnorm(n, rep(0, p - 4L), toeplitz_cov(p - 4L, 0.5))
    X <- cbind(X1, X2, X3, X4, Xadd)
    S <- sample(1:2, n, replace = TRUE, prob = c(0.4, 0.6))
    g0 <- z$alpha0 + z$beta0[1] * X[, 1] * X[, 2] /
      (X[, 1] + X[, 2] + 2) + z$beta0[2] * X[, 1]^2 * (X[, 2] + X[, 3])
    g1 <- z$alpha1 + z$beta1[1] * (X[, 2] + X[, 4]) +
      z$beta1[2] * X[, 2]^2 / exp(X[, 1] + 2)
    eps0 <- stats::rnorm(n, sd = z$sigma0)
    eps1 <- stats::rnorm(n, sd = z$sigma1)
  } else if (model == "Model4") {
    Xadd <- MASS::mvrnorm(n, rep(0, p - 2L), equi_cov(p - 2L, 0.2))
    selected <- sample(seq_len(p - 2L), floor(p / 3), replace = FALSE)
    source <- sample(1:2, length(selected), replace = TRUE)
    for (j in seq_along(selected)) {
      Xadd[, selected[j]] <- Xadd[, selected[j]] * if (source[j] == 1L) X1 else X2
    }
    X <- cbind(X1, X2, Xadd)
    S <- sample(c(1, -1), n, replace = TRUE)
    g0 <- z$alpha0 + (z$beta0[1] * X[, 1] + z$beta0[2] * X[, 2]) * S +
      z$beta0[3] * log1p(X[, 1]) * (S == 1)
    g1 <- z$alpha1 + (z$beta1[1] * X[, 1] + z$beta1[2] * X[, 2]) * S +
      z$beta1[3] * exp(X[, 2]) * (S == -1)
    eps0 <- stats::rnorm(n, sd = z$sigma0)
    eps1 <- stats::rnorm(n, sd = z$sigma1)
  } else {
    stop("Unknown DGP: ", model)
  }

  colnames(X) <- paste0("X", seq_len(ncol(X)))
  Y0 <- as.numeric(g0 + eps0)
  Y1 <- as.numeric(g1 + eps1)
  A <- draw_assignment(S, design = design, seed = assignment_seed)
  pop_tau <- population_tau(model)
  target <- pop_tau
  target_type <- "population_ate"
  list(
    model = model, design = design, A = as.integer(A), S = S, X = as.data.frame(X),
    Y0 = Y0, Y1 = Y1, Y = ifelse(A == 1L, Y1, Y0),
    g0 = as.numeric(g0), g1 = as.numeric(g1),
    target = as.numeric(target), target_type = target_type,
    population_tau = pop_tau, sample_ate = mean(Y1 - Y0),
    conditional_sample_ate = mean(g1 - g0),
    data_seed = if (is.null(data_seed)) NA_integer_ else as.integer(data_seed),
    assignment_seed = if (is.null(assignment_seed)) NA_integer_ else as.integer(assignment_seed)
  )
}

balanced_folds <- function(S, A, folds, seed) {
  if (folds < 2L) stop("folds must be at least 2")
  out <- integer(length(S))
  cells <- interaction(S, A, drop = TRUE, lex.order = TRUE)
  for (cell in levels(cells)) {
    idx <- which(cells == cell)
    if (length(idx) < folds) stop("Cell too small for balanced folds: ", cell)
    set.seed(seed_int(seed, match(cell, levels(cells))))
    idx <- sample(idx, length(idx), replace = FALSE)
    out[idx] <- rep(seq_len(folds), length.out = length(idx))
  }
  out
}

fit_rf_predict <- function(x, y, newx, seed) {
  set.seed(seed)
  fit <- randomForest::randomForest(x = as.data.frame(x), y = as.numeric(y), ntree = 50)
  as.numeric(stats::predict(fit, as.data.frame(newx)))
}

training_outcome_scale <- function(y) {
  y <- as.numeric(y)
  if (!length(y) || any(!is.finite(y))) stop("Training outcomes must be finite and nonempty")
  list(ymin = min(y), yrange = max(y) - min(y))
}

fit_nn_predict <- function(x, y, newx, seed) {
  y <- as.numeric(y)
  scale <- training_outcome_scale(y)
  if (scale$yrange <= 1e-12) return(rep(mean(y), nrow(newx)))
  train <- data.frame(y_scaled = (y - scale$ymin) / scale$yrange, as.data.frame(x))
  newdata <- as.data.frame(newx)
  set.seed(seed)
  fit <- neuralnet::neuralnet(y_scaled ~ ., data = train, hidden = 2,
                              linear.output = TRUE, lifesign = "none")
  pred <- as.numeric(stats::predict(fit, newdata = newdata)) * scale$yrange + scale$ymin
  if (length(pred) != nrow(newx) || any(!is.finite(pred))) stop("Non-finite NN prediction")
  pred
}

fit_prediction_cube <- function(train, newx, strata, seed_base) {
  nnew <- nrow(newx)
  K <- length(strata)
  rf <- array(NA_real_, dim = c(nnew, K, 2L),
              dimnames = list(NULL, as.character(strata), c("0", "1")))
  nn <- rf
  for (k in seq_along(strata)) for (a in 0:1) {
    idx <- which(train$S == strata[k] & train$A == a)
    if (length(idx) < 5L) stop("Training cell has fewer than five observations")
    cell_mean <- mean(train$Y[idx])
    rf_pred <- tryCatch(
      fit_rf_predict(train$X[idx, , drop = FALSE], train$Y[idx], newx,
                     seed_int(seed_base, k, a, 1)),
      error = function(e) rep(cell_mean, nnew)
    )
    nn_pred <- tryCatch(
      fit_nn_predict(train$X[idx, , drop = FALSE], train$Y[idx], newx,
                     seed_int(seed_base, k, a, 2)),
      error = function(e) rep(cell_mean, nnew)
    )
    rf[, k, a + 1L] <- rf_pred
    nn[, k, a + 1L] <- nn_pred
  }
  list(rf = rf, nn = nn)
}

fit_lm_prediction_cube <- function(train, newx, strata) {
  new_df <- as.data.frame(newx)
  nnew <- nrow(newx)
  K <- length(strata)
  lm <- array(NA_real_, dim = c(nnew, K, 2L),
              dimnames = list(NULL, as.character(strata), c("0", "1")))
  for (k in seq_along(strata)) for (a in 0:1) {
    idx <- which(train$S == strata[k] & train$A == a)
    if (!length(idx)) stop("Empty OLS training stratum-arm cell")
    train_df <- data.frame(y = train$Y[idx], train$X[idx, , drop = FALSE])
    fit <- stats::lm(y ~ ., data = train_df, singular.ok = TRUE)
    pred <- suppressWarnings(as.numeric(stats::predict(fit, newdata = new_df)))
    if (length(pred) != nnew || any(!is.finite(pred))) {
      stop("Non-finite OLS prediction in stratum ", strata[k], ", arm ", a)
    }
    lm[, k, a + 1L] <- pred
  }
  list(lm = lm)
}

make_oof_cube <- function(train, strata, inner_folds, seed_base) {
  folds <- balanced_folds(train$S, train$A, inner_folds, seed_int(seed_base, 91))
  n <- nrow(train$X)
  K <- length(strata)
  rf <- array(NA_real_, c(n, K, 2L), dimnames = list(NULL, as.character(strata), c("0", "1")))
  nn <- rf
  for (j in seq_len(inner_folds)) {
    val <- which(folds == j)
    tr <- which(folds != j)
    cube <- fit_prediction_cube(
      list(X = train$X[tr, , drop = FALSE], Y = train$Y[tr],
           S = train$S[tr], A = train$A[tr]),
      train$X[val, , drop = FALSE], strata, seed_int(seed_base, 100 + j)
    )
    rf[val, , ] <- cube$rf
    nn[val, , ] <- cube$nn
  }
  if (anyNA(rf) || anyNA(nn)) stop("OOF cube is incomplete")
  list(rf = rf, nn = nn)
}

target_source_index <- function(S, strata) match(as.character(S), as.character(strata))

library_predictions <- function(cube, S, strata, library, arm) {
  aidx <- as.integer(arm) + 1L
  if (library == "rfnn") {
    target <- target_source_index(S, strata)
    rows <- seq_along(S)
    return(cbind(rf = cube$rf[cbind(rows, target, aidx)],
                 nn = cube$nn[cbind(rows, target, aidx)]))
  }
  if (library == "rf") {
    target <- target_source_index(S, strata)
    rows <- seq_along(S)
    return(matrix(cube$rf[cbind(rows, target, aidx)], ncol = 1L,
                  dimnames = list(NULL, "rf")))
  }
  if (library == "nn") {
    target <- target_source_index(S, strata)
    rows <- seq_along(S)
    return(matrix(cube$nn[cbind(rows, target, aidx)], ncol = 1L,
                  dimnames = list(NULL, "nn")))
  }
  if (library == "rflin") {
    if (is.null(cube$lm)) stop("RF+lin library requires OLS predictions")
    target <- target_source_index(S, strata)
    rows <- seq_along(S)
    return(cbind(rf = cube$rf[cbind(rows, target, aidx)],
                 lm = cube$lm[cbind(rows, target, aidx)]))
  }
  if (library == "grouped_rf") return(cube$rf[, , aidx, drop = FALSE][, , 1L])
  if (library == "grouped_nn") return(cube$nn[, , aidx, drop = FALSE][, , 1L])
  stop("Unknown library: ", library)
}

calibration_proxy <- function(cube, S, strata, library) {
  cbind(library_predictions(cube, S, strata, library, 0),
        library_predictions(cube, S, strata, library, 1))
}

all_subsets <- function(p) {
  lapply(seq_len(2^p - 1L), function(mask) which(as.logical(intToBits(mask))[seq_len(p)]))
}

fit_meta_weights <- function(X, y, type = c("sl", "stack"), tol = 1e-10) {
  type <- match.arg(type)
  X <- as.matrix(X)
  y <- as.numeric(y)
  p <- ncol(X)
  if (!length(y) || nrow(X) != length(y) || !p || any(!is.finite(c(X, y)))) {
    stop("Meta-learning inputs must be finite and have compatible dimensions")
  }
  scale <- max(abs(c(X, y)))
  if (scale > 0) {
    X <- X / scale
    y <- y / scale
  }
  best_w <- if (type == "stack") rep(0, p) else rep(1 / p, p)
  best_loss <- sum((y - X %*% best_w)^2)
  for (J in all_subsets(p)) {
    XJ <- X[, J, drop = FALSE]
    if (type == "stack") {
      wJ <- as.numeric(pinv_svd(XJ)$inverse %*% y)
    } else if (length(J) == 1L) {
      wJ <- 1
    } else {
      reference <- XJ[, length(J)]
      contrasts <- sweep(XJ[, -length(J), drop = FALSE], 1L, reference)
      coefficients <- as.numeric(pinv_svd(contrasts)$inverse %*% (y - reference))
      wJ <- c(coefficients, 1 - sum(coefficients))
    }
    if (all(is.finite(wJ)) && all(wJ >= -tol)) {
      w <- rep(0, p)
      w[J] <- pmax(wJ, 0)
      if (type == "sl") {
        sw <- sum(w)
        if (sw <= 0) next
        w <- w / sw
      }
      loss <- sum((y - X %*% w)^2)
      if (is.finite(loss) && loss < best_loss) {
        best_loss <- loss
        best_w <- w
      }
    }
  }
  list(weights = best_w)
}

learn_meta <- function(train, oof_cube, strata, library, type) {
  K <- length(strata)
  pmeta <- ncol(library_predictions(oof_cube, train$S, strata, library, 0))
  weights <- array(NA_real_, c(K, 2L, pmeta),
                   dimnames = list(as.character(strata), c("0", "1"), NULL))
  for (k in seq_along(strata)) for (a in 0:1) {
    idx <- which(train$S == strata[k] & train$A == a)
    Xmeta <- library_predictions(oof_cube, train$S, strata, library, a)[idx, , drop = FALSE]
    fit <- fit_meta_weights(Xmeta, train$Y[idx], type = type)
    weights[k, a + 1L, ] <- fit$weights
  }
  weights
}

apply_meta <- function(cube, S, strata, library, arm, weights) {
  Xmeta <- library_predictions(cube, S, strata, library, arm)
  target <- target_source_index(S, strata)
  out <- numeric(length(S))
  for (k in seq_along(strata)) {
    idx <- which(target == k)
    out[idx] <- as.numeric(Xmeta[idx, , drop = FALSE] %*% weights[k, arm + 1L, ])
  }
  out
}

center_columns <- function(x) {
  x <- as.matrix(x)
  if (!nrow(x) || any(!is.finite(x))) stop("Calibration proxies must be finite and nonempty")
  centered <- sweep(x, 2L, colMeans(x))
  constant <- vapply(seq_len(ncol(x)), function(j) all(x[, j] == x[1L, j]), logical(1))
  centered[, constant] <- 0
  centered
}

build_calibration_matrix <- function(A, S, xi, strata) {
  xi <- as.matrix(xi)
  m <- length(A)
  d <- ncol(xi)
  G <- matrix(0, m, d * length(strata))
  if (!d) return(G)
  for (k in seq_along(strata)) {
    idx <- which(S == strata[k])
    pi_k <- mean(A[idx])
    centered <- center_columns(xi[idx, , drop = FALSE])
    cols <- ((k - 1L) * d + 1L):(k * d)
    G[idx, cols] <- centered * (A[idx] - pi_k)
  }
  G
}

scale_columns <- function(x) {
  scales <- apply(abs(x), 2L, max)
  scales[scales == 0] <- 1
  sweep(x, 2L, scales, "/")
}

effective_basis <- function(G, tol = SVD_TOL) {
  if (!length(G)) return(list(B = matrix(0, nrow(G), 0L)))
  scaled <- scale_columns(G)
  sg <- svd(scaled, nu = min(dim(scaled)), nv = 0)
  cutoff <- tol * max(1, if (length(sg$d)) max(sg$d) else 0)
  keep <- sg$d > cutoff
  list(B = sg$u[, keep, drop = FALSE])
}

dq_u <- function(v, tol = 1e-13, maxit = 80L) {
  u <- -sign(v) * abs(v)^(1 / 3)
  u[!is.finite(u)] <- 0
  for (iter in seq_len(maxit)) {
    f <- u - u^2 + u^3 + v
    if (max(abs(f)) < tol) break
    deriv <- 1 - 2 * u + 3 * u^2
    u <- u - f / deriv
  }
  list(u = u)
}

solve_calibration <- function(G, discrepancy = c("quadratic", "quartic"),
                              tol = BALANCE_TOL, maxit = 100L) {
  discrepancy <- match.arg(discrepancy)
  basis <- effective_basis(G)
  B <- basis$B
  if (!ncol(B)) {
    w <- rep(1, nrow(G))
    return(list(weights = w))
  }
  if (discrepancy == "quadratic") {
    lambda <- as.numeric(crossprod(B, rep(1, nrow(B))))
    w <- as.numeric(1 - B %*% lambda)
  } else {
    dual_tol <- min(1e-12, tol * 1e-4)
    lambda <- as.numeric(crossprod(B, rep(1, nrow(B))))
    converged <- FALSE
    for (iterations in seq_len(maxit)) {
      root <- dq_u(as.numeric(B %*% lambda))
      w <- 1 + root$u
      grad <- as.numeric(crossprod(B, w))
      if (max(abs(grad)) < dual_tol) {
        converged <- TRUE
        break
      }
      deriv <- 1 - 2 * root$u + 3 * root$u^2
      H <- -crossprod(B, B / deriv)
      step <- -as.numeric(pinv_svd(H)$inverse %*% grad)
      gradient_norm <- sqrt(sum(grad^2))
      accepted <- FALSE
      for (ls in 0:30) {
        candidate <- lambda + (0.5^ls) * step
        cand_root <- dq_u(as.numeric(B %*% candidate))
        cand_grad <- as.numeric(crossprod(B, 1 + cand_root$u))
        if (sqrt(sum(cand_grad^2)) < gradient_norm) {
          lambda <- candidate
          accepted <- TRUE
          break
        }
      }
      if (!accepted) stop("Quartic dual Newton line search failed")
    }
    root <- dq_u(as.numeric(B %*% lambda))
    w <- 1 + root$u
    if (!converged && max(abs(crossprod(B, w))) >= dual_tol) stop("Quartic dual Newton did not converge")
  }
  balance <- if (ncol(G)) max(abs(colMeans(scale_columns(G) * w))) else 0
  if (!is.finite(balance) || balance > tol) stop("Calibration balance residual exceeds tolerance: ", balance)
  list(weights = w)
}

fold_outcome_objects <- function(eval) {
  strata <- sort(unique(eval$S))
  m <- length(eval$Y)
  tau_k <- pi_map <- ybar0 <- ybar1 <- setNames(numeric(length(strata)), as.character(strata))
  r <- numeric(m)
  for (k in seq_along(strata)) {
    idx <- which(eval$S == strata[k])
    i0 <- idx[eval$A[idx] == 0L]
    i1 <- idx[eval$A[idx] == 1L]
    if (!length(i0) || !length(i1)) stop("Empty evaluation stratum-arm cell")
    pi_map[k] <- length(i1) / length(idx)
    ybar0[k] <- mean(eval$Y[i0])
    ybar1[k] <- mean(eval$Y[i1])
    tau_k[k] <- ybar1[k] - ybar0[k]
    r[idx] <- eval$A[idx] / pi_map[k] * (eval$Y[idx] - ybar1[k]) -
      (1 - eval$A[idx]) / (1 - pi_map[k]) * (eval$Y[idx] - ybar0[k])
  }
  list(strata = strata, tau_k = tau_k, r = r)
}

calibration_variance_fold <- function(eval, xi, strata = sort(unique(eval$S)),
                                      tol = SVD_TOL) {
  xi <- as.matrix(xi)
  m <- length(eval$Y)
  if (nrow(xi) != m) stop("Calibration variance proxy has the wrong row count")
  tau_sdim <- mean(fold_outcome_objects(eval)$tau_k[as.character(eval$S)])
  adjusted <- 0

  for (k in seq_along(strata)) {
    idx <- which(eval$S == strata[k])
    i1 <- idx[eval$A[idx] == 1L]
    i0 <- idx[eval$A[idx] == 0L]
    if (length(i1) <= 1L || length(i0) <= 1L) {
      stop("Calibration variance requires at least two observations per stratum-arm cell")
    }
    nk <- length(idx)
    pi_k <- length(i1) / nk
    pnk <- nk / m
    ybar1k <- mean(eval$Y[i1])
    ybar0k <- mean(eval$Y[i0])
    v1 <- pnk / (1 - pi_k) * stats::var(eval$Y[i0]) +
      pnk / pi_k * stats::var(eval$Y[i1])
    v2 <- pnk * (ybar1k - ybar0k - tau_sdim)^2

    centered <- center_columns(xi[idx, , drop = FALSE])
    Xi <- centered * (eval$A[idx] - pi_k)
    B <- effective_basis(Xi, tol)$B
    residual_scores <- eval$A[idx] / pi_k * (eval$Y[idx] - ybar1k) -
      (1 - eval$A[idx]) / (1 - pi_k) * (eval$Y[idx] - ybar0k)
    projected <- as.numeric(B %*% crossprod(B, residual_scores))
    v3 <- sum(projected^2) / m
    component <- as.numeric(v1 + v2 - v3)
    component_scale <- max(1, abs(v1) + abs(v2) + abs(v3))
    if (!is.finite(component) || component < -1e-10 * component_scale) {
      stop("Calibration variance component is materially negative")
    }
    component <- max(component, 0)

    dof <- ncol(B) + 1L
    if (nk <= dof) stop("Nonpositive residual degrees of freedom in a stratum")
    factor <- nk / (nk - dof)
    adjusted <- adjusted + factor * component
  }

  list(variance = adjusted / m)
}

combine_fold_variances <- function(fold_variances, fold_sizes) {
  fold_variances <- as.numeric(fold_variances)
  fold_sizes <- as.numeric(fold_sizes)
  if (!length(fold_variances) || length(fold_variances) != length(fold_sizes) ||
      any(!is.finite(fold_variances)) || any(fold_variances < 0) ||
      any(!is.finite(fold_sizes)) || any(fold_sizes <= 0)) {
    stop("Invalid fold variances or fold sizes")
  }
  weights <- fold_sizes / sum(fold_sizes)
  sum(weights^2 * fold_variances)
}

calibration_score_adjustment <- function(A, S, xi, r, strata, tol = SVD_TOL) {
  xi <- as.matrix(xi)
  adjustment <- numeric(length(A))
  for (k in seq_along(strata)) {
    idx <- which(S == strata[k])
    pi_k <- mean(A[idx])
    centered <- center_columns(xi[idx, , drop = FALSE])
    Xi <- centered * (A[idx] - pi_k)
    B <- effective_basis(Xi, tol)$B
    adjustment[idx] <- as.numeric(B %*% crossprod(B, r[idx]))
  }
  adjustment
}

calibration_fold <- function(eval, xi, discrepancy) {
  obj <- fold_outcome_objects(eval)
  G <- build_calibration_matrix(eval$A, eval$S, xi, obj$strata)
  sol <- solve_calibration(G, discrepancy)
  str_index <- match(as.character(eval$S), names(obj$tau_k))
  exact_contributions <- as.numeric(obj$tau_k[str_index]) + sol$weights * obj$r
  adjustment <- calibration_score_adjustment(eval$A, eval$S, xi, obj$r, obj$strata)
  variance <- calibration_variance_fold(eval, xi, obj$strata)
  influence_scores <- as.numeric(obj$tau_k[str_index]) + obj$r - adjustment
  scores <- influence_scores + mean(exact_contributions) - mean(influence_scores)
  list(scores = scores, estimate = mean(scores), variance = variance)
}

aipw_scores <- function(eval, m0, m1) {
  strata <- sort(unique(eval$S))
  pi_obs <- numeric(length(eval$Y))
  for (s in strata) {
    idx <- which(eval$S == s)
    pi_obs[idx] <- mean(eval$A[idx])
  }
  m1 - m0 + eval$A / pi_obs * (eval$Y - m1) -
    (1 - eval$A) / (1 - pi_obs) * (eval$Y - m0)
}

car_asymptotic_variance <- function(A, S, R, strata = sort(unique(S))) {
  A <- as.integer(A)
  R <- as.numeric(R)
  pait <- mean(A)
  if (!is.finite(pait) || pait <= 0 || pait >= 1) stop("Invalid realized treatment proportion")
  r0 <- r1 <- nk0 <- nk1 <- sr0 <- sr1 <- numeric(length(strata))
  for (k in seq_along(strata)) {
    idx <- which(S == strata[k])
    i0 <- idx[A[idx] == 0L]
    i1 <- idx[A[idx] == 1L]
    if (length(i0) < 2L || length(i1) < 2L) {
      stop("CAR variance requires at least two observations per stratum-arm cell")
    }
    nk0[k] <- length(i0)
    nk1[k] <- length(i1)
    r0[k] <- mean(R[i0])
    r1[k] <- mean(R[i1])
    sr0[k] <- stats::var(R[i0])
    sr1[k] <- stats::var(R[i1])
  }
  pn <- (nk0 + nk1) / sum(nk0 + nk1)
  r1 <- r1 - mean(R[A == 1L])
  r0 <- r0 - mean(R[A == 0L])
  out <- sum(pn * sr1) / pait + sum(pn * sr0) / (1 - pait) +
    sum(pn * (r1 - r0)^2)
  if (!is.finite(out) || out < -1e-12) stop("Invalid CAR asymptotic variance")
  out <- max(out, 0)
  out
}

estimate_car <- function(eval, m0 = NULL, m1 = NULL) {
  n <- length(eval$Y)
  if (is.null(m0)) m0 <- rep(0, n)
  if (is.null(m1)) m1 <- rep(0, n)
  m0 <- as.numeric(m0)
  m1 <- as.numeric(m1)
  if (length(m0) != n || length(m1) != n || any(!is.finite(c(m0, m1)))) {
    stop("Invalid outcome-regression predictions")
  }
  pi_obs <- numeric(n)
  strata <- sort(unique(eval$S))
  for (s in strata) {
    idx <- which(eval$S == s)
    pi_obs[idx] <- mean(eval$A[idx])
  }
  residual <- eval$Y - ((1 - pi_obs) * m1 + pi_obs * m0)
  residual_eval <- list(Y = residual, S = eval$S, A = eval$A)
  scores <- aipw_scores(eval, m0, m1)
  obj <- fold_outcome_objects(residual_eval)
  estimate <- sum(vapply(strata, function(s) {
    idx <- which(eval$S == s)
    length(idx) / n * obj$tau_k[as.character(s)]
  }, numeric(1)))
  if (!isTRUE(all.equal(mean(scores), estimate, tolerance = 1e-10))) {
    stop("AIPW score and residualized CAR estimand disagree")
  }
  asymptotic_variance <- car_asymptotic_variance(eval$A, eval$S, residual, strata)
  variance <- as.numeric(asymptotic_variance) / n
  list(scores = scores, estimate = estimate, variance = variance)
}

method_catalog <- function() {
  methods <- c(
    "cal_quadratic__rf", "cal_quartic__rf", "aipw_rf",
    "cal_quadratic__nn", "cal_quartic__nn", "aipw_nn",
    "cal_quadratic__rfnn", "cal_quartic__rfnn", "aipw_sl__rfnn", "aipw_stack__rfnn",
    "cal_quadratic__grouped_rf", "cal_quartic__grouped_rf",
    "aipw_sl__grouped_rf", "aipw_stack__grouped_rf",
    "cal_quadratic__grouped_nn", "cal_quartic__grouped_nn",
    "aipw_sl__grouped_nn", "aipw_stack__grouped_nn",
    "sdim", "lin", "cal_quadratic__rflin", "cal_quartic__rflin"
  )
  estimator_class <- c(
    "cal_quadratic", "cal_quartic", "aipw_rf",
    "cal_quadratic", "cal_quartic", "aipw_nn",
    "cal_quadratic", "cal_quartic", "aipw_sl", "aipw_stack",
    "cal_quadratic", "cal_quartic", "aipw_sl", "aipw_stack",
    "cal_quadratic", "cal_quartic", "aipw_sl", "aipw_stack",
    "sdim", "lin", "cal_quadratic", "cal_quartic"
  )
  library <- c(
    "rf", "rf", "rf", "nn", "nn", "nn",
    rep("rfnn", 4L), rep("grouped_rf", 4L), rep("grouped_nn", 4L),
    "sdim", "lin", "rflin", "rflin"
  )
  cross_fitted <- !methods %in% c("sdim", "lin")
  data.frame(method = methods, estimator_class = estimator_class, library = library,
             cross_fitted = cross_fitted, stringsAsFactors = FALSE)
}

run_replication <- function(model, rep_id, design = "SRS", n = 500L,
                            outer_folds = 2L, inner_folds = 5L) {
  require_simulation_packages()
  design <- match.arg(design, SIMULATION_DESIGNS)
  catalog <- method_catalog()
  seeds <- replication_seeds(model, n, rep_id, design)
  data_seed <- seeds$data
  assignment_seed <- seeds$assignment
  dat <- generate_dgp(model, n = n, p = 30L, design = design,
                      data_seed = data_seed, assignment_seed = assignment_seed)
  strata <- sort(unique(dat$S))
  outer <- balanced_folds(dat$S, dat$A, outer_folds,
                          seed_int(data_seed, assignment_seed, 2L))
  scores <- setNames(lapply(catalog$method, function(x) rep(NA_real_, n)), catalog$method)
  fold_variances <- setNames(
    lapply(catalog$method, function(x) rep(NA_real_, outer_folds)), catalog$method)
  full_variances <- setNames(rep(NA_real_, nrow(catalog)), catalog$method)

  full_eval <- list(X = dat$X, Y = dat$Y, S = dat$S, A = dat$A)
  sdim <- estimate_car(full_eval)
  scores$sdim <- sdim$scores
  full_variances[["sdim"]] <- sdim$variance

  full_lm <- tryCatch(fit_lm_prediction_cube(full_eval, dat$X, strata),
                      error = function(e) e)
  if (!inherits(full_lm, "error")) {
    target <- target_source_index(dat$S, strata)
    rows_index <- seq_len(n)
    lin0 <- full_lm$lm[cbind(rows_index, target, 1L)]
    lin1 <- full_lm$lm[cbind(rows_index, target, 2L)]
    lin <- estimate_car(full_eval, lin0, lin1)
    scores$lin <- lin$scores
    full_variances[["lin"]] <- lin$variance
  }

  for (f in seq_len(outer_folds)) {
    eval_idx <- which(outer == f)
    train_idx <- which(outer != f)
    train <- list(X = dat$X[train_idx, , drop = FALSE], Y = dat$Y[train_idx],
                  S = dat$S[train_idx], A = dat$A[train_idx])
    eval <- list(X = dat$X[eval_idx, , drop = FALSE], Y = dat$Y[eval_idx],
                 S = dat$S[eval_idx], A = dat$A[eval_idx])
    fold_seed <- seed_int(data_seed, assignment_seed, f, 300L)
    outer_cube <- fit_prediction_cube(train, eval$X, strata, seed_int(fold_seed, 1))
    outer_lm <- tryCatch(fit_lm_prediction_cube(train, eval$X, strata),
                         error = function(e) e)
    if (!inherits(outer_lm, "error")) outer_cube$lm <- outer_lm$lm
    oof_cube <- make_oof_cube(train, strata, inner_folds, seed_int(fold_seed, 2))

    for (lib in c("rf", "nn")) {
      mname <- paste0("aipw_", lib)
      ans <- tryCatch({
        m0 <- library_predictions(outer_cube, eval$S, strata, lib, 0)[, 1L]
        m1 <- library_predictions(outer_cube, eval$S, strata, lib, 1)[, 1L]
        estimate_car(eval, m0, m1)
      }, error = function(e) e)
      if (!inherits(ans, "error")) {
        scores[[mname]][eval_idx] <- ans$scores
        fold_variances[[mname]][f] <- ans$variance
      }
    }

    calibration_libraries <- c("rf", "nn", "rfnn", "grouped_rf", "grouped_nn", "rflin")
    for (lib in calibration_libraries) {
      if (lib == "rflin" && inherits(outer_lm, "error")) next
      xi <- calibration_proxy(outer_cube, eval$S, strata, lib)
      for (discrepancy in c("quadratic", "quartic")) {
        method <- paste0("cal_", discrepancy, "__", lib)
        ans <- tryCatch(calibration_fold(eval, xi, discrepancy), error = function(e) e)
        if (!inherits(ans, "error")) {
          scores[[method]][eval_idx] <- ans$scores
          fold_variances[[method]][f] <- ans$variance$variance
        }
      }
    }

    for (lib in c("rfnn", "grouped_rf", "grouped_nn")) {
      for (meta_type in c("sl", "stack")) {
        mname <- paste0("aipw_", meta_type, "__", lib)
        ans <- tryCatch({
          mw <- learn_meta(train, oof_cube, strata, lib, meta_type)
          m0 <- apply_meta(outer_cube, eval$S, strata, lib, 0, mw)
          m1 <- apply_meta(outer_cube, eval$S, strata, lib, 1, mw)
          estimate_car(eval, m0, m1)
        }, error = function(e) e)
        if (!inherits(ans, "error")) {
          scores[[mname]][eval_idx] <- ans$scores
          fold_variances[[mname]][f] <- ans$variance
        }
      }
    }
  }

  fold_sizes <- tabulate(outer, nbins = outer_folds)
  rows <- lapply(seq_len(nrow(catalog)), function(j) {
    method <- catalog$method[j]
    sc <- scores[[method]]
    complete <- all(is.finite(sc))
    estimate <- if (complete) mean(sc) else NA_real_
    if (!complete) {
      se <- NA_real_
    } else if (catalog$cross_fitted[j]) {
      se <- sqrt(combine_fold_variances(fold_variances[[method]], fold_sizes))
    } else {
      se <- sqrt(full_variances[[method]])
    }
    if (complete && (!is.finite(se) || se < 0)) {
      estimate <- se <- NA_real_
    }
    err <- estimate - dat$target
    covered <- abs(err) <= stats::qnorm(0.975) * se
    data.frame(
      dgp = model, design = design, replication = rep_id, method = method,
      estimator_class = catalog$estimator_class[j], library = catalog$library[j],
      cross_fitted = catalog$cross_fitted[j],
      n = n, estimate = estimate, se = se,
      target = dat$target, target_type = dat$target_type,
      population_tau = dat$population_tau, sample_ate = dat$sample_ate,
      conditional_sample_ate = dat$conditional_sample_ate,
      error = err, abs_error = abs(err), squared_error = err^2, covered = covered,
      data_seed = data_seed, assignment_seed = assignment_seed,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

empty_replication_rows <- function(model, rep_id, design, n) {
  catalog <- method_catalog()
  seeds <- replication_seeds(model, n, rep_id, design)
  data.frame(
    dgp = model, design = design, replication = rep_id, method = catalog$method,
    estimator_class = catalog$estimator_class, library = catalog$library,
    cross_fitted = catalog$cross_fitted,
    n = n, estimate = NA_real_, se = NA_real_,
    target = NA_real_, target_type = "population_ate", population_tau = NA_real_,
    sample_ate = NA_real_, conditional_sample_ate = NA_real_, error = NA_real_,
    abs_error = NA_real_, squared_error = NA_real_, covered = NA,
    data_seed = seeds$data, assignment_seed = seeds$assignment,
    stringsAsFactors = FALSE
  )
}

aggregate_replications <- function(rep_dir, dgp, design, n, expected_reps) {
  paths <- file.path(rep_dir, sprintf("rep_%03d.csv", seq_len(expected_reps)))
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing replication results: ", paste(basename(missing), collapse = ","))
  required_columns <- names(empty_replication_rows(dgp, 1L, design, n))
  raw <- do.call(rbind, lapply(paths, function(path) {
    rows <- utils::read.csv(path, stringsAsFactors = FALSE)
    if (!all(required_columns %in% names(rows))) {
      stop("Replication result schema is incomplete for ", design, "/n", n, "/", dgp)
    }
    rows[, c(required_columns, intersect("run_id", names(rows))), drop = FALSE]
  }))
  catalog <- method_catalog()
  expected_rows <- expected_reps * nrow(catalog)
  if (nrow(raw) != expected_rows) stop("Unexpected raw row count for ", dgp)
  if (any(raw$dgp != dgp) || any(raw$design != design) || any(raw$n != n)) {
    stop("Replication result scope mismatch for ", design, "/n", n, "/", dgp)
  }
  if (anyNA(raw$method) || !all(raw$method %in% catalog$method)) {
    stop("Replication result contains an unknown method")
  }
  keys <- paste(raw$design, raw$n, raw$dgp, raw$replication, raw$method)
  if (anyDuplicated(keys)) stop("Duplicate replication-method rows for ", dgp)
  summaries <- lapply(split(raw, raw$method), function(x) {
    ok <- is.finite(x$estimate) & is.finite(x$error) & is.finite(x$se) & x$se >= 0
    valid <- sum(ok)
    cp <- if (valid) mean(as.logical(x$covered[ok])) else NA_real_
    data.frame(
      dgp = dgp, design = design, n = n, method = x$method[1],
      estimator_class = x$estimator_class[1], library = x$library[1],
      cross_fitted = x$cross_fitted[1],
      reps = nrow(x), valid = valid,
      signed_bias = if (valid) mean(x$error[ok]) else NA_real_,
      absolute_bias = if (valid) abs(mean(x$error[ok])) else NA_real_,
      sd = if (valid > 1L) stats::sd(x$error[ok]) else NA_real_,
      mean_se = if (valid) mean(x$se[ok]) else NA_real_,
      rmse = if (valid) sqrt(mean(x$squared_error[ok])) else NA_real_,
      coverage = cp, coverage_mcse = if (valid) sqrt(cp * (1 - cp) / valid) else NA_real_,
      median_absolute_error = if (valid) stats::median(x$abs_error[ok]) else NA_real_,
      target_type = paste(unique(x$target_type), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  list(raw = raw, summary = summary)
}
