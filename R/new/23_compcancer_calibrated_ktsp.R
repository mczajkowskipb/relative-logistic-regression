get_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  getwd()
}
source(file.path(get_script_dir(), "07_world_protocol_utils.R"))

suppressPackageStartupMessages({
  library(jsonlite)
  library(pROC)
  library(foreign)
})

# ============================================================
# ARG PARSER
# ============================================================
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[1], fixed = TRUE)
}

dataset_id_arg   <- get_arg("dataset_id", NULL)
outer_folds      <- as.integer(get_arg("outer_folds", "5"))
outer_repeats    <- as.integer(get_arg("outer_repeats", "10"))
inner_folds      <- as.integer(get_arg("inner_folds", "5"))
seed             <- as.integer(get_arg("seed", "42"))

# comma-separated k list: --k_values=10,20,50
k_values_arg <- get_arg("k_values", "10,20,50")
k_values <- as.integer(strsplit(k_values_arg, ",", fixed = TRUE)[[1]])

# ============================================================
# CONFIG
# ============================================================
project_root <- Sys.getenv("RLR_PROJECT_ROOT", unset = normalizePath(file.path(get_script_dir(), "../..")))
comp_dir <- get_arg("comp_dir", file.path(project_root, "data", "CompCancer"))
out_root_base <- get_arg("out_root_base", file.path(project_root, "results", "supplementary", "compcancer_calibrated_ktsp"))

if (is.null(dataset_id_arg)) {
  out_root <- file.path(out_root_base, "_all")
} else {
  out_root <- file.path(out_root_base, dataset_id_arg)
}

ensure_dir(out_root)
ensure_dir(file.path(out_root, "config"))
ensure_dir(file.path(out_root, "summaries"))
ensure_dir(file.path(out_root, "worlds"))

write_json(
  list(
    comp_dir = comp_dir,
    out_root = out_root,
    dataset_id = dataset_id_arg,
    outer_folds = outer_folds,
    outer_repeats = outer_repeats,
    inner_folds = inner_folds,
    k_values = k_values,
    seed = seed,
    note = "Calibrated kTSP with proper inner-OOF calibration on outer training only"
  ),
  file.path(out_root, "config", "used_config.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

# ============================================================
# HELPERS
# ============================================================
clip_prob <- function(p) {
  pmin(pmax(p, 1e-15), 1 - 1e-15)
}

mean_or_na <- function(x) {
  if (length(x) == 0) return(NA_real_)
  mean(x, na.rm = TRUE)
}

calc_metrics_from_prob <- function(prob, y_true) {
  prob <- as.numeric(prob)
  prob <- clip_prob(prob)

  y_true <- factor(y_true)
  y01 <- as.numeric(y_true == levels(y_true)[2])

  pred <- ifelse(prob >= 0.5, levels(y_true)[2], levels(y_true)[1])
  pred <- factor(pred, levels = levels(y_true))

  auc <- as.numeric(pROC::auc(y01, prob))
  logloss <- -mean(y01 * log(prob) + (1 - y01) * log(1 - prob))
  brier <- calc_brier(prob, y01)
  ece <- calc_ece(prob, y01)

  tp <- sum(pred == levels(y_true)[2] & y_true == levels(y_true)[2])
  fp <- sum(pred == levels(y_true)[2] & y_true == levels(y_true)[1])
  fn <- sum(pred == levels(y_true)[1] & y_true == levels(y_true)[2])

  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))

  list(
    auc = auc,
    f1 = f1,
    logloss = logloss,
    brier = brier,
    ece = ece
  )
}

fit_ktsp_and_get_raw <- function(X_train, y_train, X_eval, k) {
  trainMat <- xy_to_switchbox_mat(X_train)
  evalMat  <- xy_to_switchbox_mat(X_eval)

  clf <- SWAP.KTSP.Train(
    inputMat = trainMat,
    phenoGroup = y_train,
    FilterFunc = NULL,
    krange = k
  )

  if (!is.matrix(clf$TSPs)) {
    clf$TSPs <- matrix(clf$TSPs, ncol = 2, byrow = TRUE)
  }

  realized_k <- nrow(clf$TSPs)
  if (realized_k <= 0) stop("realized_k <= 0")

  stat_train <- SWAP.KTSP.Statistics(trainMat, clf)
  stat_eval  <- SWAP.KTSP.Statistics(evalMat, clf)

  raw_train <- as.numeric(stat_train$statistics) / realized_k
  raw_eval  <- as.numeric(stat_eval$statistics) / realized_k

  # orient using TRAIN only
  y_train01 <- as.numeric(y_train == levels(y_train)[2])
  mean_pos <- mean(raw_train[y_train01 == 1])
  mean_neg <- mean(raw_train[y_train01 == 0])

  flip <- is.finite(mean_pos) && is.finite(mean_neg) && (mean_pos < mean_neg)
  if (flip) {
    raw_train <- -raw_train
    raw_eval <- -raw_eval
  }

  selected_pairs <- as.data.frame(clf$TSPs, stringsAsFactors = FALSE)
  colnames(selected_pairs) <- c("gene_left", "gene_right")
  selected_pairs$ordered_pair <- paste0(selected_pairs$gene_left, " > ", selected_pairs$gene_right)

  list(
    clf = clf,
    realized_k = realized_k,
    raw_train = raw_train,
    raw_eval = raw_eval,
    selected_pairs = selected_pairs,
    flip = flip
  )
}

fit_platt_calibrator <- function(raw, y_true) {
  raw <- as.numeric(raw)
  y_true <- factor(y_true)
  y01 <- as.numeric(y_true == levels(y_true)[2])

  fallback_predict <- function(z) clip_prob((as.numeric(z) + 1) / 2)

  if (length(unique(y01)) < 2 || length(unique(raw)) < 2) {
    return(list(
      kind = "platt",
      status = "fallback_uncalibrated",
      predict = fallback_predict,
      intercept = NA_real_,
      slope = NA_real_
    ))
  }

  fit <- tryCatch(
    suppressWarnings(glm(y01 ~ raw, family = binomial())),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(list(
      kind = "platt",
      status = "fallback_uncalibrated",
      predict = fallback_predict,
      intercept = NA_real_,
      slope = NA_real_
    ))
  }

  cc <- suppressWarnings(coef(fit))
  if (length(cc) < 2 || any(!is.finite(cc))) {
    return(list(
      kind = "platt",
      status = "fallback_uncalibrated",
      predict = fallback_predict,
      intercept = NA_real_,
      slope = NA_real_
    ))
  }

  pred_fun <- function(z) {
    z <- as.numeric(z)
    p <- suppressWarnings(predict(fit, newdata = data.frame(raw = z), type = "response"))
    clip_prob(as.numeric(p))
  }

  list(
    kind = "platt",
    status = "ok",
    predict = pred_fun,
    intercept = as.numeric(cc[1]),
    slope = as.numeric(cc[2])
  )
}

fit_isotonic_calibrator <- function(raw, y_true) {
  raw <- as.numeric(raw)
  y_true <- factor(y_true)
  y01 <- as.numeric(y_true == levels(y_true)[2])

  fallback_predict <- function(z) clip_prob((as.numeric(z) + 1) / 2)

  if (length(unique(y01)) < 2 || length(unique(raw)) < 2) {
    return(list(
      kind = "isotonic",
      status = "fallback_uncalibrated",
      predict = fallback_predict
    ))
  }

  ord <- order(raw)
  x <- raw[ord]
  y <- y01[ord]

  iso <- tryCatch(
    isoreg(x, y),
    error = function(e) NULL
  )

  if (is.null(iso)) {
    return(list(
      kind = "isotonic",
      status = "fallback_uncalibrated",
      predict = fallback_predict
    ))
  }

  # collapse duplicate x values
  yf <- iso$yf
  agg <- tapply(yf, x, mean)
  xk <- as.numeric(names(agg))
  yk <- as.numeric(agg)

  if (length(xk) < 2 || any(!is.finite(xk)) || any(!is.finite(yk))) {
    return(list(
      kind = "isotonic",
      status = "fallback_uncalibrated",
      predict = fallback_predict
    ))
  }

  pred_fun <- function(z) {
    z <- as.numeric(z)
    p <- approx(x = xk, y = yk, xout = z, method = "linear", rule = 2, ties = mean)$y
    clip_prob(as.numeric(p))
  }

  list(
    kind = "isotonic",
    status = "ok",
    predict = pred_fun
  )
}

make_metric_row <- function(world_name, dataset_id, split_id, model, calibrator,
                            requested_k, realized_k, metrics,
                            selected_pairs, calibration_train_n,
                            calibrator_status) {
  unique_genes <- length(unique(c(selected_pairs$gene_left, selected_pairs$gene_right)))
  hubs <- calc_hub_summary(selected_pairs)

  data.frame(
    world_name = world_name,
    dataset_id = dataset_id,
    split_id = split_id,
    model = model,
    calibrator = calibrator,
    requested_k = requested_k,
    realized_k = realized_k,
    calibration_train_n = calibration_train_n,
    calibrator_status = calibrator_status,
    auc = metrics$auc,
    f1 = metrics$f1,
    logloss = metrics$logloss,
    brier = metrics$brier,
    ece = metrics$ece,
    global_active_pairs = realized_k,
    unique_genes = unique_genes,
    local_active_mean = realized_k,
    local_active_median = realized_k,
    local_active_iqr = 0,
    hub_top_share = hubs$top_share,
    stringsAsFactors = FALSE
  )
}

aggregate_dataset_level <- function(df) {
  keys <- unique(df[, c("world_name", "dataset_id", "model", "calibrator"), drop = FALSE])

  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    idx <- df$world_name == key$world_name &
      df$dataset_id == key$dataset_id &
      df$model == key$model &
      df$calibrator == key$calibrator

    data.frame(
      world_name = key$world_name,
      dataset_id = key$dataset_id,
      model = key$model,
      calibrator = key$calibrator,
      requested_k = mean_or_na(df$requested_k[idx]),
      realized_k_mean = mean_or_na(df$realized_k[idx]),
      calibration_train_n_mean = mean_or_na(df$calibration_train_n[idx]),
      auc = mean_or_na(df$auc[idx]),
      f1 = mean_or_na(df$f1[idx]),
      logloss = mean_or_na(df$logloss[idx]),
      brier = mean_or_na(df$brier[idx]),
      ece = mean_or_na(df$ece[idx]),
      global_active_pairs = mean_or_na(df$global_active_pairs[idx]),
      unique_genes = mean_or_na(df$unique_genes[idx]),
      local_active_mean = mean_or_na(df$local_active_mean[idx]),
      local_active_median = mean_or_na(df$local_active_median[idx]),
      local_active_iqr = mean_or_na(df$local_active_iqr[idx]),
      hub_top_share = mean_or_na(df$hub_top_share[idx]),
      n_splits = sum(idx),
      stringsAsFactors = FALSE
    )
  }))
}

aggregate_world_level <- function(dataset_df) {
  keys <- unique(dataset_df[, c("world_name", "model", "calibrator"), drop = FALSE])

  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    idx <- dataset_df$world_name == key$world_name &
      dataset_df$model == key$model &
      dataset_df$calibrator == key$calibrator

    data.frame(
      world_name = key$world_name,
      model = key$model,
      calibrator = key$calibrator,
      requested_k = mean_or_na(dataset_df$requested_k[idx]),
      realized_k_mean = mean_or_na(dataset_df$realized_k_mean[idx]),
      calibration_train_n_mean = mean_or_na(dataset_df$calibration_train_n_mean[idx]),
      auc = mean_or_na(dataset_df$auc[idx]),
      f1 = mean_or_na(dataset_df$f1[idx]),
      logloss = mean_or_na(dataset_df$logloss[idx]),
      brier = mean_or_na(dataset_df$brier[idx]),
      ece = mean_or_na(dataset_df$ece[idx]),
      global_active_pairs = mean_or_na(dataset_df$global_active_pairs[idx]),
      unique_genes = mean_or_na(dataset_df$unique_genes[idx]),
      local_active_mean = mean_or_na(dataset_df$local_active_mean[idx]),
      local_active_median = mean_or_na(dataset_df$local_active_median[idx]),
      local_active_iqr = mean_or_na(dataset_df$local_active_iqr[idx]),
      hub_top_share = mean_or_na(dataset_df$hub_top_share[idx]),
      n_datasets = sum(idx),
      stringsAsFactors = FALSE
    )
  }))
}

# ============================================================
# MAIN
# ============================================================
files <- list.files(comp_dir, pattern = "\\.arff$", full.names = TRUE)
files <- sort(files)

if (!is.null(dataset_id_arg)) {
  files <- files[basename(files) == paste0(dataset_id_arg, ".arff")]
  if (length(files) == 0) stop("Nie znaleziono datasetu: ", dataset_id_arg)
}

all_rows <- list()
ri <- 1L
manifest_rows <- list()
mi <- 1L

for (f in files) {
  dataset_id <- sub("\\.arff$", "", basename(f))
  cat("\n==============================\n")
  cat("DATASET:", dataset_id, "\n")

  dataset_obj <- read_xy_file(f)

  outer_folds_obj <- make_repeated_stratified_folds(
    y = dataset_obj$y,
    n_folds = outer_folds,
    n_repeats = outer_repeats,
    seed = seed
  )

  ds_dir <- file.path(out_root, "worlds", "compcancer_cv", "datasets", dataset_id)
  ensure_dir(ds_dir)

  for (sp in outer_folds_obj) {
    split_dir <- file.path(ds_dir, "splits", sp$split_id)
    ensure_dir(split_dir)

    X_outer_train <- dataset_obj$X[sp$train_idx, , drop = FALSE]
    y_outer_train <- dataset_obj$y[sp$train_idx]
    X_outer_test  <- dataset_obj$X[sp$test_idx, , drop = FALSE]
    y_outer_test  <- dataset_obj$y[sp$test_idx]

    write.csv(
      data.frame(sample_id = rownames(X_outer_train)),
      file.path(split_dir, "samples_train.csv"),
      row.names = FALSE
    )
    write.csv(
      data.frame(sample_id = rownames(X_outer_test)),
      file.path(split_dir, "samples_test.csv"),
      row.names = FALSE
    )

    split_meta <- list(
      world_name = "compcancer_cv",
      dataset_id = dataset_id,
      split_id = sp$split_id,
      train_origin = "cv_train",
      test_origin = "cv_test",
      note = "Calibration fitted on inner OOF raw scores from outer training only"
    )
    write_json(split_meta, file.path(split_dir, "split_meta.json"), auto_unbox = TRUE, pretty = TRUE)

    inner_folds_obj <- make_repeated_stratified_folds(
      y = y_outer_train,
      n_folds = inner_folds,
      n_repeats = 1,
      seed = seed
    )

    for (k in k_values) {
      base_model_name <- paste0("kTSP_", k)
      model_dir <- file.path(split_dir, "baselines_calibrated", base_model_name)
      ensure_dir(model_dir)

      cat("  ->", base_model_name, "\n")

      # ------------------------------------------------------
      # 1) Inner OOF raw scores on OUTER TRAIN only
      # ------------------------------------------------------
      oof_rows <- list()
      oi <- 1L

      for (isp in inner_folds_obj) {
        Xi_train <- X_outer_train[isp$train_idx, , drop = FALSE]
        yi_train <- y_outer_train[isp$train_idx]

        Xi_valid <- X_outer_train[isp$test_idx, , drop = FALSE]
        yi_valid <- y_outer_train[isp$test_idx]

        inner_fit <- tryCatch(
          fit_ktsp_and_get_raw(Xi_train, yi_train, Xi_valid, k = k),
          error = function(e) NULL
        )

        if (!is.null(inner_fit)) {
          oof_rows[[oi]] <- data.frame(
            sample_id = rownames(Xi_valid),
            raw_score = inner_fit$raw_eval,
            y_true = as.character(yi_valid),
            stringsAsFactors = FALSE
          )
          oi <- oi + 1L
        }
      }

      if (length(oof_rows) == 0) {
        write_error_json(file.path(model_dir, "error.json"),
                         simpleError(paste("Failed to produce inner OOF scores for", base_model_name)))
        next
      }

      calibration_oof_df <- do.call(rbind, oof_rows)
      calibration_oof_df <- calibration_oof_df[!duplicated(calibration_oof_df$sample_id), , drop = FALSE]
      write.csv(calibration_oof_df, file.path(model_dir, "calibration_oof_train.csv"), row.names = FALSE)

      # ------------------------------------------------------
      # 2) Fit calibrators on INNER OOF raw scores
      # ------------------------------------------------------
      y_oof <- factor(calibration_oof_df$y_true, levels = levels(y_outer_train))
      platt_cal <- fit_platt_calibrator(calibration_oof_df$raw_score, y_oof)
      iso_cal   <- fit_isotonic_calibrator(calibration_oof_df$raw_score, y_oof)

      write_json(
        list(
          platt_status = platt_cal$status,
          platt_intercept = platt_cal$intercept,
          platt_slope = platt_cal$slope,
          isotonic_status = iso_cal$status,
          calibration_train_n = nrow(calibration_oof_df)
        ),
        file.path(model_dir, "calibrator_info.json"),
        auto_unbox = TRUE,
        pretty = TRUE
      )

      # ------------------------------------------------------
      # 3) Refit kTSP on full OUTER TRAIN
      # ------------------------------------------------------
      final_fit <- tryCatch(
        fit_ktsp_and_get_raw(X_outer_train, y_outer_train, X_outer_test, k = k),
        error = function(e) NULL
      )

      if (is.null(final_fit)) {
        write_error_json(file.path(model_dir, "error.json"),
                         simpleError(paste("Failed final outer fit for", base_model_name)))
        next
      }

      realized_k <- final_fit$realized_k
      raw_test <- final_fit$raw_eval
      selected_pairs <- final_fit$selected_pairs

      prob_uncal <- clip_prob((raw_test + 1) / 2)
      prob_platt <- platt_cal$predict(raw_test)
      prob_iso   <- iso_cal$predict(raw_test)

      pred_df <- data.frame(
        sample_id = rownames(X_outer_test),
        y_true = as.character(y_outer_test),
        raw_score = raw_test,
        prob_uncal = prob_uncal,
        prob_platt = prob_platt,
        prob_isotonic = prob_iso,
        stringsAsFactors = FALSE
      )
      write.csv(pred_df, file.path(model_dir, "predictions_test.csv"), row.names = FALSE)
      write.csv(selected_pairs, file.path(model_dir, "selected_pairs.csv"), row.names = FALSE)

      # ------------------------------------------------------
      # 4) Metrics
      # ------------------------------------------------------
      m_uncal <- calc_metrics_from_prob(prob_uncal, y_outer_test)
      m_platt <- calc_metrics_from_prob(prob_platt, y_outer_test)
      m_iso   <- calc_metrics_from_prob(prob_iso, y_outer_test)

      row_uncal <- make_metric_row(
        world_name = "compcancer_cv",
        dataset_id = dataset_id,
        split_id = sp$split_id,
        model = paste0(base_model_name, "_uncal"),
        calibrator = "uncal",
        requested_k = k,
        realized_k = realized_k,
        metrics = m_uncal,
        selected_pairs = selected_pairs,
        calibration_train_n = nrow(calibration_oof_df),
        calibrator_status = "identity"
      )

      row_platt <- make_metric_row(
        world_name = "compcancer_cv",
        dataset_id = dataset_id,
        split_id = sp$split_id,
        model = paste0(base_model_name, "_platt"),
        calibrator = "platt",
        requested_k = k,
        realized_k = realized_k,
        metrics = m_platt,
        selected_pairs = selected_pairs,
        calibration_train_n = nrow(calibration_oof_df),
        calibrator_status = platt_cal$status
      )

      row_iso <- make_metric_row(
        world_name = "compcancer_cv",
        dataset_id = dataset_id,
        split_id = sp$split_id,
        model = paste0(base_model_name, "_isotonic"),
        calibrator = "isotonic",
        requested_k = k,
        realized_k = realized_k,
        metrics = m_iso,
        selected_pairs = selected_pairs,
        calibration_train_n = nrow(calibration_oof_df),
        calibrator_status = iso_cal$status
      )

      all_rows[[ri]] <- row_uncal; ri <- ri + 1L
      all_rows[[ri]] <- row_platt; ri <- ri + 1L
      all_rows[[ri]] <- row_iso;   ri <- ri + 1L

      manifest_rows[[mi]] <- data.frame(
        world_name = "compcancer_cv",
        dataset_id = dataset_id,
        split_id = sp$split_id,
        family = "baseline_calibrated",
        model = base_model_name,
        path = model_dir,
        stringsAsFactors = FALSE
      )
      mi <- mi + 1L
    }
  }
}

manifest_df <- do.call(rbind, manifest_rows)
write.csv(manifest_df, file.path(out_root, "config", "run_manifest.csv"), row.names = FALSE)

if (length(all_rows) == 0) {
  stop("Brak poprawnych wyników calibrated kTSP.")
}

split_df <- do.call(rbind, all_rows)
write.csv(split_df, file.path(out_root, "summaries", "all_split_metrics.csv"), row.names = FALSE)

dataset_df <- aggregate_dataset_level(split_df)
write.csv(dataset_df, file.path(out_root, "summaries", "dataset_level_metrics.csv"), row.names = FALSE)

world_df <- aggregate_world_level(dataset_df)
write.csv(world_df, file.path(out_root, "summaries", "world_level_metrics.csv"), row.names = FALSE)

cat("\nDONE\n")
cat("Results saved in:\n")
cat(out_root, "\n")