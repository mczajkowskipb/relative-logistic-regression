get_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  getwd()
}
source(file.path(get_script_dir(), "07_world_protocol_utils.R"))

suppressPackageStartupMessages({
  library(jsonlite)
  library(pROC)
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

dataset_id_arg <- get_arg("dataset_id", NULL)
search_mode <- get_arg("mode", "fast")        # full | odd | fast
selection_metric <- get_arg("selection_metric", "logloss")
inner_n_folds <- as.integer(get_arg("inner_folds", "3"))
outer_n_folds <- as.integer(get_arg("outer_folds", "5"))
outer_n_repeats <- as.integer(get_arg("outer_repeats", "10"))
seed <- as.integer(get_arg("seed", "42"))
pre_top_features_arg <- get_arg("pre_top_features", NULL)

pre_top_features <- if (is.null(pre_top_features_arg) || pre_top_features_arg == "" || pre_top_features_arg == "NULL") {
  NULL
} else {
  as.integer(pre_top_features_arg)
}

# ============================================================
# CONFIG
# ============================================================
project_root <- Sys.getenv("RLR_PROJECT_ROOT", unset = normalizePath(file.path(get_script_dir(), "../..")))
comp_dir <- get_arg("comp_dir", file.path(project_root, "data", "CompCancer"))
out_root_base <- get_arg("out_root_base", file.path(project_root, "results", "sensitivity", "compcancer_adaptive_ktsp_fast"))

kmax_values <- c(10, 20, 50)

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
    kmax_values = kmax_values,
    outer_n_folds = outer_n_folds,
    outer_n_repeats = outer_n_repeats,
    inner_n_folds = inner_n_folds,
    selection_metric = selection_metric,
    search_mode = search_mode,
    pre_top_features = pre_top_features,
    seed = seed
  ),
  file.path(out_root, "config", "used_config.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

# ============================================================
# HELPERS
# ============================================================
mean_or_na <- function(x) {
  if (length(x) == 0) return(NA_real_)
  mean(x, na.rm = TRUE)
}

odd_seq <- function(kmax) {
  seq(1, kmax, by = 2)
}

coarse_odd_seq <- function(kmax) {
  if (kmax <= 10) return(odd_seq(kmax))
  if (kmax <= 20) return(odd_seq(kmax))
  unique(c(seq(1, 9, by = 2), seq(13, kmax, by = 4)))
}

order_by_metric <- function(df, metric) {
  if (metric %in% c("logloss", "brier", "ece")) {
    order(df[[metric]], df$k)
  } else if (metric %in% c("auc", "f1")) {
    order(-df[[metric]], df$k)
  } else {
    stop("Unsupported selection_metric: ", metric)
  }
}

eval_ktsp_once <- function(X_train, y_train, X_test, y_test, k) {
  trainMat <- xy_to_switchbox_mat(X_train)
  testMat  <- xy_to_switchbox_mat(X_test)

  clf <- SWAP.KTSP.Train(
    inputMat = trainMat,
    phenoGroup = y_train,
    FilterFunc = NULL,
    krange = k
  )

  if (!is.matrix(clf$TSPs)) {
    clf$TSPs <- matrix(clf$TSPs, ncol = 2, byrow = TRUE)
  }

  pred <- SWAP.KTSP.Classify(testMat, clf)
  stat_train <- SWAP.KTSP.Statistics(trainMat, clf)
  stat_test  <- SWAP.KTSP.Statistics(testMat, clf)

  realized_k <- nrow(clf$TSPs)

  prob_train <- (as.numeric(stat_train$statistics) + realized_k) / (2 * realized_k)
  prob_test  <- (as.numeric(stat_test$statistics) + realized_k) / (2 * realized_k)

  y_train01 <- as.numeric(y_train == levels(y_train)[2])
  y_test01  <- as.numeric(y_test == levels(y_test)[2])

  flip <- mean(prob_train[y_train01 == 1]) < mean(prob_train[y_train01 == 0])
  if (flip) prob_test <- 1 - prob_test

  pred <- factor(pred, levels = levels(y_test))
  eps <- 1e-15
  p <- pmin(pmax(prob_test, eps), 1 - eps)

  auc <- as.numeric(pROC::auc(y_test01, prob_test))
  logloss <- -mean(y_test01 * log(p) + (1 - y_test01) * log(1 - p))
  brier <- calc_brier(prob_test, y_test01)
  ece <- calc_ece(prob_test, y_test01)

  tp <- sum(pred == levels(y_test)[2] & y_test == levels(y_test)[2])
  fp <- sum(pred == levels(y_test)[2] & y_test == levels(y_test)[1])
  fn <- sum(pred == levels(y_test)[1] & y_test == levels(y_test)[2])

  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))

  list(
    clf = clf,
    realized_k = realized_k,
    pred = pred,
    prob = prob_test,
    auc = auc,
    f1 = f1,
    logloss = logloss,
    brier = brier,
    ece = ece
  )
}

evaluate_candidate_ks <- function(X_train, y_train, candidate_ks, inner_n_folds, seed) {
  inner_folds <- make_repeated_stratified_folds(
    y = y_train,
    n_folds = inner_n_folds,
    n_repeats = 1,
    seed = seed
  )

  rows <- list()
  ri <- 1L

  for (sp in inner_folds) {
    Xi_train <- X_train[sp$train_idx, , drop = FALSE]
    yi_train <- y_train[sp$train_idx]

    Xi_valid <- X_train[sp$test_idx, , drop = FALSE]
    yi_valid <- y_train[sp$test_idx]

    for (k in candidate_ks) {
      one <- tryCatch({
        eval_ktsp_once(Xi_train, yi_train, Xi_valid, yi_valid, k = k)
      }, error = function(e) {
        NULL
      })

      if (!is.null(one)) {
        rows[[ri]] <- data.frame(
          k = k,
          realized_k = one$realized_k,
          auc = one$auc,
          f1 = one$f1,
          logloss = one$logloss,
          brier = one$brier,
          ece = one$ece,
          stringsAsFactors = FALSE
        )
        ri <- ri + 1L
      }
    }
  }

  if (length(rows) == 0) {
    stop("Inner CV failed for all candidate k.")
  }

  inner_df <- do.call(rbind, rows)

  agg <- do.call(rbind, lapply(split(inner_df, inner_df$k), function(dd) {
    data.frame(
      k = dd$k[1],
      realized_k = mean_or_na(dd$realized_k),
      auc = mean_or_na(dd$auc),
      f1 = mean_or_na(dd$f1),
      logloss = mean_or_na(dd$logloss),
      brier = mean_or_na(dd$brier),
      ece = mean_or_na(dd$ece),
      stringsAsFactors = FALSE
    )
  }))

  agg <- agg[order_by_metric(agg, selection_metric), , drop = FALSE]

  list(
    inner_split_metrics = inner_df,
    inner_k_summary = agg
  )
}

choose_best_k_inner_cv <- function(X_train, y_train, kmax,
                                   inner_n_folds = 3,
                                   selection_metric = "logloss",
                                   seed = 42,
                                   search_mode = "fast") {
  if (search_mode == "full") {
    candidate_ks <- seq_len(kmax)

    out <- evaluate_candidate_ks(
      X_train = X_train,
      y_train = y_train,
      candidate_ks = candidate_ks,
      inner_n_folds = inner_n_folds,
      seed = seed
    )

    best_k <- out$inner_k_summary$k[1]

    return(list(
      best_k = best_k,
      candidate_ks = candidate_ks,
      inner_split_metrics = out$inner_split_metrics,
      inner_k_summary = out$inner_k_summary
    ))
  }

  if (search_mode == "odd") {
    candidate_ks <- odd_seq(kmax)

    out <- evaluate_candidate_ks(
      X_train = X_train,
      y_train = y_train,
      candidate_ks = candidate_ks,
      inner_n_folds = inner_n_folds,
      seed = seed
    )

    best_k <- out$inner_k_summary$k[1]

    return(list(
      best_k = best_k,
      candidate_ks = candidate_ks,
      inner_split_metrics = out$inner_split_metrics,
      inner_k_summary = out$inner_k_summary
    ))
  }

  if (search_mode == "fast") {
    coarse_ks <- coarse_odd_seq(kmax)

    coarse_out <- evaluate_candidate_ks(
      X_train = X_train,
      y_train = y_train,
      candidate_ks = coarse_ks,
      inner_n_folds = inner_n_folds,
      seed = seed
    )

    best_coarse <- coarse_out$inner_k_summary$k[1]
    refine_ks <- seq(max(1, best_coarse - 4), min(kmax, best_coarse + 4), by = 1)
    refine_ks <- refine_ks[refine_ks %% 2 == 1]

    final_ks <- sort(unique(c(coarse_ks, refine_ks)))

    final_out <- evaluate_candidate_ks(
      X_train = X_train,
      y_train = y_train,
      candidate_ks = final_ks,
      inner_n_folds = inner_n_folds,
      seed = seed
    )

    best_k <- final_out$inner_k_summary$k[1]

    return(list(
      best_k = best_k,
      candidate_ks = final_ks,
      inner_split_metrics = final_out$inner_split_metrics,
      inner_k_summary = final_out$inner_k_summary
    ))
  }

  stop("Unsupported search_mode: ", search_mode)
}

fit_ktsp_adaptive_train_test <- function(X_train, y_train, X_test, y_test,
                                         kmax,
                                         out_dir,
                                         inner_n_folds = 3,
                                         selection_metric = "logloss",
                                         seed = 42,
                                         search_mode = "fast",
                                         pre_top_features = NULL,
                                         extra_metrics = NULL) {
  ensure_dir(out_dir)
  set.seed(seed)

  # optional speed-up: train-only MAD prefilter
  if (!is.null(pre_top_features)) {
    pf <- prefilter_by_mad(
      X_train = X_train,
      X_test = X_test,
      top_n = pre_top_features
    )
    X_train <- pf$X_train
    X_test  <- pf$X_test
  }

  chosen <- choose_best_k_inner_cv(
    X_train = X_train,
    y_train = y_train,
    kmax = kmax,
    inner_n_folds = inner_n_folds,
    selection_metric = selection_metric,
    seed = seed,
    search_mode = search_mode
  )

  best_k <- chosen$best_k

  final_fit <- eval_ktsp_once(
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    y_test = y_test,
    k = best_k
  )

  clf <- final_fit$clf
  realized_k <- final_fit$realized_k
  pred <- final_fit$pred
  prob_test <- final_fit$prob

  selected_pairs <- as.data.frame(clf$TSPs, stringsAsFactors = FALSE)
  colnames(selected_pairs) <- c("gene_left", "gene_right")
  selected_pairs$ordered_pair <- paste0(selected_pairs$gene_left, " > ", selected_pairs$gene_right)

  unique_genes <- length(unique(c(selected_pairs$gene_left, selected_pairs$gene_right)))
  hubs <- calc_hub_summary(selected_pairs)

  local_n <- rep(realized_k, nrow(X_test))
  local_sum <- list(
    local_n = local_n,
    mean = mean(local_n),
    median = median(local_n),
    iqr = IQR(local_n),
    min = min(local_n),
    max = max(local_n),
    pct_le_1 = mean(local_n <= 1),
    pct_le_3 = mean(local_n <= 3),
    pct_le_5 = mean(local_n <= 5),
    pct_le_10 = mean(local_n <= 10)
  )

  metrics <- list(
    method_family = "baseline_adaptive",
    method = paste0("kTSP_adaptive_kmax", kmax),
    requested_kmax = kmax,
    selected_k = best_k,
    realized_k = realized_k,
    selection_metric = selection_metric,
    search_mode = search_mode,
    pre_top_features = ifelse(is.null(pre_top_features), NA, pre_top_features),
    auc = final_fit$auc,
    f1 = final_fit$f1,
    logloss = final_fit$logloss,
    brier = final_fit$brier,
    ece = final_fit$ece,
    global_active_pairs = realized_k,
    unique_genes = unique_genes,
    local_active_mean = local_sum$mean,
    local_active_median = local_sum$median,
    local_active_iqr = local_sum$iqr,
    hub_top_gene = hubs$top_gene,
    hub_top_share = hubs$top_share
  )

  if (!is.null(extra_metrics)) {
    metrics <- c(metrics, extra_metrics)
  }

  write.csv(
    data.frame(
      sample_id = rownames(X_test),
      y_true = as.character(y_test),
      y_pred = as.character(pred),
      score = prob_test,
      prob = prob_test
    ),
    file.path(out_dir, "predictions_test.csv"),
    row.names = FALSE
  )

  write.csv(selected_pairs, file.path(out_dir, "selected_pairs.csv"), row.names = FALSE)

  write.csv(
    data.frame(
      sample_id = rownames(X_test),
      local_active_n = local_n
    ),
    file.path(out_dir, "local_activity_test.csv"),
    row.names = FALSE
  )

  write.csv(chosen$inner_split_metrics, file.path(out_dir, "inner_cv_split_metrics.csv"), row.names = FALSE)
  write.csv(chosen$inner_k_summary, file.path(out_dir, "inner_cv_k_summary.csv"), row.names = FALSE)

  write_json(local_sum, file.path(out_dir, "local_activity_summary.json"), auto_unbox = TRUE, pretty = TRUE)
  write_json(metrics, file.path(out_dir, "metrics_test.json"), auto_unbox = TRUE, pretty = TRUE)

  if (nrow(hubs$counts) > 0) {
    write.csv(hubs$counts, file.path(out_dir, "hub_counts.csv"), row.names = FALSE)
  }

  saveRDS(clf, file.path(out_dir, "classifier.rds"))
  invisible(metrics)
}

aggregate_dataset_level <- function(df) {
  keys <- unique(df[, c("world_name", "dataset_id", "model"), drop = FALSE])

  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    idx <- df$world_name == key$world_name &
      df$dataset_id == key$dataset_id &
      df$model == key$model

    data.frame(
      world_name = key$world_name,
      dataset_id = key$dataset_id,
      model = key$model,
      selected_k_mean = mean_or_na(df$selected_k[idx]),
      realized_k_mean = mean_or_na(df$realized_k[idx]),
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
  keys <- unique(dataset_df[, c("world_name", "model"), drop = FALSE])

  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    idx <- dataset_df$world_name == key$world_name &
      dataset_df$model == key$model

    data.frame(
      world_name = key$world_name,
      model = key$model,
      selected_k_mean = mean_or_na(dataset_df$selected_k_mean[idx]),
      realized_k_mean = mean_or_na(dataset_df$realized_k_mean[idx]),
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

all_split_rows <- list()
ri <- 1L
manifest_rows <- list()
mi <- 1L

for (f in files) {
  dataset_id <- sub("\\.arff$", "", basename(f))
  cat("\n==============================\n")
  cat("DATASET:", dataset_id, "\n")

  dataset_obj <- read_xy_file(f)

  folds <- make_repeated_stratified_folds(
    y = dataset_obj$y,
    n_folds = outer_n_folds,
    n_repeats = outer_n_repeats,
    seed = seed
  )

  ds_dir <- file.path(out_root, "worlds", "compcancer_cv", "datasets", dataset_id)
  ensure_dir(ds_dir)

  for (sp in folds) {
    split_dir <- file.path(ds_dir, "splits", sp$split_id)
    ensure_dir(split_dir)

    X_train <- dataset_obj$X[sp$train_idx, , drop = FALSE]
    y_train <- dataset_obj$y[sp$train_idx]
    X_test  <- dataset_obj$X[sp$test_idx, , drop = FALSE]
    y_test  <- dataset_obj$y[sp$test_idx]

    write.csv(
      data.frame(sample_id = rownames(X_train)),
      file.path(split_dir, "samples_train.csv"),
      row.names = FALSE
    )
    write.csv(
      data.frame(sample_id = rownames(X_test)),
      file.path(split_dir, "samples_test.csv"),
      row.names = FALSE
    )

    split_meta <- list(
      world_name = "compcancer_cv",
      dataset_id = dataset_id,
      split_id = sp$split_id,
      train_origin = "cv_train",
      test_origin = "cv_test"
    )
    write_json(split_meta, file.path(split_dir, "split_meta.json"), auto_unbox = TRUE, pretty = TRUE)

    for (kmax in kmax_values) {
      model_name <- paste0("kTSP_adaptive_kmax", kmax)
      model_dir <- file.path(split_dir, "baselines_adaptive", model_name)
      ensure_dir(model_dir)

      cat("  -> adaptive kTSP, kmax =", kmax, " | mode =", search_mode, "\n")

      ok <- TRUE
      res <- tryCatch({
        fit_ktsp_adaptive_train_test(
          X_train = X_train,
          y_train = y_train,
          X_test = X_test,
          y_test = y_test,
          kmax = kmax,
          out_dir = model_dir,
          inner_n_folds = inner_n_folds,
          selection_metric = selection_metric,
          seed = seed,
          search_mode = search_mode,
          pre_top_features = pre_top_features,
          extra_metrics = list(
            world_name = "compcancer_cv",
            dataset_id = dataset_id,
            split_id = sp$split_id
          )
        )
      }, error = function(e) {
        ok <<- FALSE
        write_error_json(file.path(model_dir, "error.json"), e)
        NULL
      })

      manifest_rows[[mi]] <- data.frame(
        world_name = "compcancer_cv",
        dataset_id = dataset_id,
        split_id = sp$split_id,
        family = "baseline_adaptive",
        model = model_name,
        path = model_dir,
        stringsAsFactors = FALSE
      )
      mi <- mi + 1L

      if (ok && !is.null(res)) {
        all_split_rows[[ri]] <- data.frame(
          world_name = "compcancer_cv",
          dataset_id = dataset_id,
          split_id = sp$split_id,
          model = res$method,
          requested_kmax = as.numeric(res$requested_kmax),
          selected_k = as.numeric(res$selected_k),
          realized_k = as.numeric(res$realized_k),
          auc = as.numeric(res$auc),
          f1 = as.numeric(res$f1),
          logloss = as.numeric(res$logloss),
          brier = as.numeric(res$brier),
          ece = as.numeric(res$ece),
          global_active_pairs = as.numeric(res$global_active_pairs),
          unique_genes = as.numeric(res$unique_genes),
          local_active_mean = as.numeric(res$local_active_mean),
          local_active_median = as.numeric(res$local_active_median),
          local_active_iqr = as.numeric(res$local_active_iqr),
          hub_top_share = ifelse(is.null(res$hub_top_share), NA_real_, as.numeric(res$hub_top_share)),
          stringsAsFactors = FALSE
        )
        ri <- ri + 1L
      }
    }
  }
}

manifest_df <- do.call(rbind, manifest_rows)
write.csv(manifest_df, file.path(out_root, "config", "run_manifest.csv"), row.names = FALSE)

if (length(all_split_rows) == 0) {
  stop("Brak poprawnych wynikow adaptive kTSP.")
}

split_df <- do.call(rbind, all_split_rows)
write.csv(split_df, file.path(out_root, "summaries", "all_split_metrics.csv"), row.names = FALSE)

dataset_df <- aggregate_dataset_level(split_df)
write.csv(dataset_df, file.path(out_root, "summaries", "dataset_level_metrics.csv"), row.names = FALSE)

world_df <- aggregate_world_level(dataset_df)
write.csv(world_df, file.path(out_root, "summaries", "world_level_metrics.csv"), row.names = FALSE)

cat("\nDONE\n")
cat("Results saved in:\n")
cat(out_root, "\n")