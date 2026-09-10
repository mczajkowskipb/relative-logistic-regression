get_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  getwd()
}
source(file.path(get_script_dir(), "07_world_protocol_utils.R"))

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
  library(jsonlite)
  library(switchBox)
  library(foreign)
})

# ============================================================
# CONFIG
# ============================================================
project_root <- Sys.getenv("RLR_PROJECT_ROOT", unset = normalizePath(file.path(get_script_dir(), "../..")))
comp_dir <- Sys.getenv("RLR_COMPCANCER_DIR", unset = file.path(project_root, "data", "CompCancer"))
out_root <- Sys.getenv("RLR_OUT_ROOT", unset = file.path(project_root, "results", "supplementary", "compcancer_runtime_profile"))

# NULL = all datasets from comp_dir
datasets_include <- NULL
# Example:
# datasets_include <- c("alizadeh-2000-v1", "chen-2002", "gordon-2002")

split_pick <- "rep00_fold00"
seed <- 42

pool_top_n <- 1000
pool_buffer_mult <- 5

rlr_models <- list(
  list(K = 10, m = 10, d = 2, model = "RLR_K10_m10_d2_pool1000"),
  list(K = 50, m = 10, d = 2, model = "RLR_K50_m10_d2_pool1000")
)

ktsp_k_values <- c(10, 50)

# ============================================================
# SETUP
# ============================================================
ensure_dir(out_root)
ensure_dir(file.path(out_root, "logs"))
ensure_dir(file.path(out_root, "tables"))

log_path <- file.path(out_root, "logs", paste0("runtime_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cleanup_sinks <- function() {
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number() > 0) sink()
  try(close(log_con), silent = TRUE)
}
on.exit(cleanup_sinks(), add = TRUE)

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

mean_or_na <- function(x) {
  if (length(x) == 0) return(NA_real_)
  mean(x, na.rm = TRUE)
}

median_or_na <- function(x) {
  if (length(x) == 0) return(NA_real_)
  median(x, na.rm = TRUE)
}

time_expr <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  t1 <- proc.time()[["elapsed"]]
  list(value = value, sec = as.numeric(t1 - t0))
}

# ============================================================
# HELPERS
# ============================================================
pick_outer_split <- function(y, split_id = "rep00_fold00", n_folds = 5, n_repeats = 10, seed = 42) {
  folds <- make_repeated_stratified_folds(
    y = y,
    n_folds = n_folds,
    n_repeats = n_repeats,
    seed = seed
  )
  hit <- which(vapply(folds, function(x) identical(x$split_id, split_id), logical(1)))
  if (length(hit) != 1) stop("Nie znaleziono split_id: ", split_id)
  folds[[hit]]
}

benchmark_rlr_one <- function(X_train, y_train, X_test, y_test, pairs_pool, pool_time_sec, spec) {
  model_name <- spec$model
  K <- spec$K
  m <- spec$m
  d <- spec$d

  stage_rows <- list()
  sri <- 1L

  shortlist_t <- time_expr(
    apply_disjointness(pairs_pool, K = K, d = d)
  )
  shortlist <- shortlist_t$value

  encode_train_t <- time_expr(
    encode_pair_matrix(X_train, shortlist, m = m)
  )
  R_train <- encode_train_t$value

  encode_test_t <- time_expr(
    encode_pair_matrix(X_test, shortlist, m = m)
  )
  R_test <- encode_test_t$value

  y_train01 <- as.numeric(y_train == levels(y_train)[2])
  y_test01  <- as.numeric(y_test == levels(y_test)[2])

  fit_t <- time_expr(
    cv.glmnet(
      x = R_train,
      y = y_train01,
      family = "binomial",
      alpha = 0.5,
      nfolds = safe_nfolds(y_train, 5),
      type.measure = "deviance"
    )
  )
  cvfit <- fit_t$value

  infer_t <- time_expr({
    prob <- as.numeric(predict(cvfit, newx = R_test, s = "lambda.min", type = "response"))
    pred <- ifelse(prob >= 0.5, levels(y_test)[2], levels(y_test)[1])
    pred <- factor(pred, levels = levels(y_test))

    eps <- 1e-15
    p <- pmin(pmax(prob, eps), 1 - eps)

    auc <- as.numeric(pROC::auc(y_test01, prob, quiet = TRUE))
    logloss <- -mean(y_test01 * log(p) + (1 - y_test01) * log(1 - p))
    brier <- calc_brier(prob, y_test01)
    ece <- calc_ece(prob, y_test01)

    tp <- sum(pred == levels(y_test)[2] & y_test == levels(y_test)[2])
    fp <- sum(pred == levels(y_test)[2] & y_test == levels(y_test)[1])
    fn <- sum(pred == levels(y_test)[1] & y_test == levels(y_test)[2])

    precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))

    beta <- as.numeric(coef(cvfit, s = "lambda.min"))[-1]
    names(beta) <- colnames(R_train)
    active_mask <- beta != 0
    active_pairs <- shortlist[colnames(R_train) %in% names(beta)[active_mask], , drop = FALSE]
    unique_genes <- if (nrow(active_pairs) == 0) 0L else length(unique(c(active_pairs$gene_left, active_pairs$gene_right)))
    local_sum <- calc_local_activity_summary(R_test, active_mask = active_mask)
    hubs <- calc_hub_summary(active_pairs)

    list(
      auc = auc,
      f1 = f1,
      logloss = logloss,
      brier = brier,
      ece = ece,
      global_active_pairs = sum(active_mask),
      unique_genes = unique_genes,
      local_active_mean = local_sum$mean,
      hub_top_share = hubs$top_share
    )
  })
  metrics <- infer_t$value

  stage_rows[[sri]] <- data.frame(model = model_name, stage = "pair_ranking_shared", time_sec = pool_time_sec, stringsAsFactors = FALSE); sri <- sri + 1L
  stage_rows[[sri]] <- data.frame(model = model_name, stage = "shortlist_refinement", time_sec = shortlist_t$sec, stringsAsFactors = FALSE); sri <- sri + 1L
  stage_rows[[sri]] <- data.frame(model = model_name, stage = "encoding_train", time_sec = encode_train_t$sec, stringsAsFactors = FALSE); sri <- sri + 1L
  stage_rows[[sri]] <- data.frame(model = model_name, stage = "encoding_test", time_sec = encode_test_t$sec, stringsAsFactors = FALSE); sri <- sri + 1L
  stage_rows[[sri]] <- data.frame(model = model_name, stage = "final_fit_cv", time_sec = fit_t$sec, stringsAsFactors = FALSE); sri <- sri + 1L
  stage_rows[[sri]] <- data.frame(model = model_name, stage = "inference", time_sec = infer_t$sec, stringsAsFactors = FALSE); sri <- sri + 1L

  stage_df <- do.call(rbind, stage_rows)

  model_df <- data.frame(
    model = model_name,
    family = "RLR",
    requested_k = K,
    realized_k = metrics$global_active_pairs,
    auc = metrics$auc,
    f1 = metrics$f1,
    logloss = metrics$logloss,
    brier = metrics$brier,
    ece = metrics$ece,
    global_active_pairs = metrics$global_active_pairs,
    unique_genes = metrics$unique_genes,
    local_active_mean = metrics$local_active_mean,
    hub_top_share = metrics$hub_top_share,
    pair_ranking_shared_sec = pool_time_sec,
    model_specific_sec = shortlist_t$sec + encode_train_t$sec + encode_test_t$sec + fit_t$sec + infer_t$sec,
    total_with_shared_sec = pool_time_sec + shortlist_t$sec + encode_train_t$sec + encode_test_t$sec + fit_t$sec + infer_t$sec,
    stringsAsFactors = FALSE
  )

  list(stage_df = stage_df, model_df = model_df)
}

benchmark_ktsp_one <- function(X_train, y_train, X_test, y_test, k) {
  model_name <- paste0("kTSP_", k)

  train_t <- time_expr({
    trainMat <- xy_to_switchbox_mat(X_train)
    SWAP.KTSP.Train(
      inputMat = trainMat,
      phenoGroup = y_train,
      FilterFunc = NULL,
      krange = k
    )
  })
  clf <- train_t$value

  if (!is.matrix(clf$TSPs)) {
    clf$TSPs <- matrix(clf$TSPs, ncol = 2, byrow = TRUE)
  }

  infer_t <- time_expr({
    trainMat <- xy_to_switchbox_mat(X_train)
    testMat  <- xy_to_switchbox_mat(X_test)

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

    auc <- as.numeric(pROC::auc(y_test01, prob_test, quiet = TRUE))
    logloss <- -mean(y_test01 * log(p) + (1 - y_test01) * log(1 - p))
    brier <- calc_brier(prob_test, y_test01)
    ece <- calc_ece(prob_test, y_test01)

    tp <- sum(pred == levels(y_test)[2] & y_test == levels(y_test)[2])
    fp <- sum(pred == levels(y_test)[2] & y_test == levels(y_test)[1])
    fn <- sum(pred == levels(y_test)[1] & y_test == levels(y_test)[2])

    precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))

    selected_pairs <- as.data.frame(clf$TSPs, stringsAsFactors = FALSE)
    colnames(selected_pairs) <- c("gene_left", "gene_right")
    unique_genes <- length(unique(c(selected_pairs$gene_left, selected_pairs$gene_right)))
    hubs <- calc_hub_summary(selected_pairs)

    list(
      realized_k = realized_k,
      auc = auc,
      f1 = f1,
      logloss = logloss,
      brier = brier,
      ece = ece,
      global_active_pairs = realized_k,
      unique_genes = unique_genes,
      local_active_mean = realized_k,
      hub_top_share = hubs$top_share
    )
  })
  metrics <- infer_t$value

  stage_df <- rbind(
    data.frame(model = model_name, stage = "train", time_sec = train_t$sec, stringsAsFactors = FALSE),
    data.frame(model = model_name, stage = "inference", time_sec = infer_t$sec, stringsAsFactors = FALSE)
  )

  model_df <- data.frame(
    model = model_name,
    family = "kTSP",
    requested_k = k,
    realized_k = metrics$realized_k,
    auc = metrics$auc,
    f1 = metrics$f1,
    logloss = metrics$logloss,
    brier = metrics$brier,
    ece = metrics$ece,
    global_active_pairs = metrics$global_active_pairs,
    unique_genes = metrics$unique_genes,
    local_active_mean = metrics$local_active_mean,
    hub_top_share = metrics$hub_top_share,
    pair_ranking_shared_sec = NA_real_,
    model_specific_sec = train_t$sec + infer_t$sec,
    total_with_shared_sec = train_t$sec + infer_t$sec,
    stringsAsFactors = FALSE
  )

  list(stage_df = stage_df, model_df = model_df)
}

# ============================================================
# MAIN
# ============================================================
files <- list.files(comp_dir, pattern = "\\.arff$", full.names = TRUE)
files <- sort(files)

if (!is.null(datasets_include)) {
  files <- files[basename(files) %in% paste0(datasets_include, ".arff")]
}

if (length(files) == 0) stop("No datasets found.")

all_stage <- list()
all_model <- list()
si <- 1L
mi <- 1L

for (f in files) {
  dataset_id <- sub("\\.arff$", "", basename(f))
  log_msg("====================================================")
  log_msg("DATASET: ", dataset_id)

  xy <- read_xy_file(f)
  X <- xy$X
  y <- xy$y

  sp <- pick_outer_split(
    y = y,
    split_id = split_pick,
    n_folds = 5,
    n_repeats = 10,
    seed = seed
  )

  X_train <- X[sp$train_idx, , drop = FALSE]
  y_train <- y[sp$train_idx]
  X_test  <- X[sp$test_idx, , drop = FALSE]
  y_test  <- y[sp$test_idx]

  log_msg("Train N=", nrow(X_train), " | Test N=", nrow(X_test), " | P=", ncol(X_train))

  # Shared RLR pair ranking
  log_msg("Running shared pair ranking top", pool_top_n)
  pool_t <- time_expr(
    compute_top_pairs_switchbox(
      X_train = X_train,
      y_train = y_train,
      top_n = pool_top_n,
      buffer_mult = pool_buffer_mult,
      out_dir = NULL,
      save_all = FALSE
    )
  )
  pairs_pool <- pool_t$value
  pool_time_sec <- pool_t$sec
  log_msg("Shared pair ranking time: ", round(pool_time_sec, 3), " sec")

  # RLR models
  for (spec in rlr_models) {
    log_msg("Benchmarking ", spec$model)
    out <- benchmark_rlr_one(X_train, y_train, X_test, y_test, pairs_pool, pool_time_sec, spec)

    st <- out$stage_df
    st$dataset_id <- dataset_id
    st$split_id <- split_pick
    st$n_train <- nrow(X_train)
    st$n_test <- nrow(X_test)
    st$p <- ncol(X_train)
    all_stage[[si]] <- st
    si <- si + 1L

    md <- out$model_df
    md$dataset_id <- dataset_id
    md$split_id <- split_pick
    md$n_train <- nrow(X_train)
    md$n_test <- nrow(X_test)
    md$p <- ncol(X_train)
    all_model[[mi]] <- md
    mi <- mi + 1L
  }

  # kTSP baselines
  for (k in ktsp_k_values) {
    log_msg("Benchmarking kTSP_", k)
    out <- benchmark_ktsp_one(X_train, y_train, X_test, y_test, k)

    st <- out$stage_df
    st$dataset_id <- dataset_id
    st$split_id <- split_pick
    st$n_train <- nrow(X_train)
    st$n_test <- nrow(X_test)
    st$p <- ncol(X_train)
    all_stage[[si]] <- st
    si <- si + 1L

    md <- out$model_df
    md$dataset_id <- dataset_id
    md$split_id <- split_pick
    md$n_train <- nrow(X_train)
    md$n_test <- nrow(X_test)
    md$p <- ncol(X_train)
    all_model[[mi]] <- md
    mi <- mi + 1L
  }
}

stage_df <- do.call(rbind, all_stage)
model_df <- do.call(rbind, all_model)

write.csv(stage_df, file.path(out_root, "tables", "runtime_stage_level.csv"), row.names = FALSE)
write.csv(model_df, file.path(out_root, "tables", "runtime_model_level.csv"), row.names = FALSE)

# World summaries
stage_summary <- do.call(rbind, lapply(split(stage_df, list(stage_df$model, stage_df$stage), drop = TRUE), function(dd) {
  data.frame(
    model = dd$model[1],
    stage = dd$stage[1],
    mean_time_sec = mean_or_na(dd$time_sec),
    median_time_sec = median_or_na(dd$time_sec),
    min_time_sec = min(dd$time_sec, na.rm = TRUE),
    max_time_sec = max(dd$time_sec, na.rm = TRUE),
    n_datasets = nrow(dd),
    stringsAsFactors = FALSE
  )
}))
write.csv(stage_summary, file.path(out_root, "tables", "runtime_stage_summary.csv"), row.names = FALSE)

model_summary <- do.call(rbind, lapply(split(model_df, model_df$model), function(dd) {
  data.frame(
    model = dd$model[1],
    family = dd$family[1],
    requested_k = dd$requested_k[1],
    realized_k_mean = mean_or_na(dd$realized_k),
    auc = mean_or_na(dd$auc),
    f1 = mean_or_na(dd$f1),
    logloss = mean_or_na(dd$logloss),
    brier = mean_or_na(dd$brier),
    ece = mean_or_na(dd$ece),
    mean_model_specific_sec = mean_or_na(dd$model_specific_sec),
    mean_total_with_shared_sec = mean_or_na(dd$total_with_shared_sec),
    median_total_with_shared_sec = median_or_na(dd$total_with_shared_sec),
    mean_pair_ranking_shared_sec = mean_or_na(dd$pair_ranking_shared_sec),
    n_datasets = nrow(dd),
    stringsAsFactors = FALSE
  )
}))
write.csv(model_summary, file.path(out_root, "tables", "runtime_model_summary.csv"), row.names = FALSE)

# Markdown appendix draft
md_path <- file.path(out_root, "tables", "appendix_runtime_notes.md")
cat(
  "# Runtime / complexity appendix draft\n\n",
  "This appendix reports approximate wall-clock timings measured on a single representative outer split (`", split_pick, "`) per CompCancer dataset.\n\n",
  "The measurements are intended as an engineering/runtime characterization rather than a full benchmarking study. ",
  "For RLR, stage-level timings are reported for:\n",
  "- pair ranking (`top", pool_top_n, "` candidate pool),\n",
  "- shortlist refinement,\n",
  "- train encoding,\n",
  "- test encoding,\n",
  "- final penalized logistic fit,\n",
  "- inference.\n\n",
  "For classical kTSP baselines, total train and inference times are reported.\n\n",
  "Important note: for RLR, the pair-ranking stage is **shared** across multiple downstream model variants on the same split, so the total cost can be interpreted in two ways:\n",
  "- `total_with_shared_sec`: full cost when ranking is attributed to the model,\n",
  "- `model_specific_sec`: incremental downstream cost once the ranked pool is available.\n\n",
  "Files generated:\n",
  "- `runtime_stage_level.csv`\n",
  "- `runtime_model_level.csv`\n",
  "- `runtime_stage_summary.csv`\n",
  "- `runtime_model_summary.csv`\n\n",
  "In qualitative terms, the stage expected to dominate RLR runtime is pair ranking, while shortlist refinement and inference are comparatively cheap. ",
  "For kTSP, training dominates inference, but total runtime remains model-size dependent.\n",
  file = md_path
)

log_msg("DONE")
log_msg("Saved stage-level table: ", file.path(out_root, "tables", "runtime_stage_level.csv"))
log_msg("Saved model-level table: ", file.path(out_root, "tables", "runtime_model_level.csv"))
log_msg("Saved stage summary: ", file.path(out_root, "tables", "runtime_stage_summary.csv"))
log_msg("Saved model summary: ", file.path(out_root, "tables", "runtime_model_summary.csv"))
log_msg("Saved appendix draft: ", md_path)
log_msg("Log file: ", log_path)