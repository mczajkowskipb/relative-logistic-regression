get_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  getwd()
}
source(file.path(get_script_dir(), "07_world_protocol_utils.R"))
suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Użycie: Rscript 16_size_matched_trim_analysis.R <config.json>")
}

cfg <- read_json(args[1], simplifyVector = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

normalize_worlds <- function(worlds_raw) {
  if (is.null(worlds_raw)) return(list())

  if (is.data.frame(worlds_raw)) {
    return(lapply(seq_len(nrow(worlds_raw)), function(i) {
      as.list(worlds_raw[i, , drop = FALSE])
    }))
  }

  if (is.list(worlds_raw) &&
      !is.null(names(worlds_raw)) &&
      all(c("name", "protocol") %in% names(worlds_raw))) {
    return(list(worlds_raw))
  }

  if (is.list(worlds_raw)) {
    return(worlds_raw)
  }

  stop("cfg$worlds ma nieobslugiwany typ")
}

normalize_analyses <- function(analyses_raw) {
  if (is.null(analyses_raw)) return(list())

  if (is.data.frame(analyses_raw)) {
    return(lapply(seq_len(nrow(analyses_raw)), function(i) {
      as.list(analyses_raw[i, , drop = FALSE])
    }))
  }

  if (is.list(analyses_raw) &&
      !is.null(names(analyses_raw)) &&
      all(c("source_model", "target_k") %in% names(analyses_raw))) {
    return(list(analyses_raw))
  }

  if (is.list(analyses_raw)) {
    return(analyses_raw)
  }

  stop("cfg$size_match_report$analyses ma nieobslugiwany typ")
}

cfg_worlds <- normalize_worlds(cfg$worlds)

out_root <- cfg$out_root
if (is.null(out_root)) stop("Brak out_root w config.json")

manifest_path <- file.path(out_root, "config", "run_manifest.csv")
if (!file.exists(manifest_path)) stop("Brak run_manifest.csv")

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)

default_analyses <- list(
  list(source_model = "RLR_K10_m10_d2_pool1000", target_k = 10, analysis_name = "trim10"),
  list(source_model = "RLR_K50_m10_d2_pool1000", target_k = 50, analysis_name = "trim50")
)

size_match_report <- cfg$size_match_report %||% list()

analyses <- normalize_analyses(size_match_report$analyses %||% default_analyses)
worlds_to_use <- size_match_report$worlds %||% unique(manifest$world_name)

trim_dir <- file.path(out_root, "paper_size_matched")
dir.create(trim_dir, recursive = TRUE, showWarnings = FALSE)

get_world_cfg <- function(cfg_worlds, world_name) {
  idx <- which(vapply(cfg_worlds, function(x) identical(x$name, world_name), logical(1)))
  if (length(idx) == 0) stop(paste("Nie znalazlem world_name w config:", world_name))
  cfg_worlds[[idx[1]]]
}

read_samples_test <- function(split_path) {
  f <- file.path(dirname(dirname(split_path)), "samples_test.csv")
  if (!file.exists(f)) stop(paste("Brak samples_test.csv dla:", split_path))
  read.csv(f, stringsAsFactors = FALSE)$sample_id
}

read_split_meta <- function(split_path) {
  f <- file.path(dirname(dirname(split_path)), "split_meta.json")
  if (!file.exists(f)) stop(paste("Brak split_meta.json dla:", split_path))
  read_json(f, simplifyVector = TRUE)
}

reconstruct_Xy_test <- function(cfg, world_name, dataset_id, model_path) {
  world <- get_world_cfg(cfg_worlds, world_name)
  protocol <- world$protocol

  split_path <- dirname(dirname(model_path))  # .../splits/<split_id>
  test_ids <- read_samples_test(model_path)
  split_meta <- read_split_meta(model_path)

  if (protocol == "compcancer_cv") {
    arff_path <- file.path(world$comp_dir, paste0(dataset_id, ".arff"))
    xy <- read_xy_file(
      path = arff_path,
      class_col = world$class_col %||% NULL,
      sample_id_col = world$sample_id_col %||% NULL
    )
    X_test <- xy$X[test_ids, , drop = FALSE]
    y_test <- xy$y[match(test_ids, rownames(xy$X))]
    return(list(X_test = X_test, y_test = y_test))
  }

  if (protocol == "fixed_split") {
    xy <- read_xy_file(
      path = world$test_file,
      class_col = world$class_col %||% NULL,
      sample_id_col = world$sample_id_col %||% NULL
    )
    X_test <- xy$X[test_ids, , drop = FALSE]
    y_test <- xy$y[match(test_ids, rownames(xy$X))]
    return(list(X_test = X_test, y_test = y_test))
  }

  if (protocol == "loco_files") {
    test_origin <- split_meta$test_origin
    cohort_files <- world$cohort_files %||% NULL
    if (is.null(cohort_files)) {
      cohort_files <- list.files(world$cohort_dir, pattern = world$cohort_pattern %||% "\\.arff$", full.names = TRUE)
    }
    hit <- cohort_files[basename(cohort_files) == test_origin]
    if (length(hit) != 1) stop(paste("Nie mogę jednoznacznie znaleźć holdout file:", test_origin))

    xy <- read_xy_file(
      path = hit[1],
      class_col = world$class_col %||% NULL,
      sample_id_col = world$sample_id_col %||% NULL
    )
    X_test <- xy$X[test_ids, , drop = FALSE]
    y_test <- xy$y[match(test_ids, rownames(xy$X))]
    return(list(X_test = X_test, y_test = y_test))
  }

  stop(paste("Nieobsługiwany protocol:", protocol))
}

parse_pairs_from_coefficients <- function(coeff_df) {
  pp <- do.call(rbind, strsplit(coeff_df$ordered_pair, " > ", fixed = TRUE))
  coeff_df$gene_left <- pp[, 1]
  coeff_df$gene_right <- pp[, 2]
  coeff_df
}

calc_metrics_from_prob <- function(prob, y_test) {
  y_test01 <- as.numeric(y_test == levels(y_test)[2])

  pred <- ifelse(prob >= 0.5, levels(y_test)[2], levels(y_test)[1])
  pred <- factor(pred, levels = levels(y_test))

  eps <- 1e-15
  p <- pmin(pmax(prob, eps), 1 - eps)

  auc <- as.numeric(pROC::auc(y_test01, prob))
  logloss <- -mean(y_test01 * log(p) + (1 - y_test01) * log(1 - p))
  brier <- calc_brier(prob, y_test01)
  ece <- calc_ece(prob, y_test01)

  tp <- sum(pred == levels(y_test)[2] & y_test == levels(y_test)[2])
  fp <- sum(pred == levels(y_test)[2] & y_test == levels(y_test)[1])
  fn <- sum(pred == levels(y_test)[1] & y_test == levels(y_test)[2])

  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))

  list(
    auc = auc,
    f1 = f1,
    logloss = logloss,
    brier = brier,
    ece = ece,
    pred = pred
  )
}

all_rows <- list()

for (aa in analyses) {
  source_model <- aa$source_model
  target_k <- as.integer(aa$target_k)
  analysis_name <- aa$analysis_name %||% paste0("trim", target_k)

  sel <- manifest[
    manifest$family == "rlr" &
      manifest$model == source_model &
      manifest$world_name %in% worlds_to_use,
    ,
    drop = FALSE
  ]

  if (nrow(sel) == 0) next

  for (i in seq_len(nrow(sel))) {
    rr <- sel[i, ]

    coeff_path <- file.path(rr$path, "coefficients.csv")
    cvfit_path <- file.path(rr$path, "cvfit.rds")
    metrics_path <- file.path(rr$path, "metrics_test.json")

    if (!file.exists(coeff_path) || !file.exists(cvfit_path) || !file.exists(metrics_path)) next

    coeff_df <- read.csv(coeff_path, stringsAsFactors = FALSE, check.names = FALSE)
    coeff_df$beta <- suppressWarnings(as.numeric(coeff_df$beta))
    coeff_df$abs_beta <- suppressWarnings(as.numeric(coeff_df$abs_beta))
    coeff_df <- coeff_df[!is.na(coeff_df$abs_beta), , drop = FALSE]
    coeff_df <- parse_pairs_from_coefficients(coeff_df)

    nonzero <- coeff_df[coeff_df$abs_beta > 0, , drop = FALSE]
    if (nrow(nonzero) == 0) next

    nonzero <- nonzero[order(nonzero$abs_beta, decreasing = TRUE), , drop = FALSE]
    trimmed <- head(nonzero, target_k)

    mjson <- read_json(metrics_path, simplifyVector = TRUE)
    lambda_choice <- mjson$lambda_choice %||% "lambda.min"
    m_value <- as.numeric(mjson$m %||% NA)

    cvfit <- readRDS(cvfit_path)
    intercept <- as.numeric(coef(cvfit, s = lambda_choice))[1]

    test_xy <- reconstruct_Xy_test(cfg, rr$world_name, rr$dataset_id, rr$path)
    X_test <- test_xy$X_test
    y_test <- test_xy$y_test

    pairs_df <- trimmed[, c("ordered_pair", "gene_left", "gene_right"), drop = FALSE]
    R_test <- encode_pair_matrix(X_test, pairs_df, m = m_value)

    eta <- intercept + as.numeric(as.matrix(R_test) %*% trimmed$beta)
    prob <- 1 / (1 + exp(-eta))

    mm <- calc_metrics_from_prob(prob, y_test)

    local_sum <- calc_local_activity_summary(R_test, active_mask = rep(TRUE, ncol(R_test)))
    hubs <- calc_hub_summary(trimmed)

    row <- data.frame(
      world_name = rr$world_name,
      dataset_id = rr$dataset_id,
      split_id = rr$split_id,
      source_model = source_model,
      trimmed_model = paste0(source_model, "_", analysis_name),
      target_k = target_k,
      auc = mm$auc,
      f1 = mm$f1,
      logloss = mm$logloss,
      brier = mm$brier,
      ece = mm$ece,
      global_active_pairs = nrow(trimmed),
      unique_genes = length(unique(c(trimmed$gene_left, trimmed$gene_right))),
      local_active_mean = local_sum$mean,
      local_active_median = local_sum$median,
      local_active_iqr = local_sum$iqr,
      hub_top_share = hubs$top_share,
      stringsAsFactors = FALSE
    )

    all_rows[[length(all_rows) + 1]] <- row
  }
}

if (length(all_rows) == 0) stop("Brak wyników size-matched trim.")

trim_split_df <- do.call(rbind, all_rows)
write.csv(trim_split_df, file.path(trim_dir, "all_split_size_matched_metrics.csv"), row.names = FALSE)

# dataset-level means
keys <- unique(trim_split_df[, c("world_name", "dataset_id", "trimmed_model", "target_k"), drop = FALSE])

trim_dataset_df <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
  key <- keys[i, , drop = FALSE]
  idx <- trim_split_df$world_name == key$world_name &
         trim_split_df$dataset_id == key$dataset_id &
         trim_split_df$trimmed_model == key$trimmed_model

  data.frame(
    world_name = key$world_name,
    dataset_id = key$dataset_id,
    trimmed_model = key$trimmed_model,
    target_k = key$target_k,
    auc = mean(trim_split_df$auc[idx], na.rm = TRUE),
    f1 = mean(trim_split_df$f1[idx], na.rm = TRUE),
    logloss = mean(trim_split_df$logloss[idx], na.rm = TRUE),
    brier = mean(trim_split_df$brier[idx], na.rm = TRUE),
    ece = mean(trim_split_df$ece[idx], na.rm = TRUE),
    global_active_pairs = mean(trim_split_df$global_active_pairs[idx], na.rm = TRUE),
    unique_genes = mean(trim_split_df$unique_genes[idx], na.rm = TRUE),
    local_active_mean = mean(trim_split_df$local_active_mean[idx], na.rm = TRUE),
    local_active_median = mean(trim_split_df$local_active_median[idx], na.rm = TRUE),
    local_active_iqr = mean(trim_split_df$local_active_iqr[idx], na.rm = TRUE),
    hub_top_share = mean(trim_split_df$hub_top_share[idx], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

write.csv(trim_dataset_df, file.path(trim_dir, "dataset_level_size_matched_metrics.csv"), row.names = FALSE)

cat("OK: zapisano size-matched trim analysis do:\n")
cat(trim_dir, "\n")