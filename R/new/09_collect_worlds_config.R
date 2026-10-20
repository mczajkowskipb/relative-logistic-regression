suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Użycie: Rscript 09_collect_worlds_config.R <config.json>")
}

cfg <- read_json(args[1], simplifyVector = TRUE)
out_root <- cfg$out_root
if (is.null(out_root)) stop("Brak out_root w config.json")

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

read_json_scalar_df <- function(path) {
  x <- tryCatch(read_json(path, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(x)) return(NULL)

  keep <- vapply(x, function(v) length(v) == 1 && !is.list(v), logical(1))
  x2 <- x[keep]
  if (length(x2) == 0) return(NULL)

  df <- as.data.frame(x2, stringsAsFactors = FALSE, check.names = FALSE)
  df$file <- path
  df
}

rbind_fill_named <- function(df_list) {
  df_list <- df_list[!vapply(df_list, is.null, logical(1))]
  if (length(df_list) == 0) return(NULL)

  all_names <- unique(unlist(lapply(df_list, names)))

  df_list2 <- lapply(df_list, function(df) {
    miss <- setdiff(all_names, names(df))
    for (m in miss) df[[m]] <- NA
    df <- df[, all_names, drop = FALSE]
    df
  })

  do.call(rbind, df_list2)
}

extract_from_path <- function(x, pattern, default = NA_character_) {
  out <- sub(pattern, "\\1", x)
  bad <- identical(out, x)
  if (length(bad) == 1 && bad) return(rep(default, length(x)))
  out
}

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

safe_iqr <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  IQR(x, na.rm = TRUE)
}

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

# ------------------------------------------------------------
# 1. metrics_test.json
# ------------------------------------------------------------
metric_files <- list.files(
  out_root,
  pattern = "metrics_test\\.json$",
  recursive = TRUE,
  full.names = TRUE
)

metrics_list <- lapply(metric_files, read_json_scalar_df)
metrics_df <- rbind_fill_named(metrics_list)

if (is.null(metrics_df) || nrow(metrics_df) == 0) {
  stop("Nie znaleziono poprawnych metrics_test.json")
}

# metadata z path
metrics_df$world_name <- if ("world_name" %in% names(metrics_df)) {
  as.character(metrics_df$world_name)
} else {
  extract_from_path(metrics_df$file, ".*worlds[\\\\/](.*?)[\\\\/]datasets.*", NA_character_)
}

metrics_df$dataset_id <- if ("dataset_id" %in% names(metrics_df)) {
  as.character(metrics_df$dataset_id)
} else {
  extract_from_path(metrics_df$file, ".*datasets[\\\\/](.*?)[\\\\/]splits.*", NA_character_)
}

metrics_df$split_id <- if ("split_id" %in% names(metrics_df)) {
  as.character(metrics_df$split_id)
} else {
  extract_from_path(metrics_df$file, ".*splits[\\\\/](.*?)[\\\\/](baselines|rlr|rlmt).*", NA_character_)
}

metrics_df$family <- if ("method_family" %in% names(metrics_df)) {
  as.character(metrics_df$method_family)
} else {
  ifelse(grepl("[\\\\/]baselines[\\\\/]", metrics_df$file), "baseline", "RLR")
}

metrics_df$model <- if ("method" %in% names(metrics_df)) {
  as.character(metrics_df$method)
} else {
  rep(NA_character_, nrow(metrics_df))
}

path_model <- extract_from_path(
  metrics_df$file,
  ".*[\\/](TSP_1|kTSP_[0-9]+|(RLR_K[0-9]+_m[0-9]+_d[0-9]+(_pool[0-9]+)?|RLMT_K[0-9]+_m[0-9]+_d[0-9]+_pool[0-9]+_L[0-9]+)|ValueLR_MAD1000_EN|value_LR)[\\/].*",
  NA_character_
)

miss_model <- is.na(metrics_df$model) | metrics_df$model == ""
metrics_df$model[miss_model] <- path_model[miss_model]

is_rlr_path <- !is.na(path_model) & grepl("^RLR_K[0-9]+_m[0-9]+_d[0-9]+(?:_pool[0-9]+)?$", path_model)
metrics_df$model[is_rlr_path] <- path_model[is_rlr_path]

metrics_df$model[metrics_df$model %in% c("ValueLR_MAD1000_EN", "value_LR")] <- "value_LR"

# rozbij model na pola techniczne
metrics_df$pool_cfg <- NA
metrics_df$K_cfg <- NA
metrics_df$m_cfg <- NA
metrics_df$d_cfg <- NA

is_rlr_pool <- grepl("^RLR_K[0-9]+_m[0-9]+_d[0-9]+_pool[0-9]+$", metrics_df$model)
if (any(is_rlr_pool, na.rm = TRUE)) {
  mm <- regexec("^RLR_K([0-9]+)_m([0-9]+)_d([0-9]+)_pool([0-9]+)$", metrics_df$model[is_rlr_pool])
  rr <- regmatches(metrics_df$model[is_rlr_pool], mm)
  metrics_df$K_cfg[is_rlr_pool] <- sapply(rr, `[`, 2)
  metrics_df$m_cfg[is_rlr_pool] <- sapply(rr, `[`, 3)
  metrics_df$d_cfg[is_rlr_pool] <- sapply(rr, `[`, 4)
  metrics_df$pool_cfg[is_rlr_pool] <- sapply(rr, `[`, 5)
}

is_rlr_no_pool <- grepl("^RLR_K[0-9]+_m[0-9]+_d[0-9]+$", metrics_df$model)
if (any(is_rlr_no_pool, na.rm = TRUE)) {
  mm <- regexec("^RLR_K([0-9]+)_m([0-9]+)_d([0-9]+)$", metrics_df$model[is_rlr_no_pool])
  rr <- regmatches(metrics_df$model[is_rlr_no_pool], mm)
  metrics_df$K_cfg[is_rlr_no_pool] <- sapply(rr, `[`, 2)
  metrics_df$m_cfg[is_rlr_no_pool] <- sapply(rr, `[`, 3)
  metrics_df$d_cfg[is_rlr_no_pool] <- sapply(rr, `[`, 4)
}

is_ktsp <- grepl("^kTSP_[0-9]+$", metrics_df$model)
if (any(is_ktsp, na.rm = TRUE)) {
  mm <- regexec("^kTSP_([0-9]+)$", metrics_df$model[is_ktsp])
  rr <- regmatches(metrics_df$model[is_ktsp], mm)
  metrics_df$K_cfg[is_ktsp] <- sapply(rr, `[`, 2)
}

is_tsp1 <- metrics_df$model == "TSP_1"
if (any(is_tsp1, na.rm = TRUE)) {
  metrics_df$K_cfg[is_tsp1] <- "1"
}

# Konwersje numeryczne dla typowych kolumn
num_cols <- intersect(
  c("auc","f1","logloss","brier","ece",
    "global_active_pairs","unique_genes",
    "local_active_mean","local_active_median","local_active_iqr",
    "hub_top_share","requested_k","realized_k","requested_K",
    "m","d","alpha","pool","K_cfg","m_cfg","d_cfg","pool_cfg"),
  names(metrics_df)
)

for (cc in num_cols) {
  metrics_df[[cc]] <- suppressWarnings(as.numeric(metrics_df[[cc]]))
}

summaries_dir <- file.path(out_root, "summaries")
ensure_dir(summaries_dir)

write.csv(
  metrics_df,
  file.path(summaries_dir, "all_split_metrics.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 2. local_activity_test.csv
# ------------------------------------------------------------
local_csv_files <- list.files(
  out_root,
  pattern = "local_activity_test\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

local_df_list <- lapply(local_csv_files, function(f) {
  df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df)) return(NULL)

  df$file <- f
  df$world_name <- extract_from_path(f, ".*worlds[\\\\/](.*?)[\\\\/]datasets.*", NA_character_)
  df$dataset_id <- extract_from_path(f, ".*datasets[\\\\/](.*?)[\\\\/]splits.*", NA_character_)
  df$split_id <- extract_from_path(f, ".*splits[\\\\/](.*?)[\\\\/](baselines|rlr|rlmt).*", NA_character_)
  df$model <- if (grepl("[\\/]value_lr[\\/]", f)) {
    "value_LR"
  } else {
    extract_from_path(
      f,
      ".*[\\/](TSP_1|kTSP_[0-9]+|(RLR_K[0-9]+_m[0-9]+_d[0-9]+(_pool[0-9]+)?|RLMT_K[0-9]+_m[0-9]+_d[0-9]+_pool[0-9]+_L[0-9]+)|ValueLR_MAD1000_EN|value_LR)[\\/].*",
      NA_character_
    )
  }
  df
})

local_df <- rbind_fill_named(local_df_list)

if (!is.null(local_df) && nrow(local_df) > 0) {
  local_df$local_active_n <- suppressWarnings(as.numeric(local_df$local_active_n))

  write.csv(
    local_df,
    file.path(summaries_dir, "all_split_local_activity.csv"),
    row.names = FALSE
  )

  # dataset/model-level
  split_keys <- unique(local_df[, c("world_name","dataset_id","split_id","model"), drop = FALSE])

  local_split_summary <- do.call(rbind, lapply(seq_len(nrow(split_keys)), function(i) {
    key <- split_keys[i, , drop = FALSE]
    idx <- local_df$world_name == key$world_name &
           local_df$dataset_id == key$dataset_id &
           local_df$split_id == key$split_id &
           local_df$model == key$model

    z <- local_df$local_active_n[idx]

    data.frame(
      world_name = key$world_name,
      dataset_id = key$dataset_id,
      split_id = key$split_id,
      model = key$model,
      local_mean = safe_mean(z),
      local_median = safe_median(z),
      local_iqr = safe_iqr(z),
      local_min = safe_min(z),
      local_max = safe_max(z),
      stringsAsFactors = FALSE
    )
  }))

  write.csv(
    local_split_summary,
    file.path(summaries_dir, "split_level_local_activity_summary.csv"),
    row.names = FALSE
  )

  dataset_model_keys <- unique(local_split_summary[, c("world_name","dataset_id","model"), drop = FALSE])

  local_dataset_summary <- do.call(rbind, lapply(seq_len(nrow(dataset_model_keys)), function(i) {
    key <- dataset_model_keys[i, , drop = FALSE]
    idx <- local_split_summary$world_name == key$world_name &
           local_split_summary$dataset_id == key$dataset_id &
           local_split_summary$model == key$model

    data.frame(
      world_name = key$world_name,
      dataset_id = key$dataset_id,
      model = key$model,
      local_mean = safe_mean(local_split_summary$local_mean[idx]),
      local_median = safe_mean(local_split_summary$local_median[idx]),
      local_iqr = safe_mean(local_split_summary$local_iqr[idx]),
      local_min = safe_mean(local_split_summary$local_min[idx]),
      local_max = safe_mean(local_split_summary$local_max[idx]),
      stringsAsFactors = FALSE
    )
  }))

  write.csv(
    local_dataset_summary,
    file.path(summaries_dir, "dataset_level_local_activity.csv"),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 3. dataset-level metrics
# ------------------------------------------------------------
metric_cols <- intersect(
  c("auc","f1","logloss","brier","ece",
    "global_active_pairs","unique_genes",
    "local_active_mean","local_active_median","local_active_iqr",
    "hub_top_share"),
  names(metrics_df)
)

dataset_model_keys <- unique(metrics_df[, c("world_name","dataset_id","model"), drop = FALSE])

dataset_level_metrics <- do.call(rbind, lapply(seq_len(nrow(dataset_model_keys)), function(i) {
  key <- dataset_model_keys[i, , drop = FALSE]
  idx <- metrics_df$world_name == key$world_name &
         metrics_df$dataset_id == key$dataset_id &
         metrics_df$model == key$model

  row <- data.frame(
    world_name = key$world_name,
    dataset_id = key$dataset_id,
    model = key$model,
    stringsAsFactors = FALSE
  )

  for (cc in metric_cols) {
    row[[cc]] <- safe_mean(metrics_df[[cc]][idx])
  }

  # zachowaj cfg jeśli jednoznaczne
  for (cc in c("K_cfg","m_cfg","d_cfg","pool_cfg")) {
    if (cc %in% names(metrics_df)) {
      vals <- unique(na.omit(metrics_df[[cc]][idx]))
      row[[cc]] <- if (length(vals) == 1) vals[1] else NA
    }
  }

  row
}))

write.csv(
  dataset_level_metrics,
  file.path(summaries_dir, "dataset_level_metrics.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. world-level / model-level means
# ------------------------------------------------------------
world_model_keys <- unique(dataset_level_metrics[, c("world_name","model"), drop = FALSE])

world_level_metrics <- do.call(rbind, lapply(seq_len(nrow(world_model_keys)), function(i) {
  key <- world_model_keys[i, , drop = FALSE]
  idx <- dataset_level_metrics$world_name == key$world_name &
         dataset_level_metrics$model == key$model

  row <- data.frame(
    world_name = key$world_name,
    model = key$model,
    stringsAsFactors = FALSE
  )

  for (cc in metric_cols) {
    row[[cc]] <- safe_mean(dataset_level_metrics[[cc]][idx])
  }

  for (cc in c("K_cfg","m_cfg","d_cfg","pool_cfg")) {
    if (cc %in% names(dataset_level_metrics)) {
      vals <- unique(na.omit(dataset_level_metrics[[cc]][idx]))
      row[[cc]] <- if (length(vals) == 1) vals[1] else NA
    }
  }

  row
}))

write.csv(
  world_level_metrics,
  file.path(summaries_dir, "world_level_metrics.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. run manifest passthrough jeśli istnieje
# ------------------------------------------------------------
manifest_path <- file.path(out_root, "config", "run_manifest.csv")
if (file.exists(manifest_path)) {
  manifest_df <- read.csv(manifest_path, stringsAsFactors = FALSE)
  write.csv(manifest_df, file.path(summaries_dir, "run_manifest_copy.csv"), row.names = FALSE)
}

cat("OK: zapisano wyniki do:\n")
cat(summaries_dir, "\n")
cat("Liczba metrics_test.json:", nrow(metrics_df), "\n")
cat("Liczba rekordów local_activity_test:", ifelse(is.null(local_df), 0, nrow(local_df)), "\n")