get_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  getwd()
}
source(file.path(get_script_dir(), "07_world_protocol_utils.R"))

suppressPackageStartupMessages({
  library(foreign)
  library(jsonlite)
})

# ============================================================
# CONFIG
# ============================================================
project_root <- Sys.getenv("RLR_PROJECT_ROOT", unset = normalizePath(file.path(get_script_dir(), "../..")))
comp_dir <- Sys.getenv("RLR_COMPCANCER_DIR", unset = file.path(project_root, "data", "CompCancer"))
out_root <- Sys.getenv("RLR_OUT_ROOT", unset = file.path(project_root, "results", "sensitivity", "compcancer_m_ablation"))

K_values <- c(10, 20, 50)
m_values <- c(1, 5, 10, 20, 50)
d_value <- 2

n_folds <- 5
n_repeats <- 2
seed <- 42

pool_top_n <- 1000
pool_buffer_mult <- 5

# ============================================================
# SETUP
# ============================================================
ensure_dir(out_root)
ensure_dir(file.path(out_root, "logs"))
ensure_dir(file.path(out_root, "tmp"))

log_path <- file.path(
  out_root,
  "logs",
  paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

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

safe_write_csv <- function(df, path) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  write.csv(df, path, row.names = FALSE)
}

aggregate_dataset_level <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())

  keys <- unique(df[, c("dataset_id", "K", "m", "d"), drop = FALSE])

  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]

    idx <- df$dataset_id == key$dataset_id &
      df$K == key$K &
      df$m == key$m &
      df$d == key$d

    data.frame(
      dataset_id = key$dataset_id,
      K = key$K,
      m = key$m,
      d = key$d,
      auc = mean_or_na(df$auc[idx]),
      f1 = mean_or_na(df$f1[idx]),
      logloss = mean_or_na(df$logloss[idx]),
      brier = mean_or_na(df$brier[idx]),
      ece = mean_or_na(df$ece[idx]),
      global_active_pairs = mean_or_na(df$global_active_pairs[idx]),
      unique_genes = mean_or_na(df$unique_genes[idx]),
      local_active_mean = mean_or_na(df$local_active_mean[idx]),
      hub_top_share = mean_or_na(df$hub_top_share[idx]),
      n_splits = sum(idx),
      stringsAsFactors = FALSE
    )
  }))
}

aggregate_world_level <- function(df_dataset) {
  if (is.null(df_dataset) || nrow(df_dataset) == 0) return(data.frame())

  keys <- unique(df_dataset[, c("K", "m", "d"), drop = FALSE])

  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]

    idx <- df_dataset$K == key$K &
      df_dataset$m == key$m &
      df_dataset$d == key$d

    data.frame(
      K = key$K,
      m = key$m,
      d = key$d,
      auc = mean_or_na(df_dataset$auc[idx]),
      f1 = mean_or_na(df_dataset$f1[idx]),
      logloss = mean_or_na(df_dataset$logloss[idx]),
      brier = mean_or_na(df_dataset$brier[idx]),
      ece = mean_or_na(df_dataset$ece[idx]),
      global_active_pairs = mean_or_na(df_dataset$global_active_pairs[idx]),
      unique_genes = mean_or_na(df_dataset$unique_genes[idx]),
      local_active_mean = mean_or_na(df_dataset$local_active_mean[idx]),
      hub_top_share = mean_or_na(df_dataset$hub_top_share[idx]),
      n_datasets = sum(idx),
      stringsAsFactors = FALSE
    )
  }))
}

write_partial_outputs <- function(results_df, out_root) {
  if (is.null(results_df) || nrow(results_df) == 0) return(invisible(NULL))

  split_path <- file.path(out_root, "split_level_metrics.csv")
  safe_write_csv(results_df, split_path)

  dataset_df <- aggregate_dataset_level(results_df)
  safe_write_csv(dataset_df, file.path(out_root, "dataset_level_metrics.csv"))

  world_df <- aggregate_world_level(dataset_df)
  safe_write_csv(world_df, file.path(out_root, "world_level_metrics.csv"))
}

# ============================================================
# SAVE RUN CONFIG
# ============================================================
write_json(
  list(
    comp_dir = comp_dir,
    out_root = out_root,
    K_values = K_values,
    m_values = m_values,
    d_value = d_value,
    n_folds = n_folds,
    n_repeats = n_repeats,
    seed = seed,
    pool_top_n = pool_top_n,
    pool_buffer_mult = pool_buffer_mult
  ),
  file.path(out_root, "used_config.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

log_msg("Start run")
log_msg("comp_dir = ", comp_dir)
log_msg("out_root = ", out_root)
log_msg("K_values = ", paste(K_values, collapse = ", "))
log_msg("m_values = ", paste(m_values, collapse = ", "))
log_msg("d = ", d_value)
log_msg("n_folds = ", n_folds, ", n_repeats = ", n_repeats)

# ============================================================
# MAIN
# ============================================================
files <- list.files(comp_dir, pattern = "\\.arff$", full.names = TRUE)
files <- sort(files)

if (length(files) == 0) {
  stop("Nie znaleziono plików .arff w: ", comp_dir)
}

results <- list()
ri <- 1L
fit_counter <- 0L

for (f in files) {
  dataset_id <- sub("\\.arff$", "", basename(f))

  log_msg("====================================================")
  log_msg("DATASET: ", dataset_id)

  xy <- read_xy_file(f)
  X <- xy$X
  y <- xy$y

  log_msg("Samples: ", nrow(X), " | Features: ", ncol(X))

  folds <- make_repeated_stratified_folds(
    y = y,
    n_folds = n_folds,
    n_repeats = n_repeats,
    seed = seed
  )

  log_msg("Number of splits: ", length(folds))

  for (sp in folds) {
    log_msg("Split: ", sp$split_id)

    X_train <- X[sp$train_idx, , drop = FALSE]
    y_train <- y[sp$train_idx]
    X_test  <- X[sp$test_idx, , drop = FALSE]
    y_test  <- y[sp$test_idx]

    log_msg("Building TSP pool top", pool_top_n, " for split ", sp$split_id)

    pairs_pool <- tryCatch({
      compute_top_pairs_switchbox(
        X_train = X_train,
        y_train = y_train,
        top_n = pool_top_n,
        buffer_mult = pool_buffer_mult,
        out_dir = NULL,
        save_all = FALSE
      )
    }, error = function(e) {
      log_msg("ERROR while building pool for ", dataset_id, " / ", sp$split_id, ": ", e$message)
      NULL
    })

    if (is.null(pairs_pool) || nrow(pairs_pool) == 0) {
      log_msg("Skipping split due to empty/failed pool: ", sp$split_id)
      next
    }

    for (K in K_values) {
      for (m in m_values) {
        fit_counter <- fit_counter + 1L
        log_msg("Fit #", fit_counter, " | dataset=", dataset_id, " | split=", sp$split_id, " | K=", K, " | m=", m, " | d=", d_value)

        tmp_dir <- file.path(out_root, "tmp", paste0(dataset_id, "_", sp$split_id, "_K", K, "_m", m))
        ensure_dir(tmp_dir)

        res <- tryCatch({
          fit_rlr_train_test(
            X_train = X_train,
            y_train = y_train,
            X_test = X_test,
            y_test = y_test,
            pairs_pool = pairs_pool,
            K = K,
            m = m,
            d = d_value,
            out_dir = tmp_dir,
            alpha = 0.5,
            lambda_choice = "lambda.min",
            seed = seed,
            extra_metrics = list(
              dataset_id = dataset_id,
              split_id = sp$split_id
            )
          )
        }, error = function(e) {
          log_msg("ERROR in fit: ", e$message)
          NULL
        })

        if (!is.null(res)) {
          results[[ri]] <- data.frame(
            dataset_id = dataset_id,
            split_id = sp$split_id,
            K = K,
            m = m,
            d = d_value,
            auc = as.numeric(res$auc),
            f1 = as.numeric(res$f1),
            logloss = as.numeric(res$logloss),
            brier = as.numeric(res$brier),
            ece = as.numeric(res$ece),
            global_active_pairs = as.numeric(res$global_active_pairs),
            unique_genes = as.numeric(res$unique_genes),
            local_active_mean = as.numeric(res$local_active_mean),
            hub_top_share = ifelse(is.null(res$hub_top_share), NA_real_, as.numeric(res$hub_top_share)),
            stringsAsFactors = FALSE
          )
          ri <- ri + 1L
        }

        # ============================================================
        # KLUCZOWA POPRAWKA: Sprzątanie po modelu, żeby nie zapchać dysku
        # ============================================================
        unlink(tmp_dir, recursive = TRUE)

        # partial save co kilka fitów
        if ((fit_counter %% 10) == 0 && length(results) > 0) {
          log_msg("Writing partial outputs after fit #", fit_counter)
          results_df <- do.call(rbind, results)
          write_partial_outputs(results_df, out_root)
        }
      }
    }
  }
}

if (length(results) == 0) {
  stop("Brak poprawnych wyników.")
}

results_df <- do.call(rbind, results)

write_partial_outputs(results_df, out_root)

log_msg("DONE")
log_msg("Saved:")
log_msg(" - ", file.path(out_root, "split_level_metrics.csv"))
log_msg(" - ", file.path(out_root, "dataset_level_metrics.csv"))
log_msg(" - ", file.path(out_root, "world_level_metrics.csv"))
log_msg("Log file: ", log_path)