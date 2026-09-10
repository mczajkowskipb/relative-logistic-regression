get_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  getwd()
}
source(file.path(get_script_dir(), "07_world_protocol_utils.R"))
suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Użycie: Rscript 08_run_worlds_config.R <config.json>")
}

cfg_path <- args[1]
cfg <- jsonlite::read_json(cfg_path, simplifyVector = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

normalize_value <- function(v) {
  if (is.data.frame(v)) {
    return(as.list(v))
  }

  if (is.list(v)) {
    if (length(v) == 1) {
      return(normalize_value(v[[1]]))
    }

    flat <- tryCatch(unlist(v, use.names = FALSE), error = function(e) NULL)
    if (!is.null(flat)) return(flat)

    return(v)
  }

  v
}

normalize_world_row <- function(x) {
  x <- as.list(x)
  x <- lapply(x, normalize_value)
  x
}

normalize_worlds <- function(worlds_raw) {
  if (is.null(worlds_raw)) return(list())

  if (is.data.frame(worlds_raw)) {
    return(lapply(seq_len(nrow(worlds_raw)), function(i) {
      normalize_world_row(worlds_raw[i, , drop = FALSE])
    }))
  }

  if (is.list(worlds_raw) &&
      !is.null(names(worlds_raw)) &&
      all(c("name", "protocol") %in% names(worlds_raw))) {
    return(list(normalize_world_row(worlds_raw)))
  }

  if (is.list(worlds_raw)) {
    return(lapply(worlds_raw, normalize_world_row))
  }

  stop("cfg$worlds ma nieobslugiwany typ")
}

worlds_cfg <- normalize_worlds(cfg$worlds)

if (is.null(cfg$out_root)) stop("Brak cfg$out_root")
if (length(worlds_cfg) == 0) stop("Brak cfg$worlds")
if (is.null(cfg$models)) stop("Brak cfg$models")

out_root <- cfg$out_root
ensure_dir(out_root)
ensure_dir(file.path(out_root, "config"))

jsonlite::write_json(
  cfg,
  file.path(out_root, "config", "used_config.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

# ------------------------------------------------------------
# Global model config
# ------------------------------------------------------------
models_cfg <- cfg$models

run_baselines <- isTRUE(models_cfg$run_baselines %||% TRUE)
run_rlr <- isTRUE(models_cfg$run_rlr %||% TRUE)
skip_existing <- isTRUE(models_cfg$skip_existing %||% TRUE)
save_all_pairs <- isTRUE(models_cfg$save_all_pairs %||% FALSE)
include_pool_equal_K <- isTRUE(models_cfg$include_pool_equal_K %||% FALSE)

seed <- as.integer(models_cfg$seed %||% 42)
top_pool_values <- as.integer(models_cfg$top_pool_values %||% 1000)
top_pool_values <- sort(unique(top_pool_values))

baselines_k <- as.integer(models_cfg$baselines_k %||% c(1, 10, 20, 50))

lambda_choice <- models_cfg$lambda_choice %||% "lambda.min"
alpha <- as.numeric(models_cfg$alpha %||% 0.5)
pre_tsp_top_features <- models_cfg$pre_tsp_top_features %||% NULL
if (!is.null(pre_tsp_top_features)) {
  pre_tsp_top_features <- as.integer(pre_tsp_top_features)
}

# rlr grid
if (is.null(models_cfg$rlr_grid)) {
  rlr_grid <- data.frame(K = integer(0), m = integer(0), d = integer(0))
} else {
  rlr_grid <- expand.grid(
    K = as.integer(models_cfg$rlr_grid$K),
    m = as.integer(models_cfg$rlr_grid$m),
    d = as.integer(models_cfg$rlr_grid$d),
    stringsAsFactors = FALSE
  )
}

bridge_grid <- data.frame(K = integer(0), m = integer(0), d = integer(0))
if (!is.null(models_cfg$bridge_grid) && length(models_cfg$bridge_grid) > 0) {
  if (is.data.frame(models_cfg$bridge_grid)) {
    bridge_grid <- data.frame(
      K = as.integer(models_cfg$bridge_grid$K),
      m = as.integer(models_cfg$bridge_grid$m),
      d = as.integer(models_cfg$bridge_grid$d),
      stringsAsFactors = FALSE
    )
  } else if (is.list(models_cfg$bridge_grid)) {
    bridge_grid <- do.call(rbind, lapply(models_cfg$bridge_grid, function(x) {
      data.frame(
        K = as.integer(x$K),
        m = as.integer(x$m),
        d = as.integer(x$d),
        stringsAsFactors = FALSE
      )
    }))
  }
}

full_rlr_grid <- rbind(rlr_grid, bridge_grid)
if (nrow(full_rlr_grid) > 0) {
  full_rlr_grid <- unique(full_rlr_grid)
}

write.csv(full_rlr_grid, file.path(out_root, "config", "grid_rlr.csv"), row.names = FALSE)
write.csv(data.frame(k = baselines_k), file.path(out_root, "config", "grid_baselines.csv"), row.names = FALSE)

manifest_rows <- list()
manifest_i <- 1L

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
register_manifest <- function(world_name, dataset_id, split_id, family, model, path) {
  manifest_rows[[manifest_i]] <<- data.frame(
    world_name = world_name,
    dataset_id = dataset_id,
    split_id = split_id,
    family = family,
    model = model,
    path = path,
    stringsAsFactors = FALSE
  )
  manifest_i <<- manifest_i + 1L
}

get_dataset_files <- function(world) {
  protocol <- world$protocol

  if (protocol == "compcancer_cv") {
    comp_dir <- world$comp_dir %||% stop("compcancer_cv wymaga world$comp_dir")
    all_files <- list.files(comp_dir, pattern = "\\.arff$", full.names = TRUE)

    inc <- world$datasets_include %||% NULL
    exc <- world$datasets_exclude %||% NULL
    max_datasets <- world$max_datasets %||% NULL

    if (!is.null(inc)) {
      inc_full <- paste0(inc, ".arff")
      all_files <- all_files[basename(all_files) %in% inc_full]
    }
    if (!is.null(exc)) {
      exc_full <- paste0(exc, ".arff")
      all_files <- all_files[!basename(all_files) %in% exc_full]
    }

    all_files <- sort(all_files)

    if (!is.null(max_datasets)) {
      all_files <- head(all_files, as.integer(max_datasets))
    }
    return(all_files)
  }

  if (protocol == "fixed_split") {
    train_file <- world$train_file %||% stop("fixed_split wymaga train_file")
    test_file <- world$test_file %||% stop("fixed_split wymaga test_file")
    return(list(list(train_file = train_file, test_file = test_file)))
  }

  if (protocol == "loco_files") {
    cohort_files <- world$cohort_files %||% NULL
    cohort_dir <- world$cohort_dir %||% NULL
    cohort_pattern <- world$cohort_pattern %||% "\\.arff$"

    if (is.null(cohort_files)) {
      if (is.null(cohort_dir)) stop("loco_files wymaga cohort_files albo cohort_dir")
      cohort_files <- list.files(cohort_dir, pattern = cohort_pattern, full.names = TRUE)
    }

    cohort_files <- sort(cohort_files)
    return(cohort_files)
  }

  stop(paste("Nieznany protocol:", protocol))
}

make_world_splits <- function(world, dataset_obj = NULL) {
  protocol <- world$protocol

  if (protocol == "compcancer_cv") {
    n_folds <- as.integer(world$n_folds %||% 5)
    n_repeats <- as.integer(world$n_repeats %||% 10)
    max_splits <- world$max_splits_per_dataset %||% NULL

    folds <- make_repeated_stratified_folds(
      y = dataset_obj$y,
      n_folds = n_folds,
      n_repeats = n_repeats,
      seed = seed
    )

    if (!is.null(max_splits)) {
      folds <- head(folds, as.integer(max_splits))
    }
    return(folds)
  }

  if (protocol == "fixed_split") {
    split_name <- world$split_name %||% "fixed"
    return(list(make_fixed_split_spec(split_name)))
  }

  if (protocol == "loco_files") {
    cohort_files <- get_dataset_files(world)
    return(make_loco_specs_from_files(cohort_files))
  }

  stop(paste("Nieznany protocol:", protocol))
}

prepare_train_test_for_split <- function(world, split_spec, dataset_obj = NULL) {
  protocol <- world$protocol

  if (protocol == "compcancer_cv") {
    X_train <- dataset_obj$X[split_spec$train_idx, , drop = FALSE]
    y_train <- dataset_obj$y[split_spec$train_idx]
    X_test <- dataset_obj$X[split_spec$test_idx, , drop = FALSE]
    y_test <- dataset_obj$y[split_spec$test_idx]

    return(list(
      X_train = X_train,
      y_train = y_train,
      X_test = X_test,
      y_test = y_test,
      train_origin = "cv_train",
      test_origin = "cv_test"
    ))
  }

  if (protocol == "fixed_split") {
    train_xy <- read_xy_file(
      path = world$train_file,
      class_col = world$class_col %||% NULL,
      sample_id_col = world$sample_id_col %||% NULL
    )
    test_xy <- read_xy_file(
      path = world$test_file,
      class_col = world$class_col %||% NULL,
      sample_id_col = world$sample_id_col %||% NULL
    )

    aligned <- align_xy_train_test(train_xy$X, test_xy$X)

    return(list(
      X_train = aligned$X_train,
      y_train = train_xy$y,
      X_test = aligned$X_test,
      y_test = test_xy$y,
      train_origin = basename(world$train_file),
      test_origin = basename(world$test_file)
    ))
  }

  if (protocol == "loco_files") {
    train_xy_list <- lapply(split_spec$train_files, function(f) {
      read_xy_file(
        path = f,
        class_col = world$class_col %||% NULL,
        sample_id_col = world$sample_id_col %||% NULL
      )
    })

    test_xy <- read_xy_file(
      path = split_spec$holdout_file,
      class_col = world$class_col %||% NULL,
      sample_id_col = world$sample_id_col %||% NULL
    )

    train_xy <- combine_xy_list(train_xy_list)
    aligned <- align_xy_train_test(train_xy$X, test_xy$X)

    return(list(
      X_train = aligned$X_train,
      y_train = train_xy$y,
      X_test = aligned$X_test,
      y_test = test_xy$y,
      train_origin = paste(basename(split_spec$train_files), collapse = ";"),
      test_origin = basename(split_spec$holdout_file)
    ))
  }

  stop(paste("Nieznany protocol:", protocol))
}

run_one_split <- function(world_name, dataset_id, split_spec, train_test, split_dir) {
  X_train <- train_test$X_train
  y_train <- train_test$y_train
  X_test  <- train_test$X_test
  y_test  <- train_test$y_test
  X_valid <- NULL

  pf <- maybe_prefilter_split(
    X_train = X_train,
    X_test = X_test,
    X_valid = X_valid,
    pre_tsp_top_features = pre_tsp_top_features
  )

  X_train <- pf$X_train
  X_test  <- pf$X_test
  if (!is.null(pf$X_valid)) X_valid <- pf$X_valid

  if (isTRUE(pf$prefilter_applied)) {
    cat("  -> prefilter MAD top", pf$kept_n, "\n")
  }

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
    world_name = world_name,
    dataset_id = dataset_id,
    split_id = split_spec$split_id,
    train_origin = train_test$train_origin,
    test_origin = train_test$test_origin,
    prefilter_applied = isTRUE(pf$prefilter_applied),
    prefilter_kept_n = pf$kept_n
  )
  write_json(split_meta, file.path(split_dir, "split_meta.json"), auto_unbox = TRUE, pretty = TRUE)


  # ---------- pool cache ----------
  split_pool_values <- top_pool_values
  if (include_pool_equal_K && nrow(full_rlr_grid) > 0) {
    split_pool_values <- sort(unique(c(split_pool_values, full_rlr_grid$K)))
  }

  pairs_pool_cache <- list()

  for (pool_n in split_pool_values) {
    pool_key <- as.character(pool_n)
    pool_dir <- file.path(split_dir, "pairs_train", paste0("top", pool_n))
    ensure_dir(pool_dir)

    top_pairs_file <- file.path(pool_dir, paste0("top_", pool_n, "_pairs.csv"))

    if (skip_existing && file.exists(top_pairs_file)) {
      cat("  -> pairs_train top", pool_n, " (skip-existing)\n", sep = "")
      pairs_pool_cache[[pool_key]] <- read.csv(top_pairs_file, stringsAsFactors = FALSE)
    } else {
      cat("  -> pairs_train top", pool_n, "\n", sep = "")
      pairs_pool_cache[[pool_key]] <- compute_top_pairs_switchbox(
        X_train = X_train,
        y_train = y_train,
        top_n = pool_n,
        buffer_mult = 5,
        out_dir = pool_dir,
        save_all = save_all_pairs && (pool_n == max(split_pool_values))
      )
    }
  }

  # ---------- baselines ----------
  if (run_baselines) {
    base_dir <- file.path(split_dir, "baselines")
    ensure_dir(base_dir)

    for (k in baselines_k) {
      model_name <- ifelse(k == 1, "TSP_1", paste0("kTSP_", k))
      model_dir <- file.path(base_dir, model_name)
      metrics_file <- file.path(model_dir, "metrics_test.json")

      if (skip_existing && file.exists(metrics_file)) {
        cat("  -> baseline:", model_name, " (skip-existing)\n")
      } else {
        cat("  -> baseline:", model_name, "\n")
        ensure_dir(model_dir)

        extra_metrics <- list(
          world_name = world_name,
          dataset_id = dataset_id,
          split_id = split_spec$split_id
        )

        tryCatch({
          if (k == 1) {
            fit_tsp1_manual_train_test(
              X_train = X_train, y_train = y_train,
              X_test = X_test, y_test = y_test,
              pairs_pool = pairs_pool_cache[[as.character(max(split_pool_values))]],
              out_dir = model_dir,
              extra_metrics = extra_metrics
            )
          } else {
            fit_ktsp_train_test(
              X_train = X_train, y_train = y_train,
              X_test = X_test, y_test = y_test,
              k = k,
              out_dir = model_dir,
              seed = seed,
              extra_metrics = extra_metrics
            )
          }
        }, error = function(e) {
          write_error_json(file.path(model_dir, "error.json"), e)
        })
      }

      register_manifest(world_name, dataset_id, split_spec$split_id, "baseline", model_name, model_dir)
    }
  }

  # ---------- value-LR ----------
  value_lr_cfg <- cfg$value_lr
  run_value_lr <- !is.null(value_lr_cfg) && isTRUE(value_lr_cfg$enabled)

  if (run_value_lr) {
    value_dir <- file.path(split_dir, "value_lr")
    ensure_dir(value_dir)

    model_name <- if (!is.null(value_lr_cfg$model_name)) value_lr_cfg$model_name else "value_LR"
    pre_top_features <- if (!is.null(value_lr_cfg$pre_top_features)) value_lr_cfg$pre_top_features else 1000
    alpha_value <- if (!is.null(value_lr_cfg$alpha)) value_lr_cfg$alpha else 0.5
    lambda_choice <- if (!is.null(value_lr_cfg$lambda_choice)) value_lr_cfg$lambda_choice else "lambda.min"

    model_dir <- file.path(value_dir, model_name)
    metrics_file <- file.path(model_dir, "metrics_test.json")

    if (skip_existing && file.exists(metrics_file)) {
      cat("  -> value baseline:", model_name, " (skip-existing)\n")
    } else {
      cat("  -> value baseline:", model_name, "\n")
      ensure_dir(model_dir)

      extra_metrics <- list(
        world_name = world_name,
        dataset_id = dataset_id,
        split_id = split_spec$split_id
      )

      tryCatch({
        fit_value_lr_train_test(
          X_train = X_train, y_train = y_train,
          X_test = X_test, y_test = y_test,
          out_dir = model_dir,
          pre_top_features = pre_top_features,
          alpha = alpha_value,
          lambda_choice = lambda_choice,
          seed = seed,
          extra_metrics = extra_metrics
        )
      }, error = function(e) {
        write_error_json(file.path(model_dir, "error.json"), e)
      })
    }

    register_manifest(world_name, dataset_id, split_spec$split_id, "value_baseline", model_name, model_dir)
  }

  # ---------- RLR ----------
  if (run_rlr && nrow(full_rlr_grid) > 0) {
    rlr_dir <- file.path(split_dir, "rlr")
    ensure_dir(rlr_dir)

    for (i in seq_len(nrow(full_rlr_grid))) {
      g <- full_rlr_grid[i, ]

      model_pool_values <- top_pool_values
      if (include_pool_equal_K) {
        model_pool_values <- sort(unique(c(model_pool_values, g$K)))
      }

      for (pool_n in model_pool_values) {
        model_name <- paste0("RLR_K", g$K, "_m", g$m, "_d", g$d, "_pool", pool_n)
        model_dir <- file.path(rlr_dir, model_name)
        metrics_file <- file.path(model_dir, "metrics_test.json")

        if (skip_existing && file.exists(metrics_file)) {
          cat("  -> rlr:", model_name, " (skip-existing)\n")
        } else {
          cat("  -> rlr:", model_name, "\n")
          ensure_dir(model_dir)

          extra_metrics <- list(
            world_name = world_name,
            dataset_id = dataset_id,
            split_id = split_spec$split_id,
            pool = pool_n
          )

          tryCatch({
            fit_rlr_train_test(
              X_train = X_train, y_train = y_train,
              X_test = X_test, y_test = y_test,
              pairs_pool = pairs_pool_cache[[as.character(pool_n)]],
              K = g$K,
              m = g$m,
              d = g$d,
              out_dir = model_dir,
              alpha = alpha,
              lambda_choice = lambda_choice,
              seed = seed,
              extra_metrics = extra_metrics
            )
          }, error = function(e) {
            write_error_json(file.path(model_dir, "error.json"), e)
          })
        }

        register_manifest(world_name, dataset_id, split_spec$split_id, "rlr", model_name, model_dir)
      }
    }
  }
}

# ------------------------------------------------------------
# Main loop over worlds
# ------------------------------------------------------------
for (world in worlds_cfg) {
  world_name <- world$name %||% stop("Kazdy world musi miec pole 'name'")
  protocol <- world$protocol %||% stop(paste("World", world_name, "nie ma pola protocol"))

  cat("\n############################################\n")
  cat("WORLD:", world_name, " | protocol:", protocol, "\n")

  world_dir <- file.path(out_root, "worlds", world_name)
  ensure_dir(world_dir)

  if (protocol == "compcancer_cv") {
    dataset_files <- get_dataset_files(world)

    for (arff_path in dataset_files) {
      dataset_id <- sub("\\.arff$", "", basename(arff_path))
      cat("\n==============================\n")
      cat("DATASET:", dataset_id, "\n")

      dataset_obj <- read_xy_file(
        path = arff_path,
        class_col = world$class_col %||% NULL,
        sample_id_col = world$sample_id_col %||% NULL
      )

      ds_dir <- file.path(world_dir, "datasets", dataset_id)
      ensure_dir(ds_dir)

      split_specs <- make_world_splits(world, dataset_obj = dataset_obj)

      for (sp in split_specs) {
        split_dir <- file.path(ds_dir, "splits", sp$split_id)
        ensure_dir(split_dir)

        train_test <- prepare_train_test_for_split(world, sp, dataset_obj = dataset_obj)
        run_one_split(world_name, dataset_id, sp, train_test, split_dir)
      }
    }
  } else if (protocol == "fixed_split") {
    dataset_id <- world$dataset_id %||% world_name
    ds_dir <- file.path(world_dir, "datasets", dataset_id)
    ensure_dir(ds_dir)

    split_specs <- make_world_splits(world, dataset_obj = NULL)

    for (sp in split_specs) {
      split_dir <- file.path(ds_dir, "splits", sp$split_id)
      ensure_dir(split_dir)

      train_test <- prepare_train_test_for_split(world, sp, dataset_obj = NULL)
      run_one_split(world_name, dataset_id, sp, train_test, split_dir)
    }
  } else if (protocol == "loco_files") {
    dataset_id <- world$dataset_id %||% world_name
    ds_dir <- file.path(world_dir, "datasets", dataset_id)
    ensure_dir(ds_dir)

    split_specs <- make_world_splits(world, dataset_obj = NULL)

    for (sp in split_specs) {
      split_dir <- file.path(ds_dir, "splits", sp$split_id)
      ensure_dir(split_dir)

      train_test <- prepare_train_test_for_split(world, sp, dataset_obj = NULL)
      run_one_split(world_name, dataset_id, sp, train_test, split_dir)
    }
  } else {
    stop(paste("Nieobslugiwany protocol:", protocol))
  }
}

manifest <- do.call(rbind, manifest_rows)
write.csv(manifest, file.path(out_root, "config", "run_manifest.csv"), row.names = FALSE)

cat("\nGOTOWE\n")