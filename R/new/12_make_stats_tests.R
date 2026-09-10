suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Użycie: Rscript 12_make_stats_tests.R <config.json>")
}

cfg <- read_json(args[1], simplifyVector = TRUE)
out_root <- cfg$out_root
if (is.null(out_root)) stop("Brak out_root w config.json")

summ_dir <- file.path(out_root, "summaries")
stats_dir <- file.path(out_root, "paper_stats")
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

dataset_metrics_path <- file.path(summ_dir, "dataset_level_metrics.csv")
if (!file.exists(dataset_metrics_path)) stop("Brak dataset_level_metrics.csv")

df <- read.csv(dataset_metrics_path, stringsAsFactors = FALSE, check.names = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a

# ------------------------------------------------------------
# Ustawienia domyślne
# ------------------------------------------------------------
default_metrics <- c("logloss", "brier", "ece", "auc", "f1")

normalize_comparison_spec <- function(x) {
  if (is.null(x)) return(list())

  out <- NULL

  if (is.data.frame(x)) {
    out <- lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
  } else if (is.list(x)) {
    if (!is.null(names(x)) && all(c("comparison","model_a","model_b") %in% names(x))) {
      out <- list(as.list(x))
    } else {
      out <- lapply(x, as.list)
    }
  } else {
    stop("Unsupported comparison spec type")
  }

  out <- lapply(out, function(z) {
    if (is.null(z$name) && !is.null(z$comparison)) z$name <- z$comparison
    if (is.null(z$comparison) && !is.null(z$name)) z$comparison <- z$name
    z
  })

  out
}

default_comparisons <- list(
  list(comparison = "RLR10_vs_kTSP10", model_a = "RLR_K10_m10_d2_pool1000", model_b = "kTSP_10"),
  list(comparison = "RLR20_vs_kTSP20", model_a = "RLR_K20_m10_d2_pool1000", model_b = "kTSP_20"),
  list(comparison = "RLR50_vs_kTSP50", model_a = "RLR_K50_m10_d2_pool1000", model_b = "kTSP_50"),
  list(comparison = "RLR50_d2_vs_d1", model_a = "RLR_K50_m10_d2_pool1000", model_b = "RLR_K50_m10_d1_pool1000")
)

metrics_to_test <- cfg$stats_report$metrics %||% default_metrics
comparisons <- normalize_comparison_spec(cfg$stats_report$comparisons %||% default_comparisons)
paper_focus_comparisons <- normalize_comparison_spec(cfg$stats_report$paper_focus_comparisons %||% comparisons)

# na których metrykach mniejsze = lepsze
lower_is_better <- c("logloss", "brier", "ece")
higher_is_better <- c("auc", "f1")

# ------------------------------------------------------------
# Helpery
# ------------------------------------------------------------
safe_numeric <- function(x) suppressWarnings(as.numeric(x))

make_paired_table <- function(df, model_a, model_b, metric, world_name = NULL) {
  x <- df
  if (!is.null(world_name)) {
    x <- x[x$world_name == world_name, , drop = FALSE]
  }

  xa <- x[x$model == model_a, c("world_name", "dataset_id", metric), drop = FALSE]
  xb <- x[x$model == model_b, c("world_name", "dataset_id", metric), drop = FALSE]

  colnames(xa)[3] <- "value_a"
  colnames(xb)[3] <- "value_b"

  merged <- merge(xa, xb, by = c("world_name", "dataset_id"), all = FALSE)
  merged$value_a <- safe_numeric(merged$value_a)
  merged$value_b <- safe_numeric(merged$value_b)
  merged <- merged[complete.cases(merged$value_a, merged$value_b), , drop = FALSE]

  merged
}

run_one_test <- function(paired_df, metric, model_a, model_b, comparison_name, world_name = NULL) {
  n <- nrow(paired_df)

  if (n < 2) {
    return(data.frame(
      world_name = ifelse(is.null(world_name), "ALL", world_name),
      comparison = comparison_name,
      model_a = model_a,
      model_b = model_b,
      metric = metric,
      n_datasets = n,
      mean_a = NA,
      mean_b = NA,
      median_diff = NA,
      mean_diff = NA,
      p_value = NA,
      stringsAsFactors = FALSE
    ))
  }

  diff_ab <- paired_df$value_a - paired_df$value_b

  wt <- tryCatch(
    wilcox.test(paired_df$value_a, paired_df$value_b, paired = TRUE, exact = FALSE),
    error = function(e) NULL
  )

  data.frame(
    world_name = ifelse(is.null(world_name), "ALL", world_name),
    comparison = comparison_name,
    model_a = model_a,
    model_b = model_b,
    metric = metric,
    n_datasets = n,
    mean_a = mean(paired_df$value_a, na.rm = TRUE),
    mean_b = mean(paired_df$value_b, na.rm = TRUE),
    median_diff = median(diff_ab, na.rm = TRUE),
    mean_diff = mean(diff_ab, na.rm = TRUE),
    p_value = ifelse(is.null(wt), NA, wt$p.value),
    stringsAsFactors = FALSE
  )
}

add_direction_interpretation <- function(res_df) {
  res_df$direction <- NA_character_

  for (i in seq_len(nrow(res_df))) {
    m <- res_df$metric[i]
    d <- res_df$mean_diff[i]

    if (is.na(d)) {
      res_df$direction[i] <- NA_character_
    } else if (m %in% lower_is_better) {
      res_df$direction[i] <- ifelse(d < 0, "model_a_better", ifelse(d > 0, "model_b_better", "tie"))
    } else if (m %in% higher_is_better) {
      res_df$direction[i] <- ifelse(d > 0, "model_a_better", ifelse(d < 0, "model_b_better", "tie"))
    } else {
      res_df$direction[i] <- NA_character_
    }
  }

  res_df
}

# ------------------------------------------------------------
# Analiza globalna po worldach
# ------------------------------------------------------------
worlds <- unique(df$world_name)
all_results <- list()

for (w in worlds) {
  for (cmp in comparisons) {
    for (metric in metrics_to_test) {
      if (!(metric %in% colnames(df))) next

      paired_df <- make_paired_table(
        df = df,
        model_a = cmp$model_a,
        model_b = cmp$model_b,
        metric = metric,
        world_name = w
      )

      one <- run_one_test(
        paired_df = paired_df,
        metric = metric,
        model_a = cmp$model_a,
        model_b = cmp$model_b,
        comparison_name = cmp$name,
        world_name = w
      )

      all_results[[length(all_results) + 1]] <- one

      # zapis tabeli parowanej dla tej jednej analizy
      paired_out <- paired_df
      paired_out$diff_a_minus_b <- paired_out$value_a - paired_out$value_b

      out_name <- paste0("paired_", w, "_", cmp$name, "_", metric, ".csv")
      write.csv(paired_out, file.path(stats_dir, out_name), row.names = FALSE)
    }
  }
}

results_df <- do.call(rbind, all_results)
results_df <- add_direction_interpretation(results_df)

# Holm correction per world
results_df$p_holm_within_world <- NA_real_
for (w in unique(results_df$world_name)) {
  idx <- which(results_df$world_name == w & !is.na(results_df$p_value))
  if (length(idx) > 0) {
    results_df$p_holm_within_world[idx] <- p.adjust(results_df$p_value[idx], method = "holm")
  }
}

write.csv(results_df, file.path(stats_dir, "stats_tests_all.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Tabela skrócona do paperu
# ------------------------------------------------------------
paper_focus_names <- c(
  "RLR10_vs_kTSP10",
  "RLR20_vs_kTSP20",
  "RLR50_vs_kTSP50"
)

paper_focus <- results_df[
  results_df$comparison %in% paper_focus_names &
    results_df$metric %in% metrics_to_test,
  ,
  drop = FALSE
]

write.csv(paper_focus, file.path(stats_dir, "stats_tests_paper_focus.csv"), row.names = FALSE)

cat("OK: zapisano statystyki do:\n")
cat(stats_dir, "\n")