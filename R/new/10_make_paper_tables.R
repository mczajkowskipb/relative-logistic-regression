suppressPackageStartupMessages(library(jsonlite))

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Użycie: Rscript 10_make_paper_tables.R <config.json>")
}

cfg <- read_json(args[1], simplifyVector = TRUE)
out_root <- cfg$out_root
if (is.null(out_root)) stop("Brak out_root w config.json")

summ_dir <- file.path(out_root, "summaries")
tab_dir <- file.path(out_root, "paper_tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

dataset_metrics_path <- file.path(summ_dir, "dataset_level_metrics.csv")
world_metrics_path   <- file.path(summ_dir, "world_level_metrics.csv")
local_path           <- file.path(summ_dir, "dataset_level_local_activity.csv")

if (!file.exists(dataset_metrics_path)) stop("Brak dataset_level_metrics.csv")
if (!file.exists(world_metrics_path)) stop("Brak world_level_metrics.csv")

dataset_metrics <- read.csv(dataset_metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
world_metrics   <- read.csv(world_metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
local_metrics   <- if (file.exists(local_path)) {
  read.csv(local_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  NULL
}

# ------------------------------------------------------------
# Ustawienia modeli do raportu
# ------------------------------------------------------------
default_main_models <- c(
  "value_LR",
  "kTSP_10",
  "kTSP_50",
  "RLR_K10_m10_d1",
  "RLR_K50_m10_d1"
)

default_compact_models <- c(
  "kTSP_10",
  "RLR_K10_m10_d1"
)

default_large_models <- c(
  "kTSP_50",
  "RLR_K50_m10_d1"
)

paper_report <- cfg$paper_report %||% list()

main_models    <- paper_report$main_models %||% default_main_models
compact_models <- paper_report$compact_models %||% default_compact_models
large_models   <- paper_report$large_models %||% default_large_models

# ------------------------------------------------------------
# Tabela 1: world-level performance
# ------------------------------------------------------------
perf_cols <- intersect(
  c("world_name", "model", "auc", "f1", "logloss", "brier", "ece"),
  colnames(world_metrics)
)

table1 <- world_metrics[world_metrics$model %in% main_models, perf_cols, drop = FALSE]
table1 <- table1[order(table1$world_name, match(table1$model, main_models)), , drop = FALSE]

write.csv(
  table1,
  file.path(tab_dir, "table1_world_performance.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Tabela 2: world-level structure
# ------------------------------------------------------------
struct_cols <- intersect(
  c(
    "world_name", "model",
    "global_active_pairs", "unique_genes",
    "local_active_mean", "local_active_median", "local_active_iqr",
    "hub_top_share"
  ),
  colnames(world_metrics)
)

table2 <- world_metrics[world_metrics$model %in% main_models, struct_cols, drop = FALSE]
table2 <- table2[order(table2$world_name, match(table2$model, main_models)), , drop = FALSE]

write.csv(
  table2,
  file.path(tab_dir, "table2_world_structure.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Tabela 3: dataset-level comparison for selected models
# ------------------------------------------------------------
table3_cols <- intersect(
  c("world_name", "dataset_id", "model", "auc", "f1", "logloss", "brier", "ece"),
  colnames(dataset_metrics)
)

table3_models <- unique(c(compact_models, large_models))

table3 <- dataset_metrics[
  dataset_metrics$model %in% table3_models,
  table3_cols,
  drop = FALSE
]
table3 <- table3[order(table3$world_name, table3$dataset_id, match(table3$model, table3_models)), , drop = FALSE]

write.csv(
  table3,
  file.path(tab_dir, "table3_dataset_level_selected_models.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table A: full world-level metrics
# ------------------------------------------------------------
write.csv(
  world_metrics,
  file.path(tab_dir, "supp_table_world_level_metrics_full.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table B: full dataset-level metrics
# ------------------------------------------------------------
write.csv(
  dataset_metrics,
  file.path(tab_dir, "supp_table_dataset_level_metrics_full.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table C: local activity summaries
# ------------------------------------------------------------
if (!is.null(local_metrics)) {
  write.csv(
    local_metrics,
    file.path(tab_dir, "supp_table_dataset_level_local_activity.csv"),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# Supplementary Table D: ablation by m/d/pool if present
# ------------------------------------------------------------
ablation_cols <- intersect(
  c(
    "world_name", "model",
    "K_cfg", "m_cfg", "d_cfg", "pool_cfg",
    "auc", "f1", "logloss", "brier", "ece",
    "global_active_pairs", "unique_genes",
    "local_active_mean", "hub_top_share"
  ),
  colnames(world_metrics)
)

supp_ablation <- world_metrics[grepl("^RLR_", world_metrics$model), ablation_cols, drop = FALSE]

ord_cols <- intersect(c("world_name", "K_cfg", "m_cfg", "d_cfg", "pool_cfg", "model"), colnames(supp_ablation))
if (length(ord_cols) > 0 && nrow(supp_ablation) > 0) {
  supp_ablation <- supp_ablation[do.call(order, supp_ablation[ord_cols]), , drop = FALSE]
}

write.csv(
  supp_ablation,
  file.path(tab_dir, "supp_table_rlr_ablation.csv"),
  row.names = FALSE
)

cat("OK: zapisano tabele do:\n")
cat(tab_dir, "\n")