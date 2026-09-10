suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Użycie: Rscript 13_make_dataset_delta_heatmap.R <config.json>")
}

cfg <- read_json(args[1], simplifyVector = TRUE)
out_root <- cfg$out_root
if (is.null(out_root)) stop("Brak out_root w config.json")

summ_dir <- file.path(out_root, "summaries")
fig_dir  <- file.path(out_root, "paper_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dataset_metrics_path <- file.path(summ_dir, "dataset_level_metrics.csv")
if (!file.exists(dataset_metrics_path)) stop("Brak dataset_level_metrics.csv")

df <- read.csv(dataset_metrics_path, stringsAsFactors = FALSE, check.names = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a

default_comparisons <- list(
  list(name = "RLR10_vs_kTSP10", rlr = "RLR_K10_m10_d2_pool1000", baseline = "kTSP_10"),
  list(name = "RLR20_vs_kTSP20", rlr = "RLR_K20_m10_d2_pool1000", baseline = "kTSP_20"),
  list(name = "RLR50_vs_kTSP50", rlr = "RLR_K50_m10_d2_pool1000", baseline = "kTSP_50")
)

default_metrics <- c("auc", "f1", "logloss", "brier", "ece")

comparisons <- cfg$heatmap_report$comparisons %||% default_comparisons
metrics_to_plot <- cfg$heatmap_report$metrics %||% default_metrics
worlds_to_plot <- cfg$heatmap_report$worlds %||% unique(df$world_name)

lower_is_better <- c("logloss", "brier", "ece")
higher_is_better <- c("auc", "f1")

safe_num <- function(x) suppressWarnings(as.numeric(x))

make_delta_df <- function(df, rlr_model, base_model, metric, cmp_name) {
  xa <- df[df$model == rlr_model, c("world_name", "dataset_id", metric), drop = FALSE]
  xb <- df[df$model == base_model, c("world_name", "dataset_id", metric), drop = FALSE]

  if (nrow(xa) == 0 || nrow(xb) == 0) return(NULL)

  colnames(xa)[3] <- "value_rlr"
  colnames(xb)[3] <- "value_base"

  merged <- merge(xa, xb, by = c("world_name", "dataset_id"), all = FALSE)
  if (nrow(merged) == 0) return(NULL)

  merged$value_rlr <- safe_num(merged$value_rlr)
  merged$value_base <- safe_num(merged$value_base)

  if (metric %in% lower_is_better) {
    merged$delta_oriented <- merged$value_base - merged$value_rlr
  } else if (metric %in% higher_is_better) {
    merged$delta_oriented <- merged$value_rlr - merged$value_base
  } else {
    merged$delta_oriented <- merged$value_rlr - merged$value_base
  }

  merged$comparison <- cmp_name
  merged$metric <- metric
  merged
}

save_plot_both <- function(p, stem, width = 12, height = 7) {
  ggsave(filename = file.path(fig_dir, paste0(stem, ".pdf")), plot = p, width = width, height = height)
  ggsave(filename = file.path(fig_dir, paste0(stem, ".png")), plot = p, width = width, height = height, dpi = 300)
}

all_delta <- list()

for (cmp in comparisons) {
  for (metric in metrics_to_plot) {
    if (!(metric %in% colnames(df))) next
    d <- make_delta_df(df, cmp$rlr, cmp$baseline, metric, cmp$name)
    if (!is.null(d)) {
      d <- d[d$world_name %in% worlds_to_plot, , drop = FALSE]
      if (nrow(d) > 0) all_delta[[length(all_delta) + 1]] <- d
    }
  }
}

if (length(all_delta) == 0) stop("Brak danych do heatmapy.")

delta_df <- do.call(rbind, all_delta)
write.csv(delta_df, file.path(fig_dir, "dataset_delta_heatmap_data.csv"), row.names = FALSE)

# kolejność datasetów: średnia poprawa malejąco
dataset_order <- aggregate(delta_oriented ~ dataset_id, data = delta_df, FUN = mean, na.rm = TRUE)
dataset_order <- dataset_order[order(dataset_order$delta_oriented, decreasing = TRUE), , drop = FALSE]
delta_df$dataset_id <- factor(delta_df$dataset_id, levels = dataset_order$dataset_id)

delta_df$metric <- factor(delta_df$metric, levels = metrics_to_plot)
delta_df$comparison <- factor(delta_df$comparison, levels = vapply(comparisons, function(x) x$name, character(1)))

p <- ggplot(delta_df, aes(x = metric, y = dataset_id, fill = delta_oriented)) +
  geom_tile() +
  facet_grid(world_name ~ comparison, scales = "free_y", space = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    x = NULL,
    y = NULL,
    fill = "RLR improvement",
    title = "Dataset-level deltas: positive values favor RLR"
  )

save_plot_both(p, "fig_dataset_delta_heatmap", width = 13, height = 8)

cat("OK: zapisano heatmapę do:\n")
cat(fig_dir, "\n")