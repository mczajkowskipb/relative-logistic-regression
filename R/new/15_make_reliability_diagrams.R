suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Użycie: Rscript 15_make_reliability_diagrams.R <config.json>")
}

cfg <- read_json(args[1], simplifyVector = TRUE)
out_root <- cfg$out_root
if (is.null(out_root)) stop("Brak out_root w config.json")

summ_dir <- file.path(out_root, "summaries")
fig_dir  <- file.path(out_root, "paper_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a

default_models <- c(
  "kTSP_10",
  "kTSP_50",
  "RLR_K10_m10_d2_pool1000",
  "RLR_K50_m10_d2_pool1000"
)

plot_models <- cfg$reliability_report$plot_models %||% default_models
plot_worlds <- cfg$reliability_report$plot_worlds %||% NULL
n_bins <- as.integer(cfg$reliability_report$n_bins %||% 10)

pred_files <- list.files(
  out_root,
  pattern = "predictions_test\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

extract_meta <- function(path) {
  world_name <- sub(".*worlds[\\\\/](.*?)[\\\\/]datasets.*", "\\1", path)
  dataset_id <- sub(".*datasets[\\\\/](.*?)[\\\\/]splits.*", "\\1", path)
  split_id <- sub(".*splits[\\\\/](.*?)[\\\\/](baselines|rlr).*", "\\1", path)
  model <- if (grepl("[\\\\/]value_lr[\\\\/]", path)) "value_LR" else sub(".*[\\\\/](TSP_1|kTSP_[0-9]+|RLR_K[0-9]+_m[0-9]+_d[0-9]+(?:_pool[0-9]+)?|ValueLR_MAD1000_EN)[\\\\/].*", "\\1", path)
  if (identical(model, path) || is.na(model) || model == "") model <- "value_LR"

  data.frame(
    world_name = world_name,
    dataset_id = dataset_id,
    split_id = split_id,
    model = model,
    stringsAsFactors = FALSE
  )
}

pred_list <- lapply(pred_files, function(f) {
  dd <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(dd)) return(NULL)

  meta <- extract_meta(f)
  if (!(meta$model %in% plot_models)) return(NULL)
  if (!is.null(plot_worlds) && !(meta$world_name %in% plot_worlds)) return(NULL)

  dd$world_name <- meta$world_name
  dd$dataset_id <- meta$dataset_id
  dd$split_id <- meta$split_id
  dd$model <- meta$model

  dd$prob <- suppressWarnings(as.numeric(dd$prob))
  dd$y_true_num <- as.numeric(as.factor(dd$y_true)) - 1

  dd
})

pred_list <- pred_list[!vapply(pred_list, is.null, logical(1))]
if (length(pred_list) == 0) stop("Brak danych predictions_test.csv dla wybranych modeli.")

pred_df <- do.call(rbind, pred_list)
write.csv(pred_df, file.path(fig_dir, "reliability_predictions_raw.csv"), row.names = FALSE)

# fixed-width bins
breaks <- seq(0, 1, length.out = n_bins + 1)

make_bin_summary <- function(df) {
  df$bin <- cut(df$prob, breaks = breaks, include.lowest = TRUE, labels = FALSE)

  keys <- unique(df[, c("world_name", "model", "bin"), drop = FALSE])

  out <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    idx <- df$world_name == key$world_name &
           df$model == key$model &
           df$bin == key$bin

    z <- df[idx, , drop = FALSE]
    if (nrow(z) == 0) return(NULL)

    data.frame(
      world_name = key$world_name,
      model = key$model,
      bin = key$bin,
      n = nrow(z),
      mean_pred = mean(z$prob, na.rm = TRUE),
      observed = mean(z$y_true_num, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  out
}

bin_df <- make_bin_summary(pred_df)
write.csv(bin_df, file.path(fig_dir, "reliability_bin_summary.csv"), row.names = FALSE)

save_plot_both <- function(p, stem, width = 8, height = 6) {
  ggsave(filename = file.path(fig_dir, paste0(stem, ".pdf")), plot = p, width = width, height = height)
  ggsave(filename = file.path(fig_dir, paste0(stem, ".png")), plot = p, width = width, height = height, dpi = 300)
}

# main reliability line plot
p1 <- ggplot(bin_df, aes(x = mean_pred, y = observed, group = model)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_line() +
  geom_point(aes(size = n)) +
  facet_wrap(~ world_name, scales = "free") +
  theme_bw(base_size = 11) +
  labs(
    x = "Mean predicted probability",
    y = "Observed event rate",
    title = "Reliability diagrams"
  )

save_plot_both(p1, "fig_reliability_diagrams", width = 10, height = 6)

# comparison panels: kTSP vs matched RLR
pair_defs <- list(
  list(name = "k10_vs_rlr10", models = c("kTSP_10", "RLR_K10_m10_d2_pool1000")),
  list(name = "k50_vs_rlr50", models = c("kTSP_50", "RLR_K50_m10_d2_pool1000"))
)

for (pd in pair_defs) {
  subdf <- bin_df[bin_df$model %in% pd$models, , drop = FALSE]
  if (nrow(subdf) == 0) next

  p <- ggplot(subdf, aes(x = mean_pred, y = observed, group = model)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    geom_line() +
    geom_point(aes(size = n)) +
    facet_wrap(~ world_name, scales = "free") +
    theme_bw(base_size = 11) +
    labs(
      x = "Mean predicted probability",
      y = "Observed event rate",
      title = paste("Reliability:", pd$name)
    )

  save_plot_both(p, paste0("supp_reliability_", pd$name), width = 10, height = 6)
}

cat("OK: zapisano reliability diagrams do:\n")
cat(fig_dir, "\n")