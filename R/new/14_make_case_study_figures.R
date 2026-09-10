suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Użycie: Rscript 14_make_case_study_figures.R <config.json>")
}

cfg <- read_json(args[1], simplifyVector = TRUE)
out_root <- cfg$out_root
if (is.null(out_root)) stop("Brak out_root w config.json")

fig_dir <- file.path(out_root, "paper_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

manifest_path <- file.path(out_root, "config", "run_manifest.csv")
if (!file.exists(manifest_path)) stop("Brak run_manifest.csv")

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a

case_cfg <- cfg$case_study %||% list()

world_name <- case_cfg$world_name %||% stop("Brak case_study$world_name")
dataset_id <- case_cfg$dataset_id %||% stop("Brak case_study$dataset_id")
split_id   <- case_cfg$split_id %||% stop("Brak case_study$split_id")

default_models <- c(
  "kTSP_10",
  "kTSP_50",
  "RLR_K10_m10_d2_pool1000",
  "RLR_K50_m10_d2_pool1000"
)

models_to_plot <- case_cfg$models %||% default_models
top_beta_n <- as.integer(case_cfg$top_beta_n %||% 15)
top_hub_n  <- as.integer(case_cfg$top_hub_n %||% 10)

safe_read_json <- function(path) {
  tryCatch(read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}

safe_read_csv <- function(path) {
  tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
}

save_plot_both <- function(p, stem, width = 9, height = 5) {
  ggsave(filename = file.path(fig_dir, paste0(stem, ".pdf")), plot = p, width = width, height = height)
  ggsave(filename = file.path(fig_dir, paste0(stem, ".png")), plot = p, width = width, height = height, dpi = 300)
}

sel <- manifest[
  manifest$world_name == world_name &
    manifest$dataset_id == dataset_id &
    manifest$split_id == split_id &
    manifest$model %in% models_to_plot,
  ,
  drop = FALSE
]

if (nrow(sel) == 0) {
  stop("Nie znaleziono runów dla wybranego case study.")
}

# ------------------------------------------------------------
# 1. summary table for selected models
# ------------------------------------------------------------
summary_rows <- list()

for (i in seq_len(nrow(sel))) {
  mrow <- sel[i, ]
  metrics_path <- file.path(mrow$path, "metrics_test.json")
  mj <- safe_read_json(metrics_path)
  if (is.null(mj)) next

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    world_name = world_name,
    dataset_id = dataset_id,
    split_id = split_id,
    model = mrow$model,
    auc = mj$auc %||% NA,
    f1 = mj$f1 %||% NA,
    logloss = mj$logloss %||% NA,
    brier = mj$brier %||% NA,
    ece = mj$ece %||% NA,
    global_active_pairs = mj$global_active_pairs %||% NA,
    unique_genes = mj$unique_genes %||% NA,
    local_active_mean = mj$local_active_mean %||% NA,
    local_active_median = mj$local_active_median %||% NA,
    local_active_iqr = mj$local_active_iqr %||% NA,
    hub_top_share = mj$hub_top_share %||% NA,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(fig_dir, "case_metrics_summary.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 2. local activity boxplot
# ------------------------------------------------------------
local_rows <- list()

for (i in seq_len(nrow(sel))) {
  mrow <- sel[i, ]
  f <- file.path(mrow$path, "local_activity_test.csv")
  dd <- safe_read_csv(f)
  if (is.null(dd)) next
  dd$model <- mrow$model
  local_rows[[length(local_rows) + 1]] <- dd
}

if (length(local_rows) > 0) {
  local_df <- do.call(rbind, local_rows)
  write.csv(local_df, file.path(fig_dir, "case_local_activity_data.csv"), row.names = FALSE)

  p1 <- ggplot(local_df, aes(x = model, y = local_active_n)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = NULL,
      y = "Locally active relations",
      title = paste0("Local activity: ", dataset_id, " / ", split_id)
    )

  save_plot_both(p1, "case_local_activity_boxplot", width = 9, height = 5)

  p2 <- ggplot(local_df, aes(x = local_active_n)) +
    geom_histogram(bins = 20) +
    facet_wrap(~ model, scales = "free_y") +
    theme_bw(base_size = 11) +
    labs(
      x = "Locally active relations",
      y = "Count",
      title = paste0("Local activity distributions: ", dataset_id, " / ", split_id)
    )

  save_plot_both(p2, "case_local_activity_histogram", width = 11, height = 6)
}

# ------------------------------------------------------------
# 3. top beta plot for RLR models
# ------------------------------------------------------------
coef_rows <- list()

for (i in seq_len(nrow(sel))) {
  mrow <- sel[i, ]
  if (!grepl("^RLR_", mrow$model)) next

  f <- file.path(mrow$path, "coefficients.csv")
  dd <- safe_read_csv(f)
  if (is.null(dd)) next

  dd$abs_beta <- suppressWarnings(as.numeric(dd$abs_beta))
  dd$beta <- suppressWarnings(as.numeric(dd$beta))
  dd <- dd[!is.na(dd$abs_beta), , drop = FALSE]
  dd <- dd[order(dd$abs_beta, decreasing = TRUE), , drop = FALSE]
  dd <- head(dd, top_beta_n)
  dd$model <- mrow$model
  coef_rows[[length(coef_rows) + 1]] <- dd
}

if (length(coef_rows) > 0) {
  coef_df <- do.call(rbind, coef_rows)
  write.csv(coef_df, file.path(fig_dir, "case_rlr_top_betas_data.csv"), row.names = FALSE)

  coef_df$ordered_pair <- factor(coef_df$ordered_pair, levels = rev(unique(coef_df$ordered_pair)))

  p3 <- ggplot(coef_df, aes(x = ordered_pair, y = beta)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ model, scales = "free_y") +
    theme_bw(base_size = 10) +
    labs(
      x = NULL,
      y = "Coefficient",
      title = paste0("Top weighted relations: ", dataset_id, " / ", split_id)
    )

  save_plot_both(p3, "case_rlr_top_betas", width = 12, height = 8)
}

# ------------------------------------------------------------
# 4. hub counts
# ------------------------------------------------------------
hub_rows <- list()

for (i in seq_len(nrow(sel))) {
  mrow <- sel[i, ]
  f <- file.path(mrow$path, "hub_counts.csv")
  dd <- safe_read_csv(f)
  if (is.null(dd)) next

  dd$count <- suppressWarnings(as.numeric(dd$count))
  dd$share <- suppressWarnings(as.numeric(dd$share))
  dd <- dd[order(dd$count, decreasing = TRUE), , drop = FALSE]
  dd <- head(dd, top_hub_n)
  dd$model <- mrow$model
  hub_rows[[length(hub_rows) + 1]] <- dd
}

if (length(hub_rows) > 0) {
  hub_df <- do.call(rbind, hub_rows)
  write.csv(hub_df, file.path(fig_dir, "case_hub_counts_data.csv"), row.names = FALSE)

  p4 <- ggplot(hub_df, aes(x = gene, y = share)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ model, scales = "free_y") +
    theme_bw(base_size = 10) +
    labs(
      x = NULL,
      y = "Hub share",
      title = paste0("Hub concentration: ", dataset_id, " / ", split_id)
    )

  save_plot_both(p4, "case_hub_counts", width = 11, height = 7)
}

cat("OK: zapisano case-study figury do:\n")
cat(fig_dir, "\n")