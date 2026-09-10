suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Użycie: Rscript 11_make_paper_figures.R <config.json>")
}

cfg <- read_json(args[1], simplifyVector = TRUE)
out_root <- cfg$out_root
if (is.null(out_root)) stop("Brak out_root w config.json")

summ_dir <- file.path(out_root, "summaries")
fig_dir <- file.path(out_root, "paper_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dataset_metrics_path <- file.path(summ_dir, "dataset_level_metrics.csv")
if (!file.exists(dataset_metrics_path)) stop("Brak dataset_level_metrics.csv")

df <- read.csv(dataset_metrics_path, stringsAsFactors = FALSE, check.names = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a

default_plot_models <- c(
  "value_LR",
  "kTSP_10",
  "kTSP_50",
  "RLR_K10_m10_d1",
  "RLR_K50_m10_d1"
)

plot_models <- cfg$paper_report$plot_models %||% default_plot_models
plot_worlds <- cfg$paper_report$plot_worlds %||% unique(df$world_name)

df <- df[df$model %in% plot_models & df$world_name %in% plot_worlds, , drop = FALSE]
df$model <- factor(df$model, levels = plot_models)

save_plot_both <- function(p, stem, width = 9, height = 5) {
  ggsave(filename = file.path(fig_dir, paste0(stem, ".pdf")), plot = p, width = width, height = height)
  ggsave(filename = file.path(fig_dir, paste0(stem, ".png")), plot = p, width = width, height = height, dpi = 300)
}

# ------------------------------------------------------------
# Figure 1: score-quality
# ------------------------------------------------------------
score_df <- rbind(
  data.frame(world_name = df$world_name, dataset_id = df$dataset_id, model = df$model, metric = "logloss", value = df$logloss),
  data.frame(world_name = df$world_name, dataset_id = df$dataset_id, model = df$model, metric = "brier", value = df$brier),
  data.frame(world_name = df$world_name, dataset_id = df$dataset_id, model = df$model, metric = "ece", value = df$ece)
)

p1 <- ggplot(score_df, aes(x = model, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
  facet_grid(metric ~ world_name, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Value", title = "Probability-oriented metrics across datasets")

save_plot_both(p1, "fig1_score_quality", width = 12, height = 7)

# ------------------------------------------------------------
# Figure 2: discrimination
# ------------------------------------------------------------
disc_df <- rbind(
  data.frame(world_name = df$world_name, dataset_id = df$dataset_id, model = df$model, metric = "auc", value = df$auc),
  data.frame(world_name = df$world_name, dataset_id = df$dataset_id, model = df$model, metric = "f1", value = df$f1)
)

p2 <- ggplot(disc_df, aes(x = model, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
  facet_grid(metric ~ world_name, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Value", title = "Discrimination metrics across datasets")

save_plot_both(p2, "fig2_discrimination", width = 12, height = 6)

# ------------------------------------------------------------
# Figure 3: structure
# ------------------------------------------------------------
struct_df <- rbind(
  data.frame(world_name = df$world_name, dataset_id = df$dataset_id, model = df$model, metric = "global_active_pairs", value = df$global_active_pairs),
  data.frame(world_name = df$world_name, dataset_id = df$dataset_id, model = df$model, metric = "unique_genes", value = df$unique_genes),
  data.frame(world_name = df$world_name, dataset_id = df$dataset_id, model = df$model, metric = "local_active_mean", value = df$local_active_mean),
  data.frame(world_name = df$world_name, dataset_id = df$dataset_id, model = df$model, metric = "hub_top_share", value = df$hub_top_share)
)

p3 <- ggplot(struct_df, aes(x = model, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
  facet_grid(metric ~ world_name, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Value", title = "Structural summaries across datasets")

save_plot_both(p3, "fig3_structure", width = 12, height = 8)

# ------------------------------------------------------------
# Figure 4: optional ablation for RLR only
# ------------------------------------------------------------
rlr_df <- df[grepl("^RLR_", df$model), , drop = FALSE]

if (nrow(rlr_df) > 0 && all(c("K_cfg","m_cfg","d_cfg") %in% colnames(rlr_df))) {
  rlr_df$label <- paste0("K", rlr_df$K_cfg, "_m", rlr_df$m_cfg, "_d", rlr_df$d_cfg)

  abl_df <- rbind(
    data.frame(world_name = rlr_df$world_name, dataset_id = rlr_df$dataset_id, label = rlr_df$label, metric = "logloss", value = rlr_df$logloss),
    data.frame(world_name = rlr_df$world_name, dataset_id = rlr_df$dataset_id, label = rlr_df$label, metric = "local_active_mean", value = rlr_df$local_active_mean),
    data.frame(world_name = rlr_df$world_name, dataset_id = rlr_df$dataset_id, label = rlr_df$label, metric = "hub_top_share", value = rlr_df$hub_top_share)
  )

  p4 <- ggplot(abl_df, aes(x = label, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1.0) +
    facet_grid(metric ~ world_name, scales = "free_y") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "Value", title = "RLR ablation view")

  save_plot_both(p4, "supp_fig_rlr_ablation", width = 13, height = 8)
}

cat("OK: zapisano wykresy do:\n")
cat(fig_dir, "\n")