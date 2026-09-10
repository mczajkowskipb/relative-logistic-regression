get_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  getwd()
}
source(file.path(get_script_dir(), "07_world_protocol_utils.R"))

suppressPackageStartupMessages({
  library(jsonlite)
  library(foreign)
  library(ggplot2)
  library(svglite) # Dodane do obslugi SVG
})

# ============================================================
# CONFIG
# ============================================================
project_root <- Sys.getenv("RLR_PROJECT_ROOT", unset = normalizePath(file.path(get_script_dir(), "../..")))
ptb_root <- Sys.getenv("RLR_PTB_FIXED_ROOT", unset = file.path(project_root, "results", "ptb_fixed_published_state"))

train_file <- Sys.getenv("RLR_PTB_TRAIN_FILE", unset = file.path(project_root, "data", "cohort_norm", "01_train_data.arff"))

out_dir <- Sys.getenv("RLR_OUT_ROOT", unset = file.path(project_root, "results", "supplementary", "ptb_fixed_ora"))

gmt_dir <- Sys.getenv("RLR_MSIGDB_GMT_DIR", unset = file.path(project_root, "reference_outputs", "ptb_ora", "gmt"))

pre_top_features <- 2000
fdr_cutoff <- 0.05
top_n_terms <- 10

main_compare_models <- c("kTSP_20", "RLR_K20_m10_d2_pool1000")
supp_compare_models <- c("kTSP_10", "RLR_K10_m10_d2_pool1000")

# ============================================================
# SETUP
# ============================================================
ensure_dir(out_dir)
ensure_dir(file.path(out_dir, "tables"))
ensure_dir(file.path(out_dir, "figures"))
ensure_dir(file.path(out_dir, "logs"))

log_path <- file.path(out_dir, "logs", paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
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

require_file <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}

shorten_term <- function(x, max_n = 70) {
  x <- as.character(x)
  ifelse(nchar(x) <= max_n, x, paste0(substr(x, 1, max_n - 3), "..."))
}

mean_or_na <- function(x) {
  if (length(x) == 0) return(NA_real_)
  mean(x, na.rm = TRUE)
}

# ============================================================
# HELPERS: locating PTB fixed outputs
# ============================================================
find_single_subdir <- function(path) {
  dd <- list.dirs(path, recursive = FALSE, full.names = TRUE)
  if (length(dd) != 1) {
    stop("Expected exactly one subdir in: ", path, " but found: ", length(dd))
  }
  dd[1]
}

get_ptb_fixed_split_dir <- function(ptb_root) {
  worlds_dir <- file.path(ptb_root, "worlds")
  require_file(worlds_dir)

  world_dir <- file.path(worlds_dir, "ptb_norm_fixed")
  if (!dir.exists(world_dir)) {
    world_dir <- find_single_subdir(worlds_dir)
  }

  datasets_dir <- file.path(world_dir, "datasets")
  dataset_dir <- find_single_subdir(datasets_dir)

  splits_dir <- file.path(dataset_dir, "splits")
  split_dir <- find_single_subdir(splits_dir)

  split_dir
}

# ============================================================
# HELPERS: extracting gene lists
# ============================================================
extract_genes_from_selected_pairs <- function(path) {
  require_file(path)
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  needed <- c("gene_left", "gene_right")
  if (!all(needed %in% colnames(df))) {
    stop("selected_pairs.csv does not contain gene_left/gene_right: ", path)
  }

  unique(sort(c(df$gene_left, df$gene_right)))
}

extract_genes_from_active_pairs <- function(path) {
  require_file(path)
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  needed <- c("gene_left", "gene_right")
  if (!all(needed %in% colnames(df))) {
    stop("active_pairs.csv does not contain gene_left/gene_right: ", path)
  }

  unique(sort(c(df$gene_left, df$gene_right)))
}

extract_genes_from_coefficients <- function(path) {
  require_file(path)
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  if (!("ordered_pair" %in% colnames(df))) {
    stop("coefficients.csv missing ordered_pair: ", path)
  }

  beta_col <- NULL
  if ("beta" %in% colnames(df)) beta_col <- "beta"
  if (is.null(beta_col)) stop("coefficients.csv missing beta column: ", path)

  df <- df[!is.na(df[[beta_col]]) & df[[beta_col]] != 0, , drop = FALSE]
  if (nrow(df) == 0) return(character(0))

  parts <- strsplit(df$ordered_pair, " > ", fixed = TRUE)
  left <- vapply(parts, function(x) x[1], character(1))
  right <- vapply(parts, function(x) x[2], character(1))

  unique(sort(c(left, right)))
}

get_model_gene_list <- function(split_dir, model_name) {
  ktsp_path <- file.path(split_dir, "baselines", model_name, "selected_pairs.csv")
  if (file.exists(ktsp_path)) {
    return(list(
      model = model_name,
      family = "kTSP",
      genes = extract_genes_from_selected_pairs(ktsp_path),
      source_file = ktsp_path
    ))
  }

  rlr_active_path <- file.path(split_dir, "rlr", model_name, "active_pairs.csv")
  if (file.exists(rlr_active_path)) {
    return(list(
      model = model_name,
      family = "RLR",
      genes = extract_genes_from_active_pairs(rlr_active_path),
      source_file = rlr_active_path
    ))
  }

  rlr_coef_path <- file.path(split_dir, "rlr", model_name, "coefficients.csv")
  if (file.exists(rlr_coef_path)) {
    return(list(
      model = model_name,
      family = "RLR",
      genes = extract_genes_from_coefficients(rlr_coef_path),
      source_file = rlr_coef_path
    ))
  }

  stop("Could not find gene source for model: ", model_name)
}

# ============================================================
# HELPERS: background universe
# ============================================================
build_training_mad_universe <- function(train_file, top_n = 2000) {
  log_msg("Reading training file for background universe: ", train_file)
  xy <- read_xy_file(
    path = train_file,
    class_col = NULL,
    sample_id_col = NULL
  )

  pf <- prefilter_by_mad(
    X_train = xy$X,
    X_test = xy$X,
    top_n = top_n
  )

  sort(colnames(pf$X_train))
}

# ============================================================
# HELPERS: gene sets
# ============================================================
build_term2gene <- function(universe_genes) {
  log_msg("Loading local MSigDB 2026.1.Hs GMT gene sets")


  read_gmt <- function(path, source_name) {
    if (!file.exists(path)) stop("GMT file not found: ", path)

    x <- strsplit(readLines(path, warn = FALSE), "\\t", fixed = FALSE)

    rows <- lapply(x, function(v) {
      if (length(v) < 3L) return(NULL)

      genes <- v[-c(1L, 2L)]
      genes <- unique(genes[genes %in% universe_genes])

      if (length(genes) == 0L) return(NULL)

      data.frame(
        term = v[1L],
        gene = genes,
        source = source_name,
        stringsAsFactors = FALSE
      )
    })

    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0L) {
      return(data.frame(
        term = character(),
        gene = character(),
        source = character(),
        stringsAsFactors = FALSE
      ))
    }

    do.call(rbind, rows)
  }

  hallmark <- read_gmt(
    file.path(gmt_dir, "h.all.v2026.1.Hs.symbols.gmt"),
    "Hallmark"
  )

  gobp <- read_gmt(
    file.path(gmt_dir, "c5.go.bp.v2026.1.Hs.symbols.gmt"),
    "GO_BP"
  )

  reactome <- read_gmt(
    file.path(gmt_dir, "c2.cp.reactome.v2026.1.Hs.symbols.gmt"),
    "Reactome"
  )

  out <- rbind(hallmark, gobp, reactome)
  out <- out[!duplicated(out), , drop = FALSE]

  log_msg(
    "Filtered TERM2GENE: ",
    nrow(out), " rows; ",
    length(unique(out$gene)), " universe genes; ",
    length(unique(out$term[out$source == "Hallmark"])), " Hallmark terms; ",
    length(unique(out$term[out$source == "GO_BP"])), " GO BP terms; ",
    length(unique(out$term[out$source == "Reactome"])), " Reactome terms"
  )

  out
}

# ============================================================
# ORA: simple hypergeometric test
# ============================================================
run_simple_ora <- function(query_genes, universe_genes, term2gene, source_name,
                           fdr_cutoff = 0.05,
                           min_gs = 5,
                           max_gs = 500) {
  query_genes <- unique(intersect(query_genes, universe_genes))
  if (length(query_genes) == 0) return(data.frame())

  tt <- term2gene[term2gene$source == source_name, c("term", "gene"), drop = FALSE]
  if (nrow(tt) == 0) return(data.frame())

  split_terms <- split(tt$gene, tt$term)

  N <- length(universe_genes)
  K <- length(query_genes)

  rows <- list()
  ri <- 1L

  for (term in names(split_terms)) {
    geneset <- unique(intersect(split_terms[[term]], universe_genes))
    M <- length(geneset)
    if (M < min_gs || M > max_gs) next

    overlap <- sort(intersect(query_genes, geneset))
    x <- length(overlap)
    if (x == 0) next

    pval <- phyper(q = x - 1, m = M, n = N - M, k = K, lower.tail = FALSE)

    rows[[ri]] <- data.frame(
      source = source_name,
      term = term,
      query_size = K,
      universe_size = N,
      geneset_size = M,
      overlap_size = x,
      gene_ratio = x / K,
      p_value = pval,
      contributing_genes = paste(overlap, collapse = ";"),
      stringsAsFactors = FALSE
    )
    ri <- ri + 1L
  }

  if (length(rows) == 0) return(data.frame())

  df <- do.call(rbind, rows)
  df$p_adjust <- p.adjust(df$p_value, method = "BH")
  df <- df[order(df$p_adjust, df$p_value, -df$overlap_size), , drop = FALSE]
  df$term_short <- shorten_term(df$term, max_n = 80)

  sig <- df[df$p_adjust < fdr_cutoff, , drop = FALSE]
  sig
}

run_all_ora_for_model <- function(model_name, genes, universe_genes, term2gene, fdr_cutoff = 0.05) {
  sources <- unique(term2gene$source)

  all_res <- list()
  ai <- 1L

  for (src in sources) {
    rr <- run_simple_ora(
      query_genes = genes,
      universe_genes = universe_genes,
      term2gene = term2gene,
      source_name = src,
      fdr_cutoff = fdr_cutoff
    )

    if (nrow(rr) > 0) {
      rr$model <- model_name
      all_res[[ai]] <- rr
      ai <- ai + 1L
    }
  }

  if (length(all_res) == 0) return(data.frame())
  out <- do.call(rbind, all_res)
  out <- out[order(out$model, out$p_adjust, out$p_value), , drop = FALSE]
  out
}

top_n_by_model_source <- function(df, n = 10) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())

  keys <- unique(df[, c("model", "source"), drop = FALSE])

  out <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    idx <- df$model == key$model & df$source == key$source
    head(df[idx, , drop = FALSE], n)
  }))

  out
}

top_n_by_model <- function(df, n = 10) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())

  keys <- unique(df$model)

  out <- do.call(rbind, lapply(keys, function(mm) {
    dd <- df[df$model == mm, , drop = FALSE]
    dd <- dd[order(dd$p_adjust, dd$p_value), , drop = FALSE]
    head(dd, n)
  }))

  out
}

# ============================================================
# PLOTTING
# ============================================================
save_dotplot_compare <- function(enrichment_df, compare_models, out_prefix, top_n = 10) {
  df <- enrichment_df[enrichment_df$model %in% compare_models, , drop = FALSE]
  if (nrow(df) == 0) {
    log_msg("No enrichment results for plot: ", out_prefix)
    return(invisible(NULL))
  }

  df <- top_n_by_model(df, n = top_n)
  if (nrow(df) == 0) return(invisible(NULL))

  df$plot_score <- -log10(df$p_adjust)

  df$term_lab <- paste0(df$source, ": ", df$term_short)
  df <- df[order(df$model, df$plot_score), , drop = FALSE]
  df$term_model <- paste(df$model, df$term_lab, sep = "___")
  df$term_model <- factor(df$term_model, levels = unique(df$term_model))

  p <- ggplot(df, aes(x = plot_score, y = term_model, size = overlap_size, color = source)) +
    geom_point(alpha = 0.85) +
    facet_wrap(~model, scales = "free_y") +
    scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
    labs(
      x = expression(-log[10]("adjusted p-value")),
      y = NULL,
      size = "Overlap",
      color = "Source",
      title = "PTB fixed-split ORA comparison"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )

  # DODANO ZAPIS DO SVG
  ggsave(paste0(out_prefix, ".png"), device = grDevices::png, plot = p, width = 13, height = 7, dpi = 220)
  ggsave(paste0(out_prefix, ".pdf"), plot = p, width = 13, height = 7)
  ggsave(paste0(out_prefix, ".svg"), plot = p, width = 13, height = 7)
}

# ============================================================
# MAIN
# ============================================================
log_msg("Starting PTB biological case study")

split_dir <- get_ptb_fixed_split_dir(ptb_root)
log_msg("Detected PTB split dir: ", split_dir)

universe_genes <- build_training_mad_universe(train_file, top_n = pre_top_features)
write.csv(
  data.frame(gene = universe_genes, stringsAsFactors = FALSE),
  file.path(out_dir, "tables", "background_universe_mad2000.csv"),
  row.names = FALSE
)
log_msg("Universe size after MAD prefilter: ", length(universe_genes))

term2gene <- build_term2gene(universe_genes)
write.csv(term2gene, file.path(out_dir, "tables", "term2gene_universe_filtered.csv"), row.names = FALSE)

model_names <- c(
  "kTSP_10",
  "kTSP_20",
  "RLR_K10_m10_d2_pool1000",
  "RLR_K20_m10_d2_pool1000"
)

gene_sets <- lapply(model_names, function(mm) get_model_gene_list(split_dir, mm))
names(gene_sets) <- model_names

gene_list_rows <- list()
gi <- 1L
all_enrichment <- list()
ei <- 1L

for (mm in model_names) {
  info <- gene_sets[[mm]]

  genes_raw <- unique(info$genes)
  genes_in_universe <- unique(intersect(genes_raw, universe_genes))

  gene_list_rows[[gi]] <- data.frame(
    model = mm,
    family = info$family,
    n_genes_raw = length(genes_raw),
    n_genes_in_universe = length(genes_in_universe),
    source_file = info$source_file,
    stringsAsFactors = FALSE
  )
  gi <- gi + 1L

  write.csv(
    data.frame(gene = sort(genes_in_universe), stringsAsFactors = FALSE),
    file.path(out_dir, "tables", paste0("genes_", mm, ".csv")),
    row.names = FALSE
  )

  log_msg("Running ORA for ", mm, " | genes in universe = ", length(genes_in_universe))

  ora_df <- run_all_ora_for_model(
    model_name = mm,
    genes = genes_in_universe,
    universe_genes = universe_genes,
    term2gene = term2gene,
    fdr_cutoff = fdr_cutoff
  )

  if (nrow(ora_df) > 0) {
    all_enrichment[[ei]] <- ora_df
    ei <- ei + 1L
  }
}

gene_list_summary <- do.call(rbind, gene_list_rows)
write.csv(gene_list_summary, file.path(out_dir, "tables", "gene_list_summary.csv"), row.names = FALSE)

if (length(all_enrichment) == 0) {
  write.csv(data.frame(), file.path(out_dir, "tables", "enrichment_full_all.csv"), row.names = FALSE)
  write.csv(data.frame(), file.path(out_dir, "tables", "supp_table_top10_enrichment.csv"), row.names = FALSE)
  log_msg("No significant enrichment found under FDR < ", fdr_cutoff)
} else {
  enrichment_full <- do.call(rbind, all_enrichment)
  enrichment_full <- enrichment_full[order(enrichment_full$model, enrichment_full$p_adjust, enrichment_full$p_value), , drop = FALSE]
  write.csv(enrichment_full, file.path(out_dir, "tables", "enrichment_full_all.csv"), row.names = FALSE)

  top10 <- top_n_by_model_source(enrichment_full, n = top_n_terms)
  write.csv(top10, file.path(out_dir, "tables", "supp_table_top10_enrichment.csv"), row.names = FALSE)

  save_dotplot_compare(
    enrichment_df = enrichment_full,
    compare_models = main_compare_models,
    out_prefix = file.path(out_dir, "figures", "fig_ptb_fixed_ora_k20"),
    top_n = top_n_terms
  )

  save_dotplot_compare(
    enrichment_df = enrichment_full,
    compare_models = supp_compare_models,
    out_prefix = file.path(out_dir, "figures", "supp_fig_ptb_fixed_ora_k10"),
    top_n = top_n_terms
  )
}

# ============================================================
# MARKDOWN NOTES FOR PAPER
# ============================================================
md_path <- file.path(out_dir, "tables", "biological_case_study_notes.md")

cat(
  "# PTB fixed biological case study\n\n",
  "## What was done\n",
  "- Main case study: PTB fixed-split.\n",
  "- Gene lists extracted for:\n",
  "  - kTSP_10\n",
  "  - kTSP_20\n",
  "  - RLR_K10_m10_d2_pool1000\n",
  "  - RLR_K20_m10_d2_pool1000\n",
  "- For kTSP, genes were taken from `selected_pairs.csv`.\n",
  "- For RLR, genes were taken from `active_pairs.csv` (or from non-zero coefficients if needed).\n",
  "- The background universe was **not** the whole genome. It was reconstructed from the PTB training set using the same training-only MAD prefilter (`top 2000 features`) that had been used before pair scoring.\n",
  "- Over-representation analysis (ORA) was performed separately for each model using Hallmark, GO Biological Process, and Reactome gene sets.\n",
  "- Only terms with FDR < 0.05 were retained.\n\n",
  "## Intended interpretation\n",
  "This analysis should not be presented as validation of a biomarker panel. ",
  "The appropriate interpretation is that RLR does not necessarily introduce a completely different biological theme, ",
  "but may sharpen the signal toward a more coherent inflammatory / interferon / host-response core relative to classical kTSP shortlist genes.\n\n",
  "## Key output files\n",
  "- `background_universe_mad2000.csv`\n",
  "- `gene_list_summary.csv`\n",
  "- `genes_kTSP_10.csv`\n",
  "- `genes_kTSP_20.csv`\n",
  "- `genes_RLR_K10_m10_d2_pool1000.csv`\n",
  "- `genes_RLR_K20_m10_d2_pool1000.csv`\n",
  "- `enrichment_full_all.csv`\n",
  "- `supp_table_top10_enrichment.csv`\n",
  "- `fig_ptb_fixed_ora_k20.png/pdf/svg`\n",
  "- `supp_fig_ptb_fixed_ora_k10.png/pdf/svg`\n",
  file = md_path
)

log_msg("DONE")
log_msg("Outputs saved in: ", out_dir)
log_msg("Log file: ", log_path)