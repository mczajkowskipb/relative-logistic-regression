# RLMT GATE-3 audit: read-only wrt existing experiment artifacts.
# Rebuilds the exact outer/inner stratified indices from the frozen config,
# counts candidate outcome categories and inspects local class support.
# No fitting, re-ranking, test-label selection, or production code changes.
# Usage: Rscript --vanilla R/new/31_rlmt_gate3_audit.R results/rlmt_gate3_pilot_2x3/config.json

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 31_rlmt_gate3_audit.R <pilot_config.json>")
suppressPackageStartupMessages(library(jsonlite))
source("R/new/07_world_protocol_utils.R")
cfg <- jsonlite::read_json(args[[1]], simplifyVector = TRUE)

n_worlds <- if (is.data.frame(cfg$worlds)) nrow(cfg$worlds) else length(cfg$worlds)
if (n_worlds != 1L) stop("Expected exactly one world")
world <- if (is.data.frame(cfg$worlds)) as.list(cfg$worlds[1L, , drop = FALSE]) else cfg$worlds[[1L]]
datasets <- as.character(world$datasets_include[[1L]])
if (!identical(sort(datasets), c("chen-2002", "singh-2002")) ||
    !identical(as.character(world$protocol), "compcancer_cv") ||
    as.integer(world$n_folds) != 5L ||
    as.integer(world$n_repeats) != 1L ||
    as.integer(world$max_splits_per_dataset) != 3L ||
    as.integer(cfg$models$seed) != 42L ||
    as.integer(cfg$models$rlmt$inner_folds) != 3L) {
  stop("Expected original GATE-3 pilot: Chen+Singh, 5 folds, first 3, seed 42")
}
root <- cfg$out_root
status_file <- file.path(root, "status.txt")
if (!file.exists(status_file) ||
    !identical(trimws(readLines(status_file, warn = FALSE)[1]), "RUNNER_EXIT=0")) {
  stop("Pilot runner not confirmed completed successfully")
}
dest <- file.path(root, "audit_gate3_readonly_v1")
if (dir.exists(dest)) stop("STOP: Audit already exists; refusing to overwrite: ", dest)

all_candidates <- list()
fold_summary <- list()
counter <- 0L

for (dataset in datasets) {
  xy <- read_xy_file(file.path(world$comp_dir, paste0(dataset, ".arff")))
  outer <- make_repeated_stratified_folds(
    xy$y, n_folds = as.integer(world$n_folds), n_repeats = 1L,
    seed = as.integer(cfg$models$seed))
  for (j in seq_len(3L)) {
    spec <- outer[[j]]
    split <- spec$split_id
    base <- file.path(root, "worlds", world$name, "datasets",
                      dataset, "splits", split)
    l2 <- file.path(base, "rlmt", "RLMT_K10_m10_d2_pool1000_L2")
    look_file <- file.path(l2, "lookahead.csv")
    screen_file <- file.path(l2, "screen_top.csv")
    tree_file <- file.path(l2, "tree.json")
    for (p in c(look_file, screen_file, tree_file,
                file.path(base, "samples_train.csv"),
                file.path(base, "samples_test.csv"))) {
      if (!file.exists(p)) stop("Missing pilot artifact: ", p)
    }

    Xtr <- xy$X[spec$train_idx, , drop = FALSE]
    ytr <- xy$y[spec$train_idx]
    stored_train <- read.csv(file.path(base, "samples_train.csv"),
                             stringsAsFactors = FALSE)
    stored_test <- read.csv(file.path(base, "samples_test.csv"),
                            stringsAsFactors = FALSE)
    if (!identical(stored_train$sample_id, rownames(Xtr)) ||
        !identical(stored_test$sample_id, rownames(xy$X)[spec$test_idx])) {
      stop("Rebuilt outer split does not match persisted sample IDs: ", dataset, "/", split)
    }

    inner <- make_repeated_stratified_folds(
      ytr, n_folds = as.integer(cfg$models$rlmt$inner_folds),
      n_repeats = 1L, seed = as.integer(cfg$models$rlmt$seed))
    look <- read.csv(look_file, stringsAsFactors = FALSE)
    screen <- read.csv(screen_file, stringsAsFactors = FALSE)
    tree <- jsonlite::read_json(tree_file, simplifyVector = TRUE)
    if (nrow(look) != 20L || nrow(screen) != 20L ||
        !identical(as.character(look$gene_a), as.character(screen$gene_a)) ||
        !identical(as.character(look$gene_b), as.character(screen$gene_b))) {
      stop("Candidate accounting / top-B mismatch: ", dataset, "/", split)
    }
    if (!isTRUE(tree$size == 1L) || isTRUE(tree$accepted_split)) {
      stop("Pilot tree changed or includes an accepted split: ", dataset, "/", split)
    }

    for (k in seq_len(nrow(look))) {
      r <- look[k, , drop = FALSE]
      ga <- as.character(r$gene_a[[1]])
      gb <- as.character(r$gene_b[[1]])
      if (!(ga %in% names(Xtr) && gb %in% names(Xtr))) {
        stop("Candidate genes absent from raw outer training: ", ga, " / ", gb)
      }
      z <- Xtr[[ga]] > Xtr[[gb]]
      if (anyNA(z)) stop("NA in routing: ", dataset, "/", split)
      outer_class <- lapply(list(which(z), which(!z)), function(ix)
        as.integer(table(ytr[ix])))
      inner_counts <- t(vapply(inner, function(f) {
        tr <- f$train_idx
        side <- list(tr[z[tr]], tr[!z[tr]])
        class_counts <- lapply(side, function(ix) as.integer(table(ytr[ix])))
        c(min_child_n = min(lengths(side)),
          min_child_class = min(unlist(class_counts, use.names = FALSE)),
          min_child_class0 = min(vapply(class_counts, function(x) x[1L], integer(1))),
          min_child_class1 = min(vapply(class_counts, function(x) x[2L], integer(1))))
      }, numeric(4L)))
      reason <- as.character(r$rejection_reason[[1]])
      if (is.na(reason)) reason <- ""
      valid <- isTRUE(r$valid[[1]])
      group <- if (valid) "valid" else if (identical(reason, "inner_child_support")) {
        "inner_child_support"
      } else if (grepl("class has 1 or 0 observations", reason, fixed = TRUE)) {
        "glmnet_class_error"
      } else "other_fit_error"
      counter <- counter + 1L
      all_candidates[[counter]] <- data.frame(
        dataset = dataset, split_id = split, rank = k,
        gene_a = ga, gene_b = gb, valid = valid,
        reason_group = group, rejection_reason = reason,
        mean_gain = as.numeric(r$mean_gain), se_gain = as.numeric(r$se_gain),
        lower_gain = as.numeric(r$lower_gain),
        parent_oof_logloss = as.numeric(r$parent_oof_logloss),
        split_oof_logloss = as.numeric(r$split_oof_logloss),
        outer_min_child_n = min(sum(z), sum(!z)),
        outer_min_child_class = min(unlist(outer_class, use.names = FALSE)),
        inner_min_child_n = min(inner_counts[, "min_child_n"]),
        inner_min_child_class = min(inner_counts[, "min_child_class"]),
        inner_min_child_class0 = min(inner_counts[, "min_child_class0"]),
        inner_min_child_class1 = min(inner_counts[, "min_child_class1"]),
        stringsAsFactors = FALSE)
    }

    current <- do.call(rbind, tail(all_candidates, nrow(look)))
    v <- current[current$valid, , drop = FALSE]
    best_mean <- if (nrow(v)) max(v$mean_gain) else NA_real_
    best_lower <- if (nrow(v)) max(v$lower_gain) else NA_real_
    fold_summary[[length(fold_summary) + 1L]] <- data.frame(
      dataset = dataset, split_id = split,
      train_n = nrow(Xtr), train_class0 = as.integer(table(ytr)[1]),
      train_class1 = as.integer(table(ytr)[2]),
      candidates = nrow(look), valid = nrow(v),
      invalid = sum(!current$valid),
      positive_mean_gain = sum(v$mean_gain > 0),
      positive_lower_gain = sum(v$lower_gain > 0),
      best_mean_gain = best_mean, best_lower_gain = best_lower,
      accepted = isTRUE(tree$accepted_split), stringsAsFactors = FALSE)
  }
}

candidates <- do.call(rbind, all_candidates)
folds <- do.call(rbind, fold_summary)
if (nrow(candidates) != 120L || nrow(folds) != 6L ||
    sum(candidates$valid) != 96L || sum(!candidates$valid) != 24L) {
  stop("Pilot totals do not match previous artifact audit; investigate before reporting")
}

lines <- character()
emit <- function(...) {
  s <- paste0(..., collapse = "")
  cat(s, "\n", sep = "")
  lines <<- c(lines, s)
}
emit("RLMT GATE-3 READ-ONLY SCIENTIFIC AUDIT")
emit("Original experiments untouched; no fitting and no parameter changes")
emit("6 folds; 120 candidates; valid=", sum(candidates$valid),
     "; invalid=", sum(!candidates$valid),
     "; mean_gain>0=", sum(candidates$valid & candidates$mean_gain > 0, na.rm = TRUE),
     "; lower_gain>0=", sum(candidates$valid & candidates$lower_gain > 0, na.rm = TRUE))
emit("")
emit("=== Per outer fold ===")
for (i in seq_len(nrow(folds))) {
  f <- folds[i, ]
  emit(sprintf("%-11s %-12s train=%3d class=%d/%d valid=%2d/20 positive=%d best_mean=%+.6f best_lower=%+.6f",
               f$dataset, f$split_id, f$train_n, f$train_class0,
               f$train_class1, f$valid, f$positive_mean_gain,
               f$best_mean_gain, f$best_lower_gain))
}
emit("")
emit("=== Invalid causes ===")
for (group in c("glmnet_class_error", "inner_child_support", "other_fit_error")) {
  sub <- candidates[candidates$reason_group == group, , drop = FALSE]
  emit(group, " total=", nrow(sub),
       " chen=", sum(sub$dataset == "chen-2002"),
       " singh=", sum(sub$dataset == "singh-2002"))
  if (nrow(sub)) {
    emit("  inner_min_child_class: min=", min(sub$inner_min_child_class),
         " median=", stats::median(sub$inner_min_child_class),
         " max=", max(sub$inner_min_child_class),
         " <=3=", sum(sub$inner_min_child_class <= 3L),
         " <=5=", sum(sub$inner_min_child_class <= 5L))
    for (reason in unique(sub$rejection_reason))
      emit("  reason [", sum(sub$rejection_reason == reason), "]: ", reason)
  }
}
emit("")
emit("=== Inner-fold child support distribution ===")
for (group in c("valid", "glmnet_class_error", "inner_child_support", "other_fit_error")) {
  sub <- candidates[candidates$reason_group == group, , drop = FALSE]
  if (nrow(sub)) emit(sprintf(
    "%-20s n=%3d inner_min_class=[%d .. %d] median=%.1f min_inner_child_n=%d",
    group, nrow(sub), min(sub$inner_min_child_class),
    max(sub$inner_min_child_class), median(sub$inner_min_child_class),
    min(sub$inner_min_child_n)))
}
emit("")
emit("=== Best VALID candidate by 1-SE lower gain, per fold ===")
for (i in seq_len(nrow(folds))) {
  f <- folds[i, ]
  v <- candidates[candidates$dataset == f$dataset &
                  candidates$split_id == f$split_id & candidates$valid, , drop = FALSE]
  if (!nrow(v)) next
  w <- v[which.max(v$lower_gain), , drop = FALSE]
  emit(sprintf("%s %s %s > %s mean=%+.6f SE=%.6f lower=%+.6f inner_min_class=%d",
               f$dataset, f$split_id, w$gene_a, w$gene_b,
               w$mean_gain, w$se_gain, w$lower_gain, w$inner_min_child_class))
}

dir.create(dest, recursive = TRUE, showWarnings = FALSE)
write.csv(candidates, file.path(dest, "candidate_support.csv"), row.names = FALSE)
write.csv(folds, file.path(dest, "fold_summary.csv"), row.names = FALSE)
writeLines(c(lines,
             "",
             "METHODOLOGY CAUTIONS:",
             "Global pair pool is supervised and fitted on the full outer training split.",
             "Inner OOF lookahead reuses the pool and screening selected with inner validation labels.",
             "Thus inner OOF gains are for model choice, not unbiased nested estimates.",
             "The untouched outer test is separate, but with zero accepted splits no empirical improvement has been demonstrated.",
             "glmnet CV uses its own folds and may fail despite an inner child passing min_class >= 3.",
             "Do not interpret forced functional smoke outputs as scientific results."),
           file.path(dest, "audit.txt"))
emit("")
emit("AUDIT_RESULT=PASS")
emit("Output: ", dest)
