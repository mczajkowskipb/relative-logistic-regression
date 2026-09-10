suppressPackageStartupMessages({
  library(foreign)
  library(switchBox)
  library(glmnet)
  library(pROC)
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

write_error_json <- function(path, e) {
  ensure_dir(dirname(path))
  jsonlite::write_json(
    list(error = as.character(e)),
    path,
    auto_unbox = TRUE,
    pretty = TRUE
  )
}

safe_read_csv_like <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    return(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (ext %in% c("tsv", "txt")) {
    return(read.delim(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  stop(paste("Nieobslugiwane rozszerzenie tabelaryczne:", ext, "dla pliku", path))
}

read_xy_file <- function(path,
                         class_col = NULL,
                         class_col_candidates = c("class", "Class", "CLASS", "label", "Label", "y", "Y"),
                         sample_id_col = NULL) {
  ext <- tolower(tools::file_ext(path))

  dat <- NULL
  if (ext == "arff") {
    dat <- foreign::read.arff(path)
  } else if (ext %in% c("csv", "tsv", "txt")) {
    dat <- safe_read_csv_like(path)
  } else {
    stop(paste("Nieobslugiwany format pliku:", path))
  }

  if (is.null(class_col)) {
    hits <- intersect(class_col_candidates, colnames(dat))
    if (length(hits) == 0) {
      stop(paste("Nie znalazlem kolumny klasy w:", path))
    }
    if (length(hits) > 1) {
      stop(paste("Znalazlem wiele mozliwych kolumn klasy w:", path, "=>", paste(hits, collapse = ", ")))
    }
    class_col <- hits[1]
  }

  if (!(class_col %in% colnames(dat))) {
    stop(paste("Brak wskazanej kolumny klasy:", class_col, "w pliku", path))
  }

  if (!is.null(sample_id_col)) {
    if (!(sample_id_col %in% colnames(dat))) {
      stop(paste("Brak sample_id_col:", sample_id_col, "w pliku", path))
    }
    sample_ids <- as.character(dat[[sample_id_col]])
  } else if (!is.null(rownames(dat)) && !all(rownames(dat) %in% as.character(seq_len(nrow(dat))))) {
    sample_ids <- rownames(dat)
  } else {
    sample_ids <- paste0("sample_", seq_len(nrow(dat)))
  }

  y <- factor(dat[[class_col]])

  drop_cols <- c(class_col)
  if (!is.null(sample_id_col)) drop_cols <- c(drop_cols, sample_id_col)

  X <- dat[, setdiff(colnames(dat), drop_cols), drop = FALSE]

  keep_num <- vapply(X, is.numeric, logical(1))
  X <- X[, keep_num, drop = FALSE]

  if (ncol(X) == 0) {
    stop(paste("Po odrzuceniu nienumerycznych kolumn nie zostaly zadne cechy w:", path))
  }

  rownames(X) <- sample_ids

  if (length(levels(y)) != 2) {
    stop(paste("Ten runner obsluguje tylko klasyfikacje binarna. Plik:", path, "ma", length(levels(y)), "klasy."))
  }

  list(
    X = X,
    y = y,
    sample_ids = sample_ids,
    class_col = class_col,
    path = path
  )
}

prefilter_by_mad <- function(X_train, X_test, X_valid = NULL, top_n = 2000) {
  if (is.null(top_n) || ncol(X_train) <= top_n) {
    out <- list(X_train = X_train, X_test = X_test, keep = colnames(X_train))
    if (!is.null(X_valid)) out$X_valid <- X_valid
    return(out)
  }

  mad_scores <- apply(X_train, 2, mad, na.rm = TRUE)
  mad_scores[is.na(mad_scores)] <- -Inf

  keep <- names(sort(mad_scores, decreasing = TRUE))[seq_len(min(top_n, length(mad_scores)))]

  out <- list(
    X_train = X_train[, keep, drop = FALSE],
    X_test  = X_test[, keep, drop = FALSE],
    keep = keep
  )

  if (!is.null(X_valid)) {
    out$X_valid <- X_valid[, keep, drop = FALSE]
  }

  out
}

maybe_prefilter_split <- function(X_train, X_test, X_valid = NULL, pre_tsp_top_features = NULL) {
  if (is.null(pre_tsp_top_features)) {
    out <- list(
      X_train = X_train,
      X_test = X_test,
      prefilter_applied = FALSE,
      kept_n = ncol(X_train)
    )
    if (!is.null(X_valid)) out$X_valid <- X_valid
    return(out)
  }

  pf <- prefilter_by_mad(
    X_train = X_train,
    X_test = X_test,
    X_valid = X_valid,
    top_n = pre_tsp_top_features
  )

  out <- list(
    X_train = pf$X_train,
    X_test = pf$X_test,
    prefilter_applied = TRUE,
    kept_n = ncol(pf$X_train),
    keep = pf$keep
  )
  if (!is.null(pf$X_valid)) out$X_valid <- pf$X_valid
  out
}

align_xy_train_test <- function(X_train, X_test) {
  common <- intersect(colnames(X_train), colnames(X_test))
  if (length(common) == 0) {
    stop("Brak wspolnych cech miedzy train i test.")
  }

  common <- sort(common)
  list(
    X_train = X_train[, common, drop = FALSE],
    X_test = X_test[, common, drop = FALSE],
    common_features = common
  )
}

combine_xy_list <- function(xy_list) {
  if (length(xy_list) == 0) stop("Pusta lista cohort/train files do polaczenia.")

  common <- Reduce(intersect, lapply(xy_list, function(z) colnames(z$X)))
  if (length(common) == 0) {
    stop("Brak wspolnych cech miedzy plikami do polaczenia.")
  }
  common <- sort(common)

  x_blocks <- lapply(xy_list, function(z) {
    xx <- z$X[, common, drop = FALSE]
    prefix <- tools::file_path_sans_ext(basename(z$path))
    rownames(xx) <- paste0(prefix, "__", rownames(xx))
    xx
  })

  X_all <- do.call(rbind, x_blocks)
  y_all <- factor(unlist(lapply(xy_list, function(z) as.character(z$y))))

  if (nrow(X_all) != length(y_all)) {
    stop("Niezgodnosc liczby probek po polaczeniu X i y.")
  }

  list(
    X = X_all,
    y = y_all,
    sample_ids = rownames(X_all),
    common_features = common
  )
}

xy_to_switchbox_mat <- function(X) {
  mat <- t(as.matrix(X))
  storage.mode(mat) <- "numeric"
  rownames(mat) <- colnames(X)
  colnames(mat) <- rownames(X)
  mat
}

safe_nfolds <- function(y, requested = 5) {
  m <- min(table(y))
  max(2, min(as.integer(requested), as.integer(m)))
}

calc_brier <- function(prob, y01) {
  mean((prob - y01)^2)
}

calc_ece <- function(prob, y01, n_bins = 10) {
  breaks <- seq(0, 1, length.out = n_bins + 1)
  bin_id <- cut(prob, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  ece <- 0

  for (b in seq_len(n_bins)) {
    idx <- which(bin_id == b)
    if (length(idx) == 0) next
    conf <- mean(prob[idx])
    acc <- mean(y01[idx])
    ece <- ece + (length(idx) / length(prob)) * abs(acc - conf)
  }
  ece
}

make_repeated_stratified_folds <- function(y, n_folds = 5, n_repeats = 10, seed = 42) {
  set.seed(seed)
  y <- factor(y)
  idx_all <- seq_along(y)

  out <- list()
  counter <- 1L

  for (rep_id in seq_len(n_repeats)) {
    fold_assign <- integer(length(y))

    for (cls in levels(y)) {
      cls_idx <- which(y == cls)
      cls_idx <- sample(cls_idx, length(cls_idx), replace = FALSE)
      bins <- split(cls_idx, rep(seq_len(n_folds), length.out = length(cls_idx)))

      for (f in seq_len(n_folds)) {
        fold_assign[bins[[f]]] <- f
      }
    }

    for (fold_id in seq_len(n_folds)) {
      out[[counter]] <- list(
        rep = rep_id - 1L,
        fold = fold_id - 1L,
        split_id = sprintf("rep%02d_fold%02d", rep_id - 1L, fold_id - 1L),
        train_idx = idx_all[fold_assign != fold_id],
        test_idx = idx_all[fold_assign == fold_id]
      )
      counter <- counter + 1L
    }
  }

  out
}

make_fixed_split_spec <- function(name = "fixed") {
  list(
    rep = 0L,
    fold = 0L,
    split_id = name
  )
}

make_loco_specs_from_files <- function(file_paths) {
  file_paths <- sort(file_paths)
  out <- vector("list", length(file_paths))

  for (i in seq_along(file_paths)) {
    holdout <- file_paths[i]
    train_files <- file_paths[-i]
    cohort_id <- sub("\\.[^.]+$", "", basename(holdout))

    out[[i]] <- list(
      rep = 0L,
      fold = i - 1L,
      split_id = paste0("holdout_", cohort_id),
      holdout_file = holdout,
      train_files = train_files,
      cohort_id = cohort_id
    )
  }

  out
}

compute_top_pairs_switchbox <- function(X_train, y_train,
                                        top_n = 1000,
                                        buffer_mult = 5,
                                        out_dir = NULL,
                                        save_all = FALSE) {
  trainMat <- xy_to_switchbox_mat(X_train)

  scores <- SWAP.CalculateSignedScore(
    inputMat = trainMat,
    phenoGroup = y_train,
    FilterFunc = NULL
  )

  score_mat <- scores$score
  feat_names <- rownames(score_mat)

  idx <- which(upper.tri(score_mat), arr.ind = TRUE)

  all_pairs <- data.frame(
    i = idx[, 1],
    j = idx[, 2],
    feature_i = feat_names[idx[, 1]],
    feature_j = feat_names[idx[, 2]],
    signed_score = score_mat[upper.tri(score_mat)],
    stringsAsFactors = FALSE
  )
  all_pairs$abs_score <- abs(all_pairs$signed_score)
  all_pairs <- all_pairs[order(-all_pairs$abs_score, -all_pairs$signed_score), ]

  top_buffer_n <- min(nrow(all_pairs), max(top_n, top_n * buffer_mult))
  cand <- all_pairs[seq_len(top_buffer_n), , drop = FALSE]

  ranks <- t(apply(X_train, 1, rank, ties.method = "average"))
  colnames(ranks) <- colnames(X_train)
  pos <- y_train == levels(y_train)[2]
  neg <- !pos

  orient_diff <- numeric(nrow(cand))
  gene_left <- character(nrow(cand))
  gene_right <- character(nrow(cand))

  for (r in seq_len(nrow(cand))) {
    fi <- cand$feature_i[r]
    fj <- cand$feature_j[r]

    fi_gt_fj <- ranks[, fi] < ranks[, fj]
    diff_val <- mean(fi_gt_fj[pos]) - mean(fi_gt_fj[neg])

    orient_diff[r] <- diff_val
    if (diff_val >= 0) {
      gene_left[r] <- fi
      gene_right[r] <- fj
    } else {
      gene_left[r] <- fj
      gene_right[r] <- fi
    }
  }

  cand$orientation_diff <- orient_diff
  cand$gene_left <- gene_left
  cand$gene_right <- gene_right
  cand$ordered_pair <- paste0(cand$gene_left, " > ", cand$gene_right)
  cand <- cand[order(-cand$abs_score, -abs(cand$orientation_diff)), ]

  top_pairs <- cand[seq_len(min(top_n, nrow(cand))), , drop = FALSE]

  if (!is.null(out_dir)) {
    ensure_dir(out_dir)
    if (save_all) {
      write.csv(all_pairs, file.path(out_dir, "all_signed_tsp_pairs.csv"), row.names = FALSE)
      saveRDS(scores, file.path(out_dir, "switchbox_scores.rds"))
    }
    write.csv(top_pairs, file.path(out_dir, paste0("top_", top_n, "_pairs.csv")), row.names = FALSE)
  }

  top_pairs
}

apply_disjointness <- function(pairs_df, K, d = 0) {
  if (d == 0) return(head(pairs_df, K))

  selected <- list()
  gene_count <- integer(0)

  for (r in seq_len(nrow(pairs_df))) {
    g1 <- pairs_df$gene_left[r]
    g2 <- pairs_df$gene_right[r]

    c1 <- ifelse(is.na(gene_count[g1]), 0L, gene_count[g1])
    c2 <- ifelse(is.na(gene_count[g2]), 0L, gene_count[g2])

    ok <- FALSE
    if (d == 1) ok <- (c1 < 2L) && (c2 < 2L)
    if (d == 2) ok <- (c1 < 1L) && (c2 < 1L)

    if (ok) {
      selected[[length(selected) + 1L]] <- pairs_df[r, , drop = FALSE]
      gene_count[g1] <- c1 + 1L
      gene_count[g2] <- c2 + 1L
      if (length(selected) >= K) break
    }
  }

  if (length(selected) == 0) return(pairs_df[0, , drop = FALSE])
  do.call(rbind, selected)
}

encode_pair_matrix <- function(X, pairs_df, m = 10) {
  ranks <- t(apply(X, 1, rank, ties.method = "average"))
  colnames(ranks) <- colnames(X)

  R <- matrix(0, nrow = nrow(X), ncol = nrow(pairs_df))
  rownames(R) <- rownames(X)
  colnames(R) <- pairs_df$ordered_pair

  for (k in seq_len(nrow(pairs_df))) {
    g1 <- pairs_df$gene_left[k]
    g2 <- pairs_df$gene_right[k]

    delta <- ranks[, g2] - ranks[, g1]

    R[, k] <- ifelse(delta >= m, 1,
                     ifelse(delta <= -m, -1, 0))
  }

  R
}

calc_local_activity_summary <- function(R, active_mask = NULL) {
  if (!is.null(active_mask)) {
    cols <- which(active_mask)
    if (length(cols) == 0) {
      local_n <- rep(0, nrow(R))
    } else {
      local_n <- rowSums(R[, cols, drop = FALSE] != 0)
    }
  } else {
    local_n <- rowSums(R != 0)
  }

  list(
    local_n = local_n,
    mean = mean(local_n),
    median = median(local_n),
    iqr = IQR(local_n),
    min = min(local_n),
    max = max(local_n),
    pct_le_1 = mean(local_n <= 1),
    pct_le_3 = mean(local_n <= 3),
    pct_le_5 = mean(local_n <= 5),
    pct_le_10 = mean(local_n <= 10)
  )
}

calc_hub_summary <- function(pairs_df) {
  if (nrow(pairs_df) == 0) {
    return(list(top_gene = NA, top_share = NA, counts = data.frame()))
  }

  genes <- c(pairs_df$gene_left, pairs_df$gene_right)
  tab <- sort(table(genes), decreasing = TRUE)

  counts <- data.frame(
    gene = names(tab),
    count = as.integer(tab),
    share = as.numeric(tab) / sum(tab),
    stringsAsFactors = FALSE
  )

  list(
    top_gene = counts$gene[1],
    top_share = counts$share[1],
    counts = counts
  )
}

fit_rlr_train_test <- function(X_train, y_train, X_test, y_test,
                               pairs_pool,
                               K, m, d,
                               out_dir,
                               alpha = 0.5,
                               lambda_choice = c("lambda.min", "lambda.1se"),
                               seed = 42,
                               extra_metrics = NULL) {
  lambda_choice <- match.arg(lambda_choice)
  ensure_dir(out_dir)
  set.seed(seed)

  shortlist <- apply_disjointness(pairs_pool, K = K, d = d)
  if (nrow(shortlist) == 0) stop("Brak par po shortlist refinement.")

  R_train <- encode_pair_matrix(X_train, shortlist, m = m)
  R_test  <- encode_pair_matrix(X_test, shortlist, m = m)

  y_train01 <- as.numeric(y_train == levels(y_train)[2])
  y_test01  <- as.numeric(y_test == levels(y_test)[2])

  cvfit <- cv.glmnet(
    x = R_train,
    y = y_train01,
    family = "binomial",
    alpha = alpha,
    nfolds = safe_nfolds(y_train, 5),
    type.measure = "deviance"
  )

  prob <- as.numeric(predict(cvfit, newx = R_test, s = lambda_choice, type = "response"))
  pred <- ifelse(prob >= 0.5, levels(y_test)[2], levels(y_test)[1])
  pred <- factor(pred, levels = levels(y_test))

  eps <- 1e-15
  p <- pmin(pmax(prob, eps), 1 - eps)

  auc <- as.numeric(pROC::auc(y_test01, prob))
  logloss <- -mean(y_test01 * log(p) + (1 - y_test01) * log(1 - p))
  brier <- calc_brier(prob, y_test01)
  ece <- calc_ece(prob, y_test01)

  tp <- sum(pred == levels(y_test)[2] & y_test == levels(y_test)[2])
  fp <- sum(pred == levels(y_test)[2] & y_test == levels(y_test)[1])
  fn <- sum(pred == levels(y_test)[1] & y_test == levels(y_test)[2])

  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))

  beta <- as.numeric(coef(cvfit, s = lambda_choice))[-1]
  names(beta) <- colnames(R_train)
  active_mask <- beta != 0

  active_pairs <- shortlist[colnames(R_train) %in% names(beta)[active_mask], , drop = FALSE]
  unique_genes <- if (nrow(active_pairs) == 0) 0L else length(unique(c(active_pairs$gene_left, active_pairs$gene_right)))

  local_sum <- calc_local_activity_summary(R_test, active_mask = active_mask)
  hubs <- calc_hub_summary(active_pairs)

  metrics <- list(
    method_family = "RLR",
    method = paste0("RLR_K", K, "_m", m, "_d", d),
    requested_K = K,
    m = m,
    d = d,
    alpha = alpha,
    lambda_choice = lambda_choice,
    auc = auc,
    f1 = f1,
    logloss = logloss,
    brier = brier,
    ece = ece,
    global_active_pairs = sum(active_mask),
    unique_genes = unique_genes,
    local_active_mean = local_sum$mean,
    local_active_median = local_sum$median,
    local_active_iqr = local_sum$iqr,
    hub_top_gene = hubs$top_gene,
    hub_top_share = hubs$top_share
  )

  if (!is.null(extra_metrics)) {
    metrics <- c(metrics, extra_metrics)
  }

  write.csv(
    data.frame(
      sample_id = rownames(X_test),
      y_true = as.character(y_test),
      y_pred = as.character(pred),
      score = prob,
      prob = prob
    ),
    file.path(out_dir, "predictions_test.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(
      ordered_pair = names(beta),
      beta = beta,
      abs_beta = abs(beta)
    ),
    file.path(out_dir, "coefficients.csv"),
    row.names = FALSE
  )

  write.csv(active_pairs, file.path(out_dir, "active_pairs.csv"), row.names = FALSE)

  write.csv(
    data.frame(
      sample_id = rownames(X_test),
      local_active_n = local_sum$local_n
    ),
    file.path(out_dir, "local_activity_test.csv"),
    row.names = FALSE
  )

  write_json(local_sum, file.path(out_dir, "local_activity_summary.json"), auto_unbox = TRUE, pretty = TRUE)
  write_json(metrics, file.path(out_dir, "metrics_test.json"), auto_unbox = TRUE, pretty = TRUE)

  if (nrow(hubs$counts) > 0) {
    write.csv(hubs$counts, file.path(out_dir, "hub_counts.csv"), row.names = FALSE)
  }

  saveRDS(cvfit, file.path(out_dir, "cvfit.rds"))
  invisible(metrics)
}

fit_tsp1_manual_train_test <- function(X_train, y_train, X_test, y_test,
                                       pairs_pool,
                                       out_dir,
                                       extra_metrics = NULL) {
  ensure_dir(out_dir)

  pair1 <- pairs_pool[1, , drop = FALSE]
  g1 <- pair1$gene_left[1]
  g2 <- pair1$gene_right[1]

  ranks_train <- t(apply(X_train, 1, rank, ties.method = "average"))
  colnames(ranks_train) <- colnames(X_train)
  ranks_test <- t(apply(X_test, 1, rank, ties.method = "average"))
  colnames(ranks_test) <- colnames(X_test)

  train_prob <- as.numeric(ranks_train[, g1] < ranks_train[, g2])
  y_train01 <- as.numeric(y_train == levels(y_train)[2])

  flip <- mean(train_prob[y_train01 == 1]) < mean(train_prob[y_train01 == 0])

  test_prob <- as.numeric(ranks_test[, g1] < ranks_test[, g2])
  if (flip) test_prob <- 1 - test_prob

  y_test01 <- as.numeric(y_test == levels(y_test)[2])
  pred <- ifelse(test_prob >= 0.5, levels(y_test)[2], levels(y_test)[1])
  pred <- factor(pred, levels = levels(y_test))

  eps <- 1e-15
  p <- pmin(pmax(test_prob, eps), 1 - eps)

  auc <- as.numeric(pROC::auc(y_test01, test_prob))
  logloss <- -mean(y_test01 * log(p) + (1 - y_test01) * log(1 - p))
  brier <- calc_brier(test_prob, y_test01)
  ece <- calc_ece(test_prob, y_test01)

  tp <- sum(pred == levels(y_test)[2] & y_test == levels(y_test)[2])
  fp <- sum(pred == levels(y_test)[2] & y_test == levels(y_test)[1])
  fn <- sum(pred == levels(y_test)[1] & y_test == levels(y_test)[2])

  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))

  metrics <- list(
    method_family = "baseline",
    method = "TSP_1",
    requested_k = 1,
    realized_k = 1,
    score_type = "vote_derived",
    auc = auc,
    f1 = f1,
    logloss = logloss,
    brier = brier,
    ece = ece,
    global_active_pairs = 1,
    unique_genes = 2,
    local_active_mean = 1,
    local_active_median = 1,
    local_active_iqr = 0,
    hub_top_gene = NA,
    hub_top_share = NA
  )

  if (!is.null(extra_metrics)) {
    metrics <- c(metrics, extra_metrics)
  }

  write.csv(
    data.frame(
      sample_id = rownames(X_test),
      y_true = as.character(y_test),
      y_pred = as.character(pred),
      score = test_prob,
      prob = test_prob
    ),
    file.path(out_dir, "predictions_test.csv"),
    row.names = FALSE
  )

  write.csv(pair1, file.path(out_dir, "selected_pairs.csv"), row.names = FALSE)
  write_json(metrics, file.path(out_dir, "metrics_test.json"), auto_unbox = TRUE, pretty = TRUE)
  invisible(metrics)
}

fit_ktsp_train_test <- function(X_train, y_train, X_test, y_test,
                                k,
                                out_dir,
                                seed = 42,
                                extra_metrics = NULL) {
  ensure_dir(out_dir)
  set.seed(seed)

  trainMat <- xy_to_switchbox_mat(X_train)
  testMat  <- xy_to_switchbox_mat(X_test)

  clf <- SWAP.KTSP.Train(
    inputMat = trainMat,
    phenoGroup = y_train,
    FilterFunc = NULL,
    krange = k
  )

  if (!is.matrix(clf$TSPs)) {
    clf$TSPs <- matrix(clf$TSPs, ncol = 2, byrow = TRUE)
  }

  pred <- SWAP.KTSP.Classify(testMat, clf)
  stat_train <- SWAP.KTSP.Statistics(trainMat, clf)
  stat_test  <- SWAP.KTSP.Statistics(testMat, clf)

  realized_k <- nrow(clf$TSPs)

  prob_train <- (as.numeric(stat_train$statistics) + realized_k) / (2 * realized_k)
  prob_test  <- (as.numeric(stat_test$statistics) + realized_k) / (2 * realized_k)

  y_train01 <- as.numeric(y_train == levels(y_train)[2])
  y_test01  <- as.numeric(y_test == levels(y_test)[2])

  flip <- mean(prob_train[y_train01 == 1]) < mean(prob_train[y_train01 == 0])
  if (flip) prob_test <- 1 - prob_test

  pred <- factor(pred, levels = levels(y_test))
  eps <- 1e-15
  p <- pmin(pmax(prob_test, eps), 1 - eps)

  auc <- as.numeric(pROC::auc(y_test01, prob_test))
  logloss <- -mean(y_test01 * log(p) + (1 - y_test01) * log(1 - p))
  brier <- calc_brier(prob_test, y_test01)
  ece <- calc_ece(prob_test, y_test01)

  tp <- sum(pred == levels(y_test)[2] & y_test == levels(y_test)[2])
  fp <- sum(pred == levels(y_test)[2] & y_test == levels(y_test)[1])
  fn <- sum(pred == levels(y_test)[1] & y_test == levels(y_test)[2])

  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))

  selected_pairs <- as.data.frame(clf$TSPs, stringsAsFactors = FALSE)
  colnames(selected_pairs) <- c("gene_left", "gene_right")
  selected_pairs$ordered_pair <- paste0(selected_pairs$gene_left, " > ", selected_pairs$gene_right)

  unique_genes <- length(unique(c(selected_pairs$gene_left, selected_pairs$gene_right)))
  hubs <- calc_hub_summary(selected_pairs)

  local_n <- rep(realized_k, nrow(X_test))
  local_sum <- list(
    local_n = local_n,
    mean = mean(local_n),
    median = median(local_n),
    iqr = IQR(local_n),
    min = min(local_n),
    max = max(local_n),
    pct_le_1 = mean(local_n <= 1),
    pct_le_3 = mean(local_n <= 3),
    pct_le_5 = mean(local_n <= 5),
    pct_le_10 = mean(local_n <= 10)
  )

  metrics <- list(
    method_family = "baseline",
    method = paste0("kTSP_", k),
    requested_k = k,
    realized_k = realized_k,
    score_type = "vote_derived",
    auc = auc,
    f1 = f1,
    logloss = logloss,
    brier = brier,
    ece = ece,
    global_active_pairs = realized_k,
    unique_genes = unique_genes,
    local_active_mean = local_sum$mean,
    local_active_median = local_sum$median,
    local_active_iqr = local_sum$iqr,
    hub_top_gene = hubs$top_gene,
    hub_top_share = hubs$top_share
  )

  if (!is.null(extra_metrics)) {
    metrics <- c(metrics, extra_metrics)
  }

  write.csv(
    data.frame(
      sample_id = rownames(X_test),
      y_true = as.character(y_test),
      y_pred = as.character(pred),
      score = prob_test,
      prob = prob_test
    ),
    file.path(out_dir, "predictions_test.csv"),
    row.names = FALSE
  )

  write.csv(selected_pairs, file.path(out_dir, "selected_pairs.csv"), row.names = FALSE)

  write.csv(
    data.frame(
      sample_id = rownames(X_test),
      local_active_n = local_n
    ),
    file.path(out_dir, "local_activity_test.csv"),
    row.names = FALSE
  )

  write_json(local_sum, file.path(out_dir, "local_activity_summary.json"), auto_unbox = TRUE, pretty = TRUE)
  write_json(metrics, file.path(out_dir, "metrics_test.json"), auto_unbox = TRUE, pretty = TRUE)

  if (nrow(hubs$counts) > 0) {
    write.csv(hubs$counts, file.path(out_dir, "hub_counts.csv"), row.names = FALSE)
  }

  saveRDS(clf, file.path(out_dir, "classifier.rds"))
  invisible(metrics)
}
fit_value_lr_train_test <- function(X_train, y_train, X_test, y_test,
                                    out_dir,
                                    pre_top_features = 1000,
                                    alpha = 0.5,
                                    lambda_choice = "lambda.min",
                                    seed = 42,
                                    extra_metrics = NULL) {
  ensure_dir(out_dir)
  set.seed(seed)

  pf <- prefilter_by_mad(
    X_train = X_train,
    X_test = X_test,
    top_n = pre_top_features
  )

  X_train_f <- pf$X_train
  X_test_f  <- pf$X_test

  x_train <- as.matrix(X_train_f)
  x_test  <- as.matrix(X_test_f)

  storage.mode(x_train) <- "numeric"
  storage.mode(x_test)  <- "numeric"

  y_train <- factor(y_train)
  y_test  <- factor(y_test, levels = levels(y_train))

  y_train01 <- as.numeric(y_train == levels(y_train)[2])
  y_test01  <- as.numeric(y_test  == levels(y_train)[2])

  cvfit <- cv.glmnet(
    x = x_train,
    y = y_train01,
    family = "binomial",
    alpha = alpha,
    nfolds = safe_nfolds(y_train, 5),
    type.measure = "deviance"
  )

  prob <- as.numeric(predict(cvfit, newx = x_test, s = lambda_choice, type = "response"))
  pred <- ifelse(prob >= 0.5, levels(y_train)[2], levels(y_train)[1])
  pred <- factor(pred, levels = levels(y_train))

  eps <- 1e-15
  p <- pmin(pmax(prob, eps), 1 - eps)

  auc <- if (length(unique(y_test01)) < 2) NA_real_ else as.numeric(pROC::auc(y_test01, prob))
  logloss <- -mean(y_test01 * log(p) + (1 - y_test01) * log(1 - p))
  brier <- calc_brier(prob, y_test01)
  ece <- calc_ece(prob, y_test01)

  tp <- sum(pred == levels(y_train)[2] & y_test == levels(y_train)[2])
  fp <- sum(pred == levels(y_train)[2] & y_test == levels(y_train)[1])
  fn <- sum(pred == levels(y_train)[1] & y_test == levels(y_train)[2])

  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))

  beta <- as.numeric(coef(cvfit, s = lambda_choice))[-1]
  names(beta) <- colnames(x_train)
  active_mask <- beta != 0

  global_active_features <- sum(active_mask)
  local_n <- rep(global_active_features, nrow(x_test))

  metrics <- list(
    method_family = "value_baseline",
    method = "value_LR",
    pre_top_features = pre_top_features,
    alpha = alpha,
    lambda_choice = lambda_choice,
    auc = auc,
    f1 = f1,
    logloss = logloss,
    brier = brier,
    ece = ece,
    global_active_pairs = global_active_features,
    unique_genes = global_active_features,
    local_active_mean = mean(local_n),
    local_active_median = median(local_n),
    local_active_iqr = IQR(local_n),
    hub_top_gene = "",
    hub_top_share = NA_real_
  )

  if (!is.null(extra_metrics)) {
    metrics <- c(metrics, extra_metrics)
  }

  write.csv(
    data.frame(
      sample_id = rownames(X_test),
      y_true = as.character(y_test),
      y_pred = as.character(pred),
      score = prob,
      prob = prob
    ),
    file.path(out_dir, "predictions_test.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(
      feature = names(beta),
      beta = beta,
      abs_beta = abs(beta)
    ),
    file.path(out_dir, "coefficients.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(
      sample_id = rownames(X_test),
      local_active_n = local_n
    ),
    file.path(out_dir, "local_activity_test.csv"),
    row.names = FALSE
  )

  jsonlite::write_json(
    metrics,
    file.path(out_dir, "metrics_test.json"),
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  invisible(metrics)
}
