# RLMT binary prototype: isolated from the frozen RLR implementation.
# All pair search is confined to the outer-training global pool.
rlmt_defaults <- function(config = list()) {
  x <- modifyList(list(K = 10L, m = 10L, d = 2L, pool = 1000L,
                       max_depth = 2L, screen_top_B = 20L,
                       min_child_n = 16L, min_child_class = 6L,
                       alpha = 0.5, lambda_choice = "lambda.min",
                       inner_folds = 3L, seed = 42L), config)
  for (k in c("K", "m", "d", "pool", "max_depth", "screen_top_B",
              "min_child_n", "min_child_class", "inner_folds", "seed")) {
    x[[k]] <- as.integer(x[[k]])
  }
  if (x$K != 10L || x$m != 10L || x$d != 2L || x$pool != 1000L ||
      x$alpha != 0.5 || x$lambda_choice != "lambda.min") {
    stop("First RLMT gate is restricted to K10/m10/d2/pool1000/alpha0.5/lambda.min")
  }
  if (x$max_depth < 2L || x$inner_folds < 3L ||
      x$min_child_n < 6L || x$min_child_class < 3L ||
      x$screen_top_B < 1L) stop("Invalid RLMT safety parameters")
  x
}

rlmt_y01 <- function(y) as.numeric(y == levels(y)[2])

rlmt_logloss <- function(y01, prob) {
  p <- pmin(pmax(prob, 1e-15), 1 - 1e-15)
  -mean(y01 * log(p) + (1 - y01) * log(1 - p))
}

rlmt_children_ok <- function(z, y, conf) {
  if (anyNA(z)) stop("Missing values in relational routing")
  sides <- list(which(z), which(!z))
  all(vapply(sides, function(idx) {
    length(idx) >= conf$min_child_n &&
      all(table(y[idx]) >= conf$min_child_class)
  }, logical(1)))
}

# Cheap node-local empirical TSP contrast over already selected global pairs.
# No SWAP.CalculateSignedScore and no p^2 matrix in child nodes.
rerank_pool_in_node <- function(X_node, y_node, global_pool) {
  stopifnot(length(levels(y_node)) == 2L, nrow(global_pool) > 0L)
  if (!all(c("feature_i", "feature_j", "gene_left", "gene_right") %in%
           names(global_pool))) stop("Unexpected global pool schema")
  positive <- y_node == levels(y_node)[2]
  negative <- !positive
  if (!any(positive) || !any(negative)) stop("Node contains one class")
  out <- global_pool
  contrasts <- numeric(nrow(out))
  for (i in seq_len(nrow(out))) {
    a <- out$feature_i[i]; b <- out$feature_j[i]
    z <- X_node[[a]] > X_node[[b]]
    if (anyNA(z)) stop("NA in pair comparison")
    contrasts[i] <- mean(z[positive]) - mean(z[negative])
    if (contrasts[i] >= 0) {
      out$gene_left[i] <- a; out$gene_right[i] <- b
    } else {
      out$gene_left[i] <- b; out$gene_right[i] <- a
    }
  }
  out$orientation_diff <- contrasts
  out$signed_score <- contrasts
  out$abs_score <- abs(contrasts)
  out$ordered_pair <- paste0(out$gene_left, " > ", out$gene_right)
  out[order(-out$abs_score, seq_len(nrow(out))), , drop = FALSE]
}

rlmt_fit_leaf <- function(X, y, pool, conf, rerank = TRUE, seed = conf$seed) {
  if (length(levels(y)) != 2L || any(table(y) < 3L))
    stop("Insufficient leaf class support for glmnet CV")
  pairs <- if (rerank) rerank_pool_in_node(X, y, pool) else pool
  pairs <- apply_disjointness(pairs, K = conf$K, d = conf$d)
  if (nrow(pairs) == 0L) stop("Empty RLR shortlist in leaf")
  R <- encode_pair_matrix(X, pairs, m = conf$m)
  set.seed(seed)
  fit <- glmnet::cv.glmnet(
    x = R, y = rlmt_y01(y), family = "binomial",
    alpha = conf$alpha, nfolds = safe_nfolds(y, 5),
    type.measure = "deviance"
  )
  beta <- as.numeric(coef(fit, s = conf$lambda_choice))[-1]
  names(beta) <- colnames(R)
  list(fit = fit, pairs = pairs, beta = beta,
       active = pairs[beta != 0, , drop = FALSE])
}

rlmt_predict_leaf <- function(model, X, conf) {
  R <- encode_pair_matrix(X, model$pairs, m = conf$m)
  p <- as.numeric(predict(model$fit, newx = R,
                          s = conf$lambda_choice, type = "response"))
  list(prob = p, local_active_n =
         as.numeric(rowSums(R[, model$beta != 0, drop = FALSE] != 0)))
}

rlmt_parent_oof <- function(X, y, pool, conf, folds) {
  prob <- rep(NA_real_, nrow(X))
  for (j in seq_along(folds)) {
    f <- folds[[j]]
    fit <- rlmt_fit_leaf(X[f$train_idx, , drop = FALSE],
                         y[f$train_idx], pool, conf, rerank = FALSE,
                         seed = conf$seed + j)
    prob[f$test_idx] <- rlmt_predict_leaf(
      fit, X[f$test_idx, , drop = FALSE], conf)$prob
  }
  if (anyNA(prob)) stop("Parent OOF probabilities are incomplete")
  prob
}

rlmt_screen <- function(X, y, pool, residuals, conf) {
  out <- list()
  for (i in seq_len(nrow(pool))) {
    a <- pool$gene_left[i]; b <- pool$gene_right[i]
    z <- X[[a]] > X[[b]]
    if (!rlmt_children_ok(z, y, conf)) next
    nl <- sum(z); nr <- length(z) - nl
    score <- abs(mean(residuals[z]) - mean(residuals[!z])) *
      sqrt(nl * nr / length(z))
    out[[length(out) + 1L]] <- data.frame(
      pool_index = i, gene_a = a, gene_b = b,
      n_left = nl, n_right = nr, residual_screen_score = score)
  }
  if (!length(out)) return(data.frame())
  x <- do.call(rbind, out)
  head(x[order(-x$residual_screen_score, x$pool_index), , drop = FALSE],
       conf$screen_top_B)
}

rlmt_lookahead <- function(X, y, pool, parent_oof, folds, screen, conf) {
  if (!nrow(screen)) return(list(winner = NULL, audit = data.frame()))
  rows <- vector("list", nrow(screen))
  for (i in seq_len(nrow(screen))) {
    s <- screen[i, , drop = FALSE]
    z <- X[[s$gene_a]] > X[[s$gene_b]]
    child_oof <- rep(NA_real_, nrow(X))
    fold_gains <- numeric(length(folds))
    valid <- TRUE
    error_message <- ""
    for (j in seq_along(folds)) {
      f <- folds[[j]]
      tr <- f$train_idx; va <- f$test_idx
      sides <- list(tr[z[tr]], tr[!z[tr]])
      if (any(vapply(sides, function(idx) {
        length(idx) < 6L || any(table(y[idx]) < 3L)
      }, logical(1)))) {valid <- FALSE; error_message <- "inner_child_support"; break}
      result <- tryCatch({
        left <- rlmt_fit_leaf(X[sides[[1]], , drop = FALSE], y[sides[[1]]],
                              pool, conf, seed = conf$seed + j * 100L + 1L)
        right <- rlmt_fit_leaf(X[sides[[2]], , drop = FALSE], y[sides[[2]]],
                               pool, conf, seed = conf$seed + j * 100L + 2L)
        out <- numeric(length(va))
        on_left <- z[va]
        if (any(on_left)) out[on_left] <- rlmt_predict_leaf(
          left, X[va[on_left], , drop = FALSE], conf)$prob
        if (any(!on_left)) out[!on_left] <- rlmt_predict_leaf(
          right, X[va[!on_left], , drop = FALSE], conf)$prob
        out
      }, error = function(e) e)
      if (inherits(result, "error")) {
        valid <- FALSE; error_message <- conditionMessage(result); break
      }
      child_oof[va] <- result
      fold_gains[j] <- rlmt_logloss(rlmt_y01(y[va]), parent_oof[va]) -
        rlmt_logloss(rlmt_y01(y[va]), result)
    }
    mean_gain <- if (valid) mean(fold_gains) else NA_real_
    se_gain <- if (valid) stats::sd(fold_gains) / sqrt(length(folds)) else NA_real_
    rows[[i]] <- data.frame(
      pool_index = s$pool_index, gene_a = s$gene_a, gene_b = s$gene_b,
      residual_screen_score = s$residual_screen_score,
      parent_oof_logloss = rlmt_logloss(rlmt_y01(y), parent_oof),
      split_oof_logloss = if (valid) rlmt_logloss(rlmt_y01(y), child_oof) else NA_real_,
      mean_gain = mean_gain, se_gain = se_gain,
      lower_gain = mean_gain - se_gain, valid = valid,
      rejection_reason = if (valid) "" else error_message,
      stringsAsFactors = FALSE
    )
  }
  audit <- do.call(rbind, rows)
  admissible <- which(audit$valid & is.finite(audit$lower_gain) &
                        audit$mean_gain > 0 & audit$lower_gain > 0)
  if (!length(admissible)) return(list(winner = NULL, audit = audit))
  best <- admissible[order(-audit$lower_gain[admissible],
                           -audit$mean_gain[admissible],
                           audit$pool_index[admissible])][1]
  list(winner = audit[best, , drop = FALSE], audit = audit)
}

rlmt_node_row <- function(id, parent, depth, is_leaf, y, split = NULL,
                          active_count = 0L) {
  data.frame(
    node_id = id, parent_id = parent, depth = depth, is_leaf = is_leaf,
    n = length(y), class_0 = as.integer(table(y)[1]),
    class_1 = as.integer(table(y)[2]),
    split_gene_a = if (is.null(split)) NA_character_ else split$gene_a,
    split_gene_b = if (is.null(split)) NA_character_ else split$gene_b,
    residual_screen_score = if (is.null(split)) NA_real_ else split$residual_screen_score,
    parent_oof_logloss = if (is.null(split)) NA_real_ else split$parent_oof_logloss,
    split_oof_logloss = if (is.null(split)) NA_real_ else split$split_oof_logloss,
    mean_gain = if (is.null(split)) NA_real_ else split$mean_gain,
    se_gain = if (is.null(split)) NA_real_ else split$se_gain,
    left_child = if (is_leaf) NA_integer_ else 2L,
    right_child = if (is_leaf) NA_integer_ else 3L,
    leaf_globally_active_relation_count = if (is_leaf) active_count else NA_integer_
  )
}

rlmt_save_leaf <- function(model, path) {
  ensure_dir(path)
  write.csv(data.frame(ordered_pair = names(model$beta),
                       beta = model$beta, abs_beta = abs(model$beta)),
            file.path(path, "coefficients.csv"), row.names = FALSE)
  write.csv(model$active, file.path(path, "active_pairs.csv"), row.names = FALSE)
  saveRDS(model$fit, file.path(path, "cvfit.rds"))
}

rlmt_root_output <- function(X_train, y_train, X_test, y_test,
                             pool, conf, out_dir, method, extra_metrics, reason) {
  # Exact existing RLR call: parity includes RNG and all unchanged transformations.
  fit_rlr_train_test(X_train, y_train, X_test, y_test, pool,
                     K = conf$K, m = conf$m, d = conf$d, out_dir = out_dir,
                     alpha = conf$alpha, lambda_choice = conf$lambda_choice,
                     seed = conf$seed, extra_metrics = extra_metrics)
  metrics_path <- file.path(out_dir, "metrics_test.json")
  met <- jsonlite::read_json(metrics_path, simplifyVector = TRUE)
  met$method_family <- "RLMT"
  met$method <- method
  met$max_leaves <- as.integer(sub(".*_L", "", method))
  met$realized_leaves <- 1L
  jsonlite::write_json(met, metrics_path, auto_unbox = TRUE, pretty = TRUE)
  leaf_dir <- file.path(out_dir, "leaf_models", "node_1")
  ensure_dir(leaf_dir)
  for (fn in c("coefficients.csv", "active_pairs.csv", "cvfit.rds")) {
    file.copy(file.path(out_dir, fn), file.path(leaf_dir, fn), overwrite = TRUE)
  }
  pred <- read.csv(file.path(out_dir, "predictions_test.csv"))
  nodes <- rlmt_node_row(1L, NA_integer_, 0L, TRUE, y_train,
                         active_count = as.integer(met$global_active_pairs))
  write.csv(nodes, file.path(out_dir, "nodes.csv"), row.names = FALSE)
  write.csv(data.frame(sample_id = pred$sample_id, leaf_id = 1L),
            file.path(out_dir, "test_routes.csv"), row.names = FALSE)
  jsonlite::write_json(list(size = 1L, max_leaves = met$max_leaves,
                            root_id = 1L, accepted_split = FALSE,
                            stop_reason = reason,
                            nodes = list(list(node_id = 1L, is_leaf = TRUE))),
                       file.path(out_dir, "tree.json"), auto_unbox = TRUE, pretty = TRUE)
}

rlmt_run <- function(X_train, y_train, X_test, y_test,
                     global_pool, max_leaves, config, out_dir,
                     extra_metrics = list()) {
  conf <- rlmt_defaults(config)
  if (!max_leaves %in% c(1L, 2L)) stop("Prototype supports max_leaves 1 or 2 only")
  if (!identical(levels(y_train), levels(y_test)) ||
      length(levels(y_train)) != 2L) stop("Binary class levels must match")
  if (nrow(global_pool) < conf$K) stop("Global pair pool smaller than K")
  if (anyNA(X_train) || anyNA(X_test)) stop("RLMT does not support missing predictors")
  ensure_dir(out_dir)
  method <- sprintf("RLMT_K%d_m%d_d%d_pool%d_L%d",
                    conf$K, conf$m, conf$d, conf$pool, max_leaves)
  if (max_leaves == 1L) {
    rlmt_root_output(X_train, y_train, X_test, y_test, global_pool,
                     conf, out_dir, method, extra_metrics, "max_leaves=1")
    return(invisible(TRUE))
  }

  decision <- NULL
  stop_reason <- "no_admissible_split"
  if (nrow(X_train) >= 2L * conf$min_child_n &&
      all(table(y_train) >= 2L * conf$min_child_class)) {
    folds <- make_repeated_stratified_folds(
      y_train, n_folds = conf$inner_folds, n_repeats = 1L, seed = conf$seed)
    parent <- rlmt_parent_oof(X_train, y_train, global_pool, conf, folds)
    screen <- rlmt_screen(X_train, y_train, global_pool,
                          rlmt_y01(y_train) - parent, conf)
    if (nrow(screen)) write.csv(screen, file.path(out_dir, "screen_top.csv"),
                                 row.names = FALSE)
    assessed <- rlmt_lookahead(X_train, y_train, global_pool,
                                parent, folds, screen, conf)
    if (nrow(assessed$audit)) write.csv(assessed$audit,
                                       file.path(out_dir, "lookahead.csv"),
                                       row.names = FALSE)
    decision <- assessed$winner
    if (nrow(screen) == 0L) {
      stop_reason <- "no_valid_child_splits"
    } else if (is.null(decision)) {
      stop_reason <- "no_positive_1se_gain"
    }
  } else {
    stop_reason <- "insufficient_root_support"
  }

  if (is.null(decision)) {
    rlmt_root_output(X_train, y_train, X_test, y_test, global_pool,
                     conf, out_dir, method, extra_metrics, stop_reason)
    return(invisible(TRUE))
  }
  z_train <- X_train[[decision$gene_a]] > X_train[[decision$gene_b]]
  z_test <- X_test[[decision$gene_a]] > X_test[[decision$gene_b]]
  idx <- list(which(z_train), which(!z_train))
  test_idx <- list(which(z_test), which(!z_test))
  # Fit on all outer-training observations in each accepted child.
  local_models <- lapply(seq_along(idx), function(i) {
    rlmt_fit_leaf(X_train[idx[[i]], , drop = FALSE], y_train[idx[[i]]],
                  global_pool, conf, seed = conf$seed + i)
  })
  prob <- rep(NA_real_, nrow(X_test))
  local_n <- numeric(nrow(X_test))
  for (i in seq_along(test_idx)) {
    ii <- test_idx[[i]]
    if (length(ii)) {
      pp <- rlmt_predict_leaf(local_models[[i]], X_test[ii, , drop = FALSE], conf)
      prob[ii] <- pp$prob
      local_n[ii] <- pp$local_active_n
    }
    rlmt_save_leaf(local_models[[i]], file.path(out_dir, "leaf_models",
                                                paste0("node_", i + 1L)))
  }
  if (anyNA(prob) || any(!is.finite(prob)) ||
      any(prob < 0 | prob > 1)) stop("Invalid RLMT test probabilities")
  pred <- factor(ifelse(prob >= 0.5, levels(y_test)[2], levels(y_test)[1]),
                 levels = levels(y_test))
  y01 <- rlmt_y01(y_test)
  tp <- sum(pred == levels(y_test)[2] & y_test == levels(y_test)[2])
  fp <- sum(pred == levels(y_test)[2] & y_test == levels(y_test)[1])
  fn <- sum(pred == levels(y_test)[1] & y_test == levels(y_test)[2])
  precision <- if (tp + fp == 0) 0 else tp / (tp + fp)
  recall <- if (tp + fn == 0) 0 else tp / (tp + fn)
  f1 <- if (precision + recall == 0) 0 else
    2 * precision * recall / (precision + recall)
  all_active <- unique(unlist(lapply(local_models, function(x)
    x$active$ordered_pair), use.names = FALSE))
  all_genes <- unique(unlist(lapply(local_models, function(x)
    c(x$active$gene_left, x$active$gene_right)), use.names = FALSE))
  metrics <- c(list(
    method_family = "RLMT", method = method, requested_K = conf$K,
    m = conf$m, d = conf$d, alpha = conf$alpha,
    lambda_choice = conf$lambda_choice, max_leaves = max_leaves,
    realized_leaves = 2L,
    auc = as.numeric(pROC::auc(y01, prob)), f1 = f1,
    logloss = rlmt_logloss(y01, prob),
    brier = calc_brier(prob, y01), ece = calc_ece(prob, y01),
    global_active_pairs = length(all_active),
    unique_genes = length(all_genes),
    local_active_mean = mean(local_n),
    local_active_median = median(local_n),
    local_active_iqr = IQR(local_n),
    pool = conf$pool
  ), extra_metrics)
  write.csv(data.frame(sample_id = rownames(X_test),
                       y_true = as.character(y_test), y_pred = as.character(pred),
                       score = prob, prob = prob),
            file.path(out_dir, "predictions_test.csv"), row.names = FALSE)
  write.csv(data.frame(sample_id = rownames(X_test),
                       local_active_n = local_n),
            file.path(out_dir, "local_activity_test.csv"), row.names = FALSE)
  jsonlite::write_json(metrics, file.path(out_dir, "metrics_test.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  nodes <- do.call(rbind, list(
    rlmt_node_row(1L, NA_integer_, 0L, FALSE, y_train, decision),
    rlmt_node_row(2L, 1L, 1L, TRUE, y_train[idx[[1]]],
                  active_count = nrow(local_models[[1]]$active)),
    rlmt_node_row(3L, 1L, 1L, TRUE, y_train[idx[[2]]],
                  active_count = nrow(local_models[[2]]$active))
  ))
  write.csv(nodes, file.path(out_dir, "nodes.csv"), row.names = FALSE)
  write.csv(data.frame(sample_id = rownames(X_test),
                       leaf_id = ifelse(z_test, 2L, 3L)),
            file.path(out_dir, "test_routes.csv"), row.names = FALSE)
  jsonlite::write_json(list(
    size = 2L, max_leaves = 2L, root_id = 1L, accepted_split = TRUE,
    split = list(gene_a = decision$gene_a, gene_b = decision$gene_b,
                 residual_screen_score = decision$residual_screen_score,
                 parent_oof_logloss = decision$parent_oof_logloss,
                 split_oof_logloss = decision$split_oof_logloss,
                 mean_gain = decision$mean_gain, se_gain = decision$se_gain,
                 left_child = 2L, right_child = 3L),
    nodes = list(list(node_id = 1L, is_leaf = FALSE),
                 list(node_id = 2L, is_leaf = TRUE),
                 list(node_id = 3L, is_leaf = TRUE))
  ), file.path(out_dir, "tree.json"), auto_unbox = TRUE, pretty = TRUE)
  invisible(TRUE)
}
