# RLMT-DC v0: discovery/confirmation components.
# The published RLR and the frozen RLMT-v1 functions are NEVER modified here.
# Source 07_world_protocol_utils.R and 28_rlmt_utils.R before this file.

rlmt_dc_defaults <- function(config = list()) {
  conf <- rlmt_defaults(config)
  x <- modifyList(list(proposal_n = 4096L, confirm_folds = 3L,
                       discovery_oof_folds = 3L,
                       confirmation_min_n = 20L), config$dc %||% list())
  x$proposal_n <- as.integer(x$proposal_n)
  x$confirm_folds <- as.integer(x$confirm_folds)
  x$discovery_oof_folds <- as.integer(x$discovery_oof_folds)
  x$confirmation_min_n <- as.integer(x$confirmation_min_n)
  if (x$proposal_n != 4096L || x$confirm_folds != 3L ||
      x$discovery_oof_folds != 3L ||
      x$confirmation_min_n != 20L ||
      conf$min_child_n != 16L || conf$min_child_class != 6L ||
      conf$inner_folds != 3L || conf$seed != 42L) {
    stop("RLMT-DC-v0 configuration must match its frozen research contract")
  }
  list(leaf = conf, dc = x)
}

rlmt_dc_partition <- function(y, seed = 42L) {
  if (!is.factor(y) || length(levels(y)) != 2L || anyNA(y) ||
      min(table(y)) < 3L) stop("Invalid class support for DC partition")
  folds <- make_repeated_stratified_folds(
    y, n_folds = 3L, n_repeats = 1L, seed = seed)
  if (length(folds) != 3L) stop("DC partition construction failed")
  d <- as.integer(folds[[1L]]$train_idx)
  v <- as.integer(folds[[1L]]$test_idx)
  if (!length(d) || !length(v) || length(intersect(d,v)) ||
      !identical(sort(c(d,v)),seq_along(y))) {
    stop("Discovery/confirmation are not a partition")
  }
  list(discovery = d, confirmation = v)
}

# Label-agnostic and uniformly random over unordered feature-index pairs.
# Bounded O(n_pairs) memory; no p^2 candidate search / pair matrix.
# Candidate order is deterministic for a fixed feature-name vector and seed.
rlmt_dc_random_pairs <- function(feature_names, n_pairs = 4096L,
                                 seed = 42L) {
  if (!is.character(feature_names) || anyNA(feature_names) ||
      anyDuplicated(feature_names) || length(feature_names) < 2L ||
      any(!nzchar(feature_names))) stop("Invalid candidate feature names")
  n_pairs <- as.integer(n_pairs)
  if (length(n_pairs) != 1L || is.na(n_pairs) || n_pairs < 1L)
    stop("Invalid n_pairs")
  p <- length(feature_names)
  total <- as.double(p) * (p - 1) / 2
  target <- as.integer(min(n_pairs,total))
  # For tiny synthetic feature spaces, exhaust without rejection sampling.
  if (total <= target && p <= 128L) {
    idx <- t(utils::combn(seq_len(p),2L))
    return(data.frame(
      pool_index = seq_len(nrow(idx)),
      gene_a = feature_names[idx[,1L]], gene_b = feature_names[idx[,2L]],
      stringsAsFactors = FALSE))
  }

  set.seed(seed)
  seen <- new.env(hash = TRUE, parent = emptyenv())
  left <- integer(target)
  right <- integer(target)
  found <- 0L
  attempts <- 0L
  while (found < target) {
    attempts <- attempts + 1L
    if (attempts > 1000L) stop("Could not draw unique unordered pairs")
    need <- target-found
    batch <- as.integer(min(100000, max(64, 4 * need)))
    a <- sample.int(p,batch,replace=TRUE)
    b <- sample.int(p,batch,replace=TRUE)
    lo <- pmin(a,b)
    hi <- pmax(a,b)
    candidates <- which(lo != hi)
    for (j in candidates) {
      key <- paste0(lo[j], "_", hi[j])
      if (!exists(key, envir = seen, inherits = FALSE)) {
        assign(key, TRUE, envir = seen)
        found <- found+1L
        left[found] <- lo[j]
        right[found] <- hi[j]
        if (found >= target) break
      }
    }
  }
  data.frame(pool_index=seq_len(target),
             gene_a=feature_names[left],gene_b=feature_names[right],
             stringsAsFactors=FALSE)
}

# Training-only candidate feasibility and residual-heterogeneity score.
# Candidate orientations are fixed before observing any validation labels.
rlmt_dc_propose <- function(X_discovery,y_discovery,prob_oof,
                            candidate_pairs, conf) {
  if (!is.data.frame(X_discovery) ||
      !is.factor(y_discovery) ||
      length(levels(y_discovery)) != 2L ||
      nrow(X_discovery) != length(y_discovery) ||
      length(prob_oof) != nrow(X_discovery) ||
      anyNA(prob_oof) || any(!is.finite(prob_oof)) ||
      any(prob_oof < 0 | prob_oof > 1)) {
    stop("Invalid discovery data or OOF probabilities")
  }
  required <- c("pool_index","gene_a","gene_b")
  if (!all(required %in% names(candidate_pairs))) {
    stop("Candidate schema mismatch")
  }
  if (!all(candidate_pairs$gene_a %in% names(X_discovery)) ||
      !all(candidate_pairs$gene_b %in% names(X_discovery))) {
    stop("Candidate genes absent from discovery")
  }

  residual <- rlmt_y01(y_discovery) - prob_oof
  out <- vector("list",nrow(candidate_pairs))
  for (j in seq_len(nrow(candidate_pairs))) {
    a <- candidate_pairs$gene_a[[j]]
    b <- candidate_pairs$gene_b[[j]]
    z <- X_discovery[[a]] > X_discovery[[b]]
    if (anyNA(z)) stop("Missing routing value in discovery")
    side <- list(which(z),which(!z))
    size <- lengths(side)
    tab <- rbind(as.integer(table(y_discovery[side[[1L]]])),
                 as.integer(table(y_discovery[side[[2L]]])))
    feasible <- all(size >= conf$min_child_n) &&
                all(tab >= conf$min_child_class)
    score <- if (feasible) {
      abs(mean(residual[z]) - mean(residual[!z])) *
        sqrt(size[1L] * size[2L] / length(z))
    } else NA_real_
    out[[j]] <- data.frame(
      pool_index = as.integer(candidate_pairs$pool_index[j]),
      gene_a = a, gene_b = b, feasible = feasible,
      n_left=size[1L], n_right=size[2L],
      left_class0=tab[1L,1L],left_class1=tab[1L,2L],
      right_class0=tab[2L,1L],right_class1=tab[2L,2L],
      residual_screen_score = score,
      stringsAsFactors=FALSE)
  }
  audit <- do.call(rbind,out)
  eligible <- which(audit$feasible & is.finite(audit$residual_screen_score))
  nominated <- if (!length(eligible)) NULL else {
    winner <- eligible[order(-audit$residual_screen_score[eligible],
                              audit$pool_index[eligible])][1L]
    audit[winner,,drop=FALSE]
  }
  list(candidate=nominated,audit=audit,feasible_n=length(eligible))
}

rlmt_dc_confirm <- function(y_confirm,parent_prob,split_prob,
                            min_n = 20L) {
  n <- length(y_confirm)
  if (!is.factor(y_confirm) || length(levels(y_confirm)) != 2L ||
      anyNA(y_confirm) || n != length(parent_prob) ||
      n != length(split_prob) || n < 2L ||
      anyNA(parent_prob) || anyNA(split_prob) ||
      any(!is.finite(parent_prob)) || any(!is.finite(split_prob)) ||
      any(parent_prob < 0 | parent_prob > 1) ||
      any(split_prob < 0 | split_prob > 1)) {
    stop("Invalid independent confirmation observations/predictions")
  }
  y01 <- rlmt_y01(y_confirm)
  per_loss <- function(p) {
    z <- pmin(pmax(p,1e-15),1-1e-15)
    -(y01*log(z)+(1-y01)*log1p(-z))
  }
  diff <- per_loss(parent_prob) - per_loss(split_prob)
  gain <- mean(diff)
  se <- stats::sd(diff)/sqrt(n)
  lower <- gain - se
  list(
    accepted = isTRUE(n >= min_n && gain > 0 && lower > 0),
    gain = gain, se_gain = se, lower_gain = lower,
    n = n, parent_logloss = mean(per_loss(parent_prob)),
    split_logloss = mean(per_loss(split_prob)))
}

rlmt_dc_parent_and_child_on_validation <- function(
    X_discovery,y_discovery,X_confirmation,pair,discovery_pool,conf) {
  if (!is.data.frame(pair) || nrow(pair)!=1L ||
      !isTRUE(pair$feasible[[1L]])) {
    stop("Confirmation can only evaluate one nominated feasible candidate")
  }
  z <- X_discovery[[pair$gene_a[[1L]]]] >
       X_discovery[[pair$gene_b[[1L]]]]
  child <- list(which(z),which(!z))
  if (!rlmt_children_ok(z,y_discovery,conf))
    stop("Nominated discovery split violates local model support")
  parent <- rlmt_fit_leaf(X_discovery,y_discovery,discovery_pool,
                          conf,rerank=FALSE,seed=conf$seed+701L)
  left <- rlmt_fit_leaf(X_discovery[child[[1L]],,drop=FALSE],
                        y_discovery[child[[1L]]],
                        discovery_pool,conf,seed=conf$seed+702L)
  right <- rlmt_fit_leaf(X_discovery[child[[2L]],,drop=FALSE],
                         y_discovery[child[[2L]]],
                         discovery_pool,conf,seed=conf$seed+703L)
  parent_prob <- rlmt_predict_leaf(parent,X_confirmation,conf)$prob
  split_prob <- numeric(nrow(X_confirmation))
  zV <- X_confirmation[[pair$gene_a[[1L]]]] >
        X_confirmation[[pair$gene_b[[1L]]]]
  if (anyNA(zV)) stop("Missing confirmation routing value")
  if (any(zV)) {
    split_prob[zV] <- rlmt_predict_leaf(
      left,X_confirmation[zV,,drop=FALSE],conf)$prob
  }
  if (any(!zV)) {
    split_prob[!zV] <- rlmt_predict_leaf(
      right,X_confirmation[!zV,,drop=FALSE],conf)$prob
  }
  list(parent_prob=parent_prob,split_prob=split_prob,
       nominated_pair=pair[,c("gene_a","gene_b"),drop=FALSE])
}
