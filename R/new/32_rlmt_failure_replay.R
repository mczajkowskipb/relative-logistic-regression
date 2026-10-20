# Test-only replay of the three recorded non-conformable-arguments failures.
# Reads original pilot artifacts and reconstructs the same outer/inner folds.
# Only fitting for affected candidates; NO re-ranking outside the saved pool.
# Does not change RLMT methods, CV rules or the existing pilot outputs.
# Usage: Rscript --vanilla R/new/32_rlmt_failure_replay.R results/rlmt_gate3_pilot_2x3/config.json

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 32_rlmt_failure_replay.R <pilot_config.json>")
suppressPackageStartupMessages(library(jsonlite))
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")

cfg <- jsonlite::read_json(args[[1]], simplifyVector = TRUE)
n_worlds <- if (is.data.frame(cfg$worlds)) nrow(cfg$worlds) else length(cfg$worlds)
if (n_worlds != 1L) stop("Expected one pilot world")
world <- if (is.data.frame(cfg$worlds)) as.list(cfg$worlds[1, , drop=FALSE]) else cfg$worlds[[1]]
datasets <- as.character(world$datasets_include[[1L]])
if (!identical(sort(datasets), c("chen-2002", "singh-2002")) ||
    as.integer(world$n_folds) != 5L ||
    as.integer(world$max_splits_per_dataset) != 3L ||
    as.integer(cfg$models$seed) != 42L) stop("Wrong pilot configuration")

root <- cfg$out_root
audit_file <- file.path(root, "audit_gate3_readonly_v1", "candidate_support.csv")
if (!file.exists(audit_file)) stop("Missing prior GATE-3 audit: ", audit_file)
records <- read.csv(audit_file, stringsAsFactors = FALSE)
dim_errors <- records[records$reason_group == "other_fit_error", , drop = FALSE]
if (nrow(dim_errors) != 3L) stop("Expected exactly 3 previously recorded dimension errors")
if (!all(grepl("non-conformable arguments", dim_errors$rejection_reason,
               fixed=TRUE))) stop("Unexpected dimension error group")

output <- file.path(root, "audit_gate3_failure_replay_v1")
if (dir.exists(output)) stop("Refusing to overwrite prior replay: ", output)

conf <- rlmt_defaults(cfg$models$rlmt)
results <- list()
add <- function(dataset, split, candidate, fold, side, stage,
                n, nclass0, nclass1, nfeatures, cvfolds, result, detail) {
  results[[length(results)+1L]] <<- data.frame(
    dataset=dataset, split_id=split, candidate=candidate,
    inner_fold=fold, side=side, stage=stage,
    child_train_n=n, child_class0=nclass0, child_class1=nclass1,
    relation_features=nfeatures, glmnet_nfolds=cvfolds,
    result=result, detail=detail, stringsAsFactors=FALSE)
}

# One fit, with matching operations and exact seed from rlmt_fit_leaf().
fit_debug <- function(X, y, pool, seed) {
  stage <- "rerank_pool_in_node"
  info <- list(stage=stage, nfeatures=NA_integer_)
  warnings <- character(0)
  out <- withCallingHandlers(tryCatch({
    pairs <- rerank_pool_in_node(X, y, pool)
    stage <- "apply_disjointness"
    pairs <- apply_disjointness(pairs, K=conf$K, d=conf$d)
    if (!nrow(pairs)) stop("Empty shortlist")
    stage <- "encode_pair_matrix"
    R <- encode_pair_matrix(X, pairs, m=conf$m)
    stage <- "glmnet::cv.glmnet"
    set.seed(seed)
    fit <- glmnet::cv.glmnet(
      x=R, y=rlmt_y01(y), family="binomial",
      alpha=conf$alpha, nfolds=safe_nfolds(y, 5),
      type.measure="deviance")
    stage <- "coef(cvfit)"
    beta <- as.numeric(coef(fit, s=conf$lambda_choice))[-1L]
    names(beta) <- colnames(R)
    list(ok=TRUE, stage="fitted", fit=fit, pairs=pairs,
         beta=beta, nfeatures=ncol(R), warnings=warnings)
  }, error=function(e) {
    list(ok=FALSE, stage=stage, error=conditionMessage(e),
         nfeatures=if (exists("R", inherits=FALSE)) ncol(R) else NA_integer_)
  }), warning=function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  out$warnings <- warnings
  out
}

for (i in seq_len(nrow(dim_errors))) {
  case <- dim_errors[i, , drop=FALSE]
  dataset <- case$dataset
  split <- case$split_id
  xy <- read_xy_file(file.path(world$comp_dir, paste0(dataset, ".arff")))
  outer <- make_repeated_stratified_folds(
    xy$y, n_folds=5L, n_repeats=1L, seed=42L)
  spec <- outer[[as.integer(sub(".*fold", "", split)) + 1L]]
  if (!identical(spec$split_id, split)) stop("Outer split reconstruction mismatch")
  X <- xy$X[spec$train_idx, , drop=FALSE]
  y <- xy$y[spec$train_idx]
  inner <- make_repeated_stratified_folds(y,
    n_folds=conf$inner_folds, n_repeats=1L, seed=conf$seed)
  base <- file.path(root, "worlds", world$name, "datasets", dataset,
                    "splits", split)
  pool_file <- file.path(base, "pairs_train", "top1000", "top_1000_pairs.csv")
  if (!file.exists(pool_file)) stop("Missing saved global pool")
  pool <- read.csv(pool_file, stringsAsFactors=FALSE)
  ga <- case$gene_a
  gb <- case$gene_b
  if (!(ga %in% names(X) && gb %in% names(X))) stop("Missing split gene")
  z <- X[[ga]] > X[[gb]]
  candidate <- paste(ga, ">", gb)
  cat("\nREPLAY:", dataset, split, candidate, "\n")

  found_first_error <- FALSE
  for (j in seq_along(inner)) {
    f <- inner[[j]]
    train <- f$train_idx
    va <- f$test_idx
    sides <- list(train[z[train]], train[!z[train]])
    local_models <- list()
    for (side in c(1L, 2L)) {
      idx <- sides[[side]]
      classcounts <- as.integer(table(y[idx]))
      side_name <- if (side == 1L) "left" else "right"
      if (length(idx) < 6L || any(classcounts < 3L)) {
        add(dataset, split, candidate, j, side_name,
            "inner_child_support", length(idx), classcounts[1],
            classcounts[2], NA_integer_, NA_integer_,
            "INVALID_SUPPORT", "Requires at least 6 and 3 per class")
        found_first_error <- TRUE
        break
      }

      model <- fit_debug(X[idx, , drop=FALSE], y[idx], pool,
                         conf$seed + j*100L + side)
      if (!model$ok) {
        add(dataset, split, candidate, j, side_name,
            model$stage, length(idx), classcounts[1],
            classcounts[2], model$nfeatures, safe_nfolds(y[idx],5),
            "ERROR", model$error)
        cat("  FIRST ERROR inner=", j, " child=",side_name,
            " stage=",model$stage," : ",model$error,"\n",sep="")
        found_first_error <- TRUE
        break
      }
      add(dataset, split, candidate, j, side_name,
          "glmnet::cv.glmnet", length(idx), classcounts[1],
          classcounts[2], model$nfeatures, safe_nfolds(y[idx],5),
          "FIT_OK", if (length(model$warnings))
            paste(unique(model$warnings), collapse=" | ") else "")
      local_models[[side]] <- model
    }
    if (found_first_error) break

    for (side in c(1L, 2L)) {
      on_side <- if (side == 1L) z[va] else !z[va]
      sel <- va[on_side]
      if (!length(sel)) next
      stage <- "predict"
      err <- tryCatch({
        rlmt_predict_leaf(local_models[[side]], X[sel, , drop=FALSE], conf)
        NULL
      }, error=function(e) conditionMessage(e))
      if (!is.null(err)) {
        cc <- as.integer(table(y[sides[[side]]]))
        add(dataset, split, candidate, j,
            if (side == 1L) "left" else "right",
            stage, length(sides[[side]]), cc[1], cc[2],
            length(local_models[[side]]$beta),
            safe_nfolds(y[sides[[side]]],5), "ERROR", err)
        cat("  FIRST ERROR inner=",j," side=",side,
            " stage=predict : ",err,"\n",sep="")
        found_first_error <- TRUE
        break
      }
    }
    if (found_first_error) break
  }
  if (!found_first_error) cat("  REPLAY_NOT_REPRODUCED\n")
}

out <- do.call(rbind, results)
dir.create(output, recursive=TRUE, showWarnings=FALSE)
write.csv(out, file.path(output,"stage_trace.csv"), row.names=FALSE)
errors <- out[out$result == "ERROR", , drop=FALSE]
status <- if (nrow(errors) == 3L &&
              all(grepl("non-conformable arguments", errors$detail, fixed=TRUE))) {
  "THREE_DIMENSION_ERRORS_REPRODUCED"
} else "REPLAY_NEEDS_REVIEW"
writeLines(c(
  "RLMT failure replay - diagnostic only, no experimental results.",
  "Original production code and saved pilot metrics untouched.",
  paste("Observed original dimension errors:", nrow(dim_errors)),
  paste("Reproduced dimension errors:", nrow(errors)),
  paste("Status:", status),
  "See stage_trace.csv for exact failed stage and local class support."
), file.path(output,"README.txt"))
cat("\nSTAGE_REPLAY_STATUS=",status,"\n",sep="")
cat("REPLAY_ERRORS=",nrow(errors),"\n",sep="")
if (nrow(errors)) {
  for (i in seq_len(nrow(errors))) {
    r <- errors[i,]
    cat(r$dataset," ",r$split_id," ",r$stage," ",r$detail,"\n",sep="")
  }
}
cat("OUTPUT=",output,"\n",sep="")
