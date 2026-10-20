# Regression gate for the deterministic RLMT-only stratified glmnet CV fix.
# Refit only the THREE saved failures, on their unchanged local training subsets.
# Usage: Rscript --vanilla R/new/34_rlmt_stratcv_regression.R results/rlmt_gate3_pilot_2x3/config.json

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 34_rlmt_stratcv_regression.R <old pilot config.json>")
suppressPackageStartupMessages(library(jsonlite))
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")

cfg <- jsonlite::read_json(args[1], simplifyVector = TRUE)
n_worlds <- if (is.data.frame(cfg$worlds)) nrow(cfg$worlds) else length(cfg$worlds)
if (n_worlds != 1L) stop("Expected original single world")
world <- if (is.data.frame(cfg$worlds)) as.list(cfg$worlds[1, , drop=FALSE]) else cfg$worlds[[1]]
if (!identical(sort(as.character(world$datasets_include[[1]])),
               c("chen-2002", "singh-2002")) ||
    as.integer(world$n_folds) != 5L ||
    as.integer(world$max_splits_per_dataset) != 3L ||
    as.integer(cfg$models$seed) != 42L) stop("Not the fixed 2x3 pilot")

conf <- rlmt_defaults(cfg$models$rlmt)
old_root <- cfg$out_root
trace_path <- file.path(old_root, "audit_gate3_failure_replay_v1", "stage_trace.csv")
if (!file.exists(trace_path)) stop("Missing old failure-replay trace")
tr <- read.csv(trace_path, stringsAsFactors=FALSE)
fail <- tr[tr$result=="ERROR" & tr$stage=="glmnet::cv.glmnet",,drop=FALSE]
if (nrow(fail) != 3L) stop("Expected exactly 3 reproduced old failures")
regout <- file.path(old_root, "regression_stratcv_v1")
if (dir.exists(regout)) stop("Refusing to overwrite prior regression: ", regout)

# Standalone invariants on arbitrary binary labels and minimal supported leaves.
for (counts in list(c(3L,17L), c(5L,38L), c(20L,40L))) {
  y <- factor(rep(c("0", "1"), counts))
  k <- safe_nfolds(y, 5)
  a <- rlmt_stratified_foldid(y, nfolds=k, seed=42L)
  b <- rlmt_stratified_foldid(y, nfolds=k, seed=42L)
  if (!identical(a,b) || length(a)!=length(y) ||
      length(unique(a))!=k ||
      any(vapply(seq_len(k), function(i) {
        any(table(y[a==i]) < 1L) ||
        any(table(y[a!=i]) < 2L)
      }, logical(1)))) stop("Invalid foldid determinism or support")
}
bad <- try(rlmt_stratified_foldid(factor(c(rep("0",15),rep("1",2))),
                                     nfolds=3L,seed=42L), silent=TRUE)
if (!inherits(bad,"try-error")) stop("Unsafe minority support was accepted")

rows <- list()
for (i in seq_len(nrow(fail))) {
  case <- fail[i,,drop=FALSE]
  dataset <- case$dataset
  split <- case$split_id
  xy <- read_xy_file(file.path(world$comp_dir, paste0(dataset, ".arff")))
  outer <- make_repeated_stratified_folds(xy$y,
                        n_folds=5L,n_repeats=1L,seed=42L)
  fnum <- as.integer(sub("^.*fold","",split))+1L
  o <- outer[[fnum]]
  if (!identical(o$split_id,split)) stop("Outer fold changed")
  Xtr <- xy$X[o$train_idx,,drop=FALSE]
  ytr <- xy$y[o$train_idx]
  ins <- make_repeated_stratified_folds(ytr,n_folds=conf$inner_folds,
                                         n_repeats=1L,seed=conf$seed)
  j <- as.integer(case$inner_fold)
  side <- if (identical(case$side,"left")) 1L else 2L
  f <- ins[[j]]
  gb <- strsplit(case$candidate," > ",fixed=TRUE)[[1]]
  if (length(gb)!=2L) stop("Malformed gene pair")
  z <- Xtr[[gb[1L]]] > Xtr[[gb[2L]]]
  train <- f$train_idx
  idx <- if (side==1L) train[z[train]] else train[!z[train]]
  X <- Xtr[idx,,drop=FALSE]
  y <- ytr[idx]
  k <- safe_nfolds(y,5)
  seed <- conf$seed + j*100L + side
  expected_foldid <- rlmt_stratified_foldid(y, k, seed)
  pool_path <- file.path(old_root,"worlds",world$name,"datasets",
                   dataset,"splits",split,"pairs_train","top1000",
                   "top_1000_pairs.csv")
  pool <- read.csv(pool_path, stringsAsFactors=FALSE)
  warns <- character()
  fit <- withCallingHandlers(
    rlmt_fit_leaf(X,y,pool,conf,seed=seed),
    warning=function(w) {
      warns <<- c(warns,conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (!inherits(fit$fit,"cv.glmnet") ||
      !isTRUE(all.equal(as.integer(fit$fit$foldid),
                         as.integer(expected_foldid))) ||
      nrow(fit$active)>conf$K ||
      length(fit$beta)<1L) {
    stop("cv.glmnet fit/foldid mismatch: ",dataset," / ",split)
  }
  pred <- rlmt_predict_leaf(fit,X,conf)$prob
  if (length(pred)!=nrow(X) || any(!is.finite(pred)) ||
      any(pred<0 | pred>1)) stop("Invalid local predictions")
  cat("STRATCV_FIT_PASS ",dataset," ",split,
      " inner=",j," side=",case$side,
      " class=",paste(as.integer(table(y)),collapse="/"),
      " folds=",k," warnings=",length(unique(warns)),"\n",sep="")
  rows[[i]] <- data.frame(dataset=dataset,split_id=split,
                           inner_fold=j,side=case$side,
                           fit=TRUE,n=nrow(X),class0=as.integer(table(y)[1]),
                           class1=as.integer(table(y)[2]),
                           nfolds=k, warnings=length(unique(warns)),
                           stringsAsFactors=FALSE)
}

report <- do.call(rbind,rows)
if (nrow(report)!=3L || !all(report$fit)) stop("Regression not complete")
dir.create(regout,recursive=TRUE,showWarnings=FALSE)
write.csv(report,file.path(regout,"results.csv"),row.names=FALSE)
writeLines(c("Exactly three formerly failing local cv.glmnet fits now PASS.",
             "Stratified foldid checked for determinism, class support and exact glmnet use.",
             "Original GATE-3 results unchanged. Test-only fit; not scientific metrics."),
           file.path(regout,"README.txt"))
cat("STRATIFIED_CV_REGRESSION=PASS (3/3)\n")
cat("Original pilot remains unchanged; new full pilot must use distinct out_root.\n")
