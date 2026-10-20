# RLMT: diagnostic-only check of glmnet's random CV fold allocation on
# the three exactly reproduced dimension errors (GATE-3).
# No production code/benchmark outputs are changed.
# Usage: Rscript --vanilla R/new/33_rlmt_cv_foldid_probe.R results/rlmt_gate3_pilot_2x3/config.json

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 33_rlmt_cv_foldid_probe.R <pilot_config.json>")
suppressPackageStartupMessages(library(jsonlite))
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
cfg <- jsonlite::read_json(args[[1]], simplifyVector = TRUE)
n_worlds <- if (is.data.frame(cfg$worlds)) nrow(cfg$worlds) else length(cfg$worlds)
if (n_worlds != 1L) stop("Expected one world")
world <- if (is.data.frame(cfg$worlds)) as.list(cfg$worlds[1L, , drop = FALSE]) else cfg$worlds[[1L]]
datasets <- as.character(world$datasets_include[[1L]])
if (!identical(sort(datasets), c("chen-2002", "singh-2002")) ||
    as.integer(cfg$models$seed) != 42L ||
    as.integer(world$n_folds) != 5L ||
    as.integer(world$max_splits_per_dataset) != 3L) stop("Not the GATE-3 pilot")
root <- cfg$out_root
trace_file <- file.path(root, "audit_gate3_failure_replay_v1", "stage_trace.csv")
if (!file.exists(trace_file)) stop("Need previous stage-specific replay: ", trace_file)
trace <- read.csv(trace_file, stringsAsFactors = FALSE)
failed <- trace[trace$result == "ERROR" &
                  trace$stage == "glmnet::cv.glmnet", , drop = FALSE]
if (nrow(failed) != 3L) stop("Expected three cv.glmnet failures")
dest <- file.path(root, "audit_gate3_foldid_probe_v1")
if (dir.exists(dest)) stop("Refusing to overwrite prior probe: ", dest)

conf <- rlmt_defaults(cfg$models$rlmt)
result_rows <- list()
warnings_all <- list()

for (i in seq_len(nrow(failed))) {
  rec <- failed[i, , drop = FALSE]
  dataset <- rec$dataset
  split <- rec$split_id
  xy <- read_xy_file(file.path(world$comp_dir, paste0(dataset, ".arff")))
  outers <- make_repeated_stratified_folds(xy$y, n_folds=5L,
                                           n_repeats=1L, seed=42L)
  outer_id <- as.integer(sub("^.*fold", "", split)) + 1L
  outer <- outers[[outer_id]]
  if (!identical(outer$split_id, split)) stop("Reconstructed outer split differs")
  Xtr <- xy$X[outer$train_idx, , drop=FALSE]
  ytr <- xy$y[outer$train_idx]
  inner <- make_repeated_stratified_folds(ytr, n_folds=conf$inner_folds,
                                           n_repeats=1L, seed=conf$seed)
  j <- as.integer(rec$inner_fold)
  side <- if (identical(rec$side, "left")) 1L else 2L
  f <- inner[[j]]
  genes <- strsplit(rec$candidate, " > ", fixed=TRUE)[[1L]]
  if (length(genes) != 2L) stop("Bad candidate name in replay")
  z <- Xtr[[genes[1L]]] > Xtr[[genes[2L]]]
  tr <- f$train_idx
  idx <- if (side == 1L) tr[z[tr]] else tr[!z[tr]]
  X <- Xtr[idx, , drop=FALSE]
  y <- ytr[idx]
  seed <- conf$seed + j*100L + side
  class_n <- as.integer(table(y))
  if (any(class_n < 3L)) stop("Insufficient child class counts")
  pool_file <- file.path(root, "worlds", world$name, "datasets", dataset,
                         "splits", split, "pairs_train", "top1000",
                         "top_1000_pairs.csv")
  pool <- read.csv(pool_file, stringsAsFactors=FALSE)
  pairs <- apply_disjointness(rerank_pool_in_node(X, y, pool),
                              K=conf$K, d=conf$d)
  R <- encode_pair_matrix(X, pairs, m=conf$m)
  folds <- safe_nfolds(y,5)
  if (folds < 3L) stop("Cannot run glmnet CV with <3 folds")
  Y <- rlmt_y01(y)

  # Reconstruct exactly the default cv.glmnet fold assignment. The function
  # samples foldid before calling the full-data glmnet model.
  set.seed(seed)
  random_foldid <- sample(rep(seq_len(folds), length.out=length(Y)))
  random_counts <- table(factor(random_foldid, levels=seq_len(folds)), y)
  random_min_train <- min(vapply(seq_len(folds), function(k)
    min(as.integer(table(y[random_foldid != k]))), integer(1)))
  random_min_test <- min(as.integer(random_counts))

  # Deterministic stratified fold assignment within the leaf TRAINING set only.
  set.seed(seed)
  balanced_foldid <- integer(length(Y))
  for (cls in levels(y)) {
    ci <- sample(which(y == cls))
    balanced_foldid[ci] <- rep(seq_len(folds), length.out=length(ci))
  }
  balanced_counts <- table(factor(balanced_foldid,levels=seq_len(folds)), y)
  balanced_min_train <- min(vapply(seq_len(folds), function(k)
    min(as.integer(table(y[balanced_foldid != k]))), integer(1)))
  balanced_min_test <- min(as.integer(balanced_counts))

  # Test the SAME matrix with different CV partitions. Separately report
  # unstratified and stratified results; do not swallow fit warnings silently.
  tryfit <- function(explicit_balanced) {
    warns <- character(0)
    set.seed(seed)
    res <- withCallingHandlers(tryCatch({
      fit <- glmnet::cv.glmnet(
        x=R, y=Y, family="binomial", alpha=conf$alpha,
        nfolds=folds, type.measure="deviance",
        foldid=if (explicit_balanced) balanced_foldid else NULL)
      list(ok=TRUE, message="OK",
           nlambda=length(fit$lambda),
           lambda_min=fit$lambda.min)
    }, error=function(e) {
      list(ok=FALSE, message=conditionMessage(e),
           nlambda=NA_integer_, lambda_min=NA_real_)
    }), warning=function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    res$warnings <- unique(warns)
    res
  }
  old <- tryfit(FALSE)
  balanced <- tryfit(TRUE)
  showcounts <- function(x) paste(apply(x,1,function(r)
    paste(r,collapse="/")),collapse=";")

  cat("\n=== CASE ",i,": ",dataset," ",split," ",rec$candidate," ===\n",sep="")
  cat("inner_fold=",j," side=",rec$side," n=",nrow(R),
      " p=",ncol(R)," class=",paste(class_n,collapse="/"),
      " cv_folds=",folds,"\n",sep="")
  cat("default fold validation class counts: ",showcounts(random_counts),"\n",sep="")
  cat("stratified validation class counts:  ",showcounts(balanced_counts),"\n",sep="")
  cat("default:    ",if(old$ok)"FIT_OK" else paste("ERROR",old$message),"\n",sep="")
  cat("stratified: ",if(balanced$ok)"FIT_OK" else paste("ERROR",balanced$message),"\n",sep="")

  result_rows[[i]] <- data.frame(
    dataset=dataset, split_id=split, candidate=rec$candidate,
    inner_fold=j, side=rec$side, child_n=nrow(R),
    child_class0=class_n[1], child_class1=class_n[2],
    num_features=ncol(R), cv_nfolds=folds,
    default_min_valid_class=random_min_test,
    default_min_train_class=random_min_train,
    default_counts=showcounts(random_counts),
    default_ok=old$ok, default_message=old$message,
    stratified_min_valid_class=balanced_min_test,
    stratified_min_train_class=balanced_min_train,
    stratified_counts=showcounts(balanced_counts),
    stratified_ok=balanced$ok,
    stratified_message=balanced$message, stringsAsFactors=FALSE)
  warnings_all[[i]] <- c(paste0("CASE ",i," / DEFAULT"),old$warnings,
                         paste0("CASE ",i," / STRATIFIED"),balanced$warnings)
}
table_out <- do.call(rbind, result_rows)
dir.create(dest, recursive=TRUE, showWarnings=FALSE)
write.csv(table_out, file.path(dest,"foldid_comparison.csv"), row.names=FALSE)
writeLines(unlist(warnings_all), file.path(dest,"warnings.txt"))
cat("\n=== SUMMARY ===\n")
cat("Default CV errors: ",sum(!table_out$default_ok),"/3\n",sep="")
cat("Stratified CV successful: ",sum(table_out$stratified_ok),"/3\n",sep="")
if (all(!table_out$default_ok) && all(table_out$stratified_ok) &&
    all(table_out$stratified_min_valid_class >= 1L) &&
    all(table_out$stratified_min_train_class >= 2L)) {
  cat("DIAGNOSTIC_CONCLUSION=ALL_3_FIXED_BY_STRATIFIED_CV_IN_ISOLATED_TEST\n")
} else {
  cat("DIAGNOSTIC_CONCLUSION=INVESTIGATE_BEFORE_PRODUCTION_CHANGE\n")
}
cat("READ_ONLY_PROBE=COMPLETE\n")
cat("REPORT=",dest,"\n",sep="")
