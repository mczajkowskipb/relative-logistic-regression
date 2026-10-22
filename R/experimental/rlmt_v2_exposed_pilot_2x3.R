# RLMT v2 H1/H2 exact matched-paper-fold QUICK PILOT, EXPOSED data ONLY.
# The previous 19x5 paper cohorts were already examined during method
# development. This is not an independent validation or new evidence of
# generalization. NO parameter changes after seeing these results.
#
# Run from NEW benchmark-pilot worktree; one positional argument:
# absolute original frozen 285-fit REF repository.
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v2_aligned_offset.R")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("V2_PILOT_STOP: original REF path required")
REF <- normalizePath(args[[1L]],mustWork=TRUE)
OLD <- file.path(REF,"results","rlmt_dc_published_data_5fold_v1")
OUT <- "results/rlmt_v2_exposed_pilot_2x3_v1"
stopifnot(dir.exists(OUT),
          file.exists(file.path(OUT,"CODE_SHA.txt")))
if(!identical(trimws(readLines(file.path(OUT,"CODE_SHA.txt"))[[1L]]),
              system("git rev-parse HEAD",intern=TRUE)[[1L]]))
  stop("V2_PILOT_STOP: code SHA drifted")
if(!identical(trimws(readLines(file.path(OLD,"status.txt"))[[1L]]),
              "RUNNER_EXIT=0"))
  stop("V2_PILOT_STOP: original immutable models not complete")
if(file.exists(file.path(OUT,"run_manifest.csv")))
  stop("V2_PILOT_STOP: manifest exists; never overwrite")
input <- read.csv(file.path(OLD,"input_manifest.csv"),
                  stringsAsFactors=FALSE)
archive <- read.csv(file.path(OLD,"run_manifest.csv"),
                    stringsAsFactors=FALSE)
stopifnot(nrow(input)==19L,nrow(archive)==285L)
DATASETS <- c("chen-2002","gordon-2002")
FOLDS <- sprintf("rep00_fold%02d",0:2)
MODES <- c("H1","H2")
check <- function(cond,context)
  if(!isTRUE(cond))stop("V2_PILOT_STOP: ",context,call.=FALSE)
integrity_sha <- function(path) {
  h <- suppressWarnings(system2("sha256sum",
     shQuote(normalizePath(path,mustWork=TRUE)),
     stdout=TRUE,stderr=TRUE))
  check(length(h)==1L && grepl("^[0-9a-f]{64}[[:space:]]+",h[[1L]]),
        paste("hash failed",path))
  substr(h[[1L]],1L,64L)
}
results <- list()
index <- 0L
for(dataset in DATASETS) {
  data_ref <- input[input$dataset==dataset,,drop=FALSE]
  check(nrow(data_ref)==1L,paste("dataset missing",dataset))
  file <- as.character(data_ref$source_path[[1L]])
  check(identical(integrity_sha(file),
                  as.character(data_ref$sha256[[1L]])),
        paste("original input SHA differs",dataset))
  xy <- read_xy_file(file)
  rlmt_v2_check_xy(xy$X,xy$y)
  folds <- make_repeated_stratified_folds(
    xy$y,n_folds=5L,n_repeats=1L,seed=42L)
  check(length(folds)==5L &&
        identical(vapply(folds,function(f)f$split_id,character(1)),
                  sprintf("rep00_fold%02d",0:4)),
        paste("old outer fold generator differs",dataset))
  world <- as.character(data_ref$world[[1L]])
  for(f in folds[1:3]) {
    sid <- f$split_id
    check(sid %in% FOLDS,"unexpected fold ID")
    base <- archive[archive$dataset==dataset &
                    archive$world==world &
                    archive$split_id==sid &
                    archive$method=="RLR",,drop=FALSE]
    check(nrow(base)==1L && base$status[[1L]]=="OK",
          paste("frozen original RLR not found",dataset,sid))
    original_dir <- file.path(REF,as.character(base$model_dir[[1L]]))
    oldpred <- read.csv(file.path(original_dir,"predictions_test.csv"),
                        stringsAsFactors=FALSE,
                        colClasses=c(sample_id="character",
                                     y_true="character",
                                     y_pred="character"))
    oldtrain <- read.csv(file.path(dirname(original_dir),"samples_train.csv"),
                         stringsAsFactors=FALSE,
                         colClasses=c(sample_id="character"))
    tr <- f$train_idx
    te <- f$test_idx
    check(length(intersect(tr,te))==0L &&
          identical(sort(c(tr,te)),seq_len(nrow(xy$X))) &&
          identical(rownames(xy$X)[tr],oldtrain$sample_id) &&
          identical(rownames(xy$X)[te],oldpred$sample_id) &&
          identical(as.character(xy$y[te]),oldpred$y_true),
          paste("source-paired train/test IDs mismatched",dataset,sid))
    orig_prob <- as.numeric(oldpred$prob)
    check(all(is.finite(orig_prob)) && all(orig_prob>=0&orig_prob<=1),
          "frozen original RLR test probability invalid")
    Xtr <- xy$X[tr,,drop=FALSE]
    Xte <- xy$X[te,,drop=FALSE]
    ytr <- xy$y[tr]
    yte <- xy$y[te]
    aligned <- NULL
    for(mode in MODES) {
      dest <- file.path(OUT,"models",dataset,sid,mode)
      check(!dir.exists(dest),"prior mode output exists; do not overwrite")
      message("V2_PILOT_FIT=",dataset,":",sid,":",mode)
      fit <- rlmt_v2_run(
        Xtr,ytr,Xte,yte,
        out_dir=dest,mode=mode)
      preds <- read.csv(file.path(dest,"predictions_test.csv"),
                        stringsAsFactors=FALSE,
                        colClasses=c(sample_id="character",
                                     y_true="character",
                                     y_pred="character"))
      check(identical(preds$sample_id,oldpred$sample_id) &&
            identical(preds$y_true,oldpred$y_true) &&
            isTRUE(all.equal(as.numeric(preds$prob),fit$prob,
                             tolerance=1e-12)) &&
            isTRUE(all.equal(as.numeric(preds$parent_prob),
                             fit$parent,tolerance=1e-12)),
            paste("saved results mismatch",dataset,sid,mode))
      if (!fit$decision$accepted)
        check(identical(fit$prob,fit$parent),
              "rejected route not bit-exact aligned parent")
      if(is.null(aligned)) {
        aligned <- fit$parent
      } else {
        check(isTRUE(all.equal(aligned,fit$parent,tolerance=1e-12)),
              "H1/H2 different aligned parent: cannot attribute tree gain")
      }
      y01 <- rlmt_y01(yte)
      v2_ll <- rlmt_logloss(y01,fit$prob)
      aligned_ll <- rlmt_logloss(y01,fit$parent)
      orig_ll <- rlmt_logloss(y01,orig_prob)
      index <- index+1L
      results[[index]] <- data.frame(
        dataset=dataset,world=world,split_id=sid,mode=mode,
        status="OK",
        n_train=length(tr),n_test=length(te),
        accepted=isTRUE(fit$decision$accepted),
        new_logloss=v2_ll,aligned_parent_logloss=aligned_ll,
        original_RLR_logloss=orig_ll,
        delta_vs_aligned=v2_ll-aligned_ll,
        delta_vs_original=v2_ll-orig_ll,
        new_brier=mean((fit$prob-y01)^2),
        aligned_brier=mean((fit$parent-y01)^2),
        original_RLR_brier=mean((orig_prob-y01)^2),
        model_dir=dest,original_RLR_dir=as.character(base$model_dir[[1L]]),
        original_RLR_retained=TRUE,
        scientific_stage="EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY",
        stringsAsFactors=FALSE)
      write.csv(do.call(rbind,results),
                file.path(OUT,"run_manifest.csv"),row.names=FALSE)
      message("V2_PILOT_COMPLETED=",index,"/12; ",
              "mode=",mode,"; accepted=",fit$decision$accepted,
              "; delta_logloss_vs_aligned=",
              sprintf("%+.7f",v2_ll-aligned_ll))
    }
  }
}
check(index==12L,"expected 12 H1/H2 fits (2 datasets x 3 folds x 2)")
cat("RLMT_V2_PILOT_REAL_FITS=12\n")
cat("RLMT_V2_PILOT_DATASETS=2\n")
cat("RLMT_V2_PILOT_PAIRED_OUTER_FOLDS=6\n")
cat("RLMT_V2_PILOT_COMPLETE=PASS\n")
cat("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY\n")
