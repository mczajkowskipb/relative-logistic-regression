# Matched Gini/RLMT-LA biomedical pilot, pairing EXISTING H1/H2 models
# without refitting them or changing their source/results.
# EXACT old 2 datasets × 3 outer folds (6 paired folds).
# All data EXPOSED in prior development; NO independent generalization claim.
# Args: original frozen REF, already successful H1/H2 paper-pilot worktree.

source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v2_aligned_offset.R")
source("R/experimental/rlmt_v2_gini_lookahead.R")

options(warn=1L) # Make all glmnet/source warnings attributable to exact folds.

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)
  stop("FOURWAY_STOP: pass frozen original REF and completed H1/H2 pilot")
REF <- normalizePath(args[[1L]],mustWork=TRUE)
HIST <- normalizePath(args[[2L]],mustWork=TRUE)
OLD <- file.path(REF,"results","rlmt_dc_published_data_5fold_v1")
PREV <- file.path(HIST,"results","rlmt_v2_exposed_pilot_2x3_v1")
OUT <- "results/rlmt_v2_fourway_pilot_2x3_v1"
check <- function(ok,msg)
  if(!isTRUE(ok))stop("FOURWAY_STOP: ",msg,call.=FALSE)

check(dir.exists(OUT) &&
      file.exists(file.path(OUT,"CODE_SHA.txt")),
      "new results root not initialized")
SHA <- system("git rev-parse HEAD",intern=TRUE)[[1L]]
check(identical(readLines(file.path(OUT,"CODE_SHA.txt"))[[1L]],SHA),
      "code SHA mismatch")
check(identical(trimws(readLines(file.path(OLD,"status.txt"))[[1L]]),
                "RUNNER_EXIT=0"),"archived source incomplete")
check(identical(trimws(readLines(file.path(PREV,"status.txt"))[[1L]]),
                "RLMT_V2_PILOT_EXIT=0"),"H1/H2 result incomplete")
check(!file.exists(file.path(OUT,"run_manifest.csv")),
      "manifest already exists; no overwrite")
input <- read.csv(file.path(OLD,"input_manifest.csv"),
                  stringsAsFactors=FALSE)
reference <- read.csv(file.path(OLD,"run_manifest.csv"),
                      stringsAsFactors=FALSE)
oldh <- read.csv(file.path(PREV,"run_manifest.csv"),
                 stringsAsFactors=FALSE)
check(nrow(input)==19L && nrow(reference)==285L &&
      nrow(oldh)==12L,"source/pilot count differs")
check(setequal(oldh$mode,c("H1","H2")) &&
      all(oldh$status=="OK"),"H1/H2 method/status issue")

integrity <- function(file,expected) {
  r <- suppressWarnings(system2(
      "sha256sum", shQuote(normalizePath(file,mustWork=TRUE)),
      stdout=TRUE,stderr=TRUE))
  check(length(r)==1L &&
        grepl("^[a-f0-9]{64}[[:space:]]",r[[1L]]),
        "sha256sum error")
  check(identical(tolower(substr(r[[1L]],1L,64L)),
                  tolower(expected)),
        paste("archived original SHA256 drift",file))
}
readpred <- function(folder) {
  read.csv(file.path(folder,"predictions_test.csv"),
           stringsAsFactors=FALSE,
           colClasses=c(sample_id="character",
                        y_true="character",y_pred="character"))
}
keys <- function(pred)
  paste(pred$sample_id,pred$y_true,sep="\037")

DATASETS <- c("chen-2002","gordon-2002")
MODES <- c("TSPDT_GINI","RLMT_LA")
rows <- list()
n_done <- 0L

for(dataset in DATASETS) {
  entry <- input[input$dataset==dataset,,drop=FALSE]
  check(nrow(entry)==1L,paste("missing dataset",dataset))
  source_file <- as.character(entry$source_path[[1L]])
  integrity(source_file,as.character(entry$sha256[[1L]]))
  xy <- read_xy_file(source_file)
  rlmt_v2_check_xy(xy$X,xy$y)
  folds <- make_repeated_stratified_folds(
    xy$y,n_folds=5L,n_repeats=1L,seed=42L)
  check(identical(vapply(folds,function(f)f$split_id,character(1)),
                  sprintf("rep00_fold%02d",0:4)),
        "original fold scheme mismatched")
  world <- as.character(entry$world[[1L]])
  for(f in folds[1:3]) {
    split <- f$split_id
    oldmodel <- reference[
      reference$world==world &
      reference$dataset==dataset &
      reference$split_id==split &
      reference$method=="RLR",,drop=FALSE]
    check(nrow(oldmodel)==1L &&
          oldmodel$status[[1L]]=="OK",
          "missing frozen original RLR")
    oldpred <- readpred(file.path(REF,oldmodel$model_dir[[1L]]))
    tr <- as.integer(f$train_idx)
    te <- as.integer(f$test_idx)
    check(identical(sort(c(tr,te)),seq_len(nrow(xy$X))) &&
          !length(intersect(tr,te)) &&
          identical(rownames(xy$X)[te],oldpred$sample_id) &&
          identical(as.character(xy$y[te]),oldpred$y_true),
          "outer fold/outcome parity failed")
    prev <- oldh[oldh$world==world &
                 oldh$dataset==dataset &
                 oldh$split_id==split,,drop=FALSE]
    check(nrow(prev)==2L &&
          setequal(prev$mode,c("H1","H2")) &&
          all(prev$status=="OK"),
          "missing exact paired H1/H2 pilot")
    p_h1 <- readpred(file.path(HIST,
               prev$model_dir[prev$mode=="H1"]))
    p_h2 <- readpred(file.path(HIST,
               prev$model_dir[prev$mode=="H2"]))
    check(identical(keys(p_h1),keys(oldpred)) &&
          identical(keys(p_h2),keys(oldpred)) &&
          isTRUE(all.equal(p_h1$parent_prob,p_h2$parent_prob,
                           tolerance=1e-12)),
          "H1/H2 baseline predictions diverged")
    conf <- rlmt_v2_conf()$leaf
    Xtr <- xy$X[tr,,drop=FALSE]
    Xte <- xy$X[te,,drop=FALSE]
    ytr <- xy$y[tr]
    yte <- xy$y[te]
    for(mode in MODES) {
      dest <- file.path(OUT,"models",dataset,split,mode)
      check(!dir.exists(dest),"model root exists")
      message("FOURWAY_FIT=",dataset,":",split,":",mode)
      output <- rlmt_la_run(Xtr,ytr,Xte,yte,dest,method=mode)
      pred <- readpred(dest)
      check(identical(keys(pred),keys(oldpred)) &&
            isTRUE(all.equal(as.numeric(pred$prob),
                             output$prob,tolerance=1e-12)) &&
            isTRUE(all.equal(as.numeric(pred$parent_prob),
                             output$parent,tolerance=1e-12)) &&
            isTRUE(all.equal(as.numeric(pred$parent_prob),
                             as.numeric(p_h1$parent_prob),
                             tolerance=1e-12)),
            "new fit vs archival/H1/H2 predictions mismatch")
      if(identical(output$decision$status,"NO_SPLIT"))
        check(identical(output$parent,output$prob),
              "no-split not exactly aligned parent")
      yy <- rlmt_y01(yte)
      oldp <- as.numeric(oldpred$prob)
      check(all(is.finite(oldp)) &&
            all(oldp>=0 & oldp<=1),"old RLR probability invalid")
      ll <- rlmt_logloss(yy,output$prob)
      parentll <- rlmt_logloss(yy,output$parent)
      origll <- rlmt_logloss(yy,oldp)
      n_done <- n_done+1L
      rows[[n_done]] <- data.frame(
        dataset=dataset,world=world,split_id=split,mode=mode,
        status="OK",accepted=identical(output$decision$status,"ACCEPTED"),
        n_train=length(tr),n_test=length(te),
        new_logloss=ll,aligned_parent_logloss=parentll,
        original_RLR_logloss=origll,
        delta_vs_aligned=ll-parentll,delta_vs_original=ll-origll,
        new_brier=mean((output$prob-yy)^2),
        aligned_brier=mean((output$parent-yy)^2),
        original_RLR_brier=mean((oldp-yy)^2),
        model_dir=dest,
        old_RLR_dir=oldmodel$model_dir[[1L]],
        old_H1_model_dir=prev$model_dir[prev$mode=="H1"],
        old_H2_model_dir=prev$model_dir[prev$mode=="H2"],
        scientific_stage="EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY",
        stringsAsFactors=FALSE)
      write.csv(do.call(rbind,rows),
                file.path(OUT,"run_manifest.csv"),row.names=FALSE)
      message("FOURWAY_COMPLETED=",n_done,"/12; mode=",mode,
              "; accepted=",rows[[n_done]]$accepted,
              "; delta_LL_vs_aligned=",sprintf("%+.7f",ll-parentll))
    }
  }
}
check(n_done==12L && length(rows)==12L,
      "expected exactly 12 new real Gini/LA fits")
cat("RLMT_V2_FOURWAY_ADDITIONAL_FITS=12\n")
cat("RLMT_V2_FOURWAY_DATASETS=2\n")
cat("RLMT_V2_FOURWAY_PAIRED_FOLDS=6\n")
cat("RLMT_V2_FOURWAY_NEW_METHODS=TSPDT_GINI,RLMT_LA\n")
cat("RLMT_V2_FOURWAY_H1_H2_PREVIOUS_REUSED=TRUE\n")
cat("RLMT_V2_FOURWAY_RUN=PASS\n")
cat("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY\n")
