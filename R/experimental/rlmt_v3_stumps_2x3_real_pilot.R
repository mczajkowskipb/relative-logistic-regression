# New family, matched exactly 2 previously-exposed cohorts x original
# first 3 outer folds x 3 NEW method variants; 18 attempted fits.
# Train/test IDs, source hash & original RLR probability verified.
# No reuse of any heldout outcome to change algorithm parameters.
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v2_aligned_offset.R")
source("R/experimental/rlmt_v3_stumps_pruned.R")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("V3_PILOT_STOP: original REF path required")
REF <- normalizePath(args[[1L]],mustWork=TRUE)
OLD <- file.path(REF,"results","rlmt_dc_published_data_5fold_v1")
OUT <- "results/rlmt_v3_stumps_2x3_pilot_v1"
stopifnot(dir.exists(OUT),
          file.exists(file.path(OUT,"CODE_SHA.txt")))
need <- function(ok,msg)
  if(!isTRUE(ok))stop("V3_PILOT_STOP: ",msg,call.=FALSE)
need(identical(trimws(readLines(file.path(OUT,"CODE_SHA.txt"))[[1L]]),
               system("git rev-parse HEAD",intern=TRUE)[[1L]]),
     "code SHA not pinned")
need(identical(trimws(readLines(file.path(OLD,"status.txt"))[[1L]]),
               "RUNNER_EXIT=0"),"old 285 reference incomplete")
need(!file.exists(file.path(OUT,"run_manifest.csv")),
     "old pilot manifest exists")
manifest <- read.csv(file.path(OLD,"input_manifest.csv"),
                     stringsAsFactors=FALSE)
models <- read.csv(file.path(OLD,"run_manifest.csv"),
                   stringsAsFactors=FALSE)
need(nrow(manifest)==19L && nrow(models)==285L,
     "original cohorts/models not complete")
methods <- c("STUMP_TSP_BETA",
             "STUMP_RIDGE_OFFSET","TREE_D2_PRUNED_BETA")
inputs <- c("chen-2002","gordon-2002")
rows <- list()
warningrows <- list()
index <- 0L
sha <- function(f) {
  r <- system2("sha256sum",shQuote(
          normalizePath(f,mustWork=TRUE)),stdout=TRUE)
  need(length(r)==1L && grepl("^[a-f0-9]{64}[[:space:]]+",r[[1L]]),
       "could not hash source")
  substr(r[[1L]],1L,64L)
}
readpred <- function(f)
  read.csv(file.path(f,"predictions_test.csv"),
           stringsAsFactors=FALSE,
           colClasses=c(sample_id="character",
                        y_true="character",y_pred="character"))
for(dataset in inputs) {
  entry <- manifest[manifest$dataset==dataset,,drop=FALSE]
  need(nrow(entry)==1L,"missing selected dataset")
  need(identical(tolower(sha(entry$source_path[[1L]])),
                 tolower(entry$sha256[[1L]])),
       paste("source mismatch",dataset))
  xy <- read_xy_file(entry$source_path[[1L]])
  rlmt_v2_check_xy(xy$X,xy$y)
  folds <- make_repeated_stratified_folds(
    xy$y,n_folds=5L,n_repeats=1L,seed=42L)
  need(identical(vapply(folds,function(x)x$split_id,character(1)),
                 sprintf("rep00_fold%02d",0:4)),
       "fold ID generator changed")
  for(fold in folds[1:3]) {
    sid <- fold$split_id
    reference <- models[
      models$dataset==dataset &
      models$world==entry$world[[1L]] &
      models$split_id==sid &
      models$method=="RLR",,drop=FALSE]
    need(nrow(reference)==1L &&
         reference$status[[1L]]=="OK",
         "frozen original fold missing")
    olddir <- file.path(REF,reference$model_dir[[1L]])
    oldpred <- readpred(olddir)
    oldtrain <- read.csv(file.path(dirname(olddir),
                                 "samples_train.csv"),
                         stringsAsFactors=FALSE,
                         colClasses=c(sample_id="character"))
    tr <- fold$train_idx
    te <- fold$test_idx
    need(identical(sort(c(tr,te)),seq_len(nrow(xy$X))) &&
         identical(oldtrain$sample_id,rownames(xy$X)[tr]) &&
         identical(oldpred$sample_id,rownames(xy$X)[te]) &&
         identical(oldpred$y_true,as.character(xy$y[te])),
         "original frozen train/test patient pairing changed")
    XX <- xy$X[tr,,drop=FALSE]
    yy <- xy$y[tr]
    TT <- xy$X[te,,drop=FALSE]
    yt <- xy$y[te]
    old_prob <- as.numeric(oldpred$prob)
    need(all(is.finite(old_prob)) &&
         all(old_prob>=0 & old_prob<=1),
         "original frozen RLR probability invalid")
    parent <- NULL
    for(method in methods) {
      dest <- file.path(OUT,"models",dataset,sid,method)
      need(!dir.exists(dest),"model prior result dir exists")
      message("V3_PILOT_FIT=",dataset,":",sid,":",method)
      captured <- character()
      answer <- tryCatch(withCallingHandlers({
        rlmt_v3_run(XX,yy,TT,yt,dest,method)
      },warning=function(w) {
        captured <<- c(captured,conditionMessage(w))
        invokeRestart("muffleWarning")
      }),error=function(e)e)
      if(inherits(answer,"error")) {
        message("V3_PILOT_MODEL_ERROR=",conditionMessage(answer))
        stop(answer)
      }
      newp <- readpred(dest)
      need(identical(newp$sample_id,oldpred$sample_id) &&
           identical(newp$y_true,oldpred$y_true) &&
           isTRUE(all.equal(as.numeric(newp$prob),
                            answer$prob,tolerance=1e-12)) &&
           isTRUE(all.equal(as.numeric(newp$parent_prob),
                            answer$parent,tolerance=1e-12)),
           "RDS/prediction sample alignment problem")
      if(is.null(parent))parent <- answer$parent
      else need(isTRUE(all.equal(parent,answer$parent,
                                 tolerance=1e-12)),
                "different aligned global parent between new methods")
      y01 <- rlmt_y01(yt)
      index <- index+1L
      rows[[index]] <- data.frame(
        world=entry$world[[1L]],dataset=dataset,
        split_id=sid,method=method,status="OK",
        n_test=length(te),n_train=length(tr),
        realized_internal_nodes=answer$metrics$n_internal_nodes,
        chosen_depth=answer$metrics$chosen_depth,
        logloss=rlmt_logloss(y01,answer$prob),
        aligned_logloss=rlmt_logloss(y01,answer$parent),
        original_RLR_logloss=rlmt_logloss(y01,old_prob),
        delta_vs_aligned=rlmt_logloss(y01,answer$prob)-
          rlmt_logloss(y01,answer$parent),
        delta_vs_original=rlmt_logloss(y01,answer$prob)-
          rlmt_logloss(y01,old_prob),
        brier=mean((answer$prob-y01)^2),
        warning_count=length(captured),
        model_dir=dest,original_model_dir=reference$model_dir[[1L]],
        scientific_stage="EXPOSED_RETROSPECTIVE_DEV_ONLY",
        stringsAsFactors=FALSE)
      write.csv(do.call(rbind,rows),
                file.path(OUT,"run_manifest.csv"),row.names=FALSE)
      if(length(captured)) {
        tab <- as.data.frame(table(captured),stringsAsFactors=FALSE)
        names(tab) <- c("message","frequency")
        tab$dataset <- dataset
        tab$split_id <- sid
        tab$method <- method
        warningrows[[length(warningrows)+1L]] <-
          tab[,c("dataset","split_id","method",
                 "message","frequency"),drop=FALSE]
        write.csv(do.call(rbind,warningrows),
                  file.path(OUT,"warning_inventory.csv"),
                  row.names=FALSE)
      }
      message("V3_PILOT_COMPLETED=",index,"/18 method=",method,
              " nodes=",answer$metrics$n_internal_nodes,
              " delta_original=",sprintf("%+.6f",
                      rows[[index]]$delta_vs_original),
              " warnings=",length(captured))
    }
  }
}
need(index==18L,"expected all 18 v3 real fits")
cat("RLMT_V3_PILOT_REAL_FITS=18\n")
cat("RLMT_V3_PILOT_COEXISTING_RLR_ORIGINAL_READONLY=TRUE\n")
cat("RLMT_V3_PILOT_COMPLETE=PASS\n")
cat("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY\n")
