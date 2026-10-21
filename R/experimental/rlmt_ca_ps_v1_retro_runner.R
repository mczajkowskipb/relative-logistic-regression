# 95 NEW RLMT-CA-PS-v1 fits. The 285 ORIGINAL models remain read-only.
# All 19 paper cohorts have been exposed; this is retrospective development.
# Run: Rscript --vanilla R/experimental/rlmt_ca_ps_v1_retro_runner.R REF_DIR
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_ca_ps_v1.R")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("RETRO_STOP: reference repo argument required")
ref <- normalizePath(args[[1L]],mustWork=TRUE)
old <- file.path(ref,"results/rlmt_dc_published_data_5fold_v1")
out <- "results/rlmt_ca_ps_v1_retrospective_19x5"
must <- function(p,why) if(!isTRUE(p)) stop("RETRO_STOP: ",why)
must(dir.exists(out) && !file.exists(file.path(out,"run_manifest.csv")),
     "new experiment root missing or existing manifest")
must(trimws(readLines(file.path(out,"CODE_SHA.txt"))[1L])==
       system("git rev-parse HEAD",intern=TRUE)[1L],"code SHA changed")
must(trimws(readLines(file.path(old,"status.txt"))[1L])=="RUNNER_EXIT=0",
     "original reference not complete")
inputs <- read.csv(file.path(old,"input_manifest.csv"),
                   stringsAsFactors=FALSE)
orig <- read.csv(file.path(old,"run_manifest.csv"),
                 stringsAsFactors=FALSE)
must(nrow(inputs)==19L && nrow(orig)==285L &&
     all(orig$status=="OK"),"reference dimensions/status")
manifest <- list()
idx <- 0L
for(i in seq_len(nrow(inputs))) {
  row <- inputs[i,,drop=FALSE]
  ds <- as.character(row$dataset)
  world <- as.character(row$world)
  src <- as.character(row$source_path)
  digest <- suppressWarnings(system2("sha256sum",
    shQuote(normalizePath(src,mustWork=TRUE)),stdout=TRUE))
  must(length(digest)==1L &&
       substr(digest[[1L]],1L,64L)==as.character(row$sha256),
       paste("source hash differs",world,ds))
  xy <- read_xy_file(src)
  must(length(levels(xy$y))==2L &&
       !anyDuplicated(rownames(xy$X)) &&
       all(vapply(xy$X,function(x)
          is.numeric(x)&&all(is.finite(x)),logical(1))),
       paste("binary numeric cohort schema",ds))
  folds <- make_repeated_stratified_folds(
    xy$y,n_folds=5L,n_repeats=1L,seed=42L)
  must(length(folds)==5L &&
       identical(vapply(folds,function(z)z$split_id,character(1)),
                 sprintf("rep00_fold%02d",0:4)),
       paste("frozen 5-fold split mismatch",ds))
  for(f in folds) {
    key <- orig$world==world & orig$dataset==ds &
           orig$split_id==f$split_id & orig$method=="RLR"
    must(sum(key)==1L,"original RLR fold not unique")
    original_dir <- file.path(ref,as.character(orig$model_dir[key]))
    old_fold <- dirname(original_dir)
    old_train <- read.csv(file.path(old_fold,"samples_train.csv"),
                          stringsAsFactors=FALSE,
                          colClasses=c(sample_id="character"))
    old_test <- read.csv(file.path(old_fold,"samples_test.csv"),
                         stringsAsFactors=FALSE,
                         colClasses=c(sample_id="character"))
    oldpred <- read.csv(file.path(original_dir,"predictions_test.csv"),
                        stringsAsFactors=FALSE,
                        colClasses=c(sample_id="character",
                                     y_true="character",y_pred="character"))
    tr <- f$train_idx
    te <- f$test_idx
    Xtr <- xy$X[tr,,drop=FALSE]
    Xte <- xy$X[te,,drop=FALSE]
    ytr <- xy$y[tr]
    yte <- xy$y[te]
    must(identical(sort(c(tr,te)),seq_len(nrow(xy$X))) &&
         length(intersect(tr,te))==0L &&
         identical(old_train$sample_id,rownames(Xtr)) &&
         identical(old_test$sample_id,rownames(Xte)) &&
         identical(oldpred$sample_id,rownames(Xte)) &&
         identical(oldpred$y_true,as.character(yte)),
         paste("sample/truth or split mismatch",ds,f$split_id))
    dir_new <- file.path(out,"worlds",world,"datasets",ds,
                         "splits",f$split_id,"RLMT_CA_PS_v1")
    must(!dir.exists(dir_new),"refusing to overwrite new fitted fold")
    pool <- compute_top_pairs_switchbox(
      Xtr,ytr,top_n=1000L,buffer_mult=5L,
      out_dir=NULL,save_all=FALSE)
    must(nrow(pool)>=10L,"training pool too small")
    fit <- rlmt_ca_ps_v1_run(
      Xtr,ytr,Xte,yte,outer_training_pool=pool,
      out_dir=dir_new,config=list(),extra_metrics=list(
        world_name=world,dataset_id=ds,
        split_id=f$split_id,
        analysis_stage="EXPOSED_RETROSPECTIVE_DEVELOPMENT"))
    parent <- read.csv(file.path(dir_new,"rlr_reference",
                                 "predictions_test.csv"),
                       stringsAsFactors=FALSE,
                       colClasses=c(sample_id="character",
                                    y_true="character",y_pred="character"))
    must(identical(oldpred$sample_id,parent$sample_id) &&
         identical(oldpred$y_true,parent$y_true) &&
         identical(oldpred$y_pred,parent$y_pred) &&
         isTRUE(all.equal(oldpred$prob,parent$prob,
                          tolerance=1e-12)),
         paste("new parent != FROZEN original RLR",ds,f$split_id))
    pred <- read.csv(file.path(dir_new,"predictions_test.csv"),
                     stringsAsFactors=FALSE,
                     colClasses=c(sample_id="character",y_true="character"))
    must(identical(pred$sample_id,oldpred$sample_id) &&
         identical(pred$y_true,oldpred$y_true),
         paste("CA predictions don't match old heldout",ds,f$split_id))
    if(!fit$accepted_split) must(
      identical(as.numeric(pred$prob),as.numeric(parent$prob)),
      "no-split original RLR exact parity")
    idx <- idx+1L
    manifest[[idx]] <- data.frame(
      world=world,dataset=ds,split_id=f$split_id,
      method="RLMT_CA_PS_V1_K10_m10_d2_pool1000_L2",
      status="OK",n_train=length(tr),n_test=length(te),
      accepted_split=fit$accepted_split,
      realized_leaves=fit$realized_leaves,
      accuracy=fit$metrics$accuracy,
      balanced_accuracy=fit$metrics$balanced_accuracy,
      macro_f1=fit$metrics$macro_f1,
      logloss=fit$metrics$logloss,
      brier=fit$metrics$brier,
      model_dir=dir_new,
      original_RLR_model_dir=as.character(orig$model_dir[key]),
      original_parent_parity=TRUE,
      stringsAsFactors=FALSE)
    write.csv(do.call(rbind,manifest),
              file.path(out,"run_manifest.csv"),row.names=FALSE)
    message("RETRO_FOLD=",idx,"/95 ",ds," ",f$split_id,
            " accepted=",fit$accepted_split)
  }
}
must(idx==95L,"expected all 95 new matched fits")
cat("RETRO_CA_MODELS=95\n")
cat("RETRO_CA_ACCEPTED=",sum(vapply(manifest,
    function(x)isTRUE(x$accepted_split[[1L]]),logical(1))),"/95\n",sep="")
cat("RLMT_CA_PS_RETRO_19X5_RUN=PASS\n")
cat("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY\n")
