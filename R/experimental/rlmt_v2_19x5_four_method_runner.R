# RLMT-v2 all 19 previously exposed biomedical source cohorts × original
# five outer folds × four PREDECLARED algorithms = 380 attempted NEW fits.
# Old RLR/RLMT-v1/RLMT-DC models: immutable, READ ONLY, NEVER retrained.
# Identical method constants as completed paired two-cohort pilot.
# No use of outer test outcomes during training, selection or prediction.
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v2_aligned_offset.R")
source("R/experimental/rlmt_v2_gini_lookahead.R")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("V2_19X5_STOP: original frozen REF required")
REF <- normalizePath(args[[1L]],mustWork=TRUE)
OLD <- file.path(REF,"results","rlmt_dc_published_data_5fold_v1")
OUT <- "results/rlmt_v2_19x5_four_method_v1"
METHODS <- c("H1","H2","TSPDT_GINI","RLMT_LA")
check <- function(ok,context) {
  if(!isTRUE(ok))stop("V2_19X5_STOP: ",context,call.=FALSE)
}
check(dir.exists(OUT) &&
      file.exists(file.path(OUT,"CODE_SHA.txt")) &&
      identical(trimws(readLines(file.path(OUT,"CODE_SHA.txt"))[[1L]]),
                system("git rev-parse HEAD",intern=TRUE)[[1L]]),
      "new source SHA marker not matched")
check(identical(readLines(file.path(OLD,"status.txt"))[[1L]],
                "RUNNER_EXIT=0"),"original baseline incomplete")
check(!file.exists(file.path(OUT,"run_manifest.csv")),
      "prior results manifest exists, no overwriting")
original <- read.csv(file.path(OLD,"run_manifest.csv"),
                     stringsAsFactors=FALSE)
inputs <- read.csv(file.path(OLD,"input_manifest.csv"),
                   stringsAsFactors=FALSE)
check(nrow(original)==285L && nrow(inputs)==19L,
      "original fold input/model inventory")
check(!anyDuplicated(paste(original$world,original$dataset,
                           original$split_id,original$method)),
      "original duplicate model records")
hash <- function(f) {
  answer <- suppressWarnings(system2(
    "sha256sum",shQuote(normalizePath(f,mustWork=TRUE)),
    stdout=TRUE,stderr=TRUE))
  check(length(answer)==1L &&
        grepl("^[a-f0-9]{64}[[:space:]]+",answer[[1L]]),
        "source hash utility failed")
  substr(answer[[1L]],1L,64L)
}
readpred <- function(f) {
  read.csv(file.path(f,"predictions_test.csv"),
           stringsAsFactors=FALSE,
           colClasses=c(sample_id="character",
                        y_true="character",y_pred="character"))
}
full <- list()
warnings_all <- list()
completed <- 0L
failed <- 0L

for(p in seq_len(nrow(inputs))) {
  entry <- inputs[p,,drop=FALSE]
  dataset <- entry$dataset[[1L]]
  world <- entry$world[[1L]]
  src <- entry$source_path[[1L]]
  check(identical(tolower(hash(src)),tolower(entry$sha256[[1L]])),
        paste("source SHA mismatch",dataset))
  xy <- read_xy_file(src)
  rlmt_v2_check_xy(xy$X,xy$y)
  folds <- make_repeated_stratified_folds(
    xy$y,n_folds=5L,n_repeats=1L,seed=42L)
  check(identical(vapply(folds,function(x)x$split_id,character(1)),
                  sprintf("rep00_fold%02d",0:4)),
        paste("original fold grid changed",dataset))
  for(f in folds) {
    s <- f$split_id
    legacy <- original[
      original$world==world &
      original$dataset==dataset &
      original$split_id==s &
      original$method=="RLR",,drop=FALSE]
    check(nrow(legacy)==1L && legacy$status[[1L]]=="OK",
          paste("missing preserved original RLR fit",dataset,s))
    oldmodel <- file.path(REF,legacy$model_dir[[1L]])
    oldp <- readpred(oldmodel)
    oldtrain <- read.csv(
      file.path(dirname(oldmodel),"samples_train.csv"),
      stringsAsFactors=FALSE,
      colClasses=c(sample_id="character"))
    tr <- f$train_idx
    te <- f$test_idx
    check(identical(sort(c(tr,te)),seq_len(nrow(xy$X))) &&
          !length(intersect(tr,te)) &&
          identical(rownames(xy$X)[tr],oldtrain$sample_id) &&
          identical(rownames(xy$X)[te],oldp$sample_id) &&
          identical(as.character(xy$y[te]),oldp$y_true),
          paste("EXACT historical patient/fold parity mismatch",dataset,s))
    original_prob <- as.numeric(oldp$prob)
    check(length(original_prob)==length(te) &&
          all(is.finite(original_prob)) &&
          all(original_prob>=0 & original_prob<=1),
          "archived original RLR probabilities invalid")
    XX <- xy$X[tr,,drop=FALSE]
    yy <- xy$y[tr]
    TT <- xy$X[te,,drop=FALSE]
    yt <- xy$y[te]
    aligned <- NULL
    for(method in METHODS) {
      dest <- file.path(OUT,"models",dataset,s,method)
      check(!dir.exists(dest),"model already exists; never overwrite")
      message("FOURMETHOD_FIT=",p,"/19; ",dataset,":",s,":",method)
      note <- character()
      # Keep every warning message/count in a separate audit CSV.
      # Only suppress verbose console repetition after capturing warning.
      outcome <- tryCatch(
        withCallingHandlers({
          if(method %in% c("H1","H2")) {
            rlmt_v2_run(XX,yy,TT,yt,out_dir=dest,mode=method)
          } else {
            rlmt_la_run(XX,yy,TT,yt,out_dir=dest,method=method)
          }
        },warning=function(w) {
          note <<- c(note,conditionMessage(w))
          invokeRestart("muffleWarning")
        }),error=function(e)e)
      if(length(note)) {
        freq <- as.data.frame(table(note),stringsAsFactors=FALSE)
        names(freq) <- c("message","frequency")
        freq$world <- world
        freq$dataset <- dataset
        freq$split_id <- s
        freq$method <- method
        warnings_all[[length(warnings_all)+1L]] <- freq[
          ,c("world","dataset","split_id","method",
             "message","frequency"),drop=FALSE]
        write.csv(do.call(rbind,warnings_all),
                  file.path(OUT,"warning_inventory.csv"),
                  row.names=FALSE)
      }
      success <- !inherits(outcome,"error")
      if(success) {
        pred <- readpred(dest)
        check(identical(pred$sample_id,oldp$sample_id) &&
              identical(pred$y_true,oldp$y_true) &&
              isTRUE(all.equal(as.numeric(pred$prob),
                               outcome$prob,tolerance=1e-12)) &&
              isTRUE(all.equal(as.numeric(pred$parent_prob),
                               outcome$parent,tolerance=1e-12)),
              paste("model predicted sample mismatch",dataset,s,method))
        if(is.null(aligned))aligned <- outcome$parent
        else check(isTRUE(all.equal(aligned,outcome$parent,
                                    tolerance=1e-12)),
                   paste("unaligned parent procedure",dataset,s,method))
        accepted <- if(method %in% c("H1","H2"))
          isTRUE(outcome$decision$accepted) else
          identical(outcome$decision$status,"ACCEPTED")
        if(!accepted)
          check(identical(outcome$prob,outcome$parent),
                "rejected fit changes aligned parent probabilities")
        y01 <- rlmt_y01(yt)
        mdl_loss <- rlmt_logloss(y01,outcome$prob)
        par_loss <- rlmt_logloss(y01,outcome$parent)
        old_loss <- rlmt_logloss(y01,original_prob)
        row <- data.frame(
          world=world,dataset=dataset,split_id=s,mode=method,
          status="OK",error="",accepted=accepted,
          n_train=length(tr),n_test=length(te),
          train_class0=sum(yy==levels(yy)[1L]),
          train_class1=sum(yy==levels(yy)[2L]),
          new_logloss=mdl_loss,aligned_parent_logloss=par_loss,
          original_RLR_logloss=old_loss,
          delta_vs_aligned=mdl_loss-par_loss,
          delta_vs_original=mdl_loss-old_loss,
          new_brier=mean((outcome$prob-y01)^2),
          aligned_brier=mean((outcome$parent-y01)^2),
          original_RLR_brier=mean((original_prob-y01)^2),
          warning_total=length(note),
          warning_low_class=sum(grepl("fewer than 8",note)),
          model_dir=dest,original_RLR_dir=legacy$model_dir[[1L]],
          stringsAsFactors=FALSE)
      } else {
        failed <- failed+1L
        row <- data.frame(
          world=world,dataset=dataset,split_id=s,mode=method,
          status="FIT_ERROR",error=conditionMessage(outcome),
          accepted=NA,
          n_train=length(tr),n_test=length(te),
          train_class0=sum(yy==levels(yy)[1L]),
          train_class1=sum(yy==levels(yy)[2L]),
          new_logloss=NA_real_,aligned_parent_logloss=NA_real_,
          original_RLR_logloss=NA_real_,
          delta_vs_aligned=NA_real_,delta_vs_original=NA_real_,
          new_brier=NA_real_,aligned_brier=NA_real_,
          original_RLR_brier=NA_real_,
          warning_total=length(note),
          warning_low_class=sum(grepl("fewer than 8",note)),
          model_dir=dest,original_RLR_dir=legacy$model_dir[[1L]],
          stringsAsFactors=FALSE)
      }
      completed <- completed+1L
      full[[completed]] <- row
      write.csv(do.call(rbind,full),
                file.path(OUT,"run_manifest.csv"),row.names=FALSE)
      message("FOURMETHOD_COMPLETE=",completed,"/380 ",
              dataset,":",s,":",method,
              " status=",row$status,
              " accepted=",row$accepted,
              " warnings=",length(note),
              " low_class_warnings=",row$warning_low_class,
              if(success)
                paste0(" deltaLL=",sprintf("%+.7f",row$delta_vs_aligned))
              else paste0(" ERROR=",row$error))
    }
  }
}
check(completed==380L,"expected 380 attempted fold-method fits")
cat("RLMT_V2_19X5_ATTEMPTED=380\n")
cat("RLMT_V2_19X5_FAILED=",failed,"\n",sep="")
cat("RLMT_V2_19X5_ALL_FOLDS_RUN=TRUE\n")
if(failed>0L) {
  cat("RLMT_V2_19X5_RUN=INCOMPLETE_WITH_FAILURES\n")
  stop("V2_19X5_STOP: preserve failure manifest; no rerun/overwrite")
}
cat("RLMT_V2_19X5_RUN=PASS\n")
cat("PREVIOUSLY_EXPOSED_DATA=TRUE\n")
