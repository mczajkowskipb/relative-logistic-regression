# RLMT-v3B: 19 EXPOSED datasets, exact original 5 folds per dataset.
# Four methods are generated from ONE training-only parent and OOF
# to compare route-specific ridge offsets to a calibrated intercept.
# Frozen original 285 RLR/RLMT results: READ ONLY.
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v2_aligned_offset.R")
source("R/experimental/rlmt_v3_stumps_pruned.R")
source("R/experimental/rlmt_v3b_calibration_control.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("V3B_STOP: frozen original REF required")
REF <- normalizePath(args[[1L]],mustWork=TRUE)
OLD <- file.path(REF,"results","rlmt_dc_published_data_5fold_v1")
OUT <- "results/rlmt_v3b_calibration_19x5_v1"
need <- function(x,msg) if(!isTRUE(x))
  stop("V3B_STOP: ",msg,call.=FALSE)
need(dir.exists(OUT) &&
     file.exists(file.path(OUT,"CODE_SHA.txt")) &&
     identical(trimws(readLines(file.path(OUT,"CODE_SHA.txt"))[[1L]]),
               system("git rev-parse HEAD",intern=TRUE)[[1L]]),
     "source SHA drifted")
need(identical(trimws(readLines(file.path(OLD,"status.txt"))[[1L]]),
               "RUNNER_EXIT=0"),"archived 285 source models incomplete")
need(!file.exists(file.path(OUT,"run_manifest.csv")),
     "old completed/partial manifest must not be overwritten")
original <- read.csv(file.path(OLD,"run_manifest.csv"),
                     stringsAsFactors=FALSE)
inputs <- read.csv(file.path(OLD,"input_manifest.csv"),
                   stringsAsFactors=FALSE)
need(nrow(original)==285L && nrow(inputs)==19L,
     "frozen archived inventory counts diverged")
sha256 <- function(path) {
  output <- suppressWarnings(system2(
    "sha256sum",
    shQuote(normalizePath(path,mustWork=TRUE)),
    stdout=TRUE,stderr=TRUE))
  need(length(output)==1L &&
       grepl("^[a-f0-9]{64}[[:space:]]",output[[1L]]),
       paste("hash failed",path))
  substr(output[[1L]],1L,64L)
}
archival_predict <- function(dir)
  read.csv(file.path(dir,"predictions_test.csv"),
           stringsAsFactors=FALSE,
           colClasses=c(sample_id="character",
                        y_true="character",y_pred="character"))
rows <- list()
n <- 0L
warning_table <- list()
warning_rows <- 0L
for(i in seq_len(nrow(inputs))) {
  entry <- inputs[i,,drop=FALSE]
  world <- entry$world[[1L]]
  dataset <- entry$dataset[[1L]]
  source_file <- entry$source_path[[1L]]
  need(identical(tolower(sha256(source_file)),
                 tolower(entry$sha256[[1L]])),
       paste("source input hash drift",dataset))
  cohort <- read_xy_file(source_file)
  rlmt_v2_check_xy(cohort$X,cohort$y)
  folds <- make_repeated_stratified_folds(
       cohort$y,n_folds=5L,n_repeats=1L,seed=42L)
  need(identical(vapply(folds,function(f)f$split_id,character(1)),
                 sprintf("rep00_fold%02d",0:4)),
       "old fold generator changed")
  for(f in folds) {
    split <- f$split_id
    baseline <- original[
       original$world==world &
       original$dataset==dataset &
       original$split_id==split &
       original$method=="RLR",,drop=FALSE]
    need(nrow(baseline)==1L &&
         baseline$status[[1L]]=="OK",
         paste("original baseline fold missing",dataset,split))
    archivedir <- file.path(REF,baseline$model_dir[[1L]])
    archived_pred <- archival_predict(archivedir)
    archived_train <- read.csv(file.path(
         dirname(archivedir),"samples_train.csv"),
         stringsAsFactors=FALSE,
         colClasses=c(sample_id="character"))
    tr <- f$train_idx
    te <- f$test_idx
    need(!length(intersect(tr,te)) &&
         identical(sort(c(tr,te)),seq_len(nrow(cohort$X))) &&
         identical(archived_train$sample_id,
                   rownames(cohort$X)[tr]) &&
         identical(archived_pred$sample_id,
                   rownames(cohort$X)[te]) &&
         identical(archived_pred$y_true,
                   as.character(cohort$y[te])),
         paste("original heldout fold or patient mismatch",dataset,split))
    Xtr <- cohort$X[tr,,drop=FALSE]
    ytr <- cohort$y[tr]
    Xte <- cohort$X[te,,drop=FALSE]
    yte <- cohort$y[te]
    oprob <- as.numeric(archived_pred$prob)
    need(length(oprob)==length(te) &&
         all(is.finite(oprob)) &&
         all(oprob>=0&oprob<=1),
         "original frozen RLR probabilities corrupt")
    basepath <- file.path(OUT,"models",dataset,split)
    need(!dir.exists(basepath),"duplicate output or prior fit root exists")
    message("V3B_TRAINING=",dataset,":",split,
            "; ",i,"/19; ONE_SHARED_RLR_PARENT")
    warn <- character()
    fit <- withCallingHandlers(
      rlmt_v3b_fit_bundle(Xtr,ytr),
      warning=function(w) {
        warn <<- c(warn,conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    predictions <- rlmt_v3b_predict_bundle(fit,Xte)
    need(identical(names(predictions$prob),rlmt_v3b_modes),
         "missing new method probabilities")
    ensure_dir(basepath)
    saveRDS(fit,file.path(basepath,"shared_complete_bundle.rds"))
    if(length(warn)) {
      tab <- as.data.frame(table(warn),stringsAsFactors=FALSE)
      names(tab) <- c("message","frequency")
      tab$world <- world
      tab$dataset <- dataset
      tab$split_id <- split
      warning_rows <- warning_rows+1L
      warning_table[[warning_rows]] <-
         tab[,c("world","dataset","split_id",
                "message","frequency"),drop=FALSE]
      write.csv(do.call(rbind,warning_table),
        file.path(OUT,"shared_training_warning_inventory.csv"),
        row.names=FALSE)
    }
    y01 <- rlmt_y01(yte)
    parent <- predictions$parent
    aligned_loss <- rlmt_logloss(y01,parent)
    orig_loss <- rlmt_logloss(y01,oprob)
    controls <- list(GLOBAL_RIDGE_CAL=NULL,
                     STUMP_TSP_RIDGE=fit$TSP,
                     STUMP_GINI_RIDGE=fit$GINI,
                     STUMP_RESID_RIDGE=fit$RESID)
    for(method in rlmt_v3b_modes) {
      pr <- predictions$prob[[method]]
      need(length(pr)==length(te) &&
           all(is.finite(pr)) && all(pr>=0&pr<=1),
           "invalid heldout method probabilities")
      dest <- file.path(basepath,method)
      ensure_dir(dest)
      node <- controls[[method]]
      chosen <- if(is.null(node))NULL else node$candidate
      z <- if(is.null(chosen)) rep(FALSE,nrow(Xte)) else
         rlmt_v2_route(Xte,chosen)
      write.csv(data.frame(
         sample_id=rownames(Xte),y_true=as.character(yte),
         y_pred=ifelse(pr>=.5,levels(yte)[2L],levels(yte)[1L]),
         prob=pr,parent_prob=parent,
         leaf_id=if(is.null(chosen))rep(1L,nrow(Xte)) else
             ifelse(z,2L,3L)),
         file.path(dest,"predictions_test.csv"),row.names=FALSE)
      count <- as.integer(!is.null(chosen))
      effect <- 0
      if(!is.null(node)) {
        beta <- as.numeric(stats::coef(
             node$fit,s=node$lambda))
        need(length(beta)==3L,
             "expected intercept and TWO route coefficients")
        effect <- beta[[3L]]-beta[[2L]]
      }
      scored <- rlmt_logloss(y01,pr)
      n <- n+1L
      rows[[n]] <- data.frame(
          world=world,dataset=dataset,split_id=split,method=method,
          status="OK",n_train=length(tr),n_test=length(te),
          internal_nodes=count,route_effect=effect,
          gene_a=if(is.null(chosen))"" else
                 as.character(chosen$gene_a[[1L]]),
          gene_b=if(is.null(chosen))"" else
                 as.character(chosen$gene_b[[1L]]),
          calibration_offset_intercept=
               fit$calibrator$offset_intercept,
          scored_logloss=scored,
          parent_logloss=aligned_loss,
          original_RLR_logloss=orig_loss,
          delta_vs_parent=scored-aligned_loss,
          delta_vs_original=scored-orig_loss,
          brier=mean((pr-y01)^2),
          shared_parent_warnings=length(warn),
          shared_parent_low_class_warnings=
               sum(grepl("fewer than 8",warn)),
          model_dir=dest,
          original_model_dir=baseline$model_dir[[1L]],
          scientific_stage="EXPOSED_RETROSPECTIVE_DEV_ONLY",
          stringsAsFactors=FALSE)
      write.csv(do.call(rbind,rows),
          file.path(OUT,"run_manifest.csv"),row.names=FALSE)
      message("V3B_COMPLETE=",n,"/380 ",dataset,":",split,":",method,
              " nodes=",count,
              " delta_vs_original=",sprintf("%+.7f",scored-orig_loss),
              " delta_vs_parent=",sprintf("%+.7f",scored-aligned_loss))
    }
    message("V3B_FOLD_COMPLETE=",dataset,":",split,
            "; candidate_inventory=",fit$proposals$scanned_pairs,
            "; feasible=",fit$proposals$feasible_pairs,
            "; shared_warning_count=",length(warn))
  }
}
need(n==380L,"must have 380 19x5 method records")
cat("RLMT_V3B_19X5_REAL_TRAINING_BUNDLES=95\n")
cat("RLMT_V3B_19X5_METHOD_FOLD_RESULTS=380\n")
cat("RLMT_V3B_19X5_RUN=PASS\n")
cat("PUBLISHED_RLR_MODELS_REUSED_READ_ONLY=TRUE\n")
cat("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY\n")
