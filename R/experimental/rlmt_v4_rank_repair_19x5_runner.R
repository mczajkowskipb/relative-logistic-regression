# RLMT-v4: full 19 published-data EXPOSED cohorts x 5 original folds.
# Frozen original archived RLR is NEVER rewritten.
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v4_rank_repair.R")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("RLMT4_19X5_STOP: original REF required")
REF <- normalizePath(args[[1L]],mustWork=TRUE)
OLD <- file.path(REF,"results/rlmt_dc_published_data_5fold_v1")
OUT <- "results/rlmt_v4_rank_repair_19x5_v1"
assert <- function(value,reason)if(!isTRUE(value))
  stop("RLMT4_19X5_STOP: ",reason,call.=FALSE)
assert(dir.exists(OUT)&&file.exists(file.path(OUT,"CODE_SHA.txt")) &&
       !file.exists(file.path(OUT,"run_manifest.csv")),
       "new output missing or already started")
assert(identical(trimws(readLines(file.path(OUT,"CODE_SHA.txt"))[[1L]]),
  system("git rev-parse HEAD",intern=TRUE)[[1L]]),
  "code SHA mismatch")
assert(trimws(readLines(file.path(OLD,"status.txt"))[[1L]])==
       "RUNNER_EXIT=0","old original model source not complete")
orig <- read.csv(file.path(OLD,"run_manifest.csv"),
                  stringsAsFactors=FALSE)
inputs <- read.csv(file.path(OLD,"input_manifest.csv"),
                   stringsAsFactors=FALSE)
assert(nrow(orig)==285L&&nrow(inputs)==19L,
       "old reference dimensions")
sha256 <- function(path) {
 digest <- suppressWarnings(system2("sha256sum",
   shQuote(normalizePath(path,mustWork=TRUE)),stdout=TRUE))
 assert(length(digest)==1L &&
        grepl("^[a-f0-9]{64}[[:space:]]",digest[[1L]]),
        "hash command failed")
 substr(digest[[1L]],1L,64L)
}
rows <- list()
idx <- 0L
for(i in seq_len(nrow(inputs))) {
  input <- inputs[i,,drop=FALSE]
  dataset <- as.character(input$dataset[[1L]])
  assert(nrow(input)==1L,"missing original cohort input")
  assert(identical(tolower(sha256(input$source_path[[1L]])),
                   tolower(input$sha256[[1L]])),
         "source SHA mismatch")
  xy <- read_xy_file(input$source_path[[1L]])
  rlmt_v2_check_xy <- function(X,y) {
    assert(is.data.frame(X)&&is.factor(y) &&
           length(levels(y))==2L && nrow(X)==length(y),
           "invalid original cohort schema")
  }
  rlmt_v2_check_xy(xy$X,xy$y)
  folds <- make_repeated_stratified_folds(
    xy$y,n_folds=5L,n_repeats=1L,seed=42L)
  assert(identical(vapply(folds,function(f)f$split_id,character(1)),
                   sprintf("rep00_fold%02d",0:4)),
         "fold ids changed")
  for(f in folds) {
    frozen <- orig[orig$world==input$world[[1L]] &
      orig$dataset==dataset & orig$split_id==f$split_id &
      orig$method=="RLR",,drop=FALSE]
    assert(nrow(frozen)==1L &&
           frozen$status[[1L]]=="OK",
           "original frozen RLR fold missing")
    refdir <- file.path(REF,frozen$model_dir[[1L]])
    oldtrain <- read.csv(file.path(dirname(refdir),
                       "samples_train.csv"),
                         colClasses=c(sample_id="character"),
                         stringsAsFactors=FALSE)
    oldpred <- read.csv(file.path(refdir,"predictions_test.csv"),
        colClasses=c(sample_id="character",
                     y_true="character",y_pred="character"),
        stringsAsFactors=FALSE)
    train <- f$train_idx
    test <- f$test_idx
    assert(length(intersect(train,test))==0L &&
           identical(sort(c(train,test)),seq_len(nrow(xy$X))) &&
           identical(as.character(oldtrain$sample_id),
                     rownames(xy$X)[train]) &&
           identical(as.character(oldpred$sample_id),
                     rownames(xy$X)[test]) &&
           identical(as.character(oldpred$y_true),
                     as.character(xy$y[test])),
           "old patient train/test parity changed")
    X <- xy$X[train,,drop=FALSE]
    y <- xy$y[train]
    T <- xy$X[test,,drop=FALSE]
    yT <- xy$y[test]
    path <- file.path(OUT,"models",dataset,f$split_id)
    assert(!dir.exists(path),"attempted rerun/overwrite")
    message("RLMT4_FIT=",dataset,":",f$split_id,
            " original RLR + strict 3-fold OOF")
    warnings <- character()
    fit <- withCallingHandlers(rlmt4_fit(X,y),
      warning=function(w) {
        warnings <<- c(warnings,conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    prob <- rlmt4_predict(fit,T)
    oldp <- as.numeric(oldpred$prob)
    parity_err <- max(abs(prob$ORIGINAL_RLR-oldp))
    # HARD GATE. If original root changes, no experimental gain is
    # scientifically interpretable. Stop rather than relabeling parent.
    assert(is.finite(parity_err) && parity_err<=1e-9,
       paste("ORIGINAL RLR numerical parity FAIL",dataset,
             f$split_id,sprintf("max_delta=%g",parity_err)))
    dir.create(path,recursive=TRUE,showWarnings=FALSE)
    saveRDS(fit,file.path(path,"all_fitted_methods.rds"))
    y01 <- rlmt_y01(yT)
    old_ll <- rlmt_logloss(y01,oldp)
    for(method in rlmt4_modes) {
      pr <- prob[[method]]
      dest <- file.path(path,method)
      dir.create(dest,recursive=TRUE,showWarnings=FALSE)
      write.csv(data.frame(
        sample_id=rownames(T),y_true=as.character(yT),
        y_pred=ifelse(pr>=.5,levels(yT)[2L],levels(yT)[1L]),
        prob=pr,original_RLR_prob=prob$ORIGINAL_RLR),
        file.path(dest,"predictions_test.csv"),row.names=FALSE)
      idx <- idx+1L
      r <- fit$routes[[method]]
      pair <- if(is.null(r))NULL else r$pair
      metric <- rlmt_logloss(y01,pr)
      rows[[idx]] <- data.frame(
        world=input$world[[1L]],dataset=dataset,
        split_id=f$split_id,method=method,
        status="OK",n_train=length(train),n_test=length(test),
        logloss=metric,original_logloss=old_ll,
        delta_vs_original=metric-old_ll,
        brier=mean((pr-y01)^2),
        route_n=as.integer(!is.null(pair)),
        gene_a=if(is.null(pair))"" else pair$gene_a,
        gene_b=if(is.null(pair))"" else pair$gene_b,
        gamma=if(is.null(r))0 else r$coefficients$gamma,
        warning_count=length(warnings),parity_max_abs=parity_err,
        model_dir=dest,original_model_dir=frozen$model_dir[[1L]],
        scientific_stage="EXPOSED_RETROSPECTIVE_DEV_ONLY",
        stringsAsFactors=FALSE)
      write.csv(do.call(rbind,rows),file.path(OUT,"run_manifest.csv"),
                row.names=FALSE)
      message("RLMT4_19X5_COMPLETED=",idx,"/475 ",
              dataset,":",f$split_id,":",method,
              " delta_original=",sprintf("%+.7f",metric-old_ll),
              " real_route=",!is.null(pair))
    }
    message("RLMT4_FOLD_PARENT_PARITY=PASS ",
            dataset,":",f$split_id,
            " max_abs=",sprintf("%.3g",parity_err),
            " warnings=",length(warnings))
  }
}
assert(idx==475L,"not all 475 experimental predictions")
message("RLMT4_19X5_475_FOLD_METHODS=PASS")
message("RLMT4_19X5_ORIGINAL_RLR_PARITY=PASS")
message("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY")
