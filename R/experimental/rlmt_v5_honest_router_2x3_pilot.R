# RLMT-v5: six already-exposed real outer fold positions, no tuning
# or labels transferred from any original heldout patients.
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v4_rank_repair.R")
source("R/experimental/rlmt_v5_honest_router.R")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)
  stop("RLMT5_PILOT_STOP: frozen ORIGINAL_RLR REF required")
REF <- normalizePath(args[[1L]],mustWork=TRUE)
OLD <- file.path(REF,"results/rlmt_dc_published_data_5fold_v1")
OUT <- "results/rlmt_v5_honest_router_2x3_v1"
require <- function(ok,msg)
  if(!isTRUE(ok))stop("RLMT5_PILOT_STOP: ",msg,call.=FALSE)
require(dir.exists(OUT) &&
        file.exists(file.path(OUT,"CODE_SHA.txt")) &&
        !file.exists(file.path(OUT,"run_manifest.csv")),
        "output missing or prior run exists")
require(identical(trimws(readLines(file.path(OUT,"CODE_SHA.txt"))[1L]),
                  system("git rev-parse HEAD",intern=TRUE)[1L]),
        "source SHA mismatch")
require(identical(trimws(readLines(file.path(OLD,"status.txt"))[1L]),
                  "RUNNER_EXIT=0"),"original RLR archive not complete")
archived <- read.csv(file.path(OLD,"run_manifest.csv"),
       stringsAsFactors=FALSE)
input <- read.csv(file.path(OLD,"input_manifest.csv"),
       stringsAsFactors=FALSE)
require(nrow(archived)==285L && nrow(input)==19L,
        "frozen original inventory changed")
sha <- function(path) {
  d <- suppressWarnings(system2("sha256sum",
     shQuote(normalizePath(path,mustWork=TRUE)),stdout=TRUE))
  require(length(d)==1L &&
          grepl("^[a-f0-9]{64}[[:space:]]",d[[1L]]),
          "sha256sum command failed")
  substr(d[[1L]],1L,64L)
}
rows <- list()
count <- 0L
for(dataset in c("chen-2002","gordon-2002")) {
 item <- input[input$dataset==dataset,,drop=FALSE]
 require(nrow(item)==1L,paste("cohort missing",dataset))
 require(identical(tolower(sha(item$source_path[[1L]])),
           tolower(item$sha256[[1L]])),
         paste("cohort hash mismatch",dataset))
 xy <- read_xy_file(item$source_path[[1L]])
 require(is.factor(xy$y) && length(levels(xy$y))==2L,
         "input y must be binary factor")
 folds <- make_repeated_stratified_folds(
           xy$y,n_folds=5L,n_repeats=1L,seed=42L)
 require(length(folds)==5L &&
         identical(vapply(folds,function(f)f$split_id,character(1)),
               sprintf("rep00_fold%02d",0:4)),
         "frozen original seed/fold arrangement differs")
 for(f in folds[1:3]) {
   reference <- archived[archived$world==item$world[[1L]] &
     archived$dataset==dataset & archived$split_id==f$split_id &
     archived$method=="RLR",,drop=FALSE]
   require(nrow(reference)==1L &&
       reference$status[[1L]]=="OK",
       "original matched RLR missing")
   old_dir <- file.path(REF,reference$model_dir[[1L]])
   oldtrain <- read.csv(file.path(dirname(old_dir),"samples_train.csv"),
               stringsAsFactors=FALSE,colClasses=c(sample_id="character"))
   oldpred <- read.csv(file.path(old_dir,"predictions_test.csv"),
               stringsAsFactors=FALSE,
               colClasses=c(sample_id="character",y_true="character"))
   train <- f$train_idx
   test <- f$test_idx
   X <- xy$X[train,,drop=FALSE]
   y <- xy$y[train]
   T <- xy$X[test,,drop=FALSE]
   Y <- xy$y[test]
   require(length(intersect(train,test))==0L &&
         identical(sort(c(train,test)),seq_len(nrow(xy$X))) &&
         identical(as.character(oldtrain$sample_id),rownames(X)) &&
         identical(as.character(oldpred$sample_id),rownames(T)) &&
         identical(as.character(oldpred$y_true),as.character(Y)),
         "archived original patient split and labels changed")
   dest <- file.path(OUT,"models",dataset,f$split_id)
   require(!dir.exists(dest),"would overwrite previous v5 fold")
   message("RLMT5_TRAIN_START=",dataset,":",f$split_id,
       " fit 3 inner-disjoint v4 folds for routing + final original RLR")
   warning_messages <- character()
   fitted <- withCallingHandlers(rlmt5_fit(X,y),
     warning=function(w) {
       warning_messages <<- c(warning_messages,conditionMessage(w))
       invokeRestart("muffleWarning")
     })
   predicted <- rlmt5_predict(fitted,T)
   original <- as.numeric(oldpred$prob)
   parity <- max(abs(predicted$original_prob-original))
   require(is.finite(parity) && parity<1e-9,
      paste("ARCHIVED ORIGINAL RLR PARITY FAILURE",dataset,
            f$split_id,sprintf("%.12g",parity)))
   require(length(predicted$prob)==nrow(T) &&
      all(predicted$prob>=0&predicted$prob<=1),
      "v5 test probabilities invalid")
   dir.create(dest,recursive=TRUE,showWarnings=FALSE)
   saveRDS(fitted,file.path(dest,"complete_rlmt5.rds"))
   write.csv(data.frame(
      sample_id=rownames(T),y_true=as.character(Y),
      y_pred=ifelse(predicted$prob>=.5,levels(Y)[2L],levels(Y)[1L]),
      prob=predicted$prob,original_RLR_prob=predicted$original_prob,
      selected_method=predicted$model,selected_weight=predicted$weight),
      file.path(dest,"predictions_test.csv"),row.names=FALSE)
   g <- fitted$router
   evidence <- if(is.null(g$evidence)) data.frame() else g$evidence
   write.csv(evidence,
      file.path(dest,"inner_validation_evidence.csv"),
      row.names=FALSE)
   y01 <- rlmt_y01(Y)
   delta <- rlmt_logloss(y01,predicted$prob)-
       rlmt_logloss(y01,original)
   count <- count+1L
   row <- data.frame(
      world=item$world[[1L]],dataset=dataset,
      split_id=f$split_id,status="OK",
      method="RLMT5_HONEST_ROUTER",
      n_train=length(train),n_test=length(test),
      selected_method=predicted$model,
      selected_weight=predicted$weight,
      accepted=predicted$weight>0,
      route_n=as.integer(!is.null(predicted$active_pair) &&
                        predicted$weight>0),
      gene_a=if(is.null(predicted$active_pair))"" else
           predicted$active_pair$gene_a,
      gene_b=if(is.null(predicted$active_pair))"" else
           predicted$active_pair$gene_b,
      logloss=rlmt_logloss(y01,predicted$prob),
      original_logloss=rlmt_logloss(y01,original),
      delta_vs_original=delta,
      brier=mean((predicted$prob-y01)^2),
      original_brier=mean((original-y01)^2),
      max_parent_parity=parity,
      shared_warning_count=length(warning_messages),
      model_dir=dest,original_model_dir=reference$model_dir[[1L]],
      scientific_stage="EXPOSED_RETROSPECTIVE_DEV_ONLY",
      stringsAsFactors=FALSE)
   rows[[count]] <- row
   write.csv(do.call(rbind,rows),file.path(OUT,"run_manifest.csv"),
             row.names=FALSE)
   message("RLMT5_PILOT_COMPLETED=",count,"/6 ",
           dataset,":",f$split_id,
           " chosen=",predicted$model,
           " weight=",sprintf("%.3f",predicted$weight),
           " delta_vs_original=",sprintf("%+.9f",delta),
           " exact_RLR_parity=",format(parity,digits=3))
 }
}
require(count==6L,"must evaluate all six frozen pilot outer folds")
message("RLMT5_PILOT_REAL_OUTER_FOLDS=6")
message("RLMT5_PILOT_ORIGINAL_RLR_PARITY=PASS")
message("RLMT5_PILOT_ALL_TRAINING=PASS")
message("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY")
