# CRIT-v8: frozen 2 datasets x 3 original outer folds x 6 methods.
# EXPOSED RETROSPECTIVE DEVELOPMENT. No new cohort validation.
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v4_rank_repair.R")
source("R/experimental/crit_v7_original_rlr_two_leaf.R")
source("R/experimental/crit_v8_single_nomination.R")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("CRIT9_PILOT_STOP: original RLR path required")
REF <- normalizePath(args[[1L]],mustWork=TRUE)
OLD <- file.path(REF,"results/rlmt_dc_published_data_5fold_v1")
OUT <- "results/crit_v9_d_crossfit_2x3_v1"
need <- function(ok,msg)
 if(!isTRUE(ok))stop("CRIT9_PILOT_STOP: ",msg,call.=FALSE)
need(dir.exists(OUT) && file.exists(file.path(OUT,"CODE_SHA.txt")) &&
     !file.exists(file.path(OUT,"run_manifest.csv")),
     "result dir missing or result previously written")
need(identical(trimws(readLines(file.path(OUT,"CODE_SHA.txt"))[[1L]]),
               system("git rev-parse HEAD",intern=TRUE)[[1L]]),
     "worktree SHA differs from prepared runner")
need(trimws(readLines(file.path(OLD,"status.txt"))[[1L]])==
     "RUNNER_EXIT=0","original RLR archive missing or incomplete")
orig <- read.csv(file.path(OLD,"run_manifest.csv"),stringsAsFactors=FALSE)
inputs <- read.csv(file.path(OLD,"input_manifest.csv"),stringsAsFactors=FALSE)
need(nrow(orig)==285L&&nrow(inputs)==19L,
     "original manifest dimensions unexpected")
hash <- function(path) {
 res <- suppressWarnings(system2("sha256sum",
     shQuote(normalizePath(path,mustWork=TRUE)),stdout=TRUE))
 need(length(res)==1L && grepl("^[a-f0-9]{64}[[:space:]]",res[[1L]]),
      "cannot hash source dataset")
 substr(res[[1L]],1L,64L)
}
rows<-list()
i<-0L
for(dataset in c("chen-2002","gordon-2002")) {
 inp<-inputs[inputs$dataset==dataset,,drop=FALSE]
 need(nrow(inp)==1L,"dataset absent in original input inventory")
 need(tolower(hash(inp$source_path[[1L]]))==
      tolower(inp$sha256[[1L]]),"source dataset hash differs")
 xy<-read_xy_file(inp$source_path[[1L]])
 need(is.factor(xy$y) && length(levels(xy$y))==2L,
      "invalid archived label factor")
 folds<-make_repeated_stratified_folds(
       xy$y,n_folds=5L,n_repeats=1L,seed=42L)
 need(length(folds)==5L &&
      identical(vapply(folds,function(z)z$split_id,character(1)),
                sprintf("rep00_fold%02d",0:4)),
      "original fold definitions differ")
 for(f in folds[1:3]) {
   match<-orig[orig$world==inp$world[[1L]] &
             orig$dataset==dataset & orig$split_id==f$split_id &
             orig$method=="RLR",,drop=FALSE]
   need(nrow(match)==1L && match$status[[1L]]=="OK",
        "original archived fold RLR row missing")
   od<-file.path(REF,match$model_dir[[1L]])
   oldtr<-read.csv(file.path(dirname(od),"samples_train.csv"),
      stringsAsFactors=FALSE,colClasses=c(sample_id="character"))
   oldp<-read.csv(file.path(od,"predictions_test.csv"),
      stringsAsFactors=FALSE,
      colClasses=c(sample_id="character",
                   y_true="character",y_pred="character"))
   train<-f$train_idx
   test<-f$test_idx
   Xt<-xy$X[train,,drop=FALSE]
   Yt<-xy$y[train]
   Xv<-xy$X[test,,drop=FALSE]
   Yv<-xy$y[test]
   need(length(intersect(train,test))==0L &&
        identical(sort(c(train,test)),seq_len(nrow(xy$X))) &&
        identical(as.character(oldtr$sample_id),rownames(Xt)) &&
        identical(as.character(oldp$sample_id),rownames(Xv)) &&
        identical(as.character(oldp$y_true),as.character(Yv)),
        "original sample IDs/classes/split mismatch")
   dest<-file.path(OUT,"models",dataset,f$split_id)
   need(!dir.exists(dest),"would overwrite frozen fold output")
   warn<-character()
   model<-withCallingHandlers(crit9_fit(Xt,Yt),warning=function(w){
          warn<<-c(warn,conditionMessage(w))
          invokeRestart("muffleWarning")
        })
   pred<-crit9_predict(model,Xv)
   oldroot<-as.numeric(oldp$prob)
   rooterr<-max(abs(pred$probs$ORIGINAL_RLR-oldroot))
   need(is.finite(rooterr)&&rooterr<1e-9,
        paste("exact frozen ORIGINAL_RLR prediction mismatch",
              dataset,f$split_id,rooterr))
   need(model$nomination_n<=6L,
        "D-only nominee shortlist exceeded six")
   if(!is.null(model$partition))
     need(!length(intersect(model$partition$discovery,
                            model$partition$confirmation)),
          "inner training patients overlap")
   dir.create(dest,recursive=TRUE,showWarnings=FALSE)
   saveRDS(model,file.path(dest,"complete_crit9.rds"))
   y01<-rlmt_y01(Yv)
   ll0<-rlmt_logloss(y01,oldroot)
   for(method in crit9_modes) {
     prob<-as.numeric(pred$probs[[method]])
     loc<-file.path(dest,method)
     dir.create(loc,recursive=TRUE,showWarnings=FALSE)
     accepted<-if(method=="ORIGINAL_RLR")FALSE else
       isTRUE(model$accepted[[method]])
     if(!accepted)
       need(max(abs(prob-oldroot))<1e-9,
            "unaccepted control changed original predictions")
     leaf<-if(method=="CONDITIONAL_TREE" && accepted)
       pred$leaf_id else rep(0L,length(prob))
     write.csv(data.frame(
       sample_id=rownames(Xv),y_true=as.character(Yv),
       y_pred=ifelse(prob>=.5,levels(Yv)[[2L]],levels(Yv)[[1L]]),
       prob=prob,original_RLR_prob=pred$probs$ORIGINAL_RLR,
       leaf_id=leaf,gate_side=pred$leaf_id),
       file.path(loc,"predictions_test.csv"),row.names=FALSE)
     val<-if(method=="ORIGINAL_RLR" ||
             !nrow(model$validation)) NULL else
        model$validation[model$validation$method==method,,drop=FALSE]
     i<-i+1L
     records<-function(name,part) if(is.null(part)) "" else
       as.character(part[[name]][[1L]])
     rows[[i]]<-data.frame(
       world=inp$world[[1L]],dataset=dataset,split_id=f$split_id,
       method=method,status="OK",n_train=length(train),
       n_test=length(test),
       logloss=rlmt_logloss(y01,prob),
       original_logloss=ll0,
       delta_vs_original=rlmt_logloss(y01,prob)-ll0,
       brier=mean((prob-y01)^2),original_brier=mean((oldroot-y01)^2),
       accepted=as.integer(accepted),
       nomination_count=model$nomination_n,
       selected_D_rank=if(is.null(model$selected_D_rank))NA_integer_ else
          model$selected_D_rank,
       selected_D_cv_lower=if(is.null(model$selected_D_rank))NA_real_ else
          model$ranking$D_cv_lower[[model$selected_D_rank]],
       validation_available=if(is.null(val))0L else
           as.integer(val$available[[1L]]),
       confirmation_lower=if(is.null(val)) NA_real_ else
            val$validation_lower[[1L]],
       gate_gene_a=records("gene_a",model$gate),
       gate_gene_b=records("gene_b",model$gate),
       false_gene_a=records("gene_a",model$pair_false),
       false_gene_b=records("gene_b",model$pair_false),
       true_gene_a=records("gene_a",model$pair_true),
       true_gene_b=records("gene_b",model$pair_true),
       model_reason=model$reason,
       root_max_abs_err=rooterr,
       warning_count=length(warn),model_dir=loc,
       original_model_dir=match$model_dir[[1L]],
       scientific_stage="EXPOSED_RETROSPECTIVE_DEV_ONLY",
       stringsAsFactors=FALSE)
     write.csv(do.call(rbind,rows),file.path(OUT,"run_manifest.csv"),
               row.names=FALSE)
     message("CRIT9_COMPLETED=",i,"/36 ",dataset,":",f$split_id,
             ":",method," accept=",as.integer(accepted),
             " deltaLL=",sprintf("%+.9f",
                rlmt_logloss(y01,prob)-ll0))
   }
   message("CRIT9_ROOT_PARITY=PASS ",dataset,":",f$split_id,
           " max_abs=",rooterr," warnings=",length(warn),
           " nominated=",model$nomination_n,
           " selected_D_rank=",if(is.null(model$selected_D_rank))"NA" else model$selected_D_rank)
 }
}
need(i==36L,"expected 36 method-fold results")
message("CRIT9_36_METHOD_FOLDS=PASS")
message("CRIT9_EXACT_ORIGINAL_RLR_PARITY=PASS")
message("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY")
