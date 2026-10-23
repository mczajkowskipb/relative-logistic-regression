# Six genuine two-leaf CRIT7 vs frozen exact original RLR, exact parity.
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v4_rank_repair.R")
source("R/experimental/crit_v7_original_rlr_two_leaf.R")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("CRIT7_PILOT_STOP: original RLR reference needed")
REF <- normalizePath(args[[1L]],mustWork=TRUE)
OLD <- file.path(REF,"results/rlmt_dc_published_data_5fold_v1")
OUT <- "results/crit_v7_original_rlr_2x3_v1"
need <- function(ok,msg)
 if(!isTRUE(ok))stop("CRIT7_PILOT_STOP: ",msg,call.=FALSE)
need(dir.exists(OUT)&&
     file.exists(file.path(OUT,"CODE_SHA.txt")) &&
     !file.exists(file.path(OUT,"run_manifest.csv")),
     "root missing or previous result exists")
need(identical(trimws(readLines(file.path(OUT,"CODE_SHA.txt"))[[1L]]),
       system("git rev-parse HEAD",intern=TRUE)[[1L]]),
       "pinned source SHA changed")
need(trimws(readLines(file.path(OLD,"status.txt"))[[1L]])==
       "RUNNER_EXIT=0","original archive incomplete")
orig <- read.csv(file.path(OLD,"run_manifest.csv"),
                  stringsAsFactors=FALSE)
inputs <- read.csv(file.path(OLD,"input_manifest.csv"),
                   stringsAsFactors=FALSE)
need(nrow(orig)==285L&&nrow(inputs)==19L,
     "old reference manifest dimensions")
hash <- function(path) {
  txt <- suppressWarnings(system2("sha256sum",
    shQuote(normalizePath(path,mustWork=TRUE)),stdout=TRUE))
  need(length(txt)==1L &&
        grepl("^[a-f0-9]{64}[[:space:]]",txt[[1L]]),
        "sha256sum did not return SHA")
  substr(txt[[1L]],1L,64L)
}
rows <- list()
i <- 0L
for(dataset in c("chen-2002","gordon-2002")) {
  inp <- inputs[inputs$dataset==dataset,,drop=FALSE]
  need(nrow(inp)==1L,"input dataset absent")
  need(identical(tolower(hash(inp$source_path[[1L]])),
        tolower(inp$sha256[[1L]])), "original input file SHA mismatch")
  xy <- read_xy_file(inp$source_path[[1L]])
  need(is.factor(xy$y) && length(levels(xy$y))==2L,
       "invalid binary class schema")
  folds <- make_repeated_stratified_folds(
    xy$y,n_folds=5L,n_repeats=1L,seed=42L)
  need(length(folds)==5L &&
       identical(vapply(folds,function(z)z$split_id,character(1)),
        sprintf("rep00_fold%02d",0:4)),"old CV splits not replicated")
  for(f in folds[1:3]) {
    match <- orig[orig$world==inp$world[[1L]] &
       orig$dataset==dataset&orig$split_id==f$split_id&
       orig$method=="RLR",,drop=FALSE]
    need(nrow(match)==1L && match$status[[1L]]=="OK",
         "original archived RLR fold missing")
    od <- file.path(REF,match$model_dir[[1L]])
    oldtr <- read.csv(file.path(dirname(od),"samples_train.csv"),
          stringsAsFactors=FALSE,colClasses=c(sample_id="character"))
    oldp <- read.csv(file.path(od,"predictions_test.csv"),
          stringsAsFactors=FALSE,
          colClasses=c(sample_id="character",
                       y_true="character",y_pred="character"))
    train <- f$train_idx
    test <- f$test_idx
    Xt <- xy$X[train,,drop=FALSE]
    Yt <- xy$y[train]
    Xv <- xy$X[test,,drop=FALSE]
    Yv <- xy$y[test]
    need(length(intersect(train,test))==0 &&
       identical(sort(c(train,test)),seq_len(nrow(xy$X)))&&
       identical(as.character(oldtr$sample_id),rownames(Xt))&&
       identical(as.character(oldp$sample_id),rownames(Xv))&&
       identical(as.character(oldp$y_true),as.character(Yv)),
       "patient sample IDs/labels/split mismatch")
    dest <- file.path(OUT,"models",dataset,f$split_id)
    need(!dir.exists(dest),"would overwrite existing results")
    message("CRIT7_REAL_FIT=",dataset,":",f$split_id)
    warn <- character()
    model <- withCallingHandlers(crit7_fit(Xt,Yt),
      warning=function(w){
        warn <<- c(warn,conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    pred <- crit7_predict(model,Xv)
    origprob <- as.numeric(oldp$prob)
    probs <- pred$probs
    maxerr <- max(abs(probs$ORIGINAL_RLR-origprob))
    need(is.finite(maxerr) && maxerr<=1e-9,
         paste("ORIGINAL_RLR_FOLD_PARITY_FAILED",
                dataset,f$split_id,maxerr))
    dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    saveRDS(model,file.path(dest,"complete_crit7.rds"))
    y01 <- rlmt_y01(Yv)
    ref_ll <- rlmt_logloss(y01,origprob)
    for(method in crit7_modes) {
      ip <- probs[[method]]
      loc <- file.path(dest,method)
      dir.create(loc,recursive=TRUE,showWarnings=FALSE)
      write.csv(data.frame(
        sample_id=rownames(Xv),y_true=as.character(Yv),
        y_pred=ifelse(ip>=.5,levels(Yv)[2],levels(Yv)[1]),
        prob=ip,original_RLR_prob=probs$ORIGINAL_RLR,
        leaf_id=if(method=="CRIT_TREE")pred$leaf_id else
                rep(0L,length(ip))),
        file.path(loc,"predictions_test.csv"),row.names=FALSE)
      actual <- as.integer(method=="CRIT_TREE" &&
                           isTRUE(pred$split_accepted))
      gate_a <- if(is.null(model$gate)) "" else
        as.character(model$gate$gene_a[[1L]])
      gate_b <- if(is.null(model$gate)) "" else
        as.character(model$gate$gene_b[[1L]])
      false_a <- if(is.null(model$pair_false)) "" else
        as.character(model$pair_false$gene_a[[1L]])
      false_b <- if(is.null(model$pair_false)) "" else
        as.character(model$pair_false$gene_b[[1L]])
      true_a <- if(is.null(model$pair_true)) "" else
        as.character(model$pair_true$gene_a[[1L]])
      true_b <- if(is.null(model$pair_true)) "" else
        as.character(model$pair_true$gene_b[[1L]])
      i <- i+1L
      metric <- rlmt_logloss(y01,ip)
      rows[[i]] <- data.frame(
        world=inp$world[[1L]],dataset=dataset,split_id=f$split_id,
        method=method,status="OK",n_train=length(train),
        n_test=length(test),logloss=metric,
        original_logloss=ref_ll,delta_vs_original=metric-ref_ll,
        brier=mean((ip-y01)^2),
        original_brier=mean((origprob-y01)^2),
        tree_accepted=actual,
        gate_gene_a=gate_a,gate_gene_b=gate_b,
        leaf_false_gene_a=false_a,leaf_false_gene_b=false_b,
        leaf_true_gene_a=true_a,leaf_true_gene_b=true_b,
        selected_non_tree_interaction=as.integer(
          method=="GLOBAL_INTERACTION" &&
          !is.null(model$interaction_coef)),
        confirmation_lower=if(is.null(model$confirm)) NA_real_ else
            model$confirm$best_confirmed_lower,
        parent_max_abs_err=maxerr,
        warning_count=length(warn),model_dir=loc,
        original_model_dir=match$model_dir[[1L]],
        scientific_stage="EXPOSED_RETROSPECTIVE_DEV_ONLY",
        stringsAsFactors=FALSE)
      write.csv(do.call(rbind,rows),
        file.path(OUT,"run_manifest.csv"),row.names=FALSE)
      message("CRIT7_PILOT_COMPLETED=",i,"/24 ",
        dataset,":",f$split_id,":",method,
        " delta_original=",sprintf("%+.9f",metric-ref_ll),
        " actual_tree=",actual,
        " gate=",gate_a,">",gate_b)
    }
    message("CRIT7_PILOT_ROOT_PARITY=PASS ",
         dataset,":",f$split_id," max_abs=",maxerr,
         " warnings=",length(warn))
  }
}
need(i==24L,"incomplete 24 outer method-fold outcomes")
message("CRIT7_PILOT_ALL_24_METHOD_FOLDS=PASS")
message("CRIT7_PILOT_EXACT_ARCHIVED_ORIGINAL_RLR=PASS")
message("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY")
