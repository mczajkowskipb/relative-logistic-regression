# Read-only post-run audit, no fitting, no coefficient tuning.
# Rebuilds test assignments from original ARFF and the publication fold generator.
# Usage: Rscript --vanilla R/new/48_rlmt_dc_published_5fold_audit.R
source("R/new/07_world_protocol_utils.R")

root <- "results/rlmt_dc_published_data_5fold_v1"
cfg <- jsonlite::read_json(
  "configs/rlmt_dc/published_data_5fold_v1.json",simplifyVector=FALSE)
need <- function(ok,msg) if (!isTRUE(ok))
  stop("PUBLICATION_PILOT_AUDIT_STOP: ",msg,call.=FALSE)
need(identical(cfg$out_root,root),"unexpected frozen root")
need(identical(trimws(readLines(file.path(root,"CODE_SHA.txt"))[1L]),
               "EXPECTED_FIVEFOLD_RUN_SHA_PLACEHOLDER"),
     "five-fold experiment code SHA mismatch")
need(identical(trimws(readLines(file.path(root,"status.txt"))[1L]),
               "RUNNER_EXIT=0"),
     "original 57-fit runner did not finish successfully")
mfile <- file.path(root,"run_manifest.csv")
if (!file.exists(mfile)) stop("No completed run manifest yet")
tab <- read.csv(mfile,stringsAsFactors=FALSE)
inputs <- read.csv(file.path(root,"input_manifest.csv"),
                   stringsAsFactors=FALSE)
required <- c("world","dataset","split_id","n_train","n_test",
              "method","status","logloss","brier",
              "realized_leaves","accepted_split","model_dir")
need(all(required %in% names(tab)),"run manifest schema mismatch")
need(nrow(inputs)==19L &&
     nrow(tab)==285L &&
     all(tab$status %in% c("OK","ERROR")),
     "Expected 19 datasets and 57 model records (including errors)")
need(!anyDuplicated(paste(tab$world,tab$dataset,tab$split_id,
                          tab$method,sep="|")),
     "duplicate model records")
need(all(tab$split_id %in% sprintf("rep00_fold%02d",0:4)),
     "unexpected five-fold split identifiers")
need(length(unique(paste(tab$world,tab$dataset,tab$split_id)))==95L,
     "missing or duplicate original outer folds")
for (dataset in unique(paste(inputs$world,inputs$dataset,sep="/"))) {
  x <- strsplit(dataset,"/",fixed=TRUE)[[1L]]
  subset <- tab[tab$world==x[[1L]] & tab$dataset==x[[2L]],,drop=FALSE]
  need(nrow(subset)==15L &&
       identical(sort(unique(subset$split_id)),
                 sprintf("rep00_fold%02d",0:4)),
       paste("incomplete five-fold dataset",dataset))
}
need(identical(sort(unique(tab$method)),
               sort(c("RLR","RLMT_v1","RLMT_DC_v0"))),
     "missing comparator")
audit <- vector("list",length=nrow(tab))
failures <- character()
data_cache <- new.env(parent=emptyenv())
for (i in seq_len(nrow(tab))) {
  r <- tab[i,,drop=FALSE]
  key <- paste(r$world,r$dataset,r$split_id,sep="/")
  if (!exists(key,envir=data_cache,inherits=FALSE)) {
    candidates <- inputs[
      inputs$world==r$world & inputs$dataset==r$dataset,
      ,drop=FALSE]
    need(nrow(candidates)==1L,"missing input source for dataset")
    xy <- read_xy_file(candidates$source_path[[1L]])
    folds <- make_repeated_stratified_folds(
      xy$y,n_folds=5L,n_repeats=1L,seed=42L)
    fold_index <- match(r$split_id[[1L]],
                        vapply(folds,function(x)x$split_id,
                               character(1)))
    need(!is.na(fold_index),paste("missing original outer fold",key))
    sp <- folds[[fold_index]]
    assign(key,list(
      y_test=xy$y[sp$test_idx],
      test_ids=rownames(xy$X)[sp$test_idx],
      n_test=length(sp$test_idx),
      n_train=length(sp$train_idx)),
      envir=data_cache)
  }
  base <- get(key,envir=data_cache,inherits=FALSE)
  need(base$n_test==r$n_test &&
       base$n_train==r$n_train,
       paste("outer fold shape mismatch",key))
  if (r$status!="OK") {
    failures <- c(failures,sprintf("%s %s: %s",
        key,r$method,if ("error" %in% names(r)) r$error else "error"))
    audit[[i]] <- data.frame(
      world=r$world,dataset=r$dataset,method=r$method,
      split_id=r$split_id,status=r$status,logloss_recomputed=NA_real_,
      brier_recomputed=NA_real_,delta_vs_RLR=NA_real_,
      realized_leaves=NA_integer_,stringsAsFactors=FALSE)
    next
  }
  dir <- r$model_dir[[1L]]
  pred_path <- file.path(dir,"predictions_test.csv")
  met_path <- file.path(dir,"metrics_test.json")
  need(file.exists(pred_path) && file.exists(met_path),
       paste("model output missing",key,r$method))
  pred <- read.csv(pred_path,stringsAsFactors=FALSE,
                   colClasses=c(sample_id="character",
                                y_true="character",y_pred="character"))
  need(nrow(pred)==base$n_test &&
       identical(pred$sample_id,base$test_ids) &&
       identical(pred$y_true,as.character(base$y_test)),
       paste("test sample/label mismatch",key,r$method))
  need(all(is.finite(pred$prob)) &&
       all(pred$prob>=0 & pred$prob<=1),
       paste("probability invalid",key,r$method))
  pos_label <- levels(base$y_test)[2L]
  neg_label <- levels(base$y_test)[1L]
  need(identical(pred$y_pred,
                 ifelse(pred$prob>=0.5,pos_label,neg_label)),
       paste("class decision threshold mismatch",key,r$method))
  y01 <- as.numeric(base$y_test==pos_label)
  p <- pmin(pmax(pred$prob,1e-15),1-1e-15)
  ll <- -mean(y01*log(p)+(1-y01)*log1p(-p))
  bs <- mean((y01-pred$prob)^2)
  met <- jsonlite::read_json(met_path,simplifyVector=TRUE)
  need(abs(ll-as.numeric(met$logloss))<=5.1e-5 &&
       abs(bs-as.numeric(met$brier))<=5.1e-5 &&
       abs(ll-r$logloss)<=5.1e-5 &&
       abs(bs-r$brier)<=5.1e-5,
       paste("independent metrics mismatch",key,r$method))
  tree_path <- file.path(dir,"tree.json")
  tree <- if(file.exists(tree_path))
    jsonlite::read_json(tree_path,simplifyVector=TRUE) else NULL
  if (!is.null(tree)) {
    need(tree$size %in% c(1L,2L) &&
         identical(as.integer(tree$size),
                   as.integer(r$realized_leaves)) &&
         isTRUE(tree$accepted_split)==isTRUE(r$accepted_split),
         paste("tree manifest mismatch",key,r$method))
  }
  audit[[i]] <- data.frame(
    world=r$world,dataset=r$dataset,method=r$method,
    split_id=r$split_id,status=r$status,logloss_recomputed=ll,
    brier_recomputed=bs,delta_vs_RLR=NA_real_,
    realized_leaves=r$realized_leaves,
    stringsAsFactors=FALSE)
}
df <- do.call(rbind,audit)
delta <- list()
for (key in unique(paste(df$world,df$dataset,df$split_id,sep="/"))) {
  ix <- which(paste(df$world,df$dataset,df$split_id,sep="/")==key)
  sub <- df[ix,,drop=FALSE]
  need(nrow(sub)==3L && length(unique(sub$method))==3L,
       paste("unexpected per-fold method triplet",key))
  base <- sub$logloss_recomputed[sub$method=="RLR"]
  if (length(base)!=1L || !is.finite(base)) next
  df$delta_vs_RLR[ix] <- sub$logloss_recomputed-base
  # Zero-split alternatives MUST exactly reproduce frozen published RLR.
  iR <- which(tab$world==sub$world[[1L]] &
              tab$dataset==sub$dataset[[1L]] &
              tab$split_id==sub$split_id[[1L]] &
              tab$method=="RLR")
  prefr <- file.path(tab$model_dir[[iR]],"predictions_test.csv")
  if (!file.exists(prefr)) next
  rpred <- read.csv(prefr,stringsAsFactors=FALSE)
  for (method in c("RLMT_v1","RLMT_DC_v0")) {
    j <- ix[which(sub$method==method)]
    if (!length(j) || df$status[[j]]!="OK" ||
        !isTRUE(df$realized_leaves[[j]]==1L)) next
    ii <- which(tab$world==sub$world[[1L]] &
                tab$dataset==sub$dataset[[1L]] &
                tab$split_id==sub$split_id[[1L]] &
                tab$method==method)
    tpred <- read.csv(file.path(tab$model_dir[[ii]],
                                "predictions_test.csv"),
                      stringsAsFactors=FALSE)
    need(identical(rpred$sample_id,tpred$sample_id) &&
         isTRUE(all.equal(rpred$prob,tpred$prob,tolerance=1e-12)),
         paste("unsplit RLR probability parity failure",key,method))
  }
}
report <- file.path(root,"published_5fold_audit_v1")
need(!dir.exists(report),"refusing to overwrite previous audit")
dir.create(report,recursive=TRUE,showWarnings=FALSE)
write.csv(df,file.path(report,"independent_fold_metrics.csv"),
          row.names=FALSE,na="")
writeLines(c(
  "Read-only independent first-fold metrics reconstructed from published ARFF.",
  "285 model attempts, 19 fixed datasets, five folds, three methods per fold.",
  "All failures preserved; no filtering based on performance.",
  "Five folds on exposed development datasets are exploratory, not inferential.",
  "Any accepted splits belong to the historical, method-development corpus."
),file.path(report,"README.txt"))
cat("=== PUBLICATION DATA FIVE-FOLD AUDIT ===\n")
print(df[,c("dataset","method","status","logloss_recomputed",
            "delta_vs_RLR","realized_leaves")],
      row.names=FALSE,digits=6)
cat("TOTAL_MODELS=",nrow(df),"\n",sep="")
cat("AUDITED_FOLDS=",length(unique(paste(df$world,df$dataset,df$split_id))),"\n",sep="")
cat("ERROR_MODELS=",sum(df$status!="OK"),"\n",sep="")
cat("REPORT=",report,"\n",sep="")
if (length(failures)) {
  cat("=== PRESERVED MODEL FAILURES ===\n")
  cat(paste(failures,collapse="\n"),"\n")
  stop("PUBLISHED_5FOLD_AUDIT=INCOMPLETE")
}
cat("PUBLISHED_5FOLD_AUDIT=PASS\n")
