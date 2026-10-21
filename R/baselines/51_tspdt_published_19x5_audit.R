# Read-only replay/audit of 95 saved TSPDT-inspired reconstruction models.
# No R model fitting; original RLR output is READ ONLY.
source("R/new/07_world_protocol_utils.R")
source("R/baselines/tspdt_reconstruction_v1.R")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("Require absolute REFERENCE_REPO")
reference <- normalizePath(args[[1L]],mustWork=TRUE)
cfg <- jsonlite::read_json(
  "configs/tspdt/tspdt_reconstruction_19x5_v1.json",
  simplifyVector=TRUE)
out <- as.character(cfg$output_root)
base <- file.path(reference,cfg$reference_result_root)
check <- function(ok,why) if(!isTRUE(ok))
  stop("TSPDT_95FOLD_AUDIT_STOP: ",why,call.=FALSE)
check(identical(out,"results/tspdt_recon_19x5_v1"),
      "wrong result root")
check(identical(trimws(readLines(file.path(out,"status.txt"))[1L]),
                "RUNNER_EXIT=0"),
      "TSPDT run did not exit successfully")
check(identical(trimws(readLines(file.path(base,"status.txt"))[1L]),
                "RUNNER_EXIT=0"),
      "reference RLR run did not exit successfully")
report <- file.path(out,"audit_v1")
check(!dir.exists(report),"audit already exists, refuse overwrite")
manifest <- read.csv(file.path(out,"run_manifest.csv"),
                     stringsAsFactors=FALSE)
inputs <- read.csv(file.path(out,"input_manifest.csv"),
                   stringsAsFactors=FALSE)
reference_models <- read.csv(file.path(base,"run_manifest.csv"),
                             stringsAsFactors=FALSE)
check(nrow(manifest)==95L && nrow(inputs)==19L &&
      nrow(reference_models)==285L &&
      all(manifest$status=="OK"),
      "incomplete/failed 95-fit or 285-fit benchmark")
check(!anyDuplicated(paste(manifest$world,manifest$dataset,
                           manifest$split_id,sep="/")),
      "duplicate TSPDT fold entry")
records <- list()
counter <- 0L
for(i in seq_len(nrow(inputs))) {
  src <- inputs[i,,drop=FALSE]
  name <- paste(src$world,src$dataset,sep="/")
  xy <- read_xy_file(src$source_path[[1L]])
  X <- tspdt_recon_check_X(xy$X)
  y <- factor(xy$y)
  folds <- make_repeated_stratified_folds(y,5L,1L,42L)
  for(f in folds) {
    r <- manifest[manifest$world==src$world &
                  manifest$dataset==src$dataset &
                  manifest$split_id==f$split_id,,drop=FALSE]
    check(nrow(r)==1L,paste("missing fold",name,f$split_id))
    train <- f$train_idx
    test <- f$test_idx
    check(r$n_train==length(train) && r$n_test==length(test),
          paste("fold size differs",name,f$split_id))
    matchref <- reference_models[
      reference_models$world==src$world &
      reference_models$dataset==src$dataset &
      reference_models$split_id==f$split_id &
      reference_models$method=="RLR",,drop=FALSE]
    check(nrow(matchref)==1L &&
          matchref$status=="OK",
          paste("missing full RLR comparator",name,f$split_id))
    ref_path <- file.path(reference,matchref$model_dir[[1L]],
                          "predictions_test.csv")
    ref <- read.csv(ref_path,stringsAsFactors=FALSE,
                    colClasses=c(sample_id="character",
                                 y_true="character",y_pred="character"))
    saved <- r$model_dir[[1L]]
    check(dir.exists(saved) &&
          file.exists(file.path(saved,"model.rds")) &&
          file.exists(file.path(saved,"predictions_test.csv")) &&
          file.exists(file.path(saved,"nodes.csv")) &&
          file.exists(file.path(saved,"metrics_test.json")),
          paste("incomplete serialized TSPDT model",name,f$split_id))
    model <- readRDS(file.path(saved,"model.rds"))
    predictions <- read.csv(file.path(saved,"predictions_test.csv"),
                            stringsAsFactors=FALSE,
                            colClasses=c(sample_id="character",
                                         y_true="character",
                                         y_pred="character"))
    check(model$method=="TSPDT2011_INSPIRED_RECONSTRUCTION_V1" &&
          model$source_claim=="not_original_TSPDT_2011" &&
          model$n_train==length(train) &&
          identical(model$training_genes,colnames(X)) &&
          identical(model$selected_genes,
                    tspdt_recon_select(X[train,,drop=FALSE],100L)$names),
          paste("serialized model source/schema/feature-filter mismatch",
                name,f$split_id))
    replay <- tspdt_recon_predict(model,X[test,,drop=FALSE])
    check(identical(predictions$sample_id,rownames(X)[test]) &&
          identical(ref$sample_id,predictions$sample_id) &&
          identical(ref$y_true,as.character(y[test])) &&
          identical(predictions$y_true,as.character(y[test])) &&
          identical(predictions$y_pred,as.character(replay$class)) &&
          identical(as.integer(predictions$leaf_id),replay$leaf_id) &&
          isTRUE(all.equal(as.numeric(predictions$prob),
                           replay$prob_class1,tolerance=1e-12)),
          paste("test routing/replay mismatch",name,f$split_id))
    nodes <- read.csv(file.path(saved,"nodes.csv"),
                      stringsAsFactors=FALSE)
    check(nrow(nodes)==model$internal_nodes+model$leaves &&
          sum(nodes$is_leaf)==model$leaves,
          paste("node count mismatch",name,f$split_id))
    true <- as.character(y[test])
    pr <- as.character(replay$class)
    classes <- levels(y)
    p <- replay$prob_class1
    y01 <- as.numeric(true==classes[2L])
    prob <- pmin(pmax(p,1e-15),1-1e-15)
    ll <- -mean(y01*log(prob)+(1-y01)*log1p(-prob))
    bs <- mean((y01-p)^2)
    accuracy <- mean(pr==true)
    recall <- vapply(classes,function(cls) {
      ix <- which(true==cls)
      if (!length(ix)) return(NA_real_)
      mean(pr[ix]==cls)
    },numeric(1))
    bal <- if(anyNA(recall)) NA_real_ else mean(recall)
    f1 <- vapply(classes,function(cls) {
      tp <- sum(pr==cls & true==cls)
      fp <- sum(pr==cls & true!=cls)
      fn <- sum(pr!=cls & true==cls)
      if(2*tp+fp+fn==0L) return(NA_real_)
      2*tp/(2*tp+fp+fn)
    },numeric(1))
    mf1 <- if(anyNA(f1)) NA_real_ else mean(f1)
    stats <- jsonlite::read_json(file.path(saved,"metrics_test.json"),
                                 simplifyVector=TRUE)
    check(abs(ll-stats$logloss)<1e-10 &&
          abs(bs-stats$brier)<1e-10 &&
          abs(accuracy-stats$accuracy)<1e-10 &&
          (is.na(bal) && is.null(stats$balanced_accuracy) ||
           is.finite(bal) &&
           abs(bal-stats$balanced_accuracy)<1e-10) &&
          (is.na(mf1) && is.null(stats$macro_f1) ||
           is.finite(mf1) &&
           abs(mf1-stats$macro_f1)<1e-10),
          paste("stored metric mismatch",name,f$split_id))
    check(abs(ll-r$logloss)<1e-5 && abs(bs-r$brier)<1e-5,
          paste("manifest logloss/Brier mismatch",name,f$split_id))
    ref_true <- as.character(ref$y_true)
    ref_pred <- as.character(ref$y_pred)
    ref_accuracy <- mean(ref_true==ref_pred)
    count <- counter <- counter+1L
    records[[counter]] <- data.frame(
      world=src$world,dataset=src$dataset,
      split_id=f$split_id,n_test=length(test),
      rlr_accuracy=ref_accuracy,
      tspdt_accuracy=accuracy,
      delta_accuracy=accuracy-ref_accuracy,
      tspdt_balanced_accuracy=bal,
      tspdt_macro_f1=mf1,
      tspdt_logloss_laplace=ll,
      tspdt_brier_laplace=bs,
      rlr_logloss=as.numeric(matchref$logloss),
      difference_tspdt_minus_RLR_logloss=ll-as.numeric(matchref$logloss),
      internal_nodes=model$internal_nodes,leaves=model$leaves,
      stringsAsFactors=FALSE)
  }
}
check(counter==95L,"not all 95 folds independently replayed")
res <- do.call(rbind,records)
dataset <- do.call(rbind,lapply(
  split(res,paste(res$world,res$dataset,sep="/")),function(d) {
    check(nrow(d)==5L && length(unique(d$split_id))==5L,
          "incomplete five-fold dataset")
    data.frame(world=d$world[1L],dataset=d$dataset[1L],
      n_folds=5L,mean_accuracy_RLR=mean(d$rlr_accuracy),
      mean_accuracy_TSPDT=mean(d$tspdt_accuracy),
      delta_accuracy_TSPDT_minus_RLR=mean(d$delta_accuracy),
      mean_logloss_RLR=mean(d$rlr_logloss),
      mean_logloss_TSPDT_Laplace=mean(d$tspdt_logloss_laplace),
      delta_logloss_TSPDT_minus_RLR=mean(
        d$difference_tspdt_minus_RLR_logloss),
      mean_nodes=mean(d$internal_nodes),mean_leaves=mean(d$leaves),
      stringsAsFactors=FALSE)
  }))
check(nrow(dataset)==19L,"expected 19 dataset summary rows")
dir.create(report,recursive=TRUE,showWarnings=FALSE)
write.csv(res,file.path(report,"per_fold.csv"),row.names=FALSE,na="")
write.csv(dataset,file.path(report,"per_dataset.csv"),
          row.names=FALSE,na="")
writeLines(c("TSPDT-inspired reconstruction; original 95 outer folds.",
             "All RLR outputs read-only; all 95 tree model RDS files replayed.",
             "Historical data, exploratory; Laplace leaf scores not original.",
             "RLR comparison accuracy from frozen existing predictions."),
           file.path(report,"README.txt"))
print(dataset[,c("dataset","mean_accuracy_RLR",
                 "mean_accuracy_TSPDT",
                 "delta_accuracy_TSPDT_minus_RLR",
                 "mean_logloss_TSPDT_Laplace")],
      row.names=FALSE,digits=6)
cat("AUDITED_TSPDT_FOLDS=",counter,"\n",sep="")
cat("TSPDT_95FOLD_AUDIT=PASS\n")
cat("REPORT=",report,"\n",sep="")
