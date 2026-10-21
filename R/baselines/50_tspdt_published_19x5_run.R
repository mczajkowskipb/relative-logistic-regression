# TSPDT-inspired 2011 reconstruction; 95 publication-matched outer folds.
# NEW isolated model family. Reads original completed/fully audited RLMT
# output without modifying it. No BigTSP required or trained in this run.
# Run by scripts/tspdt_recon_19x5_launch.py AFTER its fail-closed preflight.
source("R/new/07_world_protocol_utils.R")
source("R/baselines/tspdt_reconstruction_v1.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Require CONFIG and absolute REFERENCE_REPO")
cfg <- jsonlite::read_json(args[[1L]],simplifyVector=TRUE)
reference <- normalizePath(args[[2L]],mustWork=TRUE)
expect <- function(v,why) if(!isTRUE(v))
  stop("TSPDT_95FOLD_STOP: ",why,call.=FALSE)
out <- as.character(cfg$output_root)
expect(identical(out,"results/tspdt_recon_19x5_v1") &&
       identical(as.character(cfg$stage),
                 "TSPDT_2011_INSPIRED_95FOLD_EXPLORATORY_V1") &&
       identical(as.character(cfg$model),
                 "TSPDT2011_INSPIRED_RECONSTRUCTION_V1"),
       "model/config identity mismatch")
expect(identical(as.integer(cfg$total_fits),95L) &&
       identical(as.integer(cfg$fold_seed),42L) &&
       identical(as.integer(cfg$top_n),100L) &&
       identical(as.integer(cfg$min_leaf),5L) &&
       identical(as.integer(cfg$max_depth),3L) &&
       identical(as.character(cfg$pruning),"wilson") &&
       identical(as.numeric(cfg$confidence),0.25) &&
       identical(as.numeric(cfg$laplace_alpha),1),
       "model hyperparameters changed")
expect(file.exists(file.path(out,"input_manifest.csv")) &&
       file.exists(file.path(out,"CODE_SHA.txt")) &&
       !file.exists(file.path(out,"run_manifest.csv")),
       "preflight manifest absent or output already processed")
sha <- trimws(readLines(file.path(out,"CODE_SHA.txt"))[1L])
expect(identical(sha,as.character(system("git rev-parse HEAD",intern=TRUE)[1L])),
       "model HEAD differs from preflight")
refout <- file.path(reference,as.character(cfg$reference_result_root))
expect(file.exists(file.path(refout,"status.txt")) &&
       identical(trimws(readLines(file.path(refout,"status.txt"))[1L]),
                 "RUNNER_EXIT=0"),
       "reference runner not complete")
expect(file.exists(file.path(refout,"published_5fold_audit_v1",
                             "independent_fold_metrics.csv")) &&
       file.exists(file.path(refout,"published_5fold_accepted_replay_v1",
                             "accepted_models.csv")) &&
       file.exists(file.path(refout,"published_5fold_firstfold_parity_v1",
                             "firstfold_prediction_parity.csv")) &&
       file.exists(file.path(refout,"published_5fold_summary_v1",
                             "descriptive_aggregates.json")),
       "reference full audit not finished")
input <- read.csv(file.path(out,"input_manifest.csv"),
                  stringsAsFactors=FALSE)
expect(nrow(input)==19L && !anyDuplicated(
       paste(input$world,input$dataset,sep="/")),
       "expected exactly 19 unique source inputs")
reftab <- read.csv(file.path(refout,"run_manifest.csv"),
                   stringsAsFactors=FALSE)
expect(nrow(reftab)==285L && all(reftab$status=="OK"),
       "reference 285-fit model manifest not clean")
expected_order <- unlist(lapply(cfg$datasets,function(d)
  paste(d$world,d$names,sep="/")),use.names=FALSE)
expect(identical(paste(input$world,input$dataset,sep="/"),expected_order),
       "source inventory/order differs from frozen config")
sha256 <- function(path) {
  lines <- suppressWarnings(system2(
    "sha256sum",shQuote(normalizePath(path,mustWork=TRUE)),
    stdout=TRUE,stderr=TRUE))
  if(length(lines)!=1L || !grepl("^[a-f0-9]{64} ",lines[1L]))
    stop("TSPDT_95FOLD_STOP: cannot verify input SHA")
  substr(lines[1L],1L,64L)
}
rows <- list()
count <- 0L
for(i in seq_len(nrow(input))) {
  row <- input[i,,drop=FALSE]
  key <- paste(row$world,row$dataset,sep="/")
  filename <- as.character(row$source_path)
  expect(identical(sha256(filename),as.character(row$sha256)),
         paste("published ARFF integrity failure",key))
  xy <- read_xy_file(filename)
  X <- tspdt_recon_check_X(xy$X)
  y <- factor(xy$y)
  splits <- make_repeated_stratified_folds(
    y,n_folds=5L,n_repeats=1L,seed=42L)
  expect(length(splits)==5L &&
         identical(vapply(splits,function(z)z$split_id,character(1)),
                   sprintf("rep00_fold%02d",0:4)),
         paste("unexpected fold generator",key))
  for(f in splits) {
    fold_id <- f$split_id
    dest <- file.path(out,"worlds",row$world,"datasets",
                      row$dataset,"splits",fold_id,
                      "TSPDT2011_INSPIRED_RECON_V1")
    expect(!dir.exists(dest),paste("refuse fold overwrite",dest))
    ixtrain <- f$train_idx
    ixtest <- f$test_idx
    expect(length(intersect(ixtrain,ixtest))==0L &&
           identical(sort(c(ixtrain,ixtest)),seq_len(nrow(X))),
           paste("outer sample partition invalid",key,fold_id))
    pair_ref <- reftab[reftab$world==row$world &
                       reftab$dataset==row$dataset &
                       reftab$split_id==fold_id &
                       reftab$method=="RLR",,drop=FALSE]
    expect(nrow(pair_ref)==1L &&
           as.integer(pair_ref$n_train)==length(ixtrain) &&
           as.integer(pair_ref$n_test)==length(ixtest),
           paste("reference outer-fold shape differs",key,fold_id))
    refbase <- dirname(pair_ref$model_dir)
    train_ids <- read.csv(file.path(refbase,"samples_train.csv"),
                          stringsAsFactors=FALSE,
                          colClasses=c(sample_id="character"))$sample_id
    test_ids <- read.csv(file.path(refbase,"samples_test.csv"),
                         stringsAsFactors=FALSE,
                         colClasses=c(sample_id="character"))$sample_id
    expect(identical(train_ids,rownames(X)[ixtrain]) &&
           identical(test_ids,rownames(X)[ixtest]),
           paste("reference outer sample IDs differ",key,fold_id))
    Xtr <- X[ixtrain,,drop=FALSE]
    Xte <- X[ixtest,,drop=FALSE]
    ytr <- y[ixtrain]
    yte <- y[ixtest]
    message("=== TSPDT-RECON ",key," ",fold_id," ===")
    dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    result <- tryCatch({
      model <- tspdt_recon_fit(
        Xtr,ytr,top_n=100L,min_leaf=5L,max_depth=3L,
        pruning="wilson",confidence=0.25,
        laplace_alpha=1)
      pred <- tspdt_recon_predict(model,Xte)
      expect(identical(levels(pred$class),levels(ytr)) &&
             length(pred$prob_class1)==length(ixtest),
             paste("test prediction schema",key,fold_id))
      label <- levels(ytr)[2L]
      true01 <- as.numeric(yte==label)
      prob <- pred$prob_class1
      margin <- pmax(pmin(prob,1-1e-15),1e-15)
      ll <- -mean(true01*log(margin)+(1-true01)*log1p(-margin))
      brier <- mean((true01-prob)^2)
      accuracy <- mean(as.character(pred$class)==as.character(yte))
      rec <- vapply(levels(ytr),function(cls) {
        yes <- which(as.character(yte)==cls)
        if(!length(yes)) return(NA_real_)
        mean(as.character(pred$class[yes])==cls)
      },numeric(1))
      balanced <- if(anyNA(rec)) NA_real_ else mean(rec)
      f1 <- vapply(levels(ytr),function(cls) {
        tp <- sum(as.character(pred$class)==cls & as.character(yte)==cls)
        fp <- sum(as.character(pred$class)==cls & as.character(yte)!=cls)
        fn <- sum(as.character(pred$class)!=cls & as.character(yte)==cls)
        denom <- 2*tp+fp+fn
        if(denom==0L) return(NA_real_)
        as.numeric(2*tp/denom)
      },numeric(1))
      macrof1 <- if(anyNA(f1)) NA_real_ else mean(f1)
      saveRDS(model,file.path(dest,"model.rds"))
      write.csv(data.frame(
        sample_id=rownames(Xte),y_true=as.character(yte),
        y_pred=as.character(pred$class),prob=prob,
        leaf_id=pred$leaf_id),
        file.path(dest,"predictions_test.csv"),row.names=FALSE)
      write.csv(tspdt_recon_nodes(model),file.path(dest,"nodes.csv"),
                row.names=FALSE)
      jsonlite::write_json(list(
        family=model$method,source_claim=model$source_claim,
        world=row$world,dataset=row$dataset,split_id=fold_id,
        n_train=length(ixtrain),n_test=length(ixtest),
        logloss=ll,brier=brier,accuracy=accuracy,
        balanced_accuracy=balanced,macro_f1=macrof1,
        internal_nodes=model$internal_nodes,leaves=model$leaves,
        selected_genes=model$selected_genes,
        probability_kind="train-only Laplace; not original TSPDT"),
        file.path(dest,"metrics_test.json"),
        auto_unbox=TRUE,pretty=TRUE,digits=12,null="null")
      list(status="OK",error="",logloss=ll,brier=brier,
           accuracy=accuracy,balanced_accuracy=balanced,
           macro_f1=macrof1,internal_nodes=model$internal_nodes,
           leaves=model$leaves)
    },error=function(e) {
      jsonlite::write_json(list(error=conditionMessage(e)),
        file.path(dest,"error.json"),auto_unbox=TRUE,pretty=TRUE)
      list(status="ERROR",error=conditionMessage(e),
           logloss=NA_real_,brier=NA_real_,accuracy=NA_real_,
           balanced_accuracy=NA_real_,macro_f1=NA_real_,
           internal_nodes=NA_integer_,leaves=NA_integer_)
    })
    count <- count+1L
    rows[[count]] <- data.frame(
      world=row$world,dataset=row$dataset,split_id=fold_id,
      method="TSPDT2011_INSPIRED_RECON_V1",
      n_train=length(ixtrain),n_test=length(ixtest),
      status=result$status,error=result$error,
      logloss=result$logloss,brier=result$brier,
      accuracy=result$accuracy,
      balanced_accuracy=result$balanced_accuracy,
      macro_f1=result$macro_f1,
      internal_nodes=result$internal_nodes,
      leaves=result$leaves,
      model_dir=dest,stringsAsFactors=FALSE)
    write.csv(do.call(rbind,rows),file.path(out,"run_manifest.csv"),
              row.names=FALSE,na="")
    message(" status=",result$status,
            " logloss=",format(result$logloss,digits=6))
  }
}
tab <- do.call(rbind,rows)
expect(nrow(tab)==95L,"wrong number of attempted fits")
cat("TSPDT_FITS=",nrow(tab),"\n",sep="")
cat("TSPDT_FIT_ERRORS=",sum(tab$status!="OK"),"\n",sep="")
if(any(tab$status!="OK"))
  stop("TSPDT_95FOLD_RUN=INCOMPLETE: retained error outputs")
cat("TSPDT_95FOLD_RUN=PASS\n")
