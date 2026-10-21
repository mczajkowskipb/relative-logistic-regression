# RLMT-DC publication-dataset exploratory pilot.
# NEW wrapper only: uses unchanged published RLR and original RLMT-v1,
# then new RLMT-DC-v0 on the same frozen outer-training/test split.
# NEVER writes paper-frozen results or historical GATE3/4 artifacts.
# Command: Rscript --vanilla R/new/44_rlmt_dc_published_pilot.R <config.json>

source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")

args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L) stop("Requires exactly one frozen config path")
config_file <- args[[1L]]
cfg <- jsonlite::read_json(config_file,simplifyVector=FALSE)
root <- as.character(cfg$out_root)
check <- function(x,msg) if (!isTRUE(x)) stop(
  "PUBLISHED_DATA_PILOT_STOP: ",msg,call.=FALSE)
check(identical(root,"results/rlmt_dc_published_data_smoke_v1"),
      "unexpected output directory")
check(identical(as.character(cfg$stage),
      "PUBLICATION_DATA_EXPLORATORY_SMOKE_V1"),
      "unexpected experiment stage")
check(as.integer(cfg$model_seed)==42L &&
      as.integer(cfg$n_folds)==5L &&
      as.integer(cfg$n_repeats)==1L &&
      as.integer(cfg$max_splits_per_dataset)==1L,
      "outer split must be the first original stratified fold, seed42")
check(file.exists(file.path(root,"input_manifest.csv")),
      "source manifest absent: launch preflight first")
check(file.exists(file.path(root,"CODE_SHA.txt")),
      "execution code SHA absent")
check(!file.exists(file.path(root,"run_manifest.csv")),
      "completed or partially run experiment must not be overwritten")
check(length(cfg$cohorts)==3L,
      "frozen three publication input shards required")

manifest <- read.csv(file.path(root,"input_manifest.csv"),
                     stringsAsFactors=FALSE,check.names=FALSE)
expected <- do.call(rbind,lapply(cfg$cohorts,function(shard) {
  data.frame(world=as.character(shard$world),
             dataset=unlist(shard$datasets,use.names=FALSE),
             stringsAsFactors=FALSE)
}))
check(nrow(manifest)==19L && nrow(expected)==19L &&
      identical(manifest$world,expected$world) &&
      identical(manifest$dataset,expected$dataset),
      "source manifest differs from complete registered datasets")
check(all(file.exists(manifest$source_path)) &&
      all(nchar(manifest$sha256)==64L) &&
      !anyDuplicated(paste(manifest$world,manifest$dataset)),
      "source manifest integrity/schema failure")

sha256 <- function(path) {
  z <- suppressWarnings(system2("sha256sum",
       shQuote(normalizePath(path,mustWork=TRUE)),
       stdout=TRUE,stderr=TRUE))
  if (length(z)!=1L ||
      !grepl("^[0-9a-f]{64}[[:space:]]+",z[[1L]])) {
    stop("Cannot obtain SHA256: ",path)
  }
  substr(z[[1L]],1L,64L)
}
if (!identical(
      as.character(trimws(readLines(file.path(root,"CODE_SHA.txt"))[[1L]])),
      as.character(system("git rev-parse HEAD",intern=TRUE)[[1L]]))) {
  stop("Code SHA changed after preflight")
}

methods <- c("RLR", "RLMT_v1", "RLMT_DC_v0")
all_results <- list()
counter <- 0L
read_metrics <- function(d) {
  p <- file.path(d,"metrics_test.json")
  check(file.exists(p),paste("missing metric artifact",p))
  z <- jsonlite::read_json(p,simplifyVector=TRUE)
  check(length(z$logloss)==1L && is.finite(z$logloss) &&
        length(z$brier)==1L && is.finite(z$brier),
        paste("invalid logloss/brier",d))
  z
}

for (i in seq_len(nrow(manifest))) {
  dataset <- as.character(manifest$dataset[[i]])
  world <- as.character(manifest$world[[i]])
  path <- as.character(manifest$source_path[[i]])
  expected_hash <- as.character(manifest$sha256[[i]])
  message(sprintf("=== WORLD %s DATASET %s (%d/19) ===",
                  world,dataset,i))
  check(identical(sha256(path),expected_hash),
        paste("input changed since preflight:",dataset))
  xy <- read_xy_file(path)
  check(identical(length(levels(xy$y)),2L) &&
        length(unique(rownames(xy$X)))==nrow(xy$X) &&
        all(vapply(xy$X,function(c)all(is.finite(c)),logical(1))),
        paste("ARFF source/schema invalid:",dataset))
  folds <- make_repeated_stratified_folds(
    xy$y,n_folds=as.integer(cfg$n_folds),
    n_repeats=as.integer(cfg$n_repeats),
    seed=as.integer(cfg$model_seed))
  check(length(folds)>=1L &&
        identical(folds[[1L]]$split_id,"rep00_fold00"),
        paste("wrong original first split generator:",dataset))

  for (sp in head(folds,as.integer(cfg$max_splits_per_dataset))) {
    base <- file.path(root,"worlds",world,"datasets",dataset,
                      "splits",sp$split_id)
    check(!dir.exists(base),
          paste("existing split output is not overwritable",base))
    dir.create(base,recursive=TRUE,showWarnings=FALSE)
    tr <- sp$train_idx
    te <- sp$test_idx
    check(length(intersect(tr,te))==0L &&
          identical(sort(c(tr,te)),seq_len(nrow(xy$X))),
          paste("overlapping split",dataset))
    Xtr <- xy$X[tr,,drop=FALSE]
    Xte <- xy$X[te,,drop=FALSE]
    ytr <- xy$y[tr]
    yte <- xy$y[te]
    write.csv(data.frame(sample_id=rownames(Xtr)),
              file.path(base,"samples_train.csv"),row.names=FALSE)
    write.csv(data.frame(sample_id=rownames(Xte)),
              file.path(base,"samples_test.csv"),row.names=FALSE)
    write.csv(data.frame(class=levels(ytr),
                         train_count=as.integer(table(ytr)),
                         test_count=as.integer(table(yte))),
              file.path(base,"class_support.csv"),row.names=FALSE)

    # The only supervised global predictor pool shared across methods
    # is computed from this outer TRAINING portion. RLMT-DC internally
    # rebuilds its discovery OOF pools inside training partitions.
    pool <- tryCatch(compute_top_pairs_switchbox(
      Xtr,ytr,top_n=1000L,buffer_mult=5L,
      out_dir=file.path(base,"pairs_train","top1000"),
      save_all=FALSE),error=function(e)e)
    if (inherits(pool,"error")) {
      stop("Cannot construct shared outer-training pool for ",
           dataset,": ",conditionMessage(pool))
    }
    check(nrow(pool)>=10L,
          paste("fewer than 10 pairs on outer training:",dataset))

    for (method in methods) {
      modir <- file.path(base,method)
      check(!dir.exists(modir),
            paste("model directory already exists",modir))
      message(sprintf("  -> %s %s",method,sp$split_id))
      extras <- list(world_name=world,dataset_id=dataset,
                     split_id=sp$split_id,pool=1000L)
      run <- tryCatch({
        if (method=="RLR") {
          dir.create(modir,recursive=TRUE,showWarnings=FALSE)
          fit_rlr_train_test(
            Xtr,ytr,Xte,yte,pairs_pool=pool,
            K=10L,m=10L,d=2L,out_dir=modir,
            alpha=0.5,lambda_choice="lambda.min",
            seed=42L,extra_metrics=extras)
        } else if (method=="RLMT_v1") {
          dir.create(modir,recursive=TRUE,showWarnings=FALSE)
          rlmt_run(Xtr,ytr,Xte,yte,pool,max_leaves=2L,
                   config=list(),out_dir=modir,
                   extra_metrics=extras)
        } else if (method=="RLMT_DC_v0") {
          rlmt_dc_run(Xtr,ytr,Xte,yte,pool,
                      config=list(),out_dir=modir,
                      extra_metrics=extras)
        }
        met <- read_metrics(modir)
        tree_file <- file.path(modir,"tree.json")
        tree <- if (file.exists(tree_file))
          jsonlite::read_json(tree_file,simplifyVector=TRUE) else NULL
        list(status="OK",error="",
             logloss=as.numeric(met$logloss),
             brier=as.numeric(met$brier),
             auc=as.numeric(met$auc),
             f1=as.numeric(met$f1),
             realized_leaves=if (is.null(tree)) 1L
                              else as.integer(tree$size),
             accepted_split=if (is.null(tree)) FALSE
                              else isTRUE(tree$accepted_split))
      },error=function(e) {
        ensure_dir(modir)
        write_error_json(file.path(modir,"error.json"),e)
        list(status="ERROR",error=conditionMessage(e),
             logloss=NA_real_,brier=NA_real_,
             auc=NA_real_,f1=NA_real_,
             realized_leaves=NA_integer_,accepted_split=NA)
      })
      counter <- counter+1L
      all_results[[counter]] <- data.frame(
        world=world,dataset=dataset,split_id=sp$split_id,
        n_train=length(tr),n_test=length(te),
        method=method, status=run$status,error=run$error,
        logloss=run$logloss,brier=run$brier,
        auc=run$auc,f1=run$f1,
        realized_leaves=run$realized_leaves,
        accepted_split=run$accepted_split,
        model_dir=modir,
        stringsAsFactors=FALSE)
      summary <- do.call(rbind,all_results)
      # Incremental durable manifest even if a later fit fails.
      write.csv(summary,file.path(root,"run_manifest.csv"),
                row.names=FALSE,na="")
      message(sprintf("     status=%s logloss=%s",
          run$status,
          if (is.finite(run$logloss))
              sprintf("%.6f",run$logloss) else "NA"))
    }
    # As an invariant, an unsplit RLMT-v1 or RLMT-DC must be identical
    # to the published RLR baseline on the same test identities.
    ref <- file.path(base,"RLR","predictions_test.csv")
    if (file.exists(ref)) {
      pred_ref <- read.csv(ref,stringsAsFactors=FALSE)
      for (method in methods[-1L]) {
        other <- file.path(base,method,"predictions_test.csv")
        if (!file.exists(other)) next
        z <- read.csv(other,stringsAsFactors=FALSE)
        row <- all_results[[which(vapply(all_results,function(rr)
            identical(rr$model_dir, file.path(base,method)),
            logical(1)))]]
        if (isTRUE(row$realized_leaves==1L)) {
          check(identical(pred_ref$sample_id,z$sample_id) &&
                identical(as.character(pred_ref$y_pred),
                          as.character(z$y_pred)) &&
                isTRUE(all.equal(pred_ref$prob,z$prob,
                                 tolerance=1e-12)),
                paste("unsplit/RLR parity failure",world,dataset,method))
        }
      }
    }
  }
}
tab <- do.call(rbind,all_results)
check(nrow(tab)==57L,"expected 19 cohorts x 3 methods")
cat("=== RLR / RLMT_v1 / RLMT_DC_v0 PUBLICATION SMOKE ===\n")
print(tab[,c("dataset","method","status","logloss",
              "realized_leaves","accepted_split")],
      row.names=FALSE,digits=6)
cat("PILOT_DATASETS=",length(unique(tab$dataset)),"\n",sep="")
cat("PILOT_METHOD_RUNS=",nrow(tab),"\n",sep="")
cat("PILOT_FAILURES=",sum(tab$status!="OK"),"\n",sep="")
if (any(tab$status!="OK")) {
  stop("PUBLISHED_DATA_PILOT=INCOMPLETE: all failures preserved in manifest")
}
cat("PUBLISHED_DATA_PILOT=PASS\n")
