# Read-only replay of every REAL accepted two-leaf publication-data pilot
# model (not the forced synthetic DEV fixture). No fitting or refitting.
# Run only after successful 45_rlmt_dc_published_audit.R.
source("R/new/07_world_protocol_utils.R")

root <- "results/rlmt_dc_published_data_smoke_v1"
report <- file.path(root,"published_accepted_replay_v1")
check <- function(ok,reason) if(!isTRUE(ok))
  stop("PUBLISHED_ACCEPTED_REPLAY_STOP: ",reason,call.=FALSE)
check(!dir.exists(report),"refuse to overwrite accepted replay")
check(file.exists(file.path(root,"published_audit_v1",
                            "independent_fold_metrics.csv")),
      "run the independent full 57-fit audit before this replay")
check(identical(trimws(readLines(file.path(root,"status.txt"))[1L]),
                "RUNNER_EXIT=0"),"original fit runner did not complete")
check(identical(trimws(readLines(file.path(root,"CODE_SHA.txt"))[1L]),
                "814fee71fe6306db1e62b2ca47b1038df82d77a8"),
      "incorrect original model code SHA")
tab <- read.csv(file.path(root,"run_manifest.csv"),
                stringsAsFactors=FALSE)
inputs <- read.csv(file.path(root,"input_manifest.csv"),
                   stringsAsFactors=FALSE)
check(nrow(tab)==57L && nrow(inputs)==19L,
      "expected complete pilot manifest")
picked <- which(!is.na(tab$accepted_split) &
                tab$accepted_split & tab$status=="OK")
records <- list()

for (j in seq_along(picked)) {
  entry <- tab[picked[[j]],,drop=FALSE]
  key <- paste(entry$world,entry$dataset,entry$method,sep="/")
  check(entry$method %in% c("RLMT_v1","RLMT_DC_v0"),
        paste("unknown accepted method",key))
  target <- inputs[inputs$world==entry$world &
                     inputs$dataset==entry$dataset,,drop=FALSE]
  check(nrow(target)==1L,paste("expected one matching ARFF",key))
  xy <- read_xy_file(target$source_path[[1L]])
  fold <- make_repeated_stratified_folds(
    xy$y,n_folds=5L,n_repeats=1L,seed=42L)[[1L]]
  check(identical(fold$split_id,entry$split_id),
        paste("incorrect frozen outer split",key))
  Xte <- xy$X[fold$test_idx,,drop=FALSE]
  yte <- xy$y[fold$test_idx]
  base <- entry$model_dir[[1L]]
  check(startsWith(base,paste0(root,"/worlds/")),
        paste("model directory outside the pilot root",key))
  required <- c("tree.json","test_routes.csv","nodes.csv",
                "predictions_test.csv","metrics_test.json")
  check(all(file.exists(file.path(base,required))),
        paste("accepted model artifacts missing",key))
  tree <- jsonlite::read_json(file.path(base,"tree.json"),
                              simplifyVector=TRUE)
  nodes <- read.csv(file.path(base,"nodes.csv"),
                    stringsAsFactors=FALSE)
  routes <- read.csv(file.path(base,"test_routes.csv"),
                     stringsAsFactors=FALSE,
                     colClasses=c(sample_id="character"))
  pred <- read.csv(file.path(base,"predictions_test.csv"),
                   stringsAsFactors=FALSE,
                   colClasses=c(sample_id="character",
                                y_true="character",
                                y_pred="character"))
  met <- jsonlite::read_json(file.path(base,"metrics_test.json"),
                             simplifyVector=TRUE)
  check(isTRUE(tree$accepted_split) &&
        as.integer(tree$size)==2L &&
        as.integer(entry$realized_leaves)==2L &&
        nrow(nodes)==3L &&
        identical(as.integer(nodes$node_id),c(1L,2L,3L)) &&
        all(as.integer(nodes$n[2:3])>=16L) &&
        all(as.integer(nodes$class_0[2:3])>=6L) &&
        all(as.integer(nodes$class_1[2:3])>=6L),
        paste("accepted tree structure/support differs",key))
  a <- as.character(tree$split$gene_a)
  b <- as.character(tree$split$gene_b)
  check(length(a)==1L && length(b)==1L &&
        a %in% names(Xte) && b %in% names(Xte),
        paste("missing split features",key))
  expected_leaf <- ifelse(Xte[[a]]>Xte[[b]],2L,3L)
  check(!anyNA(expected_leaf) &&
        identical(routes$sample_id,rownames(Xte)) &&
        identical(pred$sample_id,rownames(Xte)) &&
        identical(pred$y_true,as.character(yte)) &&
        identical(as.integer(routes$leaf_id),
                  as.integer(expected_leaf)),
        paste("route or sample identity mismatch",key))
  recovered <- rep(NA_real_,nrow(Xte))
  for (leaf in c(2L,3L)) {
    fold_idx <- which(expected_leaf==leaf)
    d <- file.path(base,"leaf_models",
                   paste0("node_",leaf))
    check(all(file.exists(file.path(
      d,c("cvfit.rds","coefficients.csv","active_pairs.csv")))),
      paste("serialized leaf artifacts missing",key,leaf))
    coefficients <- read.csv(file.path(d,"coefficients.csv"),
                             stringsAsFactors=FALSE)
    active <- read.csv(file.path(d,"active_pairs.csv"),
                       stringsAsFactors=FALSE)
    check(nrow(coefficients)>0L &&
          all(c("ordered_pair","beta") %in% names(coefficients)) &&
          nrow(active)==sum(coefficients$beta!=0),
          paste("active/coefficients mismatch",key,leaf))
    if (!length(fold_idx)) next
    fragments <- strsplit(
      coefficients$ordered_pair," > ",fixed=TRUE)
    check(all(lengths(fragments)==2L),
          paste("unrecognized leaf pair labels",key,leaf))
    pairs <- data.frame(
      gene_left=vapply(fragments,function(x)x[[1L]],character(1)),
      gene_right=vapply(fragments,function(x)x[[2L]],character(1)),
      ordered_pair=coefficients$ordered_pair,
      stringsAsFactors=FALSE)
    mat <- encode_pair_matrix(
      Xte[fold_idx,,drop=FALSE],pairs,m=10L)
    fit <- readRDS(file.path(d,"cvfit.rds"))
    recovered[fold_idx] <- as.numeric(
      predict(fit,newx=mat,s="lambda.min",type="response"))
  }
  check(!anyNA(recovered) && all(is.finite(recovered)) &&
        isTRUE(all.equal(recovered,pred$prob,tolerance=1e-11)),
        paste("serialized leaf model replay mismatch",key))
  y01 <- as.numeric(yte==levels(yte)[2L])
  clipped <- pmin(pmax(recovered,1e-15),1-1e-15)
  ll <- -mean(y01*log(clipped)+(1-y01)*log1p(-clipped))
  check(abs(ll-as.numeric(met$logloss))<=5.1e-5,
        paste("replayed logloss mismatch",key))
  if (entry$method=="RLMT_DC_v0") {
    decision <- jsonlite::read_json(
      file.path(base,"dc","decision.json"),
      simplifyVector=TRUE)
    check(isTRUE(decision$confirmation$accepted) &&
          decision$status=="ACCEPTED" &&
          decision$confirmation$gain>0 &&
          decision$confirmation$lower_gain>0 &&
          isTRUE(all.equal(
            as.numeric(decision$confirmation$lower_gain),
            as.numeric(tree$split$confirmation_lower_gain),
            tolerance=1e-9)),
          paste("disjoint confirmation audit mismatch",key))
  } else {
    check(tree$split$mean_gain>0 &&
          tree$split$mean_gain-tree$split$se_gain>0,
          paste("v1 split violates 1-SE",key))
  }
  records[[j]] <- data.frame(
    dataset=entry$dataset,method=entry$method,
    n_test=nrow(Xte),n_left=sum(expected_leaf==2L),
    n_right=sum(expected_leaf==3L),
    saved_logloss=as.numeric(met$logloss),
    replayed_logloss=ll,stringsAsFactors=FALSE)
}
dir.create(report,recursive=TRUE,showWarnings=FALSE)
if (length(records)) {
  summary <- do.call(rbind,records)
  write.csv(summary,file.path(report,"accepted_models.csv"),
            row.names=FALSE)
  print(summary,row.names=FALSE,digits=8)
} else {
  write.csv(data.frame(dataset=character(),method=character()),
            file.path(report,"accepted_models.csv"),row.names=FALSE)
}
writeLines(c(
  "Read-only replay for all actually accepted two-leaf models, not forced DEV.",
  "Original published ARFF -> fixed first fold -> route -> cvfit.rds",
  "-> pair matrix -> probability -> independently recomputed test logloss.",
  "No parameter changes, fit calls, new models or outer-test optimization."
),file.path(report,"README.txt"))
cat("ACCEPTED_MODEL_COUNT=",length(picked),"\n",sep="")
cat("PUBLISHED_ACCEPTED_REPLAY=PASS\n")
cat("REPORT=",report,"\n",sep="")
