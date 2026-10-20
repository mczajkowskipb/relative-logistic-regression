# GATE-4 Gordon end-of-run scientific / artifact check. No new fitting.
# Usage: Rscript --vanilla R/new/38_rlmt_gate4_check.R
suppressPackageStartupMessages(library(jsonlite))
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")

root <- "results/rlmt_gate4_gordon_frozen_5fold"
cfg <- jsonlite::read_json("configs/rlmt/gate4_gordon_frozen.json",
                            simplifyVector=TRUE)
require_true <- function(ok,msg) if (!isTRUE(ok)) stop("GATE4 AUDIT FAIL: ",msg)
require_true(identical(cfg$out_root,root),"unexpected experiment output path")
require_true(length(cfg$worlds$datasets_include)==1L &&
             identical(as.character(cfg$worlds$datasets_include[1]),"gordon-2002") &&
             as.integer(cfg$worlds$n_folds)==5L &&
             as.integer(cfg$worlds$n_repeats)==1L &&
             as.integer(cfg$worlds$max_splits_per_dataset)==5L,
             "not the locked Gordon five-fold protocol")
require_true(file.exists(file.path(root,"status.txt")) &&
             identical(trimws(readLines(file.path(root,"status.txt"))[1]),
                       "RUNNER_EXIT=0"),"runner status is not successful")
require_true(file.exists(file.path(root,"CODE_SHA.txt")),
             "missing recorded code SHA")
expected_config <- normalizePath("configs/rlmt/gate4_gordon_frozen.json")
run_config <- file.path(root,"config","used_config.json")
require_true(file.exists(run_config),"missing used_config.json")
stored <- jsonlite::read_json(run_config,simplifyVector=TRUE)
expected <- jsonlite::read_json(expected_config,simplifyVector=TRUE)
# jsonlite::write_json(auto_unbox=TRUE) can turn singleton JSON arrays
# into scalars in the persisted config. Check semantic run controls instead.
flatten <- function(x) unlist(x,use.names=FALSE)
for (key in c("out_root")) {
  require_true(identical(as.character(flatten(stored[[key]])),
                         as.character(flatten(expected[[key]]))),
               paste("used_config mismatch in",key))
}
for (key in c("seed","run_baselines","run_rlr","run_rlmt",
              "top_pool_values","alpha","lambda_choice","rlr_grid")) {
  require_true(isTRUE(all.equal(flatten(stored$models[[key]]),
                               flatten(expected$models[[key]]))),
               paste("used_config model mismatch in",key))
}
for (key in names(expected$models$rlmt)) {
  require_true(isTRUE(all.equal(flatten(stored$models$rlmt[[key]]),
                               flatten(expected$models$rlmt[[key]]))),
               paste("used_config RLMT mismatch in",key))
}
w_exp <- if(is.data.frame(expected$worlds))
  as.list(expected$worlds[1L,,drop=FALSE]) else expected$worlds[[1L]]
w_act <- if(is.data.frame(stored$worlds))
  as.list(stored$worlds[1L,,drop=FALSE]) else stored$worlds[[1L]]
for (key in c("name","protocol","comp_dir","datasets_include",
              "n_folds","n_repeats","max_splits_per_dataset")) {
  require_true(isTRUE(all.equal(flatten(w_act[[key]]),
                               flatten(w_exp[[key]]))),
               paste("used_config world mismatch in",key))
}
xy <- read_xy_file(file.path(cfg$worlds$comp_dir,"gordon-2002.arff"))
require_true(nrow(xy$X)==181L && ncol(xy$X)==1626L &&
             identical(as.integer(table(xy$y)),c(150L,31L)) &&
             identical(levels(xy$y),c("AD","MPM")),
             "Gordon source schema or class counts changed")
# Derive the expected digest from the publication manifest instead of
# maintaining a second, error-prone hardcoded copy.
manifest_lines <- readLines("checksums/DATA_SHA256_EXPECTED.txt", warn=FALSE)
manifest_matches <- grep(
  "^[[:xdigit:]]{64}[[:space:]]+data/CompCancer/gordon-2002[.]arff$",
  trimws(manifest_lines), value=TRUE)
require_true(length(manifest_matches)==1L,
             "missing or ambiguous publication Gordon checksum")
expected_sha <- tolower(sub("^([[:xdigit:]]{64}).*$","\\1",
                            manifest_matches[[1L]]))
require_true(identical(expected_sha,"42bb0900f7a22851244c27243b626514917c0674d1a3c63dbf77bcab43361bab"),
             "publication manifest SHA is not the frozen Gordon input")
sha_output <- system2("sha256sum",
                     shQuote(file.path(cfg$worlds$comp_dir,
                                       "gordon-2002.arff")),
                     stdout=TRUE, stderr=TRUE)
require_true(length(sha_output)==1L &&
             grepl("^[[:xdigit:]]{64}[[:space:]]+",sha_output[1L]),
             "sha256sum failed on Gordon input")
sha <- tolower(sub("^([[:xdigit:]]{64}).*$","\\1",sha_output[[1L]]))
require_true(identical(sha,expected_sha),
             "source SHA mismatch against publication manifest")

pdirs <- list(
  rlr=c("rlr","RLR_K10_m10_d2_pool1000"),
  l1=c("rlmt","RLMT_K10_m10_d2_pool1000_L1"),
  l2=c("rlmt","RLMT_K10_m10_d2_pool1000_L2"))
load_preds <- function(d) read.csv(file.path(d,"predictions_test.csv"),
                                   stringsAsFactors=FALSE,
                                   colClasses=c(sample_id="character",
                                                y_true="character",
                                                y_pred="character"))
read_metric <- function(d) jsonlite::read_json(file.path(d,"metrics_test.json"),
                                              simplifyVector=TRUE)
close_num <- function(a,b,tol=1e-12) {
  length(a)==length(b) && !anyNA(a) && !anyNA(b) &&
    all(is.finite(a)) && all(is.finite(b)) &&
    isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol))
}
metric_from_prob <- function(y01,p) {
  clipped <- pmin(pmax(p,1e-15),1-1e-15)
  list(logloss=-mean(y01*log(clipped)+(1-y01)*log1p(-clipped)),
       brier=mean((p-y01)^2))
}
validate_metric <- function(calc,actual,label) {
  require_true(length(actual)==1L && is.finite(actual) &&
               abs(calc-actual)<=5.1e-5,
               paste("metric mismatch",label,
                     sprintf("computed %.9f json %.9f",calc,actual)))
}
folds <- make_repeated_stratified_folds(xy$y,n_folds=5L,n_repeats=1L,seed=42L)
results <- list()
invalid_reasons <- character()

for (i in seq_along(folds)) {
  f <- folds[[i]]
  fold <- f$split_id
  base <- file.path(root,"worlds","compcancer_cv","datasets",
                    "gordon-2002","splits",fold)
  dirs <- lapply(pdirs,function(e) file.path(base,e[1],e[2]))
  for (method in names(dirs)) {
    d <- dirs[[method]]
    for (name in c("metrics_test.json","predictions_test.csv")) {
      require_true(file.exists(file.path(d,name)),
                   paste(fold,method,"missing",name))
    }
    require_true(!file.exists(file.path(d,"error.json")),
                 paste(fold,method,"error.json exists"))
  }
  p <- lapply(dirs,load_preds)
  m <- lapply(dirs,read_metric)
  xte <- xy$X[f$test_idx,,drop=FALSE]
  yte <- xy$y[f$test_idx]
  y01 <- as.numeric(yte=="MPM")
  ids <- rownames(xte)

  for (method in names(p)) {
    z <- p[[method]]
    require_true(nrow(z)==length(ids) &&
                 identical(z$sample_id,ids) &&
                 identical(z$y_true,as.character(yte)) &&
                 all(is.finite(z$prob)) &&
                 all(z$prob>=0 & z$prob<=1) &&
                 identical(z$y_pred,
                           ifelse(z$prob>=0.5,"MPM","AD")),
                 paste(fold,method,"prediction sample/class/threshold mismatch"))
    computed <- metric_from_prob(y01,z$prob)
    validate_metric(computed$logloss,m[[method]]$logloss,
                    paste(fold,method,"logloss"))
    validate_metric(computed$brier,m[[method]]$brier,
                    paste(fold,method,"brier"))
  }
  require_true(identical(p$rlr$y_pred,p$l1$y_pred) &&
               close_num(p$rlr$prob,p$l1$prob) &&
               close_num(p$rlr$score,p$l1$score),
               paste(fold,"RLR/L1 prediction parity"))
  for (metric in c("auc","f1","logloss","brier","ece",
                   "global_active_pairs","unique_genes")) {
    require_true(close_num(as.numeric(m$rlr[[metric]]),
                           as.numeric(m$l1[[metric]])),
                 paste(fold,"RLR/L1 metric parity",metric))
  }
  for (file in c("coefficients.csv","active_pairs.csv")) {
    a <- read.csv(file.path(dirs$rlr,file),stringsAsFactors=FALSE)
    b <- read.csv(file.path(dirs$l1,file),stringsAsFactors=FALSE)
    require_true(identical(a,b),paste(fold,"RLR/L1 parity",file))
  }

  ldir <- dirs$l2
  for (file in c("tree.json","nodes.csv","test_routes.csv"))
    require_true(file.exists(file.path(ldir,file)),
                 paste(fold,"missing",file))
  tree <- jsonlite::read_json(file.path(ldir,"tree.json"),
                              simplifyVector=TRUE)
  nodes <- read.csv(file.path(ldir,"nodes.csv"))
  routes <- read.csv(file.path(ldir,"test_routes.csv"),
                     stringsAsFactors=FALSE,
                     colClasses=c(sample_id="character"))
  realized <- as.integer(tree$size)
  accepted <- isTRUE(tree$accepted_split)
  require_true(realized %in% c(1L,2L) &&
               accepted==(realized==2L) &&
               as.integer(m$l2$realized_leaves)==realized &&
               sum(nodes$is_leaf)==realized &&
               identical(routes$sample_id,ids) &&
               nrow(routes)==length(ids) &&
               !anyDuplicated(routes$sample_id),
               paste(fold,"tree/nodes/routes structure mismatch"))

  screen_path <- file.path(ldir,"screen_top.csv")
  audit_path <- file.path(ldir,"lookahead.csv")
  scr <- if(file.exists(screen_path)) read.csv(screen_path,
                         stringsAsFactors=FALSE) else NULL
  look <- if(file.exists(audit_path)) read.csv(audit_path,
                         stringsAsFactors=FALSE) else NULL
  require_true(is.null(scr)==is.null(look),
               paste(fold,"screen/lookahead artifacts mismatch"))
  valid <- 0L
  positive <- 0L
  stable <- 0L
  if (!is.null(look)) {
    require_true(nrow(look)==nrow(scr) && nrow(look)<=20L &&
                 identical(as.character(look$gene_a),as.character(scr$gene_a)) &&
                 identical(as.character(look$gene_b),as.character(scr$gene_b)),
                 paste(fold,"lookahead candidate accounting mismatch"))
    ok <- look[which(!is.na(look$valid)&look$valid),,drop=FALSE]
    valid <- nrow(ok)
    if (valid>0) {
      positive <- sum(ok$mean_gain>0)
      stable <- sum(ok$lower_gain>0 & ok$mean_gain>0)
      require_true(all(is.finite(ok$mean_gain)) &&
                   all(is.finite(ok$se_gain)) &&
                   all(is.finite(ok$lower_gain)),
                   paste(fold,"nonfinite lookahead gain"))
    }
    bad <- look[which(is.na(look$valid) | !look$valid),,drop=FALSE]
    if (nrow(bad)) invalid_reasons <- c(invalid_reasons,
                                        as.character(bad$rejection_reason))
  }
  if (!accepted) {
    require_true(realized==1L && nrow(nodes)==1L &&
                 all(routes$leaf_id==1L) &&
                 close_num(p$rlr$prob,p$l2$prob),
                 paste(fold,"stopped L2 differs from frozen RLR"))
  } else {
    require_true(!is.null(look) && stable>0L &&
                 nrow(nodes)==3L && all(routes$leaf_id %in% c(2L,3L)),
                 paste(fold,"accepted tree without valid 1-SE evidence"))
    a <- as.character(tree$split$gene_a)
    b <- as.character(tree$split$gene_b)
    require_true(a %in% names(xte) && b %in% names(xte) &&
                 identical(as.integer(routes$leaf_id),
                           as.integer(ifelse(xte[[a]]>xte[[b]],2L,3L))),
                 paste(fold,"relational routing mismatch"))
    winner <- look[which(look$valid & look$gene_a==a &
                         look$gene_b==b & look$lower_gain>0),,drop=FALSE]
    require_true(nrow(winner)==1L &&
                 tree$split$mean_gain-tree$split$se_gain>0,
                 paste(fold,"accepted split fails 1-SE"))
    for (leaf in c(2L,3L)) {
      ix <- which(routes$leaf_id==leaf)
      ld <- file.path(ldir,"leaf_models",paste0("node_",leaf))
      for (fn in c("cvfit.rds","coefficients.csv","active_pairs.csv"))
        require_true(file.exists(file.path(ld,fn)),
                     paste(fold,"leaf",leaf,"missing",fn))
      coeff <- read.csv(file.path(ld,"coefficients.csv"),
                        stringsAsFactors=FALSE)
      active <- read.csv(file.path(ld,"active_pairs.csv"),
                         stringsAsFactors=FALSE)
      require_true(nrow(coeff)>=1L && nrow(active)==sum(coeff$beta!=0),
                   paste(fold,"leaf coefficient/activity mismatch"))
      if (length(ix)) {
        pairparts <- strsplit(coeff$ordered_pair," > ",fixed=TRUE)
        require_true(all(lengths(pairparts)==2L),
                     paste(fold,"unrecognized stored pair labels"))
        pairs <- data.frame(
          gene_left=vapply(pairparts,function(x)x[1L],character(1)),
          gene_right=vapply(pairparts,function(x)x[2L],character(1)),
          ordered_pair=coeff$ordered_pair,stringsAsFactors=FALSE)
        R <- encode_pair_matrix(xte[ix,,drop=FALSE],pairs,
                                 m=as.integer(cfg$models$rlmt$m))
        fit <- readRDS(file.path(ld,"cvfit.rds"))
        recovered <- as.numeric(predict(
          fit,newx=R,s=cfg$models$rlmt$lambda_choice,type="response"))
        require_true(close_num(recovered,p$l2$prob[ix]),
                     paste(fold,"leaf",leaf,"model replay mismatch"))
      }
    }
  }

  results[[i]] <- data.frame(
    fold=fold,n_test=length(ids), n_test_AD=sum(yte=="AD"),
    n_test_MPM=sum(yte=="MPM"),
    screened=if(is.null(look))0L else nrow(look),
    valid=valid,invalid=if(is.null(look))0L else nrow(look)-valid,
    positive_mean=positive,positive_1SE=stable,
    accepted_split=accepted,realized_leaves=realized,
    rlr_logloss=as.numeric(m$rlr$logloss),
    l2_logloss=as.numeric(m$l2$logloss),
    delta_logloss=as.numeric(m$l2$logloss)-as.numeric(m$rlr$logloss),
    rlr_brier=as.numeric(m$rlr$brier),
    l2_brier=as.numeric(m$l2$brier),
    delta_brier=as.numeric(m$l2$brier)-as.numeric(m$rlr$brier),
    rlr_auc=as.numeric(m$rlr$auc), l2_auc=as.numeric(m$l2$auc),
    rlr_f1=as.numeric(m$rlr$f1), l2_f1=as.numeric(m$l2$f1),
    stringsAsFactors=FALSE)
}
tab <- do.call(rbind,results)
require_true(nrow(tab)==5L && sum(tab$n_test)==181L,
             "Expected exactly 5 folds covering 181 test samples")
cat("=== GATE-4 GORDON 5-FOLD RESULT ===\n")
print(tab[,c("fold","n_test","valid","invalid","positive_mean",
              "positive_1SE","accepted_split","realized_leaves",
              "rlr_logloss","l2_logloss","delta_logloss")],
      row.names=FALSE,digits=6)
cat("\n=== DESCRIPTIVE AGGREGATES (FOLDS NOT INDEPENDENT) ===\n")
cat("Valid candidate assessments: ",sum(tab$valid),"/",
    sum(tab$screened),"\n",sep="")
cat("Accepted two-leaf trees: ",sum(tab$accepted_split),"/5\n",sep="")
cat(sprintf("Mean fold RLR logloss: %.6f\n",mean(tab$rlr_logloss)))
cat(sprintf("Mean fold L2 logloss: %.6f\n",mean(tab$l2_logloss)))
cat(sprintf("Mean paired delta (L2-RLR): %+.6f\n",
            mean(tab$delta_logloss)))
cat(sprintf("Mean paired Brier delta: %+.6f\n",mean(tab$delta_brier)))
if (length(invalid_reasons)) {
  cat("Invalid lookahead reasons:\n")
  print(table(invalid_reasons))
}
report <- file.path(root,"gate4_audit_v1")
require_true(!dir.exists(report),"refusing to overwrite GATE4 audit")
dir.create(report,recursive=TRUE,showWarnings=FALSE)
write.csv(tab,file.path(report,"fold_summary.csv"),row.names=FALSE)
writeLines(c(
  "Validated five Gordon outer folds against original ARFF class labels.",
  "RLR/L1 parity verified; test logloss/Brier independently recomputed.",
  "For accepted L2 splits: 1-SE evidence, route and serialized leaf predictions checked.",
  "Five folds are dependent; this is a historical benchmark, not external validation.",
  paste("CODE_SHA:",trimws(readLines(file.path(root,"CODE_SHA.txt"))[1L])),
  "No fitting or post-hoc tuning performed during this audit."
),file.path(report,"README.txt"))
cat("GATE4_GORDON_AUDIT=PASS\n")
cat("REPORT=",report,"\n",sep="")
