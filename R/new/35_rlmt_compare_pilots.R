# Compare the original 2x3 RLMT pilot with the stratified-CV rerun.
# Correctly reconstruct original binary class levels from each ARFF.
# No training, changes to core, or modification of existing model artifacts.
# Usage: Rscript --vanilla R/new/35_rlmt_compare_pilots.R

suppressPackageStartupMessages(library(jsonlite))
source("R/new/07_world_protocol_utils.R")

old_root <- "results/rlmt_gate3_pilot_2x3"
new_root <- "results/rlmt_gate3_stratcv_2x3"
roots <- list(old = old_root, new = new_root)

for (label in names(roots)) {
  root <- roots[[label]]
  for (f in c("config.json", "status.txt")) {
    if (!file.exists(file.path(root, f)))
      stop("Missing ", label, " pilot artifact: ", f)
  }
  status <- trimws(readLines(file.path(root, "status.txt"), warn = FALSE))
  if (!identical(status, "RUNNER_EXIT=0")) {
    stop(label, " pilot not confirmed completed: ", paste(status, collapse = " / "))
  }
}

# Ensure that only the output directory differs between pilot configurations.
cfg_old_raw <- jsonlite::read_json(file.path(old_root, "config.json"),
                                   simplifyVector = FALSE)
cfg_new_raw <- jsonlite::read_json(file.path(new_root, "config.json"),
                                   simplifyVector = FALSE)
cfg_old_raw$out_root <- cfg_new_raw$out_root
if (!identical(cfg_old_raw, cfg_new_raw))
  stop("Pilot configurations differ in fields other than out_root")

cfg <- jsonlite::read_json(file.path(new_root, "config.json"),
                           simplifyVector = TRUE)
num_worlds <- if (is.data.frame(cfg$worlds)) nrow(cfg$worlds) else length(cfg$worlds)
if (num_worlds != 1L) stop("Expected exactly one world")
world <- if (is.data.frame(cfg$worlds)) as.list(cfg$worlds[1L, , drop = FALSE]) else cfg$worlds[[1L]]
datasets <- as.character(world$datasets_include[[1L]])
if (!identical(sort(datasets), c("chen-2002", "singh-2002")) ||
    as.integer(world$n_folds) != 5L ||
    as.integer(world$n_repeats) != 1L ||
    as.integer(world$max_splits_per_dataset) != 3L ||
    as.integer(cfg$models$seed) != 42L)
  stop("Unexpected original pilot design")

expected_code_sha <- "7b9f716c932a65fdd359c1f5e83facae13b3ff73"
sha_file <- file.path(new_root, "CODE_SHA.txt")
if (!file.exists(sha_file) ||
    !identical(trimws(readLines(sha_file, warn = FALSE)[1L]),
               expected_code_sha))
  stop("New pilot is not identified as the tested stratified CV commit")

models <- list(
  rlr = c("rlr", "RLR_K10_m10_d2_pool1000"),
  l1 = c("rlmt", "RLMT_K10_m10_d2_pool1000_L1"),
  l2 = c("rlmt", "RLMT_K10_m10_d2_pool1000_L2")
)
isclose <- function(x, y, tol = 1e-12) {
  is.numeric(x) && is.numeric(y) &&
    length(x) == length(y) &&
    !anyNA(x) && !anyNA(y) &&
    all(is.finite(x)) && all(is.finite(y)) &&
    all(abs(x - y) <= tol + tol * pmax(abs(x), abs(y)))
}
require_true <- function(x, msg) if (!isTRUE(x)) stop(msg)
read_metrics <- function(d) jsonlite::read_json(
  file.path(d, "metrics_test.json"), simplifyVector = TRUE)
read_csv <- function(d, name) {
  file <- file.path(d, name)
  if (identical(name, "predictions_test.csv")) {
    # read.csv otherwise coerces numeric-looking class labels to numbers.
    return(read.csv(file, stringsAsFactors = FALSE,
                    colClasses = c(sample_id = "character",
                                   y_true = "character",
                                   y_pred = "character")))
  }
  read.csv(file, stringsAsFactors = FALSE)
}
close_json <- function(recomputed, stored, name) {
  require_true(length(stored) == 1L && is.finite(stored) &&
               abs(recomputed - stored) <= 5.1e-5,
               paste("JSON rounding or metric mismatch:", name,
                     sprintf(" recalculated=%.10f stored=%.10f",recomputed,stored)))
}
rows <- list()
out_i <- 0L
invalid_rows <- list()
inv_i <- 0L
positive_label_map <- list()

for (dataset in datasets) {
  rawpath <- file.path(world$comp_dir, paste0(dataset, ".arff"))
  require_true(file.exists(rawpath), paste("Missing raw CompCancer:",rawpath))
  xy <- read_xy_file(rawpath)
  require_true(length(levels(xy$y)) == 2L,
               paste("Expected binary classes:", dataset))
  negative_label <- levels(xy$y)[1L]
  positive_label <- levels(xy$y)[2L]
  positive_label_map[[dataset]] <- c(negative = negative_label,
                                    positive = positive_label)
  cat(sprintf("CLASS_LEVELS %s negative=%s positive=%s\n",
              dataset, negative_label, positive_label))
  outer <- make_repeated_stratified_folds(xy$y, n_folds=5L,
                                           n_repeats=1L, seed=42L)

  for (j in seq_len(3L)) {
    fold <- outer[[j]]
    expected_ids <- rownames(xy$X)[fold$test_idx]
    true_labels <- as.character(xy$y[fold$test_idx])
    y01 <- as.numeric(xy$y[fold$test_idx] == positive_label)
    Xtest <- xy$X[fold$test_idx, , drop=FALSE]
    old_rlr <- NULL

    for (version in names(roots)) {
      root <- roots[[version]]
      base <- file.path(root,"worlds",world$name,"datasets",
                        dataset,"splits",fold$split_id)
      dirs <- lapply(models, function(m)
        file.path(base, m[[1L]],m[[2L]]))
      for (name in names(dirs)) {
        dir <- dirs[[name]]
        for (fn in c("metrics_test.json","predictions_test.csv"))
          require_true(file.exists(file.path(dir,fn)),
                       paste("Missing",version,dataset,fold$split_id,name,fn))
        require_true(!file.exists(file.path(dir,"error.json")),
                     paste("Model ERROR",version,dataset,fold$split_id,name))
      }
      pred <- lapply(dirs,function(d)read_csv(d,"predictions_test.csv"))
      metrics <- lapply(dirs,read_metrics)
      for (name in names(pred)) {
        p <- pred[[name]]
        required <- c("sample_id","y_true","y_pred","prob")
        require_true(all(required %in% names(p)),
                     paste("Prediction CSV schema",version,dataset,fold$split_id,name))
        require_true(nrow(p) == length(expected_ids) &&
                     identical(p$sample_id,expected_ids) &&
                     identical(p$y_true,true_labels) &&
                     all(p$y_pred %in% levels(xy$y)) &&
                     all(is.finite(p$prob)) &&
                     all(p$prob>=0 & p$prob<=1),
                     paste("Wrong outer test labels, IDs or probabilities",
                           version,dataset,fold$split_id,name))
        predicted_labels <- ifelse(p$prob>=0.5,positive_label,negative_label)
        require_true(identical(p$y_pred,predicted_labels),
                     paste("Prediction threshold/class mismatch",
                           version,dataset,fold$split_id,name))
        pr <- pmin(pmax(p$prob,1e-15),1-1e-15)
        ll <- -mean(y01*log(pr)+(1-y01)*log1p(-pr))
        brier <- mean((p$prob-y01)^2)
        close_json(ll,metrics[[name]]$logloss,
                   paste(version,dataset,fold$split_id,name,"logloss"))
        close_json(brier,metrics[[name]]$brier,
                   paste(version,dataset,fold$split_id,name,"brier"))
      }
      require_true(identical(pred$rlr$y_pred,pred$l1$y_pred) &&
                   isclose(pred$rlr$prob,pred$l1$prob) &&
                   identical(pred$rlr$sample_id,pred$l1$sample_id),
                   paste("RLR vs L1 parity mismatch",version,dataset,fold$split_id))
      for (key in c("auc","f1","logloss","brier","ece",
                    "global_active_pairs","unique_genes")) {
        require_true(isclose(as.numeric(metrics$rlr[[key]]),
                             as.numeric(metrics$l1[[key]])),
                     paste("RLR/L1 metric mismatch",version,dataset,fold$split_id,key))
      }

      if (version == "old") old_rlr <- pred$rlr
      if (version == "new") {
        require_true(identical(old_rlr$y_pred,pred$rlr$y_pred) &&
                     isclose(old_rlr$prob,pred$rlr$prob),
                     paste("Frozen RLR changed between pilots",dataset,fold$split_id))
      }

      treefile <- file.path(dirs$l2,"tree.json")
      nodesfile <- file.path(dirs$l2,"nodes.csv")
      routesfile <- file.path(dirs$l2,"test_routes.csv")
      require_true(all(file.exists(c(treefile,nodesfile,routesfile))),
                   paste("Missing tree/node/routes",version,dataset,fold$split_id))
      tree <- jsonlite::read_json(treefile,simplifyVector=TRUE)
      nodes <- read.csv(nodesfile)
      routes <- read.csv(routesfile,stringsAsFactors=FALSE)
      accepted <- isTRUE(tree$accepted_split)
      leaves <- as.integer(tree$size)
      require_true(leaves %in% c(1L,2L) &&
                   identical(accepted,leaves==2L) &&
                   nrow(routes)==length(expected_ids) &&
                   identical(routes$sample_id,expected_ids) &&
                   !anyDuplicated(routes$sample_id) &&
                   as.integer(metrics$l2$realized_leaves)==leaves &&
                   sum(nodes$is_leaf)==leaves,
                   paste("Tree/routes mismatch",version,dataset,fold$split_id))

      lookfile <- file.path(dirs$l2,"lookahead.csv")
      scrfile <- file.path(dirs$l2,"screen_top.csv")
      scr <- if (file.exists(scrfile)) read.csv(scrfile,
                          stringsAsFactors=FALSE) else NULL
      look <- if (file.exists(lookfile)) read.csv(lookfile,
                          stringsAsFactors=FALSE) else NULL
      if (xor(is.null(scr),is.null(look)))
        stop("Only one of screening/lookahead exists: ",
             version," ",dataset," ",fold$split_id)
      n_candidates <- if (is.null(look)) 0L else nrow(look)
      valid <- 0L
      positive <- 0L
      stable <- 0L
      best_lower <- NA_real_
      if (!is.null(look)) {
        require_true(n_candidates==nrow(scr) &&
                     n_candidates<=20L &&
                     identical(as.character(scr$gene_a),as.character(look$gene_a)) &&
                     identical(as.character(scr$gene_b),as.character(look$gene_b)),
                     paste("Screen/lookahead misaligned",version,dataset,fold$split_id))
        v <- look[which(!is.na(look$valid) & look$valid),,drop=FALSE]
        valid <- nrow(v)
        if (valid) {
          require_true(all(is.finite(v$mean_gain)) &&
                       all(is.finite(v$lower_gain)),
                       "Non-finite gain for valid candidate")
          positive <- sum(v$mean_gain>0)
          stable <- sum(v$lower_gain>0 & v$mean_gain>0)
          best_lower <- max(v$lower_gain)
        }
        invalid <- look[which(is.na(look$valid) | !look$valid),,drop=FALSE]
        if (nrow(invalid)) {
          invalid$dataset <- dataset
          invalid$split_id <- fold$split_id
          invalid$version <- version
          inv_i <- inv_i + 1L
          invalid_rows[[inv_i]] <- invalid[,c("dataset","split_id","version",
                                               "gene_a","gene_b","rejection_reason")]
        }
      }

      if (accepted) {
        require_true(!is.null(look) && stable>0L &&
                     leaves==2L && nrow(nodes)==3L &&
                     all(routes$leaf_id %in% c(2L,3L)),
                     paste("Accepted split invalid",version,dataset,fold$split_id))
        a <- as.character(tree$split$gene_a)
        b <- as.character(tree$split$gene_b)
        require_true(a %in% colnames(Xtest) && b %in% colnames(Xtest),
                     "Split genes missing")
        actual_routes <- ifelse(Xtest[[a]]>Xtest[[b]],2L,3L)
        require_true(identical(as.integer(routes$leaf_id),
                               as.integer(actual_routes)),
                     paste("Accepted relational routing mismatch",
                           version,dataset,fold$split_id))
        require_true(is.finite(tree$split$mean_gain - tree$split$se_gain) &&
                     tree$split$mean_gain-tree$split$se_gain>0,
                     paste("Accepted split violating 1-SE",
                           version,dataset,fold$split_id))
        winner_match <- look$gene_a==a & look$gene_b==b &
          look$valid & is.finite(look$lower_gain) & look$lower_gain>0
        require_true(any(winner_match),
                     paste("Winning split not in passing audit",
                           version,dataset,fold$split_id))
        for (leaf in c("node_2","node_3")) {
          leafdir <- file.path(dirs$l2,"leaf_models",leaf)
          require_true(all(file.exists(file.path(leafdir,
                      c("cvfit.rds","coefficients.csv","active_pairs.csv")))),
                       paste("Missing accepted leaf model",leaf))
        }
      } else {
        require_true(leaves==1L && nrow(nodes)==1L &&
                     all(routes$leaf_id==1L) &&
                     isclose(pred$l2$prob,pred$rlr$prob),
                     paste("Rejected tree differs from RLR",
                           version,dataset,fold$split_id))
      }

      out_i <- out_i + 1L
      rows[[out_i]] <- data.frame(
        dataset=dataset, split_id=fold$split_id, version=version,
        positive_class=positive_label, test_n=length(y01),
        candidates=n_candidates, valid=valid, invalid=n_candidates-valid,
        mean_gain_positive=positive, lower_gain_positive=stable,
        best_lower_gain=best_lower, accepted=accepted,
        leaves=leaves, l2_logloss=as.numeric(metrics$l2$logloss),
        rlr_logloss=as.numeric(metrics$rlr$logloss),
        l2_auc=as.numeric(metrics$l2$auc),
        rlr_auc=as.numeric(metrics$rlr$auc),
        stringsAsFactors=FALSE)
    }
  }
}
tab <- do.call(rbind,rows)
require_true(nrow(tab)==12L && sum(tab$version=="old")==6L &&
             sum(tab$version=="new")==6L,
             "Expected exactly 12 pilot evaluations")
old <- tab[tab$version=="old",,drop=FALSE]
new <- tab[tab$version=="new",,drop=FALSE]
require_true(sum(old$candidates)==120L && sum(old$valid)==96L &&
             sum(old$accepted)==0L,
             "Old pilot failed historical 120/96/zero-split check")

cat("\n=== FOLD-BY-FOLD COMPARISON ===\n")
print(tab[,c("dataset","split_id","version","candidates","valid",
              "mean_gain_positive","lower_gain_positive","leaves",
              "l2_logloss")],row.names=FALSE,digits=6)

cat("\n=== OLD / STRATIFIED-CV SUMMARY ===\n")
for (version in c("old","new")) {
  d <- tab[tab$version==version,,drop=FALSE]
  cat(sprintf("%-3s complete=%d candidates=%d valid=%d invalid=%d positive_mean=%d positive_1SE=%d accepted=%d\n",
              version,nrow(d),sum(d$candidates),sum(d$valid),
              sum(d$invalid),sum(d$mean_gain_positive),
              sum(d$lower_gain_positive),sum(d$accepted)))
}
if (length(invalid_rows)) {
  iv <- do.call(rbind,invalid_rows)
  cat("\n=== INVALID REASONS ===\n")
  for (version in c("old","new")) {
    v <- iv[iv$version==version,,drop=FALSE]
    cat("VERSION:",version,"\n")
    reasons <- sort(table(v$rejection_reason),decreasing=TRUE)
    for (i in seq_along(reasons))
      cat(sprintf("  %3d %s\n",as.integer(reasons[i]),names(reasons)[i]))
  }
} else iv <- data.frame()

# Output is separate from both model result trees; never overwrite.
reportdir <- file.path(new_root,"comparison_old_vs_stratcv_v1")
if (dir.exists(reportdir)) stop("Refusing to overwrite: ",reportdir)
dir.create(reportdir,recursive=TRUE,showWarnings=FALSE)
write.csv(tab,file.path(reportdir,"fold_comparison.csv"),row.names=FALSE)
if (nrow(iv)) write.csv(iv,file.path(reportdir,"invalid_candidates.csv"),
                        row.names=FALSE)
writeLines(c(
  "Validated BOTH pilot configurations and all 12 completed evaluations.",
  "TRUE classes reconstructed from foreign::read.arff / factor levels.",
  "Outer-test IDs, labels, probabilities, logloss, Brier and L1 parity checked.",
  "RLR probabilities must match old versus new pilot exactly within 1e-12.",
  "Selected split must pass mean_gain - SE > 0; otherwise L2 must match RLR.",
  "Pilot comparison is technical, not a claim of predictive superiority."
), file.path(reportdir,"README.txt"))

cat("PILOT_COMPARISON_AUDIT=PASS\n")
cat("REPORT=",reportdir,"\n",sep="")
