# Read-only reason analysis for Gordon GATE4: no screened split candidates.
# No fitting, no parameter changes, no outer-test predictions or test labels.
# Usage: Rscript --vanilla R/new/39_rlmt_gate4_screen_support.R
suppressPackageStartupMessages(library(jsonlite))
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")

root <- "results/rlmt_gate4_gordon_frozen_5fold"
out <- file.path(root, "gate4_screen_support_audit_v1")
cfg <- jsonlite::read_json("configs/rlmt/gate4_gordon_frozen.json",
                            simplifyVector = TRUE)
check <- function(ok, what) if (!isTRUE(ok))
  stop("GATE4 SUPPORT AUDIT FAIL: ", what, call. = FALSE)
check(!dir.exists(out), "refuse to overwrite prior support audit")
check(identical(cfg$out_root, root), "configuration root")
check(identical(trimws(readLines(file.path(root,"status.txt"))[1]),
                "RUNNER_EXIT=0"), "unfinished runner")
check(file.exists(file.path(root,"gate4_audit_v1","fold_summary.csv")),
      "independent preceding GATE4 result audit missing")
check(identical(trimws(readLines(file.path(root,"CODE_SHA.txt"))[1]),
                "193fab667acfd4ab542ebc51c96cb3a3df31ddba"),
      "the evaluated frozen model commit differs")
conf <- rlmt_defaults(cfg$models$rlmt)
world <- as.list(cfg$worlds[1L,,drop=FALSE])
check(identical(as.character(world$datasets_include[[1L]]),
                "gordon-2002"), "unexpected dataset")
check(as.integer(world$n_folds)==5L &&
      as.integer(world$max_splits_per_dataset)==5L &&
      as.integer(cfg$models$seed)==42L, "unexpected fold definition")
raw <- read_xy_file(file.path(world$comp_dir,"gordon-2002.arff"))
check(nrow(raw$X)==181L && ncol(raw$X)==1626L &&
      identical(levels(raw$y),c("AD","MPM")) &&
      identical(as.integer(table(raw$y)),c(150L,31L)),
      "original data differ")
folds <- make_repeated_stratified_folds(raw$y,5L,1L,42L)

detail <- list()
summary <- list()
count <- 0L

for (foldno in seq_len(5L)) {
  f <- folds[[foldno]]
  base <- file.path(root, "worlds","compcancer_cv","datasets",
                    "gordon-2002","splits",f$split_id)
  l2 <- file.path(base,"rlmt","RLMT_K10_m10_d2_pool1000_L2")
  poolfile <- file.path(base,"pairs_train","top1000","top_1000_pairs.csv")
  pool <- read.csv(poolfile,stringsAsFactors=FALSE)
  check(nrow(pool)==1000L, paste(f$split_id,"pool is not exactly 1000"))
  check(all(c("gene_left","gene_right") %in% names(pool)),
        paste(f$split_id,"invalid saved pair pool columns"))
  check(!file.exists(file.path(l2,"screen_top.csv")) &&
        !file.exists(file.path(l2,"lookahead.csv")),
        paste(f$split_id,"unexpected screening artifacts"))
  tree <- jsonlite::read_json(file.path(l2,"tree.json"),simplifyVector=TRUE)
  check(isFALSE(tree$accepted_split) && tree$size==1L &&
        identical(tree$stop_reason,"no_valid_child_splits"),
        paste(f$split_id,"stopping reason not due to candidate feasibility"))
  X <- raw$X[f$train_idx,,drop=FALSE]
  y <- raw$y[f$train_idx]
  check(identical(read.csv(file.path(base,"samples_train.csv"),
                            colClasses="character")$sample_id,
                  rownames(X)),
        paste(f$split_id,"outer-training sample IDs differ"))
  check(identical(read.csv(file.path(base,"samples_test.csv"),
                            colClasses="character")$sample_id,
                  rownames(raw$X)[f$test_idx]),
        paste(f$split_id,"outer-test sample IDs differ"))
  counts <- as.integer(table(y))
  check(nrow(X)>=2L*conf$min_child_n &&
        all(counts>=2L*conf$min_child_class),
        paste(f$split_id,"root support is insufficient"))
  check(all(pool$gene_left %in% names(X)) &&
        all(pool$gene_right %in% names(X)),
        paste(f$split_id,"pool references missing genes"))

  foldrows <- vector("list",nrow(pool))
  for (j in seq_len(nrow(pool))) {
    z <- X[[pool$gene_left[j]]] > X[[pool$gene_right[j]]]
    check(!anyNA(z),paste(f$split_id,j,"NA in routing"))
    sides <- list(which(z),which(!z))
    sizes <- lengths(sides)
    ctab <- rbind(as.integer(table(y[sides[[1L]]])),
                  as.integer(table(y[sides[[2L]]])))
    colnames(ctab) <- levels(y)
    min_n <- min(sizes)
    min_class <- min(ctab)
    n_gate <- all(sizes >= conf$min_child_n)
    class_gate <- all(ctab >= conf$min_child_class)
    full_gate <- rlmt_children_ok(z,y,conf)
    check(identical(full_gate,n_gate && class_gate),
          paste(f$split_id,j,"audit differs from production gate"))
    foldrows[[j]] <- data.frame(
      split_id=f$split_id,pool_rank=j,
      gene_a=pool$gene_left[j],gene_b=pool$gene_right[j],
      left_n=sizes[1L],right_n=sizes[2L],
      left_AD=ctab[1L,1L],left_MPM=ctab[1L,2L],
      right_AD=ctab[2L,1L],right_MPM=ctab[2L,2L],
      min_child_n=min_n,min_child_class=min_class,
      size_gate=n_gate,class_gate=class_gate,eligible=full_gate,
      stringsAsFactors=FALSE)
  }
  df <- do.call(rbind,foldrows)
  count <- count + 1L
  detail[[count]] <- df
  summary[[count]] <- data.frame(
    split_id=f$split_id,
    training_n=nrow(X),
    training_AD=counts[1L],training_MPM=counts[2L],
    pool_count=nrow(pool),
    passed_size_only=sum(df$size_gate),
    passed_class_only=sum(df$class_gate),
    passed_both=sum(df$eligible),
    failed_size_only=sum(!df$size_gate & df$class_gate),
    failed_class_only=sum(df$size_gate & !df$class_gate),
    failed_both=sum(!df$size_gate & !df$class_gate),
    max_min_child_n=max(df$min_child_n),
    max_min_child_class=max(df$min_child_class),
    max_min_MPM_in_both_sides=max(pmin(df$left_MPM,df$right_MPM)),
    median_min_MPM_in_both_sides=median(pmin(df$left_MPM,df$right_MPM)),
    stringsAsFactors=FALSE)
}
all <- do.call(rbind,detail)
summ <- do.call(rbind,summary)
check(nrow(summ)==5L && nrow(all)==5000L, "unexpected total candidate count")
check(all(!all$eligible) && all(summ$passed_both==0L),
      "recorded no-split outcome contradicts saved pair pool")
cat("=== GORDON GATE4: ROOT-POOL SPLIT SUPPORT ===\n")
print(summ,row.names=FALSE,digits=5)
cat("\n=== TOTAL (SAVED 5x1000 PAIRS, NO REFIT) ===\n")
cat("Pool pairs inspected: ",nrow(all),"\n",sep="")
cat("Size gate passed: ",sum(all$size_gate),"\n",sep="")
cat("Class support gate passed: ",sum(all$class_gate),"\n",sep="")
cat("Both gates passed: ",sum(all$eligible),"\n",sep="")
cat("Size-only failures: ",sum(!all$size_gate & all$class_gate),"\n",sep="")
cat("Class-only failures: ",sum(all$size_gate & !all$class_gate),"\n",sep="")
cat("Both failures: ",sum(!all$size_gate & !all$class_gate),"\n",sep="")
cat("MAX_MIN_CHILD_CLASS_ANY_PAIR=",max(all$min_child_class),"\n",sep="")
cat("MIN_CHILD_CLASS_REQUIRED=",conf$min_child_class,"\n",sep="")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(all,file.path(out,"candidate_support.csv"),row.names=FALSE)
write.csv(summ,file.path(out,"fold_support_summary.csv"),row.names=FALSE)
writeLines(c(
  "Read-only deterministic reconstruction of GATE4 Gordon no-split screening.",
  "All 5000 persisted outer-training top1000 pairs inspected.",
  "Original holdout predictions, metrics, datasets and methodology unchanged.",
  "The screen gate requires >=16 observations per child and >=6 of each class per child.",
  "No what-if relaxation or test-label-based model tuning performed.",
  "See candidate_support.csv for class composition of both potential children."
),file.path(out,"README.txt"))
cat("GATE4_NO_SPLIT_CAUSE_AUDIT=PASS\n")
cat("REPORT=",out,"\n",sep="")
