# RLMT: evaluate the ONLY genuinely accepted split in the stratified 2x3 pilot.
# READ ONLY for model outputs. Reconstructs original Singh classes and outer fold.
# No training, hyperparameter search, test-based model selection, or edits to paper.
# Run: Rscript --vanilla R/new/36_rlmt_accepted_split_audit.R
suppressPackageStartupMessages(library(jsonlite))
source("R/new/07_world_protocol_utils.R")

root <- "results/rlmt_gate3_stratcv_2x3"
cfg <- jsonlite::read_json(file.path(root,"config.json"),simplifyVector=TRUE)
world <- if (is.data.frame(cfg$worlds)) as.list(cfg$worlds[1L,,drop=FALSE]) else cfg$worlds[[1L]]
stopifnot(as.integer(cfg$models$seed)==42L,
          as.integer(world$n_folds)==5L,
          as.integer(world$max_splits_per_dataset)==3L)
dataset <- "singh-2002"
fold_id <- "rep00_fold00"
base <- file.path(root,"worlds",world$name,"datasets",dataset,
                  "splits",fold_id)
rdir <- file.path(base,"rlr","RLR_K10_m10_d2_pool1000")
ldir <- file.path(base,"rlmt","RLMT_K10_m10_d2_pool1000_L2")
dest <- file.path(root,"accepted_singh_fold00_audit_v1")
if (dir.exists(dest)) stop("Refusing to overwrite audit: ",dest)
read_preds <- function(dir) {
  p <- file.path(dir,"predictions_test.csv")
  if (!file.exists(p)) stop("Missing predictions ",p)
  read.csv(p,stringsAsFactors=FALSE,
    colClasses=c(sample_id="character",y_true="character",y_pred="character"))
}
require <- function(ok,msg) if (!isTRUE(ok)) stop("AUDIT FAILURE: ",msg)
rlr <- read_preds(rdir)
l2 <- read_preds(ldir)
tree <- jsonlite::read_json(file.path(ldir,"tree.json"),simplifyVector=TRUE)
nodes <- read.csv(file.path(ldir,"nodes.csv"),stringsAsFactors=FALSE)
routes <- read.csv(file.path(ldir,"test_routes.csv"),stringsAsFactors=FALSE,
                   colClasses=c(sample_id="character"))
look <- read.csv(file.path(ldir,"lookahead.csv"),stringsAsFactors=FALSE)
met_r <- jsonlite::read_json(file.path(rdir,"metrics_test.json"),simplifyVector=TRUE)
met_l <- jsonlite::read_json(file.path(ldir,"metrics_test.json"),simplifyVector=TRUE)
require(isTRUE(tree$accepted_split) && tree$size==2L &&
        nrow(nodes)==3L && met_l$realized_leaves==2L,
        "expected single accepted 2-leaf tree")
xy <- read_xy_file(file.path(world$comp_dir,paste0(dataset,".arff")))
outer <- make_repeated_stratified_folds(xy$y,5L,1L,42L)
fold <- outer[[1L]]
require(identical(fold$split_id,fold_id),"outer fold id")
Xtr <- xy$X[fold$train_idx,,drop=FALSE]
ytr <- xy$y[fold$train_idx]
Xte <- xy$X[fold$test_idx,,drop=FALSE]
yte <- xy$y[fold$test_idx]
test_ids <- rownames(Xte)
require(nrow(rlr)==length(test_ids) && nrow(l2)==length(test_ids) &&
        identical(rlr$sample_id,test_ids) &&
        identical(l2$sample_id,test_ids) &&
        identical(routes$sample_id,test_ids) &&
        identical(as.character(yte),l2$y_true) &&
        identical(as.character(yte),rlr$y_true),
        "holdout sample identity, order or class labels")
require(all(c("prob","score") %in% colnames(l2)) &&
        all(is.finite(l2$prob)) && all(l2$prob>=0 & l2$prob<=1) &&
        all(is.finite(rlr$prob)) && all(rlr$prob>=0 & rlr$prob<=1),
        "invalid predictions")
positive <- levels(xy$y)[2L]
negative <- levels(xy$y)[1L]
y01 <- as.numeric(yte==positive)
require(identical(as.character(l2$y_pred),
        ifelse(l2$prob>=0.5,positive,negative)),"class threshold")
a <- as.character(tree$split$gene_a)
b <- as.character(tree$split$gene_b)
require(a %in% colnames(Xte) && b %in% colnames(Xte),
        "absent split gene")
zt <- Xtr[[a]] > Xtr[[b]]
zv <- Xte[[a]] > Xte[[b]]
expected_route <- ifelse(zv,2L,3L)
require(identical(as.integer(routes$leaf_id),as.integer(expected_route)),
        "route mismatch")
require(identical(as.integer(nodes$n[2:3]),
                  as.integer(c(sum(zt),sum(!zt)))),
        "training leaf count discrepancy")
sel <- look[which(look$valid & look$gene_a==a & look$gene_b==b),,drop=FALSE]
require(nrow(sel)==1L && is.finite(sel$lower_gain) &&
        sel$lower_gain>0 &&
        abs((sel$mean_gain-sel$se_gain)-sel$lower_gain)<1e-10,
        "winner missing or fails 1-SE")
valid <- look[which(look$valid),,drop=FALSE]
winner <- valid[which.max(valid$lower_gain),,drop=FALSE]
require(identical(as.character(winner$gene_a),a) &&
        identical(as.character(winner$gene_b),b),
        "accepted candidate not the best admissible lower bound")

safe_prob <- function(x) pmin(pmax(as.numeric(x),1e-15),1-1e-15)
per_ll <- function(y,p) {p<-safe_prob(p);-(y*log(p)+(1-y)*log1p(-p))}
lr <- per_ll(y01,rlr$prob)
ll <- per_ll(y01,l2$prob)
br <- (rlr$prob-y01)^2
bl <- (l2$prob-y01)^2
close_json <- function(x,v,name) {
  require(is.finite(v) && abs(x-v)<=5.1e-5,
          paste("exported metric mismatch:",name))
}
close_json(mean(lr),met_r$logloss,"RLR logloss")
close_json(mean(ll),met_l$logloss,"L2 logloss")
close_json(mean(br),met_r$brier,"RLR Brier")
close_json(mean(bl),met_l$brier,"L2 Brier")

result <- list()
for (id in c(2L,3L)) {
  tix <- which(expected_route==id)
  trix <- which(if(id==2L)zt else !zt)
  nt <- table(factor(ytr[trix],levels=levels(xy$y)))
  nv <- table(factor(yte[tix],levels=levels(xy$y)))
  fdir <- file.path(ldir,"leaf_models",paste0("node_",id))
  coeff_path <- file.path(fdir,"coefficients.csv")
  active_path <- file.path(fdir,"active_pairs.csv")
  fit_path <- file.path(fdir,"cvfit.rds")
  require(all(file.exists(c(coeff_path,active_path,fit_path))),
          "missing persisted leaf model")
  co <- read.csv(coeff_path,stringsAsFactors=FALSE)
  ac <- read.csv(active_path,stringsAsFactors=FALSE)
  require(nrow(ac)==sum(co$beta!=0) &&
          nrow(ac)==nodes$leaf_globally_active_relation_count[id],
          "leaf active pair mismatch")
  if (length(tix)) {
    relation_parts <- strsplit(co$ordered_pair," > ",fixed=TRUE)
    require(all(lengths(relation_parts)==2L),
            "unparseable stored relation")
    pairs <- data.frame(
      gene_left=vapply(relation_parts,function(x)x[1L],character(1)),
      gene_right=vapply(relation_parts,function(x)x[2L],character(1)),
      ordered_pair=co$ordered_pair,stringsAsFactors=FALSE)
    design <- encode_pair_matrix(Xte[tix,,drop=FALSE],pairs,
                                 m=as.integer(cfg$models$rlmt$m))
    fit <- readRDS(fit_path)
    restored <- as.numeric(predict(fit,newx=design,
                                   s=cfg$models$rlmt$lambda_choice,
                                   type="response"))
    require(length(restored)==length(tix) &&
            isTRUE(all.equal(restored,as.numeric(l2$prob[tix]),
                             tolerance=1e-12)),
            "leaf serialized probability replay mismatch")
  }
  result[[length(result)+1L]] <- data.frame(
    leaf=id,train_n=length(trix),
    train_negative=as.integer(nt[1]),train_positive=as.integer(nt[2]),
    test_n=length(tix),test_negative=as.integer(nv[1]),
    test_positive=as.integer(nv[2]),
    active_relations=nrow(ac),
    rlr_logloss=if(length(tix))mean(lr[tix])else NA_real_,
    l2_logloss=if(length(tix))mean(ll[tix])else NA_real_,
    delta_logloss=if(length(tix))mean(ll[tix]-lr[tix])else NA_real_,
    rlr_brier=if(length(tix))mean(br[tix])else NA_real_,
    l2_brier=if(length(tix))mean(bl[tix])else NA_real_,
    delta_brier=if(length(tix))mean(bl[tix]-br[tix])else NA_real_,
    worse_than_rlr=sum(ll[tix]>lr[tix]+1e-12),
    better_than_rlr=sum(ll[tix]<lr[tix]-1e-12),
    max_sample_delta=if(length(tix))max(ll[tix]-lr[tix])else NA_real_,
    stringsAsFactors=FALSE)
}
by_leaf <- do.call(rbind,result)
require(sum(by_leaf$test_n)==length(y01) &&
        abs(weighted.mean(by_leaf$l2_logloss,by_leaf$test_n)-mean(ll))<1e-12,
        "weighted leaf logloss must reproduce total logloss")
summary <- c(
  "RLMT accepted Singh split: test-only diagnostic (not a new experiment)",
  paste("Dataset:",dataset,"fold:",fold_id),
  paste("Class order negative/positive:",negative,"/",positive),
  paste("Split:",a,">",b),
  sprintf("OOF parent loss %.8f; OOF child loss %.8f; pooled gain %.8f",
          sel$parent_oof_logloss,sel$split_oof_logloss,
          sel$parent_oof_logloss-sel$split_oof_logloss),
  sprintf("CV fold mean gain %.8f; SE %.8f; 1-SE lower bound %.8f",
          sel$mean_gain,sel$se_gain,sel$lower_gain),
  sprintf("Outer-test RLR logloss %.8f; L2 logloss %.8f; delta %.8f",
          mean(lr),mean(ll),mean(ll)-mean(lr)),
  sprintf("Outer-test RLR Brier %.8f; L2 Brier %.8f; delta %.8f",
          mean(br),mean(bl),mean(bl)-mean(br)),
  sprintf("Worse on %d/%d outer-test samples",
          sum(ll>lr+1e-12),length(y01)),
  "No retraining. No modification of the frozen experiments."
)
cat(paste(summary,collapse="\n"),"\n")
cat("\nPER-LEAF HOLDOUT DECOMPOSITION\n")
print(by_leaf,row.names=FALSE,digits=7)
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
write.csv(by_leaf,file.path(dest,"leaf_decomposition.csv"),row.names=FALSE)
writeLines(c(summary,"","Interpretation caution:",
             "This is one accepted split on one outer fold.",
             "The OOF candidate selection reuses the same outer-training data;",
             "it is not an independent nested estimate of selection generalization.",
             "Outer-test labels are used only in this post-hoc diagnostic."),
           file.path(dest,"audit.txt"))
cat("ACCEPTED_SPLIT_DIAGNOSTIC=PASS\n")
cat("REPORT=",dest,"\n",sep="")
