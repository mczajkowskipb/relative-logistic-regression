# Purely SYNTHETIC regression checks for 2011-inspired tree + BigTSP adapter.
# Run: Rscript --vanilla R/baselines/test_tspdt_reconstruction_v1.R
source("R/baselines/tspdt_reconstruction_v1.R")
source("R/baselines/bigtsp_pair_cart_v1.R")

assert <- function(ok,label) {
  if (!isTRUE(ok))
    stop("TSPDT_RECON_TEST_FAIL: ",label,call.=FALSE)
  cat("PASS ",label,"\n",sep="")
}

# A class-informative ROOT relation, followed by two different
# node-specific corrections; tests recursive node-local pair search.
n <- 120L
g <- rep(c(TRUE,FALSE),each=60L)
labels <- rep("N",n)
labels[1:54] <- "PR"
labels[61:66] <- "PR"
y <- factor(labels,levels=c("N","PR"))
X <- cbind(
  a=ifelse(g,0.3,1.7),
  b=rep(1,n),
  c=ifelse(seq_len(n) %in% 55:60,0.8,1.2),
  d=rep(1,n),
  e=ifelse(seq_len(n) %in% 61:66,0.8,1.2),
  f=rep(1,n))
rownames(X) <- sprintf("synthetic_%03d",seq_len(n))

fit <- tspdt_recon_fit(X,y,top_n=6L,min_leaf=3L,
                        max_depth=3L,pruning="none")
fit2 <- tspdt_recon_fit(X,y,top_n=6L,min_leaf=3L,
                         max_depth=3L,pruning="none")
pred <- tspdt_recon_predict(fit,X)
assert(identical(fit,fit2),"deterministic_tree_fit")
assert(fit$method=="TSPDT2011_INSPIRED_RECONSTRUCTION_V1" &&
       fit$source_claim=="not_original_TSPDT_2011",
       "explicit_reconstruction_label")
assert(identical(fit$tree$split$gene_a,"a") &&
       identical(fit$tree$split$gene_b,"b"),
       "root_matches_classic_conditional_TSP")
assert(fit$internal_nodes>=3L && fit$leaves>=4L,
       "node_local_recursive_growth")
assert(!is.null(fit$tree$left$split) &&
       !is.null(fit$tree$right$split) &&
       paste(fit$tree$left$split$gene_a,
             fit$tree$left$split$gene_b)!="a b" &&
       paste(fit$tree$right$split$gene_a,
             fit$tree$right$split$gene_b)!="a b",
       "child_TSP_tests_selected_afresh")
assert(identical(as.character(pred$class),as.character(y)),
       "synthetic_interpretability_rules_recover_labels")
assert(all(is.finite(pred$prob_class1)) &&
       all(pred$prob_class1>0 & pred$prob_class1<1) &&
       all(!is.na(pred$leaf_id)),
       "train_only_Laplace_probability_finite")

nodes <- tspdt_recon_nodes(fit)
assert(nrow(nodes)==fit$internal_nodes+fit$leaves &&
       sum(nodes$is_leaf)==fit$leaves,
       "nodes_export_consistency")

# Saved objects route identically on held-out, unlabeled X.
tmp <- tempfile(fileext=".rds")
saveRDS(fit,tmp)
replay <- tspdt_recon_predict(readRDS(tmp),X)
assert(identical(pred,replay),"serialized_tree_prediction_parity")
unlink(tmp)

# Tie convention: strict <, so equality routes right.
Xtie <- X[1:12,,drop=FALSE]
Xtie[,c("a","b")] <- 1
Xforced_right <- Xtie
Xforced_right[,"a"] <- 2
Xforced_right[,"b"] <- 1
test_prediction <- tspdt_recon_predict(fit,Xtie)
right_prediction <- tspdt_recon_predict(fit,Xforced_right)
assert(identical(test_prediction$leaf_id,right_prediction$leaf_id),
       "ties_route_false_right_at_root")

# Same training, modified test-only X must never refit the model.
set.seed(814L)
X_test <- X[1:20,,drop=FALSE]+matrix(rnorm(20L*6L,sd=0.01),20L,6L)
colnames(X_test) <- colnames(X)
pred_test1 <- tspdt_recon_predict(fit,X_test)
pred_test2 <- tspdt_recon_predict(fit,X_test)
assert(identical(pred_test1,pred_test2),
       "test_predictions_deterministic_without_test_labels")
assert(!"y_test" %in% names(formals(tspdt_recon_predict)) &&
       !"y_test" %in% names(formals(bigtsp_pair_cart_predict)),
       "prediction_api_never_consumes_test_labels")

# Class-label reversal changes the orientation of positive class,
# but cannot be used to leak outer holdout labels.
y_reversed <- factor(ifelse(y=="PR","N","PR"),levels=c("N","PR"))
fit_reversed <- tspdt_recon_fit(X,y_reversed,top_n=6L,
                                min_leaf=3L,max_depth=3L,
                                pruning="none")
assert(identical(fit$selected_genes,fit_reversed$selected_genes),
       "train_only_MAD_filter_label_independent")

# Optional training-only Wilson pruning is deterministic.
pruned <- tspdt_recon_fit(X,y,top_n=6L,min_leaf=3L,
                          max_depth=3L,pruning="wilson")
assert(pruned$internal_nodes<=fit$internal_nodes &&
       pruned$leaves<=fit$leaves,
       "Wilson_postpruning_never_expands_tree")
assert(all(is.finite(tspdt_recon_predict(pruned,X)$prob_class1)),
       "pruned_tree_predicts")

# Bounded feature selection test. 300 input genes -> 100 selected;
# no dense global pair matrix, including when X_test differs.
set.seed(901L)
Xwide <- matrix(rnorm(30L*300L),nrow=30L)
colnames(Xwide) <- sprintf("g_%03d",seq_len(300L))
fitwide <- tspdt_recon_fit(
  Xwide,factor(rep(c("N","PR"),each=15L),levels=c("N","PR")),
  top_n=100L,min_leaf=5L,max_depth=1L,pruning="none")
assert(length(fitwide$selected_genes)==100L &&
       identical(length(fitwide$training_genes),300L),
       "train_only_100_gene_guard")

# Invalid leakage-/schema-sensitive inputs fail closed.
err <- tryCatch(
  tspdt_recon_predict(fit,X[,rev(colnames(X)),drop=FALSE]),
  error=function(e)e)
assert(inherits(err,"error"),"schema_mismatch_rejected")
err <- tryCatch(
  tspdt_recon_fit(X,y,top_n=1001L),
  error=function(e)e)
assert(inherits(err,"error"),"unbounded_gene_cap_rejected")

if (requireNamespace("BigTSP",quietly=TRUE) &&
    as.character(utils::packageVersion("BigTSP"))=="1.0" &&
    requireNamespace("tree",quietly=TRUE)) {
  bt <- bigtsp_pair_cart_fit(X,y,top_n=6L)
  bp <- bigtsp_pair_cart_predict(bt,X_test)
  assert(length(bp$class)==nrow(X_test) && !anyNA(bp$class) &&
         all(is.na(bp$prob_class1)) &&
         bt$source_claim=="not_original_TSPDT_2011",
         "optional_archived_BigTSP_pair_CART_smoke")
  cat("BIGTSP_PAIR_CART_OPTIONAL=PASS\n")
} else {
  cat("BIGTSP_PAIR_CART_OPTIONAL=SKIP_NOT_INSTALLED\n")
}
cat("TSPDT_RECON_SYNTHETIC=PASS\n")
