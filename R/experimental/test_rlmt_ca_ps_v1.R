# Synthetic-only unit/integration test of RLMT-CA-PS-v1.
# Run from a separate worktree:
#   Rscript --vanilla R/experimental/test_rlmt_ca_ps_v1.R
#
# IMPORTANT: no biomedical ARFF, no original model fit. Training callbacks
# are deterministic TEST DOUBLES; verifies workflow invariants, not real
# biomedical efficacy. Production deps are exercised in a separate smoke.
stopifnot(requireNamespace("jsonlite",quietly=TRUE))
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/experimental/rlmt_ca_ps_v1.R")

assert <- function(ok,message) {
  if (!isTRUE(ok)) stop("RLMT_CA_PS_TEST_FAIL: ",message,call.=FALSE)
  cat("PASS ",message,"\n",sep="")
}

# Pure numeric properties of the class-support penalty and convex bridge.
assert(identical(rlmt_ca_ps_node_balance(10L,10L),1),
       "class-balanced_10_10_has_index_1")
assert(abs(rlmt_ca_ps_node_balance(9L,1L)-0.36)<1e-14,
       "class_skew_is_penalized")
assert(rlmt_ca_ps_node_balance(0L,0L)==0 &&
       rlmt_ca_ps_node_balance(10L,0L)==0,
       "empty_or_pure_branch_is_not_fatal")
assert(rlmt_ca_ps_node_weight(10L,10L,10L) >
       rlmt_ca_ps_node_weight(19L,1L,10L),
       "fixed_size_balance_prefers_reliable_support")
assert(rlmt_ca_ps_node_weight(50L,50L,10L)>
       rlmt_ca_ps_node_weight(10L,10L,10L),
       "larger_balanced_leaf_gets_more_local_weight")
assert(rlmt_ca_ps_node_weight(6L,6L,10L)<1 &&
       rlmt_ca_ps_node_weight(30L,30L,10L)<1,
       "no_unrestricted_local_extrapolation")
p0 <- c(0.01,0.5,0.99,0)
p1 <- c(0.99,0.9,0.01,1)
route <- c(TRUE,FALSE,TRUE,FALSE)
mixed <- rlmt_ca_ps_probs(p0,p1,route,0.2,0.7)
assert(isTRUE(all.equal(mixed,
       p0+ifelse(route,0.2,0.7)*(p1-p0),tolerance=0)),
       "convex_shrink_rule_exact")
assert(all(mixed>=pmin(p0,p1) & mixed<=pmax(p0,p1)),
       "mixture_never_extrapolates_probability")
assert(identical(rlmt_ca_ps_probs(p0,p1,route,0,0),p0),
       "zero_shrink_is_exact_parent")
assert(!"y_test" %in% names(formals(rlmt_ca_ps_propose)) &&
       !"y_test" %in% names(formals(rlmt_ca_ps_probs)),
       "split_screen_and_mixture_no_outer_test_labels")
error <- tryCatch(rlmt_ca_ps_probs(p0,p1,c(TRUE,NA,TRUE,FALSE),.2,.2),
                  error=function(e)e)
assert(inherits(error,"error"),"NA_test_routing_fails_closed")

# Deterministic held-out-free candidate screen. Both children 20/20.
n <- 120L
ys <- factor(ifelse(seq_len(n)%%2L==1L,"P","N"),
             levels=c("N","P"))
Xs <- data.frame(
  g1=ifelse(seq_len(n)%%4L %in% c(0L,1L),2,0),
  g2=rep(1,n),
  signal=ifelse(ys=="P",2,0),
  g4=seq_len(n)/n)
rownames(Xs) <- sprintf("syn_%03d",seq_len(n))
one <- data.frame(pool_index=1L,gene_a="g1",gene_b="g2")
conf <- list(min_child_n=16L,min_child_class=6L,K=10L)
cs <- rlmt_ca_ps_propose(Xs,ys,rep(.5,n),one,conf)
assert(nrow(cs$audit)==1L &&
       isTRUE(cs$audit$feasible[[1L]]) &&
       cs$audit$class_balance_left[[1L]]==1 &&
       cs$audit$class_balance_right[[1L]]==1,
       "balanced_candidate_not_rejected")
assert(identical(cs$candidate$gene_a[[1L]],"g1"),
       "deterministic_candidate_nomination")
support <- rlmt_ca_ps_support(Xs$g1>Xs$g2,ys,K=10L)
assert(identical(support$counts_true,c(30L,30L)) &&
       identical(support$counts_false,c(30L,30L)),
       "class_counts_in_both_children")

# Use actual frozen confirm mathematics, not a test double.
zV <- ys[81:120]
confirmed <- rlmt_dc_confirm(
  zV,rep(.5,40L),
  rlmt_ca_ps_probs(
    rep(.5,40L),ifelse(zV=="P",.95,.05),
    Xs$g1[81:120]>Xs$g2[81:120],
    .6,.6),min_n=20L)
assert(isTRUE(confirmed$accepted) &&
       confirmed$gain>0 && confirmed$lower_gain>0,
       "original_confirm_rule_accepts_true_training_only_gain")
failed <- rlmt_dc_confirm(
  zV,rep(.5,40L),rep(.5,40L),min_n=20L)
assert(!isTRUE(failed$accepted),
       "zero_confirmation_gain_never_accepted")

# End-to-end workflow tests use explicit synthetic training callbacks.
# These doubles forbid any inference from biological datasets and exercise
# exactly the file/status/fallback interface of a real experiment.
# 07_world_protocol_utils.R is intentionally NOT sourced here because it
# loads real dependencies. Stub its pure filesystem helper and null coalescer,
# both of which are required by the authentic DC defaults/run functions.
`%||%` <- function(a,b) if (is.null(a)) b else a
ensure_dir <- function(path)
  dir.create(path,recursive=TRUE,showWarnings=FALSE)
assert(is.function(ensure_dir) && is.function(`%||%`),
       "synthetic_07_helpers_explicitly_mocked")
rlmt_defaults <- function(config=list()) list(
  K=10L,m=10L,d=2L,pool=1000L,alpha=.5,lambda_choice="lambda.min",
  min_child_n=16L,min_child_class=6L,inner_folds=3L,seed=42L)
rlmt_dc_partition <- function(y,seed=42L) {
  assert(length(y)==120L,"synthetic_training_size_contract")
  list(discovery=seq_len(80L),confirmation=81:120)
}
make_repeated_stratified_folds <- function(y,n_folds=3L,
                                           n_repeats=1L,seed=42L) {
  list(list(train_idx=seq_len(length(y)-10L),
            test_idx=(length(y)-9L):length(y)))
}
compute_top_pairs_switchbox <- function(X_train,y_train,top_n,
                                        buffer_mult,out_dir,save_all) {
  data.frame(gene_left=paste0("g",seq(4L,22L,by=2L)),
             gene_right=paste0("g",seq(5L,23L,by=2L)),
             ordered_pair=paste0("g",seq(4L,22L,by=2L),">g",
                                 seq(5L,23L,by=2L)),
             stringsAsFactors=FALSE)
}
rlmt_dc_parent_oof_strict <- function(X,y,conf,folds)
  rep(.5,nrow(X))
rlmt_dc_random_pairs <- function(feature_names,n_pairs,seed)
  one
rlmt_dc_parent_and_child_on_validation <- function(
  X_discovery,y_discovery,X_confirmation,pair,discovery_pool,conf) {
  list(parent_prob=rep(.5,nrow(X_confirmation)),
       split_prob=ifelse(X_confirmation$signal>1,.95,.05),
       nominated_pair=pair[,c("gene_a","gene_b"),drop=FALSE])
}
fit_rlr_train_test <- function(X_train,y_train,X_test,y_test,pairs_pool,
                               K,m,d,out_dir,alpha,lambda_choice,
                               seed,extra_metrics) {
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  p <- rep(.5,nrow(X_test))
  write.csv(data.frame(sample_id=rownames(X_test),
                       y_true=as.character(y_test),
                       y_pred=rep(levels(y_test)[2L],length(y_test)),
                       prob=p,score=p),
            file.path(out_dir,"predictions_test.csv"),
            row.names=FALSE)
  saveRDS(list(stub=TRUE),file.path(out_dir,"cvfit.rds"))
  invisible(TRUE)
}
rlmt_fit_leaf <- function(X,y,pool,conf,rerank=TRUE,seed=42L) {
  list(fit=list(stub=TRUE),pairs=pool,
       active=pool[1,,drop=FALSE],beta=1,
       stub_signal=TRUE)
}
rlmt_predict_leaf <- function(model,X,conf)
  list(prob=ifelse(X$signal>1,.95,.05),
       local_active_n=rep(1,nrow(X)))
rlmt_save_leaf <- function(model,path) {
  dir.create(path,recursive=TRUE,showWarnings=FALSE)
  write.csv(data.frame(gene_left="g4",gene_right="g5"),
            file.path(path,"active_pairs.csv"),row.names=FALSE)
  write.csv(data.frame(ordered_pair="g4>g5",beta=1,abs_beta=1),
            file.path(path,"coefficients.csv"),row.names=FALSE)
}
rlmt_dc_write_json <- function(value,path)
  jsonlite::write_json(value,path,auto_unbox=TRUE,
                       pretty=TRUE,digits=12,null="null")
apply_disjointness <- function(pairs_df,K,d)
  pairs_df[seq_len(min(K,nrow(pairs_df))),,drop=FALSE]

# Synthetic feature schema with independent ROUTE group and SIGNAL class.
makeX <- function(indices,y) {
  out <- data.frame(
    g1=ifelse(indices%%4L %in% c(0L,1L),2,0),
    g2=rep(1,length(indices)),
    signal=ifelse(y=="P",2,0))
  for (j in 4L:26L) out[[paste0("g",j)]] <-
    sin(indices*(0.07+j/100))+j/100
  rownames(out) <- sprintf("synth_%03d",indices)
  out
}
X_train <- makeX(1:120,ys)
y_train <- ys
y_test <- factor(rep(c("N","P"),length.out=20L),
                 levels=levels(ys))
X_test <- makeX(201:220,y_test)
pool <- compute_top_pairs_switchbox(X_train,y_train,1000L,5L,NULL,FALSE)
base <- tempfile(pattern="rlmt_ca_ps_synthetic_")
dir.create(base)
good <- rlmt_ca_ps_v1_run(
  X_train,y_train,X_test,y_test,pool,
  out_dir=file.path(base,"accepted"))
assert(isTRUE(good$accepted_split) &&
       good$realized_leaves==2L,
       "synthetic_independent_gate_accepts_new_split")
assert(all(good$prob>0 & good$prob<1) &&
       any(abs(good$prob-.5)>1e-3),
       "accepted_leaf_shrink_changes_prob_but_is_conservative")
component <- read.csv(file.path(base,"accepted",
                                "prediction_components_test.csv"))
assert(isTRUE(all.equal(component$mixed_prob,good$prob,
                        tolerance=1e-12)) &&
       all(component$shrink_weight<1),
       "per_prediction_shrink_weights_export_and_replay")
assert(file.exists(file.path(base,"accepted","leaf_models",
                             "node_2","model_complete.rds")) &&
       file.exists(file.path(base,"accepted","leaf_models",
                             "node_3","model_complete.rds")),
       "complete_local_models_are_persisted")

# Out-of-sample TRUE labels may change evaluation, NEVER fitted
# split/routing/weights/conditioned probabilities.
permuted <- factor(rev(as.character(y_test)),levels=levels(y_test))
changed_y <- rlmt_ca_ps_v1_run(
  X_train,y_train,X_test,permuted,pool,
  out_dir=file.path(base,"changed_holdout_labels"))
assert(identical(good$prob,changed_y$prob) &&
       identical(good$decision$confirmation,
                 changed_y$decision$confirmation),
       "test_label_permutation_does_not_change_predictions_or_gate")

# No split: accepted=FALSE, original parent prediction reproduced BIT EXACT.
smallx <- X_train[1:36,,drop=FALSE]
smally <- y_train[1:36]
root <- rlmt_ca_ps_v1_run(
  smallx,smally,X_test,y_test,pool,
  out_dir=file.path(base,"no_split"))
parity <- read.csv(file.path(base,"no_split","rlr_reference",
                             "predictions_test.csv"))
assert(!root$accepted_split && root$realized_leaves==1L &&
       identical(root$prob,as.numeric(parity$prob)) &&
       all(root$prob==.5),
       "no_split_is_bit_exact_original_RLR_parent")
unlink(base,recursive=TRUE,force=TRUE)
cat("RLMT_CA_PS_SYNTHETIC_UNIT_AND_MOCK_INTEGRATION=PASS\n")
cat("NO_OMICS_DATA_TOUCHED=TRUE\n")
