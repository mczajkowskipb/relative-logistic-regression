# Real SWAP/glmnet RLMT-v2 H1/H2 synthetic fit and exact RDS replay.
# This test deliberately mocks ONLY confirmation$accepted to force
# the real fitted two-leaf route. A forced gate is NOT scientific proof.
# Reuses no old paper ARFF, no old patient cohorts and no test-driven tuning.
for(pkg in c("foreign","switchBox","glmnet","pROC","jsonlite")) {
  if(!requireNamespace(pkg,quietly=TRUE))
    stop("V2_TEST_STOP: missing installed package ",pkg)
}
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v2_aligned_offset.R")

verify <- function(ok,label) {
  if(!isTRUE(ok))stop("V2_REAL_SYN_TEST_FAIL: ",label)
  cat("V2_PASS ",label,"\n",sep="")
}
syn <- function(n,seed,prefix) {
  set.seed(seed)
  p <- 42L
  mat <- matrix(rnorm(n*p),nrow=n,ncol=p)
  colnames(mat) <- sprintf("synthetic_gene_%03d",seq_len(p))
  z <- mat[,1L]>mat[,2L]
  eta <- -.15+1.0*(mat[,3L]-mat[,4L])+
    ifelse(z,1,-1)*.85*(mat[,5L]-mat[,6L])+
    ifelse(z,.6,-.6)
  y <- factor(ifelse(runif(n)<plogis(eta),"P","N"),
              levels=c("N","P"))
  X <- as.data.frame(mat,check.names=FALSE)
  rownames(X) <- sprintf("%s_%03d",prefix,seq_len(n))
  list(X=X,y=y)
}
train <- syn(150L,29010L,"train")
test <- syn(48L,29011L,"holdout")
verify(all(table(train$y)>=36L) && all(table(test$y)>=5L),
       "synthetic_train_and_holdout_class_support")

conf <- rlmt_v2_conf()$leaf
small <- train$X[1:90,,drop=FALSE]
ysmall <- train$y[1:90]
aligned <- rlmt_v2_fit_parent(small,ysmall,conf)
p_direct <- rlmt_v2_predict_parent(aligned,small,conf)
verify(identical(aligned$fitting_procedure,"rlmt_v2_fit_parent") &&
       identical(aligned$identity,"RLR_ALIGNED_V2") &&
       length(aligned$foldid)==nrow(small) &&
       length(unique(aligned$foldid))==5L &&
       all(p_direct>0 & p_direct<1),
       "aligned_parent_same_explicit_stratified_five_fold_cv")
f <- tempfile(pattern="v2_aligned_parent_")
saveRDS(aligned,f)
verify(isTRUE(all.equal(
  p_direct,rlmt_v2_predict_parent(readRDS(f),small,conf),
  tolerance=0)),
  "original_aligned_parent_RDS_exact_replay")
unlink(f)
oof <- rlmt_v2_parent_oof(small,ysmall,conf)
coverage <- integer(nrow(small))
for(f in oof$folds) {
  verify(length(intersect(f$train_idx,f$test_idx))==0L,
         "OOF_fold_no_train_validation_overlap")
  coverage[f$test_idx] <- coverage[f$test_idx]+1L
}
verify(all(coverage==1L) && all(oof$prob>0&oof$prob<1),
       "training_only_parent_OOF_exact_coverage")
verify(abs(rlmt_v2_logit(.5))<1e-14 &&
       max(abs(stats::plogis(rlmt_v2_logit(c(.1,.9)))-
               c(.1,.9)))<1e-14,
       "fixed_epsilon_logit_link_consistency")
# Exactly zero offset correction reconstructs parent; no extra
# intercept is permitted in H2 by construction.
verify(all.equal(stats::plogis(
  rlmt_v2_logit(p_direct)+rep(0,length(p_direct))),
  p_direct,tolerance=1e-12)==TRUE,
  "zero_correction_is_aligned_parent")
for(mode in c("H1","H2")) {
  # Natural, never-forced split gating: correct irrespective of acceptance.
  base <- tempfile(pattern=paste0("rlmt_v2_",mode,"_natural_"))
  outcome <- rlmt_v2_run(train$X,train$y,test$X,test$y,
                         out_dir=base,mode=mode)
  raw <- read.csv(file.path(base,"predictions_test.csv"),
                  stringsAsFactors=FALSE)
  verify(length(outcome$prob)==nrow(test$X) &&
         isTRUE(all.equal(as.numeric(raw$prob),outcome$prob,
                          tolerance=1e-12)) &&
         identical(outcome$decision$parent_identity,"RLR_ALIGNED_V2"),
         paste0(mode,"_natural_real_fit_output"))
  if(!outcome$decision$accepted) {
    verify(identical(outcome$prob,outcome$parent) &&
           file.exists(file.path(base,"complete_aligned_parent.rds")),
           paste0(mode,"_no_split_exact_aligned_parent"))
  } else {
    verify(file.exists(file.path(base,"complete_fitted_tree.rds")),
           paste0(mode,"_natural_accepted_two_leaf_model"))
  }
  unlink(base,recursive=TRUE,force=TRUE)
}
original_confirm <- rlmt_dc_confirm
forced_calls <- 0L
rlmt_dc_confirm <- function(y_confirm,parent_prob,split_prob,min_n=20L) {
  actual <- original_confirm(y_confirm,parent_prob,split_prob,min_n)
  verify(is.finite(actual$gain) &&
         is.finite(actual$lower_gain) &&
         actual$n==length(y_confirm),
         "test_force_only_after_real_confirmation_computation")
  forced_calls <<- forced_calls+1L
  actual$accepted <- TRUE
  actual
}
for(mode in c("H1","H2")) {
  work <- tempfile(pattern=paste0("rlmt_v2_",mode,"_forced_"))
  a <- rlmt_v2_run(train$X,train$y,test$X,test$y,
                    out_dir=work,mode=mode)
  verify(isTRUE(a$decision$accepted) &&
         a$decision$reason=="independent_train_only_confirmation_accepted" &&
         file.exists(file.path(work,"complete_fitted_tree.rds")),
         paste0(mode,"_real_two_leaf_fit_test_forced"))
  readback <- readRDS(file.path(work,"complete_fitted_tree.rds"))
  p_replayed <- rlmt_v2_predict_tree(readback,test$X,conf)
  verify(isTRUE(all.equal(a$prob,p_replayed$prob,tolerance=1e-12)) &&
         isTRUE(all.equal(a$parent,p_replayed$parent,tolerance=1e-12)) &&
         all(p_replayed$route %in% c(2L,3L)),
         paste0(mode,"_parent_and_two_leaves_full_RDS_replay"))
  if(mode=="H2") {
    verify(identical(readback$leaves[[1L]]$intercept,FALSE) &&
           identical(readback$leaves[[2L]]$lambda,1) &&
           identical(readback$leaves[[1L]]$alpha,0) &&
           file.exists(file.path(work,"parent_oof_training.rds")),
           "H2_only_fixed_ridge_no_intercept_and_saved_OOF")
    offsets <- readRDS(file.path(work,"parent_oof_training.rds"))
    verify(all(offsets$prob>0&offsets$prob<1) &&
           length(offsets$prob)==nrow(train$X) &&
           length(offsets$fold_ids)==3L,
           "H2_full_train_OOF_offsets_auditable")
  }
  # Outer-test labels are not part of the learner. Permute them,
  # rerun the entire training-only fit/nomination/confirmation.
  permuted <- factor(rev(as.character(test$y)),
                     levels=levels(test$y))
  work2 <- tempfile(pattern=paste0("rlmt_v2_",mode,"_perturbed_"))
  b <- rlmt_v2_run(train$X,train$y,test$X,permuted,
                   out_dir=work2,mode=mode)
  verify(identical(a$decision$nominated_pair,
                    b$decision$nominated_pair) &&
         identical(a$decision$confirmation,
                    b$decision$confirmation) &&
         isTRUE(all.equal(a$prob,b$prob,tolerance=0)),
         paste0(mode,"_outer_label_permutation_no_fit_effect"))
  unlink(work,recursive=TRUE,force=TRUE)
  unlink(work2,recursive=TRUE,force=TRUE)
}
verify(forced_calls==4L,"exactly_four_forced_test_confirmation_calls")
cat("RLMT_V2_REAL_SYNTHETIC_H1_H2_GATE=PASS\n")
cat("H1_FIT_PROCEDURE_ALIGNED=TRUE\n")
cat("H2_FIXED_LAMBDA_OFFSET_OOF=TRUE\n")
cat("TEST_ONLY_FORCED_CONFIRMATION_DISCLOSED=TRUE\n")
cat("NO_BIOMEDICAL_DATA=TRUE\n")
