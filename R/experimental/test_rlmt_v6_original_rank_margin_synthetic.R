# RLMT-v6 REAL synthetic switchBox/glmnet + exact manuscript RLR
# parity (not merely dummy algebra). Original paper reference untouched.
for(pkg in c("switchBox","glmnet","jsonlite","pROC"))
 if(!requireNamespace(pkg,quietly=TRUE))
   stop("RLMT6_SYNTH_STOP: missing ",pkg)
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v4_rank_repair.R")
source("R/experimental/rlmt_v6_original_rank_margin.R")
check <- function(ok,msg) {
 if(!isTRUE(ok))stop("RLMT6_SYNTH_STOP: ",msg,call.=FALSE)
 message("PASS ",msg)
}
example <- function(n,seed,label) {
 set.seed(seed)
 p <- 45L
 x <- matrix(rnorm(n*p),ncol=p)
 colnames(x) <- sprintf("gene_%02d",seq_len(p))
 eta <- .7*(x[,1L]>x[,2L])-.9*(x[,3L]>x[,4L])+
     .9*(x[,5L]-x[,6L])+.3*(x[,7L]>x[,8L])
 y <- factor(ifelse(runif(n)<plogis(eta),"P","N"),
    levels=c("N","P"))
 X <- as.data.frame(x,check.names=FALSE)
 rownames(X) <- sprintf("%s_%04d",label,seq_len(n))
 list(X=X,y=y)
}
train <- example(175L,1454L,"train")
test <- example(52L,1455L,"test")
check(min(table(train$y))>=40L,"two classes in real synthetic")
fit <- rlmt6_fit(train$X,train$y)
prediction <- rlmt6_predict(fit,test$X)
check(identical(names(prediction),rlmt6_modes) &&
  all(vapply(prediction,function(z)
     length(z)==nrow(test$X) && all(z>0 & z<1),logical(1))),
  "all five actual RLR-anchored model predictions")
check(identical(fit$parent$model_identity,
  "FROZEN_ORIGINAL_RLR_K10_M10_D2"),"exact original root identity")
check(!fit$root_only &&
  length(fit$oof_partition)==3L &&
  length(fit$oof_training_parent_probs)==nrow(train$X),
  "strict original parent 3-fold OOF")
check(all(vapply(fit$oof_partition,function(z)
  length(intersect(z$train,z$heldout))==0L,logical(1))),
  "OOF parent train/heldout disjoint")
check(length(fit$feature_screen_scores)<=10 &&
  fit$selected_original_pair>=1L &&
  fit$selected_original_pair<=10L,
  "select correction only among <=10 original RLR pairs")
check(all(vapply(fit$coefficients,function(z)
  all(is.finite(c(z$beta,z$gamma))),logical(1))),
  "all low-dimensional gamma coefficients finite")
check(identical(names(train$X),fit$transform$training_genes) &&
  all(abs(fit$transform$conditional_means)<=1),
  "rank distance transformer trained without heldout labels")
cal <- fit$coefficients$GLOBAL_CAL
testp <- stats::plogis(rlmt6_logit(prediction$ORIGINAL_RLR)+
                    pmax(-2,pmin(2,cal$beta)))
check(isTRUE(all.equal(prediction$GLOBAL_CAL,
                       as.numeric(testp),tolerance=1e-14)),
  "explicit original-RLR global calibration control parity")

# Verify independent *original paper function* on real synthetic input
pool <- compute_top_pairs_switchbox(train$X,train$y,
  top_n=1000L,buffer_mult=5L,out_dir=NULL,save_all=FALSE)
original_dir <- tempfile("rlmt6_original_real_")
fit_rlr_train_test(train$X,train$y,test$X,test$y,
  pairs_pool=pool,K=10L,m=10L,d=2L,
  out_dir=original_dir,
  alpha=.5,lambda_choice="lambda.min",seed=42L)
reference <- read.csv(file.path(original_dir,"predictions_test.csv"),
                      stringsAsFactors=FALSE)
check(identical(as.character(reference$sample_id),
                rownames(test$X)) &&
    isTRUE(all.equal(as.numeric(reference$prob),
                     prediction$ORIGINAL_RLR,tolerance=1e-12)),
    "full original manuscript model output probabilities")
unlink(original_dir,recursive=TRUE)

# Relational invariance under strictly monotonic affine transform of
# all within-patient features: pair encodings and rank differences same.
transformed <- as.data.frame(
 lapply(test$X,function(x)2*x+10),check.names=FALSE)
names(transformed) <- names(test$X)
rownames(transformed) <- rownames(test$X)
same <- rlmt6_predict(fit,transformed)
check(all(vapply(rlmt6_modes,function(m)
  isTRUE(all.equal(prediction[[m]],same[[m]],tolerance=1e-12)),
  logical(1))),"positive affine within-sample order invariance")

save <- tempfile(fileext=".rds")
saveRDS(fit,save)
resaved <- rlmt6_predict(readRDS(save),test$X)
unlink(save)
check(identical(prediction,resaved),"bit-identical full RDS replay")
message("RLMT6_ORIGINAL_RLR_REAL_SYNTHETIC_PARITY=PASS")
message("BIOMEDICAL_DATA_USED=FALSE")
