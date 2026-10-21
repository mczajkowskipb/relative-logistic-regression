# End-to-end REAL algorithm / FORCE-ACCEPTED BRANCH on GENERATED DATA.
#
# A real SWAP, CV-glmnet, RLR parent, two fitted RLR leaf models and
# all final saved model RDS are exercised, with independent replay.
# ONLY the final Boolean confirmation decision is deliberately forced
# TRUE in this synthetic testing harness; this is NOT evidence that
# RLMT-CA-PS achieves a gain or would accept the split scientifically.
# Synthetic holdout labels are not read by ANY training-side callback.
for (pkg in c("glmnet","pROC","jsonlite","switchBox")) {
  if (!requireNamespace(pkg,quietly=TRUE))
    stop("RLMT_CA_PS_FORCED_STOP: missing ",pkg)
}
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_ca_ps_v1.R")

make_data <- function(n,prefix,seed) {
  set.seed(seed)
  p <- 42L
  mat <- matrix(rnorm(n*p),nrow=n,ncol=p)
  colnames(mat) <- sprintf("syn_%03d",seq_len(p))
  route <- mat[,1L]>mat[,2L]
  eta <- -.05 + 1.05*(mat[,3L]-mat[,4L]) +
    ifelse(route,1,-1)*1.15*(mat[,5L]-mat[,6L]) +
    ifelse(route,.65,-.65)
  y <- factor(ifelse(runif(n)<plogis(eta),"P","N"),levels=c("N","P"))
  X <- as.data.frame(mat,check.names=FALSE)
  rownames(X) <- sprintf("%s_%03d",prefix,seq_len(n))
  list(X=X,y=y)
}
training <- make_data(240L,"train",4107L)
heldout <- make_data(72L,"test",4108L)
stopifnot(all(table(training$y)>=60L),all(table(heldout$y)>=10L))
pool <- compute_top_pairs_switchbox(
  training$X,training$y,top_n=1000L,buffer_mult=5L,
  out_dir=NULL,save_all=FALSE)
stopifnot(nrow(pool)>=10L)

# The one intentional test double. Real split nomination, routing,
# D-trained and full-T re-fitted parent/leaf models, and probability
# shrinkage remain completely unmodified.
original_confirmation <- rlmt_dc_confirm
tracker <- new.env(parent=emptyenv())
tracker$calls <- 0L
rlmt_dc_confirm <- function(y_confirm,parent_prob,split_prob,min_n=20L) {
  r <- original_confirmation(y_confirm,parent_prob,split_prob,min_n)
  stopifnot(length(y_confirm)>=min_n,
            is.finite(r$gain),
            is.finite(r$se_gain),
            is.finite(r$lower_gain),
            identical(r$n,length(y_confirm)))
  tracker$calls <- tracker$calls+1L
  tracker$original <- r
  r$accepted <- TRUE
  r
}
out <- tempfile(pattern="rlmt_ca_ps_real_forced_")
run <- rlmt_ca_ps_v1_run(
  training$X,training$y,heldout$X,heldout$y,
  outer_training_pool=pool,out_dir=out)
stopifnot(tracker$calls==1L,
          isTRUE(run$accepted_split),run$realized_leaves==2L,
          run$decision$status=="ACCEPTED",
          all(run$prob>=0&run$prob<=1))
root <- file.path(out,"rlr_reference")
cmp <- read.csv(file.path(out,"prediction_components_test.csv"),
                stringsAsFactors=FALSE)
prediction <- read.csv(file.path(out,"predictions_test.csv"),
                       stringsAsFactors=FALSE)
stopifnot(nrow(cmp)==nrow(heldout$X),
          identical(as.character(cmp$sample_id),rownames(heldout$X)),
          isTRUE(all.equal(cmp$mixed_prob,run$prob,tolerance=1e-12)),
          isTRUE(all.equal(prediction$prob,run$prob,tolerance=1e-12)),
          all(cmp$shrink_weight>0 & cmp$shrink_weight<1))

# INDEPENDENT replay of original RLR model from its saved RDS/pair list:
root_fit <- readRDS(file.path(root,"cvfit.rds"))
root_pairs <- readRDS(file.path(root,"selected_pair_shortlist.rds"))
root_X <- encode_pair_matrix(heldout$X,root_pairs,m=10L)
root_prob <- as.numeric(predict(root_fit,newx=root_X,
                                s="lambda.min",type="response"))
stopifnot(isTRUE(all.equal(root_prob,cmp$parent_prob,tolerance=1e-12)))

# INDEPENDENT replay of BOTH complete saved fitted leaf models.
split <- run$decision$nominated_pair
z <- heldout$X[[split$gene_a]]>heldout$X[[split$gene_b]]
z_train <- training$X[[split$gene_a]]>training$X[[split$gene_b]]
support <- rlmt_ca_ps_support(z_train,training$y,10L)
replayed_child <- numeric(nrow(heldout$X))
for (i in 1:2) {
  ix <- which(if (i==1L) z else !z)
  model <- readRDS(file.path(out,"leaf_models",
                             sprintf("node_%s",i+1L),
                             "model_complete.rds"))
  stopifnot(is.list(model),inherits(model$fit,"cv.glmnet"))
  if (length(ix)) {
    replayed_child[ix] <- rlmt_predict_leaf(
      model,heldout$X[ix,,drop=FALSE],rlmt_defaults())$prob
  }
}
replayed <- rlmt_ca_ps_probs(root_prob,replayed_child,z,
                             support$weight_true,support$weight_false)
stopifnot(isTRUE(all.equal(replayed_child,cmp$child_prob,
                           tolerance=1e-12)),
          isTRUE(all.equal(replayed,run$prob,tolerance=1e-12)),
          isTRUE(all.equal(
            cmp$shrink_weight,
            ifelse(z,support$weight_true,support$weight_false),
            tolerance=1e-12)))

# Change ONLY the synthetic TEST labels. Under identical initial random
# seed and training data the entire fitted estimator must be invariant.
flipped <- factor(rev(as.character(heldout$y)),
                  levels=levels(heldout$y))
run2 <- rlmt_ca_ps_v1_run(
  training$X,training$y,heldout$X,flipped,
  outer_training_pool=pool,
  out_dir=paste0(out,"_test_labels_permuted"))
stopifnot(isTRUE(run2$accepted_split),
          isTRUE(all.equal(run$prob,run2$prob,tolerance=1e-12)),
          identical(run$decision$nominated_pair,
                    run2$decision$nominated_pair),
          tracker$calls==2L)

cat("RLMT_CA_PS_REAL_FORCED_ACCEPTED_BRANCH=PASS\n")
cat("FORCED_CONFIRMATION_BOOL_ONLY=TRUE\n")
cat("ORIGINAL_CONFIRMATION_GAIN=",
    format(tracker$original$gain,digits=8),"\n",sep="")
cat("REAL_PARENT_AND_TWO_LEAF_RDS_REPLAY=PASS\n")
cat("HOLDOUT_LABEL_PERMUTATION_INVARIANCE=PASS\n")
cat("NO_BIOMEDICAL_DATA=TRUE\n")
unlink(out,recursive=TRUE,force=TRUE)
unlink(paste0(out,"_test_labels_permuted"),recursive=TRUE,force=TRUE)
