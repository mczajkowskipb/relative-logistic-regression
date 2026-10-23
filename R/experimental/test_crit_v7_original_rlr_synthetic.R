# Genuine ORIGINAL-RLR / switchBox + 2-leaf CRIT-v7 synthetic tests.
# Synthetic success on a strong planted interaction is a feasibility
# study, NOT a heldout biological claim.
for(pkg in c("switchBox","glmnet","pROC","jsonlite"))
 if(!requireNamespace(pkg,quietly=TRUE))
   stop("CRIT7_SYNTH_STOP: missing R package ",pkg)
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v4_rank_repair.R")
source("R/experimental/crit_v7_original_rlr_two_leaf.R")
chk<-function(q,msg) {
 if(!isTRUE(q))stop("CRIT7_SYNTH_STOP: ",msg,call.=FALSE)
 message("PASS ",msg)
}
make_xy<-function(n,seed,prefix,mode) {
 set.seed(seed)
 m<-matrix(stats::rnorm(n*10),nrow=n,ncol=10)
 colnames(m)<-paste0("gene",sprintf("%02d",seq_len(10)))
 gate<-ifelse(m[,1]>m[,2],1,-1)
 rel<-ifelse(m[,3]>m[,4],1,-1)
 eta<-switch(mode,interaction=2.7*gate*rel,
             global=1.2*rel,noise=rep(0,n))
 y<-factor(ifelse(stats::runif(n)<stats::plogis(eta),"P","N"),
            levels=c("N","P"))
 X<-as.data.frame(m,check.names=FALSE)
 rownames(X)<-paste0(prefix,seq_len(n))
 list(X=X,y=y)
}
train<-make_xy(210L,1491L,"t_","interaction")
test<-make_xy(82L,1492L,"v_","interaction")
chk(min(table(train$y))>=50L && min(table(test$y))>=15L,
    "synthetic interactions have both classes")
fit<-crit7_fit(train$X,train$y)
pred<-crit7_predict(fit,test$X)
chk(identical(names(pred$probs),crit7_modes),
    "exact four predeclared methods returned")
chk(identical(fit$original$model_identity,
    "FROZEN_ORIGINAL_RLR_K10_M10_D2"),
    "genuine original paper RLR root fitted")
chk(length(intersect(fit$partition$train,
                     fit$partition$confirmation))==0L,
    "disjoint train-only discovery and confirmation patients")
chk(!is.null(fit$confirm) &&
    fit$confirm$evaluated_candidate_count<=crit7_cfg()$shortlist,
    "D/V candidate budget respected")

# Source-paper wrapper parity, not an independent surrogate.
pool<-compute_top_pairs_switchbox(train$X,train$y,top_n=1000L,
         buffer_mult=5L,out_dir=NULL,save_all=FALSE)
old<-tempfile("crit7_original_rlr_")
fit_rlr_train_test(train$X,train$y,test$X,test$y,
     pairs_pool=pool,K=10L,m=10L,d=2L,out_dir=old,
     alpha=.5,lambda_choice="lambda.min",seed=42L)
archive<-read.csv(file.path(old,"predictions_test.csv"),
                  stringsAsFactors=FALSE)
chk(identical(as.character(archive$sample_id),
              rownames(test$X)) &&
    max(abs(archive$prob-pred$probs$ORIGINAL_RLR))<1e-9,
    "exact published RLR patient-level synthetic parity")
unlink(old,recursive=TRUE)
chk(all(vapply(pred$probs,function(x)
             length(x)==nrow(test$X)&&
             all(is.finite(x))&&all(x>=0&x<=1),logical(1))),
    "all heldout probabilities valid")
if(pred$split_accepted) {
 chk(length(unique(pred$leaf_id))==2L &&
     !is.null(fit$gate)&&!is.null(fit$pair_false)&&
     !is.null(fit$pair_true)&&length(fit$tree_coef)==4L,
     "accepted true two-leaf conditional rank models")
} else {
 chk(identical(pred$probs$CRIT_TREE,
               pred$probs$ORIGINAL_RLR),
     "rejected tree EXACT original RLR probabilities")
}
chk(!is.null(fit$interaction_coef) &&
     length(fit$interaction_coef)==4L,
     "same candidate non-tree interaction control fitted")

# Strictly increasing per-patient affine expression transformations
# preserve all within-patient comparison decisions.
X2<-as.data.frame(lapply(test$X,function(a)1.7*a+13),
                   check.names=FALSE)
names(X2)<-names(test$X)
rownames(X2)<-rownames(test$X)
other<-crit7_predict(fit,X2)
chk(all(vapply(crit7_modes,function(m)
  isTRUE(all.equal(pred$probs[[m]],other$probs[[m]],
                   tolerance=1e-12)),logical(1))),
    "affine rank preservation across root/tree/interaction")

archive_model<-tempfile(fileext=".rds")
saveRDS(fit,archive_model)
replay<-crit7_predict(readRDS(archive_model),test$X)
unlink(archive_model)
chk(identical(pred,replay),"bit-exact RDS replay")

# Positive mechanism test: on a NEW planted XOR interaction,
# a prior constant-probability model is unable to express the signal.
# D and V are disjoint; the actual tree must be accepted on V.
xorD<-make_xy(320L,33015L,"xor_d_","interaction")
xorV<-make_xy(180L,33016L,"xor_v_","interaction")
strong<-crit7_select(xorD$X,xorD$y,rep(.5,nrow(xorD$X)),
        xorV$X,xorV$y,rep(.5,nrow(xorV$X)),
        fit$original,crit7_cfg())
chk(isTRUE(strong$accepted) &&
      !is.null(strong$selected) &&
      strong$best_confirmed_lower>crit7_cfg()$confirmation_min_gain,
    "planted disjoint-gene XOR gives a confirmed REAL rank split")
chk(strong$selected$rel_false != strong$selected$gate &&
      strong$selected$rel_true != strong$selected$gate,
    "two separately chosen leaf predictors not equal to gate")

# Zero corrected-residual evidence must not force a model tree.
# Fit on artificial perfect parent OOF scores and perfect V scores;
# zero of the D residual scores then no positive heldout benefit.
toy<-make_xy(150L,1451L,"toy_","global")
part<-make_repeated_stratified_folds(toy$y,3L,1L,seed=777L)[[1L]]
dd<-part$train_idx
vv<-part$test_idx
real_parent<-rlmt4_parent(toy$X[dd,,drop=FALSE],toy$y[dd])
perfectD<-as.numeric(toy$y[dd]==levels(toy$y)[2])
perfectV<-as.numeric(toy$y[vv]==levels(toy$y)[2])
selection<-crit7_select(toy$X[dd,,drop=FALSE],toy$y[dd],
         perfectD,toy$X[vv,,drop=FALSE],toy$y[vv],
         perfectV,real_parent)
chk(!selection$accepted,
  "perfect parent predictor refuses unnecessary split")
message("CRIT7_REAL_SYNTHETIC_ORIGINAL_ROOT_AND_TREE=PASS")
message("BIOMEDICAL_DATA_USED=FALSE")
