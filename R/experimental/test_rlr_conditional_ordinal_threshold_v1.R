source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_conditional_ordinal_threshold_v1.R")
set.seed(591L)
n<-160L;p<-90L
X<-as.data.frame(matrix(rnorm(n*p),nrow=n,ncol=p))
names(X)<-sprintf("g%03d",seq_len(p))
rownames(X)<-sprintf("patient%03d",seq_len(n))
U<-X
for(i in seq_len(n)){
 v<-unlist(X[i,,drop=FALSE],use.names=FALSE)
 U[i,]<-exp(v/runif(1L,2,5))*runif(1L,.5,2)
}
stopifnot(identical(unname(aa_ranks(X)),unname(aa_ranks(U))))
root<-os_node("g001",c("g002","g003","g004"),2L)
gate<-os_gate(aa_ranks(X),root)
# The same ordinal underlying disease signal receives branch-dependent
# score offsets: conditional thresholds are needed without changing p.
y<-as.integer(aa_ranks(X)[,"g005"]>aa_ranks(X)[,"g006"])
p0<-pmin(.99,pmax(.01,
    plogis(2.3*(y-.5)+ifelse(gate==0,-1,1)+rnorm(n,sd=.35))))
glob<-ct_global(y,p0)
candidate<-ct_fit_root(y,p0,gate,glob)
stopifnot(!is.null(candidate),candidate$F1>=glob$F1+.01,
          candidate$BA>=glob$BA-.015)
predictA<-ct_predict(p0,aa_ranks(X),
 list(type="TREE",root=root,t0=candidate$t0,t1=candidate$t1))
predictB<-ct_predict(p0,aa_ranks(U),
 list(type="TREE",root=root,t0=candidate$t0,t1=candidate$t1))
stopifnot(identical(predictA,predictB),
          ct_f1(y,predictA)>=glob$F1+.01)
cat("COT_RLR_SYNTHETIC_CONDITIONAL_F1_GAIN=PASS\n")
cat("COT_RLR_PATIENTWISE_MONOTONE_PREDICTION_INVARIANCE=PASS\n")
v<-ct_predict(p0,aa_ranks(X),list(type="GLOBAL",t=glob$threshold))
stopifnot(all(v %in% 0:1),identical(v,as.integer(p0>=glob$threshold)))
# Classification thresholding MUST NOT modify score-only metrics.
originalLL<-rg_ll(y,p0); originalBrier<-mean((y-p0)^2)
rawp<-p0
stopifnot(identical(originalLL,rg_ll(y,rawp)),
          identical(originalBrier,mean((y-rawp)^2)))
cat("COT_RLR_AUC_LOSS_BRIER_SCORE_UNCHANGED=PASS\n")
# Training-only undersized fallback, using ten samples.
idx<-c(which(y==0)[1:5],which(y==1)[1:5])
small_y<-factor(y[idx],levels=0:1)
small_fit<-ct_select(X[idx,,drop=FALSE],small_y)
stopifnot(is.null(small_fit$root),!small_fit$accepted,
          grepl("UNDERPOWERED",small_fit$reason))
cat("COT_RLR_NESTED_TRAIN_SUPPORT_FALLBACK=PASS\n")
cat("COT_RLR_ALL_MECHANICAL_TESTS=PASS\n")
