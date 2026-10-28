source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_weighted_pairpair_tree_core_v4.R")
set.seed(20261028L)
n<-110L;p<-72L
X<-as.data.frame(matrix(rnorm(n*p),nrow=n,ncol=p))
names(X)<-sprintf("g%03d",seq_len(p))
rownames(X)<-sprintf("patient%03d",seq_len(n))
R<-rg_ranks(X)
y<-factor(as.integer(
  (R[,"g002"]-R[,"g001"]>10) +
  (R[,"g004"]-R[,"g003"]>10)>0),levels=0:1)
stopifnot(min(table(y))>15L)
pairs<-data.frame(gene_left=sprintf("g%03d",seq(1,15,2)),
                  gene_right=sprintf("g%03d",seq(2,16,2)),
                  stringsAsFactors=FALSE)
pairs$ordered_pair<-paste(pairs$gene_left,">",pairs$gene_right)
z<-encode_pair_matrix(X,pairs,m=10L)
cv<-glmnet::cv.glmnet(z,as.integer(y==levels(y)[2L]),
     family="binomial",alpha=.5,nfolds=5L,type.measure="deviance")
parent<-list(selected=pairs,fit=cv)
teacher<-pmin(.99,pmax(.01,.20+.60*as.numeric(y)))
for(mode in c("BETA","SIGN")){
 model<-wt_fit(X,y,parent,teacher,coeff_mode=mode,alpha=.25,
               leaf_alpha=0,max_depth=2L,max_internal=3L,
               min_leaf=5L)
 result<-wt_predict(model,X)
 stopifnot(model$n_internal<=3L,model$n_unique_pairs<=6L,
           model$n_unique_genes<=12L,all(result$n_evaluated<=4L),
           all(result$n_evaluated%%2L==0L),
           all(result$prob>0&result$prob<1),
           all(is.finite(result$prob)),
           !("fit" %in% names(model)),
           !("teacher" %in% names(model)))
 Y<-X
 for(i in seq_len(n)){
   row<-unlist(X[i,,drop=FALSE],use.names=FALSE)
   Y[i,]<-exp(row/runif(1,2,4))*runif(1,.4,2)
 }
 pred<-wt_predict(model,Y)
 stopifnot(identical(result$prob,pred$prob),
           identical(result$n_evaluated,pred$n_evaluated),
           identical(result$path,pred$path))
}
cat("RLRWT_ORIGINAL_RLR_BETA_TWO_PAIR_NODE_AND_SIGN_ABLATION=PASS\n")
cat("RLRWT_AUTONOMOUS_MAX_FOUR_PAIR_COMPARISONS_PER_PATIENT=PASS\n")
cat("RLRWT_STRICT_PATIENTWISE_MONOTONE_RANK_INVARIANCE=PASS\n")
cat("RLRWT_MECHANICAL_TESTS=PASS\n")
