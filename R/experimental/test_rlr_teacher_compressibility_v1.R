source("R/experimental/rlr_teacher_compressibility_core_v1.R")
set.seed(20261029L)
n<-105L;p<-32L
D<-as.data.frame(matrix(rnorm(n*p),nrow=n,ncol=p))
colnames(D)<-paste0("g",seq_len(p))
rownames(D)<-paste0("sample",seq_len(n))
pairs<-data.frame(gene_left=paste0("g",seq(1,15,2)),
                  gene_right=paste0("g",seq(2,16,2)))
A<-rrv5_matrix(D,pairs)
y<-factor(as.integer(A[,1]==1L | A[,2]==1L),levels=0:1)
stopifnot(min(table(y))>=10L)
soft<-plogis(-1.3+3*A[,3]-2.5*A[,4]+.2*rnorm(n))
F<-lapply(c("TEACHER","LABEL"),function(goal)
  lapply(0:3,function(k)rc_fit(D,y,pairs,soft,cap=k,
                                min_leaf=6L,goal=goal)))
for(family in F){
  scores<-vapply(family,function(z)z$train_loss,numeric(1))
  stopifnot(all(diff(scores)<=1e-8))
  for(fit in family){
    pr<-rc_predict(fit,D)
    stopifnot(fit$internal<=fit$cap,fit$n_unique_pairs<=3L,
              fit$n_unique_genes<=6L,max(pr$checks)<=2L,
              !("p_oof" %in% names(fit)))
    Y<-D
    for(i in seq_len(n)){
      row<-as.numeric(unlist(D[i,,drop=FALSE],use.names=FALSE))
      Y[i,]<-exp(row/runif(1L,2,4))*runif(1L,.6,1.5)
    }
    qq<-rc_predict(fit,Y)
    stopifnot(identical(qq$prob_teacher,pr$prob_teacher),
              identical(qq$prob_label,pr$prob_label),
              identical(qq$checks,pr$checks),
              identical(qq$path,pr$path))
  }
}
# Exhaustively enumerate every root and both children independently
# to check TRAIN-ONLY TEACHER SSE optimum for every capacity.
ell<-qlogis(soft)
sse<-function(ix)sum((ell[ix]-mean(ell[ix]))^2)
splitloss<-function(ix,j){
 if(j==0L)return(sse(ix))
 lo<-ix[A[ix,j]==0L];hi<-ix[A[ix,j]==1L]
 if(min(length(lo),length(hi))<6L)return(Inf)
 sse(lo)+sse(hi)
}
for(k in 0:3){
 best<-sse(seq_len(n))
 if(k>0L)for(j in seq_len(ncol(A))){
  lo<-which(A[,j]==0L);hi<-which(A[,j]==1L)
  if(min(length(lo),length(hi))<6L)next
  for(a in 0:ncol(A))for(b in 0:ncol(A)){
   if(a==j||b==j||1L+as.integer(a>0)+as.integer(b>0)>k)next
   best<-min(best,splitloss(lo,a)+splitloss(hi,b))
  }
 }
 stopifnot(abs(best-F[[1]][[k+1]]$train_loss)<1e-7)
}
cat("RLR_CAPACITY_GLOBAL_EXHAUSTIVE_BRUTE_FORCE=PASS\n")
cat("RLR_CAPACITY_TRAINING_FRONTIER_MONOTONE=PASS\n")
cat("RLR_CAPACITY_STANDALONE_MONOTONE_RANK_PREDICTION=PASS\n")
cat("RLR_CAPACITY_UNIT_TEST=PASS\n")
