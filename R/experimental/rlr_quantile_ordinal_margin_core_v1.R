# QRLR v1 — training-only pair-specific ORDINAL dead zone.
# Preserve original switchBox top1000, d=2, K10/20, elastic-net alpha=.5,
# lambda.min. Rather than fixed m=10 for every pair, choose per-pair
# m_j = max(10, empirical type1 q-quantile of abs(rank_j-rank_i) in T).
# Entirely class-agnostic threshold learnt from T; select q=0/.25/.5 fixed
# before any inner V. No raw expression cutoffs / ratios / out-of-range
# calibration. Original RLR q=0 MUST match exactly.
source("R/experimental/rlr_rank_gate_core_v1.R")
qr_need<-rg_need

qr_threshold <- function(X,pairs,q=0.25){
 qr_need(is.data.frame(X) && nrow(pairs)>0L &&
  length(q)==1L && is.finite(q) && q>=0 && q<=.75,
  "invalid original T ranks/pair pool/quantile")
 gaps<-rg_gaps(rg_ranks(X),pairs)
 if(q==0) return(rep(10, ncol(gaps)))
 threshold<-vapply(seq_len(ncol(gaps)),function(j){
   max(10,as.numeric(stats::quantile(abs(gaps[,j]),probs=q,
                                      type=1,names=FALSE)))
 },numeric(1))
 qr_need(all(is.finite(threshold)) && all(threshold>=10),
          "invalid nonclass rank thresholds")
 threshold
}
qr_encode<-function(X,pairs,margin){
 qr_need(nrow(pairs)==length(margin)&&all(margin>=10),
         "missing rank thresholds per pair")
 g<-rg_gaps(rg_ranks(X),pairs)
 cutoff<-sweep(abs(g),2L,margin,">=")
 zz<-sign(g)*cutoff
 storage.mode(zz)<-"double"
 dimnames(zz)<-list(rownames(X),pairs$ordered_pair)
 zz
}
qr_parent <- function(T,yT,E,pool,K=20L,q=.25){
 qr_need(identical(names(T),names(E))&&nrow(T)==length(yT),
         "schema or train size mismatch")
 pairs<-apply_disjointness(pool,K=K,d=2L)
 qr_need(nrow(pairs)>1L,"pair panel too small")
 margin<-qr_threshold(T,pairs,q)
 zt<-qr_encode(T,pairs,margin)
 ze<-qr_encode(E,pairs,margin)
 set.seed(42L)
 cv<-glmnet::cv.glmnet(x=zt,y=as.numeric(yT==levels(yT)[2L]),
   family="binomial",alpha=.5,nfolds=safe_nfolds(yT,5L),
   type.measure="deviance")
 co<-as.numeric(stats::coef(cv,s="lambda.min"))[-1L]
 qr_need(length(co)==ncol(ze),"coefficient alignment")
 active<-abs(co)>1e-9
 votes<-if(any(active))rowSums(ze[,active,drop=FALSE]!=0L)
        else rep(0L,nrow(E))
 tt<-if(any(active))rowSums(zt[,active,drop=FALSE]!=0L)
        else rep(0L,nrow(T))
 list(test=as.numeric(predict(cv,newx=ze,s="lambda.min",type="response")),
  fit=cv,selected=pairs,margin=margin,quantile=q,
  coeff=co,global_active=sum(active),
  test_active=as.integer(votes),train_active=as.integer(tt),
  n_tested=sum(active),
  max_test_rank=ncol(E))
}
