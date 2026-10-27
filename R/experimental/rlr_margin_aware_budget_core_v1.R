# MBRLR v1: rank-margin-aware disjoint pair ordering, explicitly aligned
# with ORIGINAL frozen RLR signed ternary feature representation m10.
#
# Outer T-only supervised re-order of the original top1000 SWITCHBOX pool;
# original RLR logistic and d=2 global gene-disjointness unchanged.
# NO test-time raw gene values, ratios, sample-range calibration or E labels.
source("R/experimental/rlr_rank_gate_core_v1.R")

mb_need<-rg_need
mb_rank_margin_scores<-function(X,y,pairs,m=10L){
  mb_need(is.factor(y)&&length(levels(y))==2L &&
          nrow(X)==length(y) && nrow(pairs)>0L,
          "invalid training set or binary labels")
  r<-rg_ranks(X)
  gap<-rg_gaps(r,pairs)
  z<-ifelse(gap>=m,1,ifelse(gap<= -m,-1,0))
  ix0<-which(y==levels(y)[1L])
  ix1<-which(y==levels(y)[2L])
  mb_need(length(ix0)>0L && length(ix1)>0L,"one class absent")
  mu0<-colMeans(z[ix0,,drop=FALSE])
  mu1<-colMeans(z[ix1,,drop=FALSE])
  ac0<-colMeans(abs(z[ix0,,drop=FALSE]))
  ac1<-colMeans(abs(z[ix1,,drop=FALSE]))
  # Direct empirical class difference in EXACT feature used by RLR,
  # unlike TSP's binary sign scoring. Mild support penalty discourages
  # pairs whose feature is near-zero in either class.
  diff<-abs(mu1-mu0)
  score<-diff * (0.25 + 0.75*pmin(ac0,ac1))
  data.frame(score=score,margin_diff=diff,
             active_class0=ac0,active_class1=ac1,
             pct_zero=1-colMeans(abs(z)),
             stringsAsFactors=FALSE)
}
mb_reorder<-function(X,y,pool,m=10L){
  mb_need(nrow(pool)>1L,"pair pool insufficient")
  stats<-mb_rank_margin_scores(X,y,pool,m)
  order<-order(-stats$score,-stats$margin_diff,
               seq_len(nrow(pool)))
  re<-pool[order,,drop=FALSE]
  re$margin_aware_score<-stats$score[order]
  re$margin_aware_diff<-stats$margin_diff[order]
  re$margin_zero_pct<-stats$pct_zero[order]
  rownames(re)<-NULL
  re
}
mb_active<-function(fit,X){
  pairs<-fit$selected
  mat<-encode_pair_matrix(X,pairs,m=10L)
  co<-as.numeric(stats::coef(fit$fit,s="lambda.min"))[-1L]
  mb_need(length(co)==ncol(mat),"glmnet coefficients and pairs mismatch")
  selected<-abs(co)>1e-9
  loc<-if(any(selected))
    rowSums(mat[,selected,drop=FALSE]!=0L) else rep(0L,nrow(X))
  list(nominal=nrow(pairs),globally_nonzero=sum(selected),
       locally_active=as.integer(loc),
       tested=rep(sum(selected),nrow(X)),
       n_genes=if(any(selected))
          length(unique(c(pairs$gene_left[selected],
                          pairs$gene_right[selected]))) else 0L)
}
mb_budget<-function(x,max_k=50L){
  k<-as.integer(round(x))
  min(max(2L,k),max_k)
}
