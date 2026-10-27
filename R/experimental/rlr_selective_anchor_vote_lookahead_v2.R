# SAVS v2 — train-only ROOT LOOKAHEAD: evaluate a small set of
# ordinal anchor gates on THEIR FINAL sparse branch-vote predictions
# rather than fitting a good root-only classifier and hoping branch voting
# benefits. No E outcomes, ratios, raw expression thresholds, leaf glmnet.
source("R/experimental/rlr_selective_anchor_vote_stump_v1.R")
sv2_need<-rg_need
sv2_roots<-function(D,y,pool,max_hubs=5L,max_roots=12L){
 R<-rg_ranks(D)
 hubs<-aa_hubs(pool,R,max_hubs=max_hubs)
 seen<-new.env(hash=TRUE,parent=emptyenv())
 candidate<-list()
 for(g in hubs)for(k in c(3L,5L)){
  anc<-tryCatch(aa_anchors(R,y,pool,g,k),error=function(e)NULL)
  if(is.null(anc))next
  for(t in seq_len(k)){
   node<-os_node(g,anc,t)
   split<-os_split_quality(R,y,node,
       leaf_min=max(3L,as.integer(floor(.08*nrow(D)))),
       each_class_min=0L)
   if(is.null(split))next
   z<-os_gate(R,node)
   key<-paste(z,collapse="")
   if(exists(key,seen,inherits=FALSE))next
   assign(key,TRUE,seen)
   candidate[[length(candidate)+1L]]<-list(
      root=node,criterion=split$loss+.0025*(k-1L),
      nbranch=split$branch_counts)
  }
 }
 if(!length(candidate))return(list())
 o<-order(vapply(candidate,function(c)c$criterion,numeric(1)),
           seq_along(candidate))
 head(candidate[o],max_roots)
}
sv2_metrics<-function(y,p){
 yy<-sv_truth(y);h<-as.integer(p>=.5)
 f1<-function(cl){
  tp<-sum(yy==cl&h==cl);fp<-sum(yy!=cl&h==cl)
  fn<-sum(yy==cl&h!=cl)
  if(2*tp+fp+fn==0L)0 else 2*tp/(2*tp+fp+fn)
 }
 ba<-if(min(table(factor(yy,levels=0:1)))>0)
  .5*(mean(h[yy==0L]==0L)+mean(h[yy==1L]==1L))
 else NA_real_
 list(macroF1=mean(c(f1(0L),f1(1L))),positiveF1=f1(1L),
      BA=ba,LL=rg_ll(yy,p))
}
sv2_cv_evaluate<-function(D,y,pool,root,mode,folds){
 p<-rep(NA_real_,nrow(D))
 for(split in folds){
  tr<-split$train_idx;va<-split$test_idx
  yy<-factor(as.character(y[tr]),levels=levels(y))
  tree<-sv_fit(D[tr,,drop=FALSE],yy,pool=pool,
               root=root,mode=mode)
  p[va]<-sv_predict(tree,D[va,,drop=FALSE])$p
 }
 sv2_need(all(is.finite(p))&&all(p>0&p<1),
          "nonfinite OOF probability for a fixed root")
 sv2_metrics(y,p)
}
sv2_train_only_select<-function(D,y,pool,mode="RLR10",max_roots=12L){
 fallback<-list(accepted=FALSE,root=NULL,
   proposed=0L,evaluated=0L,mode=mode,
   cv_macroF1=NA_real_,cv_BA=NA_real_,cv_LL=NA_real_,
   cv_global_macroF1=NA_real_,
   reason="NOT_ENOUGH_TRAINING_INFORMATION")
 if(min(table(y))<2L)return(fallback)
 folds<-make_repeated_stratified_folds(y,2L,1L,731L)
 candidates<-sv2_roots(D,y,pool,max_roots=max_roots)
 fallback$proposed<-length(candidates)
 if(!length(candidates))return(fallback)
 baseline<-rep(NA_real_,nrow(D))
 for(sp in folds) {
  tr<-sp$train_idx;te<-sp$test_idx
  yy<-factor(as.character(y[tr]),levels=levels(y))
  m<-sv_global4(D[tr,,drop=FALSE],yy,pool,
                  sv_m(ncol(D),mode))
  baseline[te]<-sv_predict(m,D[te,,drop=FALSE])$p
 }
 global_metrics<-sv2_metrics(y,baseline)
 fallback$cv_global_macroF1<-global_metrics$macroF1
 results<-lapply(candidates,function(c){
  v<-sv2_cv_evaluate(D,y,pool,c$root,mode,folds)
  list(root=c$root,macroF1=v$macroF1,BA=v$BA,LL=v$LL,
       criterion=v$macroF1+.10*v$BA-.02*v$LL -
         .001*(length(c$root$anchors)-3L))
 })
 fallback$evaluated<-length(results)
 o<-order(-vapply(results,function(z)z$criterion,numeric(1)),
         seq_along(results))
 best<-results[[o[1L]]]
 # Nested OOF D selection only. Global four-rule alternative remains
 # mandatory; do not force a tree just for tree novelty.
 fallback$accepted<-isTRUE(
    best$macroF1>=global_metrics$macroF1+.01 &&
    best$BA>=global_metrics$BA-.02)
 fallback$root<-best$root
 fallback$cv_macroF1<-best$macroF1
 fallback$cv_BA<-best$BA
 fallback$cv_LL<-best$LL
 fallback$reason<-"FIXED_TWO_FOLD_D_RANK_LOOKAHEAD_VS_GLOBAL_FOUR"
 fallback
}
sv2_final_fit<-function(D,y,pool,select){
 if(isTRUE(select$accepted)){
  sv_fit(D,y,pool=pool,root=select$root,mode=select$mode)
 }else{
  sv_global4(D,y,pool,sv_m(ncol(D),select$mode))
 }
}
