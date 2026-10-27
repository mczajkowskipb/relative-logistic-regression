# COT-RLR V1: Conditional Ordinal Threshold stump atop original RLR-50.
# Threshold-only decision layer; original RLR continuous probabilities,
# ROC/AUC, Brier and logloss remain exactly unchanged.
# All learned decisions use within-patient rank-based gate and RLR score.
source("R/experimental/rlr_selective_anchor_gated_v1.R")

ct_truth<-sa_truth
ct_need<-rg_need
ct_f1<-function(y,pred) {
  tp<-sum(y==1L & pred==1L)
  fp<-sum(y==0L & pred==1L)
  fn<-sum(y==1L & pred==0L)
  if(2*tp+fp+fn==0L)0 else 2*tp/(2*tp+fp+fn)
}
ct_macro_f1<-function(y,pred) {
  f<-function(cl) {
    tp<-sum(y==cl & pred==cl)
    fp<-sum(y!=cl & pred==cl)
    fn<-sum(y==cl & pred!=cl)
    if(2*tp+fp+fn==0L)0 else 2*tp/(2*tp+fp+fn)
  }
  mean(c(f(0L),f(1L)))
}
ct_ba<-function(y,pred){
  if(!all(c(0L,1L) %in% y))return(NA_real_)
  .5*(mean(pred[y==0L]==0L)+mean(pred[y==1L]==1L))
}
ct_grid<-round(seq(.10,.90,by=.05),digits=8L)

ct_global<-function(y,p,grid=ct_grid) {
  candidates<-lapply(grid,function(t) {
    h<-as.integer(p>=t)
    list(threshold=t,F1=ct_f1(y,h),BA=ct_ba(y,h),
         M=ct_macro_f1(y,h))
  })
  # Guard against an apparently high F1 obtained by predicting almost
  # every sample positive. The ordinary .5 decision is always feasible.
  base_BA<-ct_ba(y,as.integer(p>=.5))
  valid<-which(vapply(candidates,function(t)t$BA>=base_BA-.02,logical(1)))
  ct_need(length(valid)>0L,"no valid global threshold")
  # Fixed tie-breaking: maximize positive F1 at non-inferior BA, then
  # prefer BA and a threshold close to .5. Only D OOF evidence is used.
  o<-order(-vapply(candidates[valid],function(t)t$F1,numeric(1)),
           -vapply(candidates[valid],function(t)t$BA,numeric(1)),
           abs(grid[valid]-.5),grid[valid])
  candidates[[valid[o[1L]]]]
}

ct_confusions<-function(y,p,ids,grid) {
  vapply(grid,function(t) {
    h<-as.integer(p[ids]>=t); yy<-y[ids]
    c(TP=sum(h==1L&yy==1L),
      FP=sum(h==1L&yy==0L),
      FN=sum(h==0L&yy==1L),
      TN=sum(h==0L&yy==0L))
  },numeric(4L))
}
ct_cm_metrics<-function(cm) {
  tp<-cm[1L];fp<-cm[2L];fn<-cm[3L];tn<-cm[4L]
  f1<-if(2*tp+fp+fn==0)0 else 2*tp/(2*tp+fp+fn)
  ba<-if(min(tp+fn,tn+fp)<=0)NA_real_ else
    .5*(tp/(tp+fn)+tn/(tn+fp))
  list(F1=f1,BA=ba)
}

# Fit both thresholds conditional on same rank gate.
# Requiring lower F1 gains on OOF D and no BA loss prevents
# a degenerate all-positive classifier.
ct_fit_root<-function(y,p,gate,global,grid=ct_grid) {
  i0<-which(gate==0L);i1<-which(gate==1L)
  if(min(length(i0),length(i1))<9L)return(NULL)
  A<-ct_confusions(y,p,i0,grid)
  B<-ct_confusions(y,p,i1,grid)
  best<-NULL
  for(i in seq_along(grid))for(j in seq_along(grid)){
    cm<-A[,i]+B[,j]
    q<-ct_cm_metrics(cm)
    if(!is.finite(q$BA)||q$BA < global$BA-.015)next
    if(q$F1 < global$F1+.01)next
    candidate<-list(t0=grid[i],t1=grid[j],F1=q$F1,
        BA=q$BA,gain=q$F1-global$F1)
    if(is.null(best)||
       candidate$F1>best$F1+1e-10||
       (abs(candidate$F1-best$F1)<1e-10 &&
        (candidate$BA>best$BA+1e-10 ||
         (abs(candidate$BA-best$BA)<1e-10 &&
          (abs(candidate$t0-.5)+abs(candidate$t1-.5) <
            abs(best$t0-.5)+abs(best$t1-.5)-1e-10)))))
      best<-candidate
  }
  if(!is.null(best) && abs(best$t0-best$t1)<1e-10)return(NULL)
  best
}

ct_predict<-function(p,R,model) {
  if(is.null(model))return(as.integer(p>=.5))
  if(model$type=="GLOBAL")return(as.integer(p>=model$t))
  ct_need(model$type=="TREE","invalid CT decision model")
  g<-os_gate(R,model$root)
  as.integer(p>=ifelse(g==0L,model$t0,model$t1))
}

ct_select<-function(T,y) {
  sp<-rg_inner(y)
  D<-T[sp$D,,drop=FALSE];V<-T[sp$V,,drop=FALSE]
  yd<-factor(as.character(y[sp$D]),levels=levels(y))
  yv<-factor(as.character(y[sp$V]),levels=levels(y))
  yy<-ct_truth(yd); vv<-ct_truth(yv)
  fallback<-list(accepted=FALSE,root=NULL,t0=NA_real_,t1=NA_real_,
      global_t=.5,proposed=0L,candidates=0L,D_F1_gain=NA_real_,
      V_F1_gain=NA_real_,V_BA_gain=NA_real_,reason="NO_VALID_ROOT")
  # A strictly original, train-only crossfitted K50 score is necessary to
  # avoid choosing thresholds from apparent in-sample RLR probabilities.
  # Guard because original glmnet internal 5-fold cannot fit <4/class.
  if(safe_nfolds(yd,5L)<4L)
    return(modifyList(fallback,list(reason="UNDERPOWERED_D_FOR_RLR")))
  plans<-make_repeated_stratified_folds(yd,3L,1L,731L)
  if(any(vapply(plans,function(q){
     yy2<-factor(as.character(yd[q$train_idx]),levels=levels(yd))
     safe_nfolds(yy2,5L)<4L
  },logical(1))))
    return(modifyList(fallback,list(reason="UNDERPOWERED_OOF_SUBTRAIN")))
  pool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
  oof<-sa_parent_oof(D,yd,K=50L,seed=731L)
  global<-ct_global(yy,oof)
  fallback$global_t<-global$threshold
  # All candidate anchors depend ONLY on D and its labels and ranks.
  proposals<-sa_candidate_roots(D,yd,pool,max_proposals=240L)
  fallback$proposed<-length(proposals)
  ranked<-list()
  for(z in proposals){
    best<-ct_fit_root(yy,oof,z$g,global)
    if(is.null(best))next
    ranked[[length(ranked)+1L]]<-list(
      root=z$root,t0=best$t0,t1=best$t1,
      DF1=best$F1,DBA=best$BA,
      gain=best$gain)
  }
  fallback$candidates<-length(ranked)
  if(length(ranked)==0L)
    return(modifyList(fallback,list(reason="NO_D_OUTPERFORMING_THRESHOLD")))
  order_idx<-order(
      -vapply(ranked,function(z)z$gain,numeric(1)),
      -vapply(ranked,function(z)z$DBA,numeric(1)))
  shortlist<-ranked[head(order_idx,3L)]
  parent<-rg_parent(D,yd,V,pool,K=50L)
  global_V<-as.integer(parent$test>=global$threshold)
  orig_V<-as.integer(parent$test>=.5)
  fVglob<-ct_f1(vv,global_V)
  fVorig<-ct_f1(vv,orig_V)
  bVglob<-ct_ba(vv,global_V)
  best<-NULL
  for(z in shortlist){
    pred<-as.integer(parent$test>=
         ifelse(os_gate(aa_ranks(V),z$root)==0L,z$t0,z$t1))
    f1<-ct_f1(vv,pred)
    ba<-ct_ba(vv,pred)
    vgain<-f1-max(fVglob,fVorig)
    if(is.null(best)||vgain>best$vgain) {
      best<-list(z=z,vgain=vgain,bgain=ba-bVglob)
    }
  }
  # This acceptance gate is deliberately strict: root must outperform
  # BOTH ordinary 0.5 RLR and a GLOBAL D-OOF-tuned decision threshold.
  accepted<-best$vgain>=.02 && best$bgain>=-.02
  list(accepted=accepted,root=best$z$root,
    t0=best$z$t0,t1=best$z$t1,global_t=global$threshold,
    proposed=length(proposals),candidates=length(ranked),
    D_F1_gain=best$z$gain,V_F1_gain=best$vgain,
    V_BA_gain=best$bgain,reason="FIXED_D_OOF_V_THRESHOLD_GATE")
}

ct_evaluate<-function(T,yT,E,yE) {
  sel<-ct_select(T,yT)
  poolT<-compute_top_pairs_switchbox(T,yT,top_n=1000L,buffer_mult=5L)
  original<-rg_parent(T,yT,E,poolT,K=50L)
  p<-original$test
  baseline<-ct_predict(p,aa_ranks(E),NULL)
  global<-ct_predict(p,aa_ranks(E),list(type="GLOBAL",t=sel$global_t))
  candidate<-if(is.null(sel$root))global else
     ct_predict(p,aa_ranks(E),list(
       type="TREE",root=sel$root,t0=sel$t0,t1=sel$t1))
  gated<-if(isTRUE(sel$accepted))candidate else global
  list(probabilities=p,labels=list(
        RLR50_NATIVE_05=baseline,
        RLR50_DOOF_GLOBAL_THRESHOLD=global,
        RLR50_CONDITIONAL_TREE_UNGATED=candidate,
        RLR50_CONDITIONAL_TREE_V_GATED=gated),
       selection=sel,truth=ct_truth(yE))
}
