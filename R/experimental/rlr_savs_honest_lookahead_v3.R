# SAVS V3. Honest rank-conditional shallow-tree nomination.
# Key V2 flaw corrected: V2 shortlisted tree root AND switchBox genes from
# ALL D before using its own OOF CV. The CV ranking was therefore partly
# supervised by its validation labels. V3 chooses both genes and root only
# on S (strict subset of D), chooses structure on H (D-internal), and
# evaluates on V untouched by either operation. Original outer E unread.
source("R/experimental/rlr_selective_anchor_vote_lookahead_v2.R")

sv3_split<-function(y,seed=991L) {
  parts<-make_repeated_stratified_folds(y,3L,1L,seed)
  a<-parts[[1L]]
  stopifnot(length(a$train_idx)>0L,length(a$test_idx)>0L,
    !length(intersect(a$train_idx,a$test_idx)))
  list(S=a$train_idx,H=a$test_idx)
}
sv3_select <- function(D,y,mode="RLR10",max_roots=12L) {
  ids<-sv3_split(y)
  S<-D[ids$S,,drop=FALSE]
  H<-D[ids$H,,drop=FALSE]
  ys<-factor(as.character(y[ids$S]),levels=levels(y))
  yh<-factor(as.character(y[ids$H]),levels=levels(y))
  fallback<-list(accepted=FALSE,root=NULL,mode=mode,
     n_S=nrow(S),n_H=nrow(H),n_proposals=0L,
     H_global_F1=NA_real_,H_best_tree_F1=NA_real_,
     H_global_BA=NA_real_,H_best_tree_BA=NA_real_,
     reason="NO_SUPPORT")
  if(min(table(ys))<2L || min(table(yh))<1L)
    return(fallback)
  # TRAIN-ONLY source. No labels or rank features from H are used here.
  poolS<-compute_top_pairs_switchbox(S,ys,top_n=1000L,buffer_mult=5L)
  roots<-sv2_roots(S,ys,poolS,max_hubs=8L,max_roots=max_roots)
  fallback$n_proposals<-length(roots)
  ctrl<-sv_global4(S,ys,poolS,sv_m(ncol(D),mode))
  base<-sv2_metrics(yh,sv_predict(ctrl,H)$p)
  fallback$H_global_F1<-base$macroF1
  fallback$H_global_BA<-base$BA
  if(length(roots)==0L){
    fallback$reason<-"NO_VALID_S_ONLY_ROOT"
    return(fallback)
  }
  evaluate<-lapply(roots,function(z){
    trained<-sv_fit(S,ys,pool=poolS,root=z$root,mode=mode)
    v<-sv2_metrics(yh,sv_predict(trained,H)$p)
    list(root=z$root,positiveF1=v$positiveF1,
         macroF1=v$macroF1,BA=v$BA,LL=v$LL)
  })
  ix<-order(
   -vapply(evaluate,function(z)z$macroF1,numeric(1)),
   -vapply(evaluate,function(z)z$BA,numeric(1)),
    vapply(evaluate,function(z)z$LL,numeric(1)),
    vapply(evaluate,function(z)length(z$root$anchors),integer(1)),
    seq_along(evaluate))
  best<-evaluate[[ix[1L]]]
  fallback$H_best_tree_F1<-best$macroF1
  fallback$H_best_tree_BA<-best$BA
  # Prevent false wins from tiny H partitions: require both a measurable
  # Macro-F1 advantage and no class-balanced accuracy regression.
  # The control must be trained on the SAME S and evaluated on SAME H.
  fallback$accepted<-isTRUE(best$macroF1>=base$macroF1+.01 &&
      best$BA>=base$BA-.02)
  fallback$root<-best$root
  fallback$reason<-if(fallback$accepted)
        "HONEST_S_H_TREE_ACCEPTED" else "HONEST_S_H_GLOBAL_FALLBACK"
  fallback
}
sv3_refit<-function(D,y,poolD,decision){
  if(isTRUE(decision$accepted)) {
    sv_fit(D,y,pool=poolD,root=decision$root,mode=decision$mode)
  }else {
    sv_global4(D,y,poolD,sv_m(ncol(D),decision$mode))
  }
}
