# SA-RLR V1: Selective Anchor-Gated Relative Logistic Regression.
# Prototype scientific algorithm; NOT evaluated on original heldout E cohorts.
# Original RLR source/lock untouched. New and original predictors use only
# within-patient gene ranks, never raw values or expression ratios.
source("R/experimental/rlr_anchor_lookahead_child_rlr_v1.R")

sa_need <- rg_need
sa_candidate_roots <- function(D,y,pool,max_hubs=8L,
                               max_swap_alternatives=4L,
                               max_proposals=240L) {
  R<-aa_ranks(D);stats<-aa_rank_stats(R,y)
  hubs<-aa_hubs(pool,R,max_hubs)
  candidates<-list()
  seen<-new.env(parent=emptyenv(),hash=TRUE)
  for(g in hubs)for(k in c(3L,5L)) {
    main<-tryCatch(aa_anchors(R,y,pool,g,k),error=function(e)NULL)
    if(is.null(main))next
    neighbors<-unique(c(
      as.character(pool$gene_left[pool$gene_right==g]),
      as.character(pool$gene_right[pool$gene_left==g])))
    universe<-setdiff(unique(c(neighbors,
      names(sort(stats$shift+0.15*stats$spread))[seq_len(
        min(length(stats$shift),30L))])),g)
    alternatives<-head(setdiff(universe,main),max_swap_alternatives)
    panels<-list(main)
    # One-swap lookahead on reference gene IDs; still polynomial and capped.
    for(pos in seq_len(k))for(a in alternatives) {
      trial<-main
      trial[pos]<-a
      if(length(unique(trial))==k)panels[[length(panels)+1L]]<-trial
    }
    for(panel in panels)for(threshold in seq_len(k)) {
      root<-os_node(g,panel,threshold)
      groups<-os_gate(R,root)
      idx0<-which(groups==0L);idx1<-which(groups==1L)
      if(min(length(idx0),length(idx1))<10L)next
      # A single eligible leaf is enough; previous lookahead mistakenly
      # required both child RLR models to be learnable.
      truth<-sa_truth(y)
      possible<-vapply(list(idx0,idx1),function(idx)
         length(idx)>=12L &&
         min(tabulate(truth[idx]+1L,nbins=2L))>=3L,logical(1))
      if(!any(possible))next
      # Do not spend expensive model fits on equivalent patient partitions.
      key<-paste(groups,collapse="")
      if(exists(key,seen,inherits=FALSE))next
      assign(key,TRUE,seen)
      candidates[[length(candidates)+1L]]<-list(
         root=root,g=groups,valid_branches=possible,
         n_branch=c(length(idx0),length(idx1)))
      if(length(candidates)>=max_proposals)return(candidates)
    }
  }
  candidates
}
sa_truth <- function(y) as.integer(y==levels(y)[2L])

sa_parent_oof <- function(D,y,K=50L,seed=731L) {
  plans<-make_repeated_stratified_folds(y,3L,1L,seed)
  p<-rep(NA_real_,nrow(D))
  for(sp in plans){
    tr<-sp$train_idx;te<-sp$test_idx
    yy<-factor(as.character(y[tr]),levels=levels(y))
    pool<-compute_top_pairs_switchbox(D[tr,,drop=FALSE],yy,
                 top_n=1000L,buffer_mult=5L)
    fit<-rg_parent(D[tr,,drop=FALSE],yy,D[te,,drop=FALSE],pool,K=K)
    p[te]<-fit$test
  }
  sa_need(all(is.finite(p))&&all(p>0&p<1),"OOF missing")
  p
}

sa_screen_roots <- function(D,y,parent_oof,pool,
                            max_top=6L,max_proposals=240L){
  proposals<-sa_candidate_roots(D,y,pool,
             max_proposals=max_proposals)
  if(length(proposals)==0L)
    return(list(roots=list(),scores=numeric(),proposed=0L))
  shortlist<-pool[seq_len(min(nrow(pool),180L)),,drop=FALSE]
  Z<-encode_pair_matrix(D,shortlist,m=10L)
  yy<-sa_truth(y)
  residual<-yy-parent_oof
  weights<-pmax(1e-5,parent_oof*(1-parent_oof))
  # Global correction proxy is subtracted. A branch-specific feature must
  # show incremental improvement beyond one global main-effect coefficient.
  fullG<-colSums(Z*residual)
  fullH<-colSums((Z^2)*weights)
  baseline<-fullG^2/(2*(fullH+5))
  score<-vapply(proposals,function(entry) {
    branch_score<-numeric()
    for(branch in 0:1) {
      if(!entry$valid_branches[[branch+1L]])next
      id<-which(entry$g==branch)
      G<-colSums(Z[id,,drop=FALSE]*residual[id])
      H<-colSums((Z[id,,drop=FALSE]^2)*weights[id])
      candidate<-G^2/(2*(H+5)) -
                 baseline*(length(id)/nrow(D))
      branch_score<-c(branch_score,max(candidate,0))
    }
    if(length(branch_score)==0L)return(-Inf)
    max(branch_score)-0.001*(length(entry$root$anchors)-1L)
  },numeric(1))
  o<-order(-score,seq_along(score))
  top<-head(o,max_top)
  list(roots=proposals[top],scores=score[top],
       proposed=length(proposals))
}

sa_fit <- function(D,y,original_D,root,fixed_pairs=NULL) {
  g<-os_gate(aa_ranks(D),root)
  child<-vector("list",2L)
  for(b in 0:1) {
    idx<-which(g==b)
    yy<-factor(as.character(y[idx]),levels=levels(y))
    viable<-length(idx)>=12L && min(table(yy))>=3L
    if(!viable || (!is.null(fixed_pairs) &&
                    is.null(fixed_pairs[[b+1L]]))){
      child[[b+1L]]<-list(pairs=NULL,fit=NULL,active=0L,
              state=if(viable)"frozen_inactive" else "not_viable")
      next
    }
    if(is.null(fixed_pairs)){
      candidate<-la_fit_child(D[idx,,drop=FALSE],yy,
                         original_D[idx],select_pairs=TRUE)
    } else {
      candidate<-la_fit_child(D[idx,,drop=FALSE],yy,
         original_D[idx],pairs=fixed_pairs[[b+1L]],
         select_pairs=FALSE)
    }
    if(is.null(candidate$fit)||candidate$active==0L){
      candidate$fit<-NULL
      candidate$active<-0L
    }
    child[[b+1L]]<-candidate
  }
  list(root=root,children=child,
       pairs=lapply(child,function(x)if(is.null(x$fit))NULL else x$pairs),
       n_fitted=sum(vapply(child,function(x)!is.null(x$fit),logical(1))),
       n_active=sum(vapply(child,function(x)x$active,integer(1))))
}

sa_predict <- function(E,original_E,fit) {
  g<-os_gate(aa_ranks(E),fit$root)
  p<-original_E
  for(b in 0:1) {
    idx<-which(g==b)
    if(length(idx)==0L||is.null(fit$children[[b+1L]]$fit))next
    p[idx]<-la_child_predict(E[idx,,drop=FALSE],
                         original_E[idx],fit$children[[b+1L]])
  }
  p
}

# D proposes one-swap reference sets; V ranks <=6 candidates with actual
# child RLR fits and accepts only if both parent/global comparators improve.
sa_select <- function(T,yT) {
  sp<-rg_inner(yT)
  D<-T[sp$D,,drop=FALSE];V<-T[sp$V,,drop=FALSE]
  yd<-factor(as.character(yT[sp$D]),levels=levels(yT))
  yv<-factor(as.character(yT[sp$V]),levels=levels(yT))
  # cv.glmnet rejects nfolds <=3. Original RLR on full outer T can
  # still be fitted; therefore insufficient nested D training support
  # must return an explicit SAFE fallback, not crash or invent tree labels.
  small_fallback<-function(reason)list(root=NULL,pairs=NULL,
      accepted=FALSE,considered=0L,proposed=0L,
      working=0L,active=0L,v_gain_parent=NA_real_,
      v_gain_samepairs=NA_real_,v_BA=NA_real_,reason=reason)
  if(safe_nfolds(yd,5L)<4L)
    return(small_fallback("D_CLASS_TOO_SMALL_FOR_GLMNET_4FOLD"))
  oofplans<-make_repeated_stratified_folds(yd,3L,1L,731L)
  if(any(vapply(oofplans,function(q) {
    yy<-factor(as.character(yd[q$train_idx]),levels=levels(yd))
    safe_nfolds(yy,5L)<4L
  },logical(1))))
    return(small_fallback("D_OOF_SUBTRAIN_TOO_SMALL_FOR_GLMNET_4FOLD"))
  pool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
  parent<-rg_parent(D,yd,V,pool,K=50L)
  parent_oof<-sa_parent_oof(D,yd,K=50L)
  ranked<-sa_screen_roots(D,yd,parent_oof,pool,
                            max_top=6L,max_proposals=240L)
  fallback<-list(root=NULL,pairs=NULL,accepted=FALSE,
      considered=length(ranked$roots),proposed=ranked$proposed,
      working=0L,active=0L,v_gain_parent=NA_real_,
      v_gain_samepairs=NA_real_,v_BA=NA_real_,reason="no_viable_root")
  if(length(ranked$roots)==0L)return(fallback)
  yy<-sa_truth(yv)
  evaluated<-list()
  for(entry in ranked$roots) {
    fit<-sa_fit(D,yd,parent$train,entry$root)
    if(fit$n_fitted<1L ||fit$n_active==0L)next
    tree<-sa_predict(V,parent$test,fit)
    shared<-la_global_control(D,yd,parent$train,
                  V,parent$test,fit)
    ll<-rg_ll(yy,tree)
    evaluated[[length(evaluated)+1L]]<-list(
         fit=fit,ll=ll,global=shared,tree=tree)
  }
  if(length(evaluated)==0L)return(fallback)
  best<-evaluated[[which.min(vapply(evaluated,
                   function(x)x$ll,numeric(1)))]]
  gain0<-rg_ll(yy,parent$test)-best$ll
  gain1<-rg_ll(yy,best$global)-best$ll
  ba<-rg_ba(yy,best$tree)-rg_ba(yy,parent$test)
  accept<-gain0>=.005 && gain1>=.002 && ba>=-.02
  list(root=best$fit$root,pairs=best$fit$pairs,
       accepted=accept,working=best$fit$n_fitted,
       active=best$fit$n_active,
       proposed=ranked$proposed,considered=length(ranked$roots),
       v_gain_parent=gain0,v_gain_samepairs=gain1,
       v_BA=ba,reason="D_OOF_V_CHECK")
}

sa_evaluate <- function(T,yT,E,yE) {
  choice<-sa_select(T,yT)
  pool<-compute_top_pairs_switchbox(T,yT,top_n=1000L,buffer_mult=5L)
  parent<-rg_parent(T,yT,E,pool,K=50L)
  raw<-parent$test;shared<-parent$test
  n_fit<-n_active<-0L
  if(!is.null(choice$root)) {
    final<-sa_fit(T,yT,parent$train,choice$root,
                 fixed_pairs=choice$pairs)
    n_fit<-final$n_fitted;n_active<-final$n_active
    raw<-sa_predict(E,parent$test,final)
    shared<-la_global_control(T,yT,parent$train,
                      E,parent$test,final)
  }
  actual<-if(isTRUE(choice$accepted))raw else parent$test
  list(probs=list(ORIGINAL_RLR50=parent$test,
      SA_RLR_UNGATED_DIAGNOSTIC=raw,
      SA_RLR_SAMEPAIRS_GLOBAL=shared,
      SA_RLR_V_GATED=actual),
       choice=choice,fit_branches=n_fit,active=n_active,
       truth=sa_truth(yE))
}
