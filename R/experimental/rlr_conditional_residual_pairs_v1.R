# CRP-RLR V1: CONDITIONAL RESIDUAL PAIR SELECTION for local gene experts.
# Unlike local switchBox, candidates are ranked on D-only 3-fold original
# K50 RLR OOF residuals, so selected relations target *unexplained* signal.
# Original RLR code, original gene universe and E scores remain untouched.
source("R/experimental/rlr_selective_anchor_gated_v1.R")

crp_need<-rg_need
crp_select_pairs<-function(D,y,p_oof,gate,branch,candidates,Z=NULL,
                           K=8L,penalty=5) {
  crp_need(length(p_oof)==nrow(D) && length(gate)==nrow(D) &&
           branch %in% c(0L,1L), "invalid D OOF branch dimensions")
  if(is.null(Z))Z<-encode_pair_matrix(D,candidates,m=10L)
  crp_need(ncol(Z)==nrow(candidates),"misaligned candidate Z")
  idx<-which(gate==branch)
  yy<-as.integer(y==levels(y)[2L])
  if(length(idx)<12L || min(tabulate(yy[idx]+1L,nbins=2L))<3L)
    return(candidates[0,,drop=FALSE])
  e<-yy-p_oof
  w<-pmax(1e-5,p_oof*(1-p_oof))
  Zi<-Z[idx,,drop=FALSE]
  Gl<-colSums(Zi*e[idx])
  Hl<-colSums((Zi^2)*w[idx])
  G<-colSums(Z*e)
  H<-colSums((Z^2)*w)
  # Difference from the approximate benefit of a common global slope.
  # All sample- and class-level quantities here are D-only and
  # the p_oof vector was obtained by OOF original RLR50 fits on D.
  gain<-(Gl^2)/(2*(Hl+penalty))-
     length(idx)/nrow(D)*(G^2)/(2*(H+penalty))
  gain[!is.finite(gain)|Hl<=1e-9]<- -Inf
  if(all(!is.finite(gain)) || max(gain)<=0)
    return(candidates[0,,drop=FALSE])
  take<-order(-gain,seq_along(gain))
  used<-character();chosen<-integer()
  for(j in take){
    if(!is.finite(gain[j])||gain[j]<=0)break
    genes<-c(as.character(candidates$gene_left[j]),
             as.character(candidates$gene_right[j]))
    # Gene-disjoint within the local expert; overlap between alternative
    # tree branches is allowed, and no raw abundance condition is used.
    if(any(genes %in% used))next
    used<-c(used,genes)
    chosen<-c(chosen,j)
    if(length(chosen)>=K)break
  }
  z<-candidates[chosen,,drop=FALSE]
  if(length(chosen)>0){
    z$residual_score<-gain[chosen]
  }
  rownames(z)<-NULL
  z
}

crp_candidate_pool<-function(D,pool,max_pool=650L,max_extra=450L) {
  # Existing helper uses ONLY within-patient RANK variation for extra
  # candidate pairs; the original source gene/probe universe is frozen.
  rg_pairs(D,pool,max_pool=max_pool,max_extra=max_extra)
}

crp_fit<-function(D,y,p_train,p_oof,root,pool,
                    lambda=.10,allow_one_branch=TRUE){
  gate<-os_gate(aa_ranks(D),root)
  allpairs<-crp_candidate_pool(D,pool)
  Z<-encode_pair_matrix(D,allpairs,m=10L)
  children<-vector("list",2L);pairs<-vector("list",2L)
  for(b in 0:1) {
    ids<-which(gate==b)
    yy<-factor(as.character(y[ids]),levels=levels(y))
    chosen<-crp_select_pairs(D,y,p_oof,gate,as.integer(b),
                             allpairs,Z,K=8L,penalty=5)
    pairs[[b+1L]]<-if(nrow(chosen))chosen else NULL
    if(nrow(chosen)<2L){
      children[[b+1L]]<-list(pairs=NULL,fit=NULL,active=0L,
                             state="no_positive_residual_pairs")
      pairs[[b+1L]]<-NULL
      next
    }
    # Same fixed penalized local logistic architecture as previous SA-RLR.
    fit<-la_fit_child(D[ids,,drop=FALSE],yy,p_train[ids],
           pairs=chosen,select_pairs=FALSE,lambda=lambda)
    if(is.null(fit$fit)||fit$active==0L){
      fit$fit<-NULL;fit$active<-0L;pairs[[b+1L]]<-NULL
    }
    children[[b+1L]]<-fit
  }
  active<-sum(vapply(children,function(z)z$active,integer(1)))
  fitted<-sum(vapply(children,function(z)!is.null(z$fit),logical(1)))
  list(root=root,children=children,pairs=pairs,n_fitted=fitted,
       n_active=active,n_candidates=nrow(allpairs),
       n_candidate_branch_pairs=vapply(pairs,function(z)
          if(is.null(z))0L else nrow(z),integer(1)))
}

crp_eval_V<-function(D,yd,V,yv) {
  sourcepool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
  parent<-rg_parent(D,yd,V,sourcepool,K=50L)
  oof<-sa_parent_oof(D,yd,K=50L,seed=731L)
  roots<-sa_screen_roots(D,yd,oof,sourcepool,
                         max_top=6L,max_proposals=240L)
  yy<-as.integer(yv==levels(yv)[2L])
  listout<-list()
  for(i in seq_along(roots$roots)){
    root<-roots$roots[[i]]$root
    old<-sa_fit(D,yd,parent$train,root)
    new<-crp_fit(D,yd,parent$train,oof,root,sourcepool,lambda=.10)
    ps<-sa_predict(V,parent$test,old)
    pn<-sa_predict(V,parent$test,new)
    ng<-la_global_control(D,yd,parent$train,V,parent$test,new)
    listout[[length(listout)+1L]]<-list(
       root=root,rank_score=roots$scores[i],
       old=old,new=new,
       probs=list(ORIGINAL_RLR50=parent$test,
                  CLASS_TSP_CHILD_RLR=ps,
                  RESIDUAL_SCORED_CHILD_RLR=pn,
                  RESIDUAL_SELECTED_GLOBAL_NONTREE=ng),
       truth=yy)
  }
  list(roots=listout,n_proposed=roots$proposed,yy=yy,
       parent_V=parent$test)
}
