# D-only original RLR K50 OOF parent offset versus full-D in-sample offset,
# for EXACTLY THE SAME, ONCE-NOMINATED residual-selected conditional pairs.
# Entire experiment uses D/V only, never original paper heldout outer E.
source("R/experimental/rlr_conditional_residual_pairs_v1.R")

cfo_nominate <- function(D,yd,p_oof,root,pool) {
  gate<-os_gate(aa_ranks(D),root)
  candidates<-crp_candidate_pool(D,pool)
  matrixZ<-encode_pair_matrix(D,candidates,m=10L)
  lapply(0:1,function(b) {
    p<-crp_select_pairs(D,yd,p_oof,gate,as.integer(b),
                        candidates,matrixZ,K=8L,penalty=5)
    if(nrow(p)<2L)NULL else p
  })
}

cfo_fit_frozen <- function(D,yd,offset_D,root,nominations) {
  gate<-os_gate(aa_ranks(D),root)
  crp_need(length(nominations)==2L, "expected exactly two child panels")
  children<-vector("list",2L)
  pairs_active<-vector("list",2L)
  for(b in 0:1){
    idx<-which(gate==b)
    yy<-factor(as.character(yd[idx]),levels=levels(yd))
    panel<-nominations[[b+1L]]
    viable<-length(idx)>=12L && min(table(yy))>=3L
    if(is.null(panel)||!viable){
      children[[b+1L]]<-list(
         pairs=NULL,fit=NULL,active=0L,state="NO_NOMINATED_PAIRS_OR_SUPPORT")
      next
    }
    fitted<-la_fit_child(D[idx,,drop=FALSE],yy,offset_D[idx],
                         pairs=panel,select_pairs=FALSE,lambda=.10)
    if(!is.null(fitted$fit) && fitted$active>0L)
      pairs_active[[b+1L]]<-panel
    else {
      fitted$fit<-NULL
      fitted$active<-0L
    }
    children[[b+1L]]<-fitted
  }
  list(root=root,children=children,pairs=pairs_active,
       nominated_pairs=nominations,
       n_active=sum(vapply(children,function(x)x$active,integer(1))),
       n_fitted=sum(vapply(children,function(x)!is.null(x$fit),logical(1))))
}

cfo_run <- function(D,yd,V,yv){
  pool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
  parent<-rg_parent(D,yd,V,pool,K=50L)
  p_oof<-sa_parent_oof(D,yd,K=50L,seed=731L)
  roots<-sa_screen_roots(D,yd,p_oof,pool,
                        max_top=6L,max_proposals=240L)
  y01<-as.integer(yv==levels(yv)[2L])
  result<-list()
  for(i in seq_along(roots$roots)){
    root<-roots$roots[[i]]$root
    # Nominating once, prior to coefficient fitting, is required. One
    # offset can zero a coefficient while the other keeps it active:
    # that is a FINDING, not a change in supervised feature nomination.
    nominations<-cfo_nominate(D,yd,p_oof,root,pool)
    A<-cfo_fit_frozen(D,yd,parent$train,root,nominations)
    B<-cfo_fit_frozen(D,yd,p_oof,root,nominations)
    crp_need(identical(A$nominated_pairs,B$nominated_pairs),
              "same D-only nomination violated")
    pA<-sa_predict(V,parent$test,A)
    pB<-sa_predict(V,parent$test,B)
    # For no-tree controls use the SAME nominated union regardless of
    # whether the local fitted coefficient was shrunken exactly to zero.
    AA<-A
    BB<-B
    AA$pairs<-nominations
    BB$pairs<-nominations
    pGA<-la_global_control(D,yd,parent$train,V,parent$test,AA)
    pGB<-la_global_control(D,yd,p_oof,V,parent$test,BB)
    result[[i]]<-list(root=root,
        D_in_sample_parent_active=A$n_active,
        D_crossfit_parent_active=B$n_active,
        D_same_pairs=identical(A$nominated_pairs,B$nominated_pairs),
        D_nomination_count=sum(vapply(nominations,function(z)
                                  if(is.null(z))0L else nrow(z),integer(1))),
        prob=list(
          ORIGINAL_RLR50=parent$test,
          RESIDUAL_LOCAL_INSAMPLE_OFFSET=pA,
          RESIDUAL_LOCAL_CROSSFIT_OFFSET=pB,
          SAMEPAIRS_GLOBAL_INSAMPLE_OFFSET=pGA,
          SAMEPAIRS_GLOBAL_CROSSFIT_OFFSET=pGB),
        y=y01)
  }
  list(roots=result,proposed=roots$proposed,
       n_D=nrow(D),n_V=nrow(V))
}
