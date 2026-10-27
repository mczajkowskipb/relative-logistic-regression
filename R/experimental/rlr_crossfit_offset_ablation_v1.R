# D-only repeated original RLR OOF parent offset vs full-D in-sample offset
# for identical residual-selected conditional child gene-pair candidates.
# Scientific ablation: no modification to heldout outer E or original RLR50.
source("R/experimental/rlr_conditional_residual_pairs_v1.R")

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
    # crp_fit computes exactly the SAME residual-selected pairs for both
    # models; only the offset used to fit their local coefficients differs.
    A<-crp_fit(D,yd,parent$train,p_oof,root,pool,lambda=.10)
    B<-crp_fit(D,yd,p_oof,p_oof,root,pool,lambda=.10)
    crp_need(identical(A$pairs,B$pairs),
             "crossfit offset changed pair nomination or ordering")
    pA<-sa_predict(V,parent$test,A)
    pB<-sa_predict(V,parent$test,B)
    pGA<-la_global_control(D,yd,parent$train,V,parent$test,A)
    pGB<-la_global_control(D,yd,p_oof,V,parent$test,B)
    result[[i]]<-list(root=root,
        D_in_sample_parent_active=A$n_active,
        D_crossfit_parent_active=B$n_active,
        D_same_pairs=identical(A$pairs,B$pairs),
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
