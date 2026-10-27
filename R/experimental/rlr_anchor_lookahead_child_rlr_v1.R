# Anchor Look-Ahead RLR Model Stump V1.
# A TRUE split tree: one ordinal anchor root, two fitted branch-specific
# relational logistic RLR corrections with SHARED original RLR offset.
# No raw-value thresholds / expression ratios. No outer E outcomes in search.
source("R/experimental/rlr_ordinal_anchor_score_tree_v1.R")

la_need <- rg_need
la_truth <- os_truth
la_logit <- function(p)qlogis(pmin(1-1e-10,pmax(1e-10,p)))
la_eps <- 1e-10
la_lambda <- 0.10

la_parent_oof <- function(D,y,folds=3L,seed=712L) {
  plans<-make_repeated_stratified_folds(y,n_folds=folds,
                                       n_repeats=1L,seed=seed)
  result<-rep(NA_real_,nrow(D))
  for(i in seq_along(plans)){
    j<-plans[[i]]
    dx<-D[j$train_idx,,drop=FALSE]
    vy<-D[j$test_idx,,drop=FALSE]
    yy<-factor(as.character(y[j$train_idx]),levels=levels(y))
    pp<-compute_top_pairs_switchbox(dx,yy,top_n=1000L,
                                    buffer_mult=5L)
    original<-rg_parent(dx,yy,vy,pp,K=20L)
    result[j$test_idx]<-original$test
  }
  la_need(all(is.finite(result)) && all(result>0&result<1),
          "3-fold D-only parent out-of-fold probabilities missing")
  result
}

la_roots <- function(D,y,pool,parent_oof,max_hubs=12L,
                     top_budget=6L,min_branch=12L,min_class=3L){
  R<-aa_ranks(D)
  candidates<-list()
  for(target in aa_hubs(pool,R,max_hubs=max_hubs)){
    for(k in c(3L,5L)){
      anchors<-tryCatch(aa_anchors(R,y,pool,target,k),
                        error=function(e)NULL)
      if(is.null(anchors))next
      for(threshold in seq_len(k)){
        node<-os_node(target,anchors,threshold)
        g<-os_gate(R,node)
        idx0<-which(g==0L); idx1<-which(g==1L)
        if(min(length(idx0),length(idx1))<min_branch)next
        truth<-la_truth(y)
        if(min(tabulate(truth[idx0]+1L,nbins=2),
               tabulate(truth[idx1]+1L,nbins=2))<min_class)next
        candidates[[length(candidates)+1L]]<-list(
           node=node,group=g,counts=c(length(idx0),length(idx1)))
      }
    }
  }
  if(length(candidates)==0L)return(list(candidates=list(),
                                         cheap_scores=numeric(),
                                         n_proposed=0L))
  # Cheap OOF residual score for different branch-specific RLR coefficients.
  # From a fixed D-only pair pool: no locally fitted conditional experts yet.
  probe<-pool[seq_len(min(160L,nrow(pool))),,drop=FALSE]
  Z<-encode_pair_matrix(D,probe,m=10L)
  truth<-la_truth(y)
  error<-truth-parent_oof
  weight<-pmax(1e-5,parent_oof*(1-parent_oof))
  z2<-Z*Z
  G<-colSums(Z * error)
  H<-colSums(z2 * weight)
  ordinary<-G^2/(2*(H+1))
  screen<-function(rows){
    gg<-colSums(Z[rows,,drop=FALSE]*error[rows])
    hh<-colSums(z2[rows,,drop=FALSE]*weight[rows])
    gg^2/(2*(hh+1))
  }
  scores<-vapply(candidates,function(c){
    il<-which(c$group==0L);ir<-which(c$group==1L)
    contextual<-screen(il)+screen(ir)-ordinary
    # Cut off negative due to ridge denominator approximations.
    max(0,max(contextual,na.rm=TRUE)) -
           .0008*(length(c$node$anchors)-1L)
  },numeric(1))
  pick<-order(-scores,
     vapply(candidates,function(c)c$node$target,character(1)),
     vapply(candidates,function(c)c$node$threshold,integer(1)))
  chosen<-head(pick,top_budget)
  list(candidates=candidates[chosen],cheap_scores=scores[chosen],
       n_proposed=length(candidates))
}

la_fit_child <- function(D,y,base_prob,pairs=NULL,
                         select_pairs=TRUE,lambda=la_lambda){
  if(length(y)<12L||min(table(y))<3L)
    return(list(pairs=NULL,fit=NULL,active=0L,state="branch_support"))
  if(select_pairs){
    pool<-tryCatch(compute_top_pairs_switchbox(
      D,y,top_n=600L,buffer_mult=3L),error=function(e)NULL)
    if(is.null(pool))
      return(list(pairs=NULL,fit=NULL,active=0L,state="no_local_pool"))
    pairs<-apply_disjointness(pool,K=8L,d=1L)
  }
  if(is.null(pairs)||nrow(pairs)<2L)
    return(list(pairs=NULL,fit=NULL,active=0L,state="no_local_pairs"))
  R<-encode_pair_matrix(D,pairs,m=10L)
  if(sum(apply(R,2L,stats::var)>0)<1L)
    return(list(pairs=pairs,fit=NULL,active=0L,state="constant_pairs"))
  fit<-tryCatch(
    suppressWarnings(glmnet::glmnet(x=R,y=la_truth(y),
      family="binomial",alpha=.5,lambda=lambda,
      offset=la_logit(base_prob),intercept=FALSE,
      standardize=FALSE)),
    error=function(e)NULL)
  if(is.null(fit))
    return(list(pairs=pairs,fit=NULL,active=0L,state="glmnet_failure"))
  beta<-as.numeric(stats::coef(fit,s=lambda))[-1]
  list(pairs=pairs,fit=fit,active=sum(abs(beta)>1e-9),state="fitted")
}

la_child_predict <- function(X,parent_prob,fit_obj,lambda=la_lambda) {
  if(is.null(fit_obj$fit))return(parent_prob)
  R<-encode_pair_matrix(X,fit_obj$pairs,m=10L)
  p<-tryCatch(as.numeric(predict(fit_obj$fit,newx=R,
       newoffset=la_logit(parent_prob),s=lambda,type="response")),
       error=function(e)NULL)
  if(is.null(p)||any(!is.finite(p)))return(parent_prob)
  p
}

la_fit_structure <- function(D,y,base_D,root,pairs_fixed=NULL){
  R<-aa_ranks(D); groups<-os_gate(R,root)
  children<-vector("list",2L)
  for(b in 0:1){
    idx<-which(groups==b)
    yy<-factor(as.character(y[idx]),levels=levels(y))
    pp<-if(is.null(pairs_fixed))NULL else pairs_fixed[[b+1L]]
    children[[b+1L]]<-la_fit_child(D[idx,,drop=FALSE],yy,
          base_D[idx],pairs=pp,select_pairs=is.null(pairs_fixed))
  }
  list(root=root,children=children,
    active=sum(vapply(children,function(z)z$active,integer(1))),
    working=sum(vapply(children,function(z)!is.null(z$fit),logical(1))),
    pairs=lapply(children,function(z)z$pairs))
}

la_predict_tree <- function(E,base_E,model){
  g<-os_gate(aa_ranks(E),model$root)
  prediction<-base_E
  for(b in 0:1) {
    idx<-which(g==b)
    if(length(idx)==0L)next
    prediction[idx]<-la_child_predict(E[idx,,drop=FALSE],
         base_E[idx],model$children[[b+1L]])
  }
  as.numeric(prediction)
}

# Same branch-specific pair *union* but no branch-specific coefficients.
# Includes ordinal gate itself to avoid misleading "tree vs no information".
la_global_control <- function(D,y,base_D,V,base_V,model){
  local<-Filter(Negate(is.null),model$pairs)
  if(length(local)==0L)return(base_V)
  union<-do.call(rbind,local)
  key<-paste(union$gene_left,union$gene_right,sep="|")
  union<-union[!duplicated(key),,drop=FALSE]
  RD<-encode_pair_matrix(D,union,m=10L)
  RV<-encode_pair_matrix(V,union,m=10L)
  gate_D<-os_gate(aa_ranks(D),model$root)
  gate_V<-os_gate(aa_ranks(V),model$root)
  x<-cbind(root_gate=gate_D,RD)
  xx<-cbind(root_gate=gate_V,RV)
  if(ncol(x)<2L)x<-cbind(x,inert=0); 
  if(ncol(xx)<ncol(x))xx<-cbind(xx,inert=0)
  if(!any(apply(x,2,stats::var)>0))return(base_V)
  fit<-tryCatch(suppressWarnings(glmnet::glmnet(x=x,
        y=la_truth(y),family="binomial",alpha=.5,lambda=la_lambda,
        intercept=FALSE,standardize=FALSE,
        offset=la_logit(base_D))),error=function(e)NULL)
  if(is.null(fit))return(base_V)
  ans<-tryCatch(as.numeric(predict(fit,newx=xx,s=la_lambda,
         type="response",newoffset=la_logit(base_V))),
         error=function(e)NULL)
  if(is.null(ans)||any(!is.finite(ans)))base_V else ans
}

la_gate_control <- function(D,y,base_D,V,base_V,root){
  gd<-os_gate(aa_ranks(D),root)
  gv<-os_gate(aa_ranks(V),root)
  design<-cbind(gd,inert=rep(0,nrow(D)))
  vdesign<-cbind(gv,inert=rep(0,nrow(V)))
  fit<-tryCatch(suppressWarnings(glmnet::glmnet(
       x=design,y=la_truth(y),family="binomial",alpha=0,
       lambda=la_lambda,intercept=FALSE,standardize=FALSE,
       offset=la_logit(base_D))),error=function(e)NULL)
  if(is.null(fit))return(base_V)
  tryCatch(as.numeric(predict(fit,newx=vdesign,s=la_lambda,
       type="response",newoffset=la_logit(base_V))),
       error=function(e)base_V)
}

la_val_candidate <- function(D,y,base_D,V,base_V,root) {
  model<-la_fit_structure(D,y,base_D,root)
  tree<-la_predict_tree(V,base_V,model)
  global<-la_global_control(D,y,base_D,V,base_V,model)
  gate<-la_gate_control(D,y,base_D,V,base_V,root)
  list(root=root,structure=model,p_tree=tree,p_global=global,
       p_gate=gate,working=model$working,active=model$active)
}

# Inner D contains all proposal and local-pair selection. Inner V chooses
# from at most SIX root designs and independently evaluates equal-feature
# non-tree and gate-only controls.
la_select <- function(T,yT){
  fold<-rg_inner(yT)
  D<-T[fold$D,,drop=FALSE];V<-T[fold$V,,drop=FALSE]
  yd<-factor(as.character(yT[fold$D]),levels=levels(yT))
  yv<-factor(as.character(yT[fold$V]),levels=levels(yT))
  pool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
  parent<-rg_parent(D,yd,V,pool,K=20L)
  yy<-la_truth(yv)
  oof<-la_parent_oof(D,yd,folds=3L,seed=712L)
  roots<-la_roots(D,yd,pool,oof,max_hubs=12L,top_budget=6L)
  greedy<-os_choose_node(aa_ranks(D),yd,pool,
           leaf_min=12L,each_class_min=3L,max_targets=12L)
  fallback<-list(accepted=FALSE,model=NULL,root=NULL,pairs=NULL,
       v_gain_parent=NA_real_,v_gain_global=NA_real_,
       v_gain_gate=NA_real_,v_gain_greedy=NA_real_,
       v_ba_delta=NA_real_,working=0L,active=0L,
       proposed=roots$n_proposed,considered=length(roots$candidates),
       greedy_root=if(is.null(greedy))NULL else greedy$node,
       reason="no_valid_tree")
  cand<-lapply(roots$candidates,function(c){
     la_val_candidate(D,yd,parent$train,V,parent$test,c$node)
  })
  cand<-Filter(function(c)c$working==2L && c$active>0L,cand)
  greedyV<-if(!is.null(greedy))
     la_val_candidate(D,yd,parent$train,V,parent$test,greedy$node)
     else NULL
  if(length(cand)==0L)return(fallback)
  losses<-vapply(cand,function(c)rg_ll(yy,c$p_tree),numeric(1))
  best<-cand[[which.min(losses)]]
  g0<-rg_ll(yy,parent$test)-rg_ll(yy,best$p_tree)
  g1<-rg_ll(yy,best$p_global)-rg_ll(yy,best$p_tree)
  g2<-rg_ll(yy,best$p_gate)-rg_ll(yy,best$p_tree)
  gg<-if(is.null(greedyV))NA_real_ else
      rg_ll(yy,greedyV$p_tree)-rg_ll(yy,best$p_tree)
  bad<-rg_ba(yy,best$p_tree)-rg_ba(yy,parent$test)
  accept<-g0>=.008 && g1>=.004 && g2>=.004 &&
          bad>=-.02
  list(accepted=accept,model=best$structure,root=best$root,
       pairs=best$structure$pairs,
       v_gain_parent=g0,v_gain_global=g1,v_gain_gate=g2,
       v_gain_greedy=gg,v_ba_delta=bad,
       working=best$working,active=best$active,
       proposed=roots$n_proposed,considered=length(roots$candidates),
       greedy_root=if(is.null(greedy))NULL else greedy$node,
       reason="V_GATE")
}

la_evaluate <- function(T,yT,E,yE) {
  sel<-la_select(T,yT)
  pool<-compute_top_pairs_switchbox(T,yT,top_n=1000L,buffer_mult=5L)
  parent<-rg_parent(T,yT,E,pool,K=20L)
  original<-parent$test
  p_tree<-p_global<-p_gate<-p_greedy<-original
  actual_active<-0L
  if(!is.null(sel$root)){
    # Freeze structure and D-only nominated pairs. Refit coefficients on T.
    fitted<-la_fit_structure(T,yT,parent$train,sel$root,
                              pairs_fixed=sel$pairs)
    p_tree<-la_predict_tree(E,original,fitted)
    p_global<-la_global_control(T,yT,parent$train,E,original,fitted)
    p_gate<-la_gate_control(T,yT,parent$train,E,original,sel$root)
    actual_active<-fitted$active
  }
  # Greedy-root comparator uses the identical child-RLR fitter; it has an
  # independent D-only root, and new T-only local pair nominations.
  if(!is.null(sel$greedy_root)){
    greedyfit<-la_fit_structure(T,yT,parent$train,sel$greedy_root)
    p_greedy<-la_predict_tree(E,original,greedyfit)
  }
  gated<-if(isTRUE(sel$accepted))p_tree else original
  list(probs=list(ORIGINAL_RLR=original,
       GREEDY_ROOT_CHILD_RLR=p_greedy,
       LOOKAHEAD_CHILD_RLR=p_tree,
       LOOKAHEAD_SAME_PAIRS_GLOBAL=p_global,
       LOOKAHEAD_GATE_ONLY=p_gate,
       V_GATED_LOOKAHEAD_CHILD_RLR=gated),
       selection=sel,truth=la_truth(yE),realized_active=actual_active)
}
