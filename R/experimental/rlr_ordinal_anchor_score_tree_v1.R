# Ordinal Anchor-Score Decision Tree (OAS-DT) — V1.
# Each split uses COUNT of within-patient target>reference comparisons, not
# raw expression, expression ratios, or cohort-specific value thresholds.
# Original RLR is sourced unmodified from R/new/07_world_protocol_utils.R.
source("R/experimental/rlr_anchor_majority_tree_v1.R")

os_need <- rg_need
os_truth <- aa_y

os_score <- function(R, node) {
  os_need(node$target %in% colnames(R) &&
          all(node$anchors %in% colnames(R)) &&
          !node$target %in% node$anchors &&
          !anyDuplicated(node$anchors) &&
          length(node$anchors) %in% c(3L,5L),
          "ordinal anchor node has invalid gene schema")
  # Exact per-patient integer count from 0 to k.
  comparisons <- vapply(node$anchors, function(a)
       as.integer(R[,node$target] > R[,a]), integer(nrow(R)))
  if(is.null(dim(comparisons)))comparisons <- matrix(comparisons,ncol=1L)
  as.integer(rowSums(comparisons))
}

os_gate <- function(R,node) as.integer(os_score(R,node)>=node$threshold)

os_node <- function(target,anchors,threshold) {
  stopifnot(length(anchors) %in% c(3L,5L),
            threshold %in% seq_len(length(anchors)))
  list(target=as.character(target),anchors=as.character(anchors),
       threshold=as.integer(threshold))
}

os_prior <- function(y01, fallback=.5, strength=8) {
  if(length(y01)==0L) return(fallback)
  (sum(y01)+strength*fallback)/(length(y01)+strength)
}

os_posterior <- function(y01, group, prior, strength=8) {
  p0<-os_prior(y01[group==0L],prior,strength)
  p1<-os_prior(y01[group==1L],prior,strength)
  c(p0,p1)
}

os_split_quality <- function(R,y,node,leaf_min=8L,
                             each_class_min=1L,prior_strength=8) {
  z<-os_gate(R,node)
  idx0<-which(z==0L);idx1<-which(z==1L)
  if(min(length(idx0),length(idx1))<leaf_min) return(NULL)
  yy<-os_truth(y)
  if(each_class_min>0L) {
    if(min(tabulate(yy[idx0]+1L,2L),
           tabulate(yy[idx1]+1L,2L))<each_class_min)return(NULL)
  }
  prior<-mean(yy)
  prob<-os_posterior(yy,z,prior,prior_strength)
  loss<-rg_ll(yy,prob[z+1L])
  base<-rg_ll(yy,rep(os_prior(yy,prior,prior_strength),length(yy)))
  list(loss=loss,gain=base-loss,branch_counts=c(length(idx0),length(idx1)),
       branch_prob=prob)
}

# Node optimization uses a *frozen* list of target candidates and one fixed
# training-only star-panel formation; only threshold t in 1..k is searched.
# Complexity penalty is per anchor comparison, with no test-E access.
os_choose_node <- function(R,y,pool,leaf_min=8L,each_class_min=1L,
                           max_targets=14L) {
  targets<-aa_hubs(pool,R,max_hubs=max_targets)
  opts<-list()
  for(target in targets)for(k in c(3L,5L)) {
    anchors<-tryCatch(aa_anchors(R,y,pool,target,k),
                        error=function(e)NULL)
    if(is.null(anchors))next
    for(t in seq_len(k)){
      node<-os_node(target,anchors,t)
      q<-os_split_quality(R,y,node,leaf_min,each_class_min)
      if(is.null(q))next
      # Parsimonious node induction, penalizing larger anchor panels.
      objective<-q$loss+0.0025*(k-1L)
      opts[[length(opts)+1L]]<-list(node=node,quality=q,
                                   objective=objective)
    }
  }
  if(length(opts)==0L)return(NULL)
  win<-which.min(vapply(opts,function(o)o$objective,numeric(1)))
  opts[[win]]
}

os_build_tree <- function(D,yD,root_pool,depth=2L) {
  os_need(depth %in% 1:2,"supported depth is one or two")
  ranks<-aa_ranks(D)
  root<-os_choose_node(ranks,yD,root_pool,leaf_min=10L,
                         each_class_min=2L,max_targets=14L)
  if(is.null(root))return(NULL)
  children<-list(NULL,NULL)
  diagnostics<-list()
  if(depth==2L){
    branches<-os_gate(ranks,root$node)
    for(branch in 0:1){
      idx<-which(branches==branch)
      yy<-factor(as.character(yD[idx]),levels=levels(yD))
      if(length(idx)<22L||min(table(yy))<5L) {
        diagnostics[[branch+1L]]<-list(state="insufficient_node_support",
                                        size=length(idx))
        next
      }
      local_pool<-tryCatch(compute_top_pairs_switchbox(
        D[idx,,drop=FALSE],yy,top_n=600L,buffer_mult=4L),
        error=function(e)NULL)
      if(is.null(local_pool)){
        diagnostics[[branch+1L]]<-list(state="no_local_candidates")
        next
      }
      candidate<-os_choose_node(ranks[idx,,drop=FALSE],yy,local_pool,
                  leaf_min=6L,each_class_min=1L,max_targets=10L)
      if(is.null(candidate)){
        diagnostics[[branch+1L]]<-list(state="no_valid_local_split")
        next
      }
      # Prespecified lower D gain than old binary anchor tree's .025.
      # D only; presence is NOT a V/E success claim.
      if(candidate$quality$gain < .004){
        diagnostics[[branch+1L]]<-list(state="insufficient_D_gain",
                        gain=candidate$quality$gain)
        next
      }
      children[[branch+1L]]<-candidate$node
      diagnostics[[branch+1L]]<-list(state="child_proposed_D_only",
                         gain=candidate$quality$gain,size=length(idx))
    }
  }
  list(root=root$node,children=children,diagnostics=diagnostics,
       root_D_gain=root$quality$gain,
       root_D_counts=root$quality$branch_counts)
}

# Leaf routes 1=left unsplit; 2=right unsplit; 3,4 left children;
# 5,6 right children; fixed semantics across train and predict.
os_route <- function(R,tree,depth=2L) {
  b<-os_gate(R,tree$root)
  ids<-b+1L
  if(depth==2L)for(branch in 0:1){
    child<-tree$children[[branch+1L]]
    if(is.null(child))next
    ii<-which(b==branch)
    if(length(ii)==0L)next
    c<-os_gate(R[ii,,drop=FALSE],child)
    ids[ii]<-3L+2L*branch+c
  }
  as.integer(ids)
}

os_fit_leaf <- function(R,y,tree,depth=2L,prior_strength=8) {
  groups<-os_route(R,tree,depth)
  truth<-os_truth(y)
  prior<-mean(truth)
  probs<-vapply(1:6,function(g)
   os_prior(truth[groups==g],fallback=prior,
            strength=prior_strength),numeric(1))
  list(tree=tree,depth=depth,probs=probs)
}

os_predict_leaf <- function(model,R) {
  ids<-os_route(R,model$tree,model$depth)
  as.numeric(model$probs[ids])
}

# A single selected anchor compared with the same target, a proper
# positive-probability one-pair comparator. The pair is selected by TRAIN
# logloss and then heldout E is untouched.
os_single_fit <- function(R,y,target,anchors,strength=8) {
  truth<-os_truth(y)
  prior<-mean(truth)
  opt<-list(loss=Inf)
  for(a in anchors) {
    g<-as.integer(R[,target]>R[,a])
    q<-os_posterior(truth,g,prior,strength)
    loss<-rg_ll(truth,q[g+1L])
    if(loss<opt$loss)opt<-list(target=target,anchor=a,prob=q,loss=loss)
  }
  opt
}
os_single_predict <- function(R,model) {
  g<-as.integer(R[,model$target]>R[,model$anchor])
  as.numeric(model$prob[g+1L])
}

# A non-tree logistic comparator sees exactly the ordinal count score of
# every selected node globally; unlike a tree, it has no split-dependent
# coefficients/decision-path routing.
os_count_design <- function(R,tree,depth=2L) {
  nodes<-list(tree$root)
  if(depth==2L)nodes<-c(nodes,Filter(Negate(is.null),tree$children))
  ans<-do.call(cbind,lapply(nodes,function(nd)os_score(R,nd)))
  if(is.null(dim(ans))) ans<-matrix(ans,ncol=1L)
  colnames(ans)<-paste0("ordinal_count_",seq_len(ncol(ans)))
  ans
}
os_global_fit <- function(R,y,tree,depth=2L) {
  x<-os_count_design(R,tree,depth)
  fit<-glmnet::glmnet(x=x,y=os_truth(y),family="binomial",
      alpha=0,lambda=.10,standardize=FALSE)
  list(fit=fit,tree=tree,depth=depth)
}
os_global_predict <- function(R,model){
  as.numeric(predict(model$fit,
    newx=os_count_design(R,model$tree,model$depth),
    s=.10,type="response"))
}

os_inner <- function(T,yT) {
  halves<-rg_inner(yT)
  D<-T[halves$D,,drop=FALSE]; V<-T[halves$V,,drop=FALSE]
  yd<-factor(as.character(yT[halves$D]),levels=levels(yT))
  yv<-factor(as.character(yT[halves$V]),levels=levels(yT))
  pool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
  parent<-rg_parent(D,yd,V,pool,K=20L)
  tree<-os_build_tree(D,yd,pool,depth=2L)
  fallback<-list(tree=NULL,mode="NONE",accepted=FALSE,
       v_gain_rlr=NA_real_,v_gain_single=NA_real_,
       v_gain_global=NA_real_,v_BA_diff=NA_real_,
       proposed_children=0L,root_gene="",k=0L,
       threshold=NA_integer_,anchor_names="",reason="no_root")
  if(is.null(tree))return(fallback)
  rD<-aa_ranks(D);rV<-aa_ranks(V)
  ps<-os_predict_leaf(os_fit_leaf(rD,yd,tree,1L),rV)
  pd<-os_predict_leaf(os_fit_leaf(rD,yd,tree,2L),rV)
  single<-os_single_fit(rD,yd,tree$root$target,tree$root$anchors)
  p_single<-os_single_predict(rV,single)
  ctrl<-os_global_fit(rD,yd,tree,depth=2L)
  p_global<-os_global_predict(rV,ctrl)
  yy<-os_truth(yv)
  # At least 0.005 improvement on V vs stump for true extra depth.
  mode<-if(any(!vapply(tree$children,is.null,logical(1))) &&
             rg_ll(yy,pd) <= rg_ll(yy,ps)-.005) "DEPTH2" else "STUMP"
  predicted<-if(mode=="DEPTH2")pd else ps
  gain0<-rg_ll(yy,parent$test)-rg_ll(yy,predicted)
  gain1<-rg_ll(yy,p_single)-rg_ll(yy,predicted)
  gain2<-rg_ll(yy,p_global)-rg_ll(yy,predicted)
  ba_diff<-rg_ba(yy,predicted)-rg_ba(yy,parent$test)
  # Both baseline and single-reference wins must be real in V,
  # otherwise report diagnostic tree but return original RLR to E.
  accepted<-gain0>=.005 && gain1>=.002 && ba_diff>=-.02
  list(tree=tree,mode=mode,accepted=accepted,
       v_gain_rlr=gain0,v_gain_single=gain1,
       v_gain_global=gain2,v_BA_diff=ba_diff,
       proposed_children=sum(!vapply(tree$children,is.null,logical(1))),
       root_gene=tree$root$target,k=length(tree$root$anchors),
       threshold=tree$root$threshold,
       anchor_names=paste(tree$root$anchors,collapse="|"),
       reason="fixed_D_V_gate")
}

os_evaluate <- function(T,yT,E,yE) {
  choice<-os_inner(T,yT)
  poolT<-compute_top_pairs_switchbox(T,yT,top_n=1000L,buffer_mult=5L)
  parent<-rg_parent(T,yT,E,poolT,K=20L)
  original<-parent$test
  stump<-deep<-single<-global<-original
  realized_children<-0L
  if(!is.null(choice$tree)){
    train_r<-aa_ranks(T);test_r<-aa_ranks(E)
    stump<-os_predict_leaf(os_fit_leaf(train_r,yT,choice$tree,1L),
                           test_r)
    deep<-os_predict_leaf(os_fit_leaf(train_r,yT,choice$tree,2L),
                          test_r)
    one<-os_single_fit(train_r,yT,choice$tree$root$target,
                      choice$tree$root$anchors)
    single<-os_single_predict(test_r,one)
    global<-os_global_predict(test_r,
               os_global_fit(train_r,yT,choice$tree,depth=2L))
    realized_children<-choice$proposed_children
  }
  selected<-if(choice$mode=="DEPTH2")deep else
            if(choice$mode=="STUMP")stump else original
  gated<-if(isTRUE(choice$accepted))selected else original
  list(probs=list(ORIGINAL_RLR=original,ANCHOR_SCORE_STUMP=stump,
        ANCHOR_SCORE_TREE_D2=deep,
        SAME_ROOT_BEST_SINGLE_PAIR=single,
        SAME_COUNT_RIDGE=global,
        V_SELECTED_ANCHOR_SCORE_TREE=selected,
        RLR_V_GATED_ANCHOR_SCORE_TREE=gated),
    selection=choice,truth=os_truth(yE),children=realized_children)
}
