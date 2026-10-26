# Anchor-majority decision stump and optional depth-2 top-down tree.
# Prediction uses only ordinal comparisons within a single patient.
# The frozen original RLR remains unchanged.
source("R/experimental/rlr_rank_gate_core_v1.R")

aa_need <- rg_need
aa_y <- function(y) as.integer(y==levels(y)[2L])
aa_ranks <- rg_ranks

aa_vote <- function(R, target, anchors) {
  aa_need(length(anchors)%%2L==1L && length(anchors)>=1L,
          "an odd, nonempty number of anchors required")
  aa_need(target %in% colnames(R) &&
            all(anchors %in% colnames(R)) &&
            !target %in% anchors && !anyDuplicated(anchors),
          "anchor/target schema invalid")
  comp <- vapply(anchors, function(a)
    as.integer(R[,target] > R[,a]), integer(nrow(R)))
  if(is.null(dim(comp))) comp <- matrix(comp,ncol=1L)
  as.integer(rowSums(comp) > length(anchors)/2)
}

aa_node <- function(target, anchors) {
  list(target=target,anchors=anchors)
}
aa_call <- function(R,node) aa_vote(R,node$target,node$anchors)
aa_leaf_p <- function(y,alpha=1.5) {
  (sum(y)+alpha)/(length(y)+2*alpha)
}
aa_probvec <- function(y,indices) aa_leaf_p(y[indices])

aa_rank_stats <- function(R,y) {
  yp<-aa_y(y)==1L; yn<-!yp
  p<-ncol(R)
  delta<-abs(colMeans(R[yp,,drop=FALSE])-colMeans(R[yn,,drop=FALSE]))/p
  spread<-apply(R,2L,stats::mad)/p
  list(shift=delta,spread=spread)
}

# Pair hub frequency from original supervised top1000; do not enforce d2 here:
# a STAR intentionally reuses its central gene.
aa_hubs <- function(pool, R, max_hubs=14L) {
  hits<-table(c(as.character(pool$gene_left),
                as.character(pool$gene_right)))
  hits<-sort(hits,decreasing=TRUE)
  hubs<-intersect(names(hits),colnames(R))
  unique(head(hubs,max_hubs))
}

aa_anchors <- function(R,y,pool,target,k) {
  p<-ncol(R)
  stats<-aa_rank_stats(R,y)
  ispos<-aa_y(y)==1L
  near<-unique(c(as.character(pool$gene_left[pool$gene_right==target]),
                 as.character(pool$gene_right[pool$gene_left==target])))
  near<-intersect(near,colnames(R))
  # A second, unsupervised train-only candidate source: reference genes with
  # low rank shift and reasonable rank proximity to the target's class midpoint.
  tmid<-mean(c(median(R[ispos,target]),median(R[!ispos,target])))
  distance<-abs(colMeans(R)-tmid)/p
  anchor_rank<-order(stats$shift+.2*stats$spread+.2*distance,
                     names(stats$shift))
  extra<-colnames(R)[head(anchor_rank,70L)]
  candidates<-setdiff(unique(c(near,extra)),target)
  aa_need(length(candidates)>=k,"insufficient anchor genes")
  u<-vapply(candidates,function(a)
     as.integer(R[,target]>R[,a]),integer(nrow(R)))
  if(is.null(dim(u)))u<-matrix(u,ncol=1L)
  signal<-abs(colMeans(u[ispos,,drop=FALSE])-
              colMeans(u[!ispos,,drop=FALSE]))
  # Stable anchor is defined only in the rank representation;
  # 'stable expression' in raw units is NEVER assumed.
  anchor_score<-signal-0.30*stats$shift[candidates]-
                       0.10*stats$spread[candidates]-
                       0.05*distance[candidates]
  orderidx<-order(-anchor_score,candidates)
  chosen<-character()
  for(ii in orderidx){
    cand<-candidates[ii]
    if(length(chosen)==0L){
      chosen<-cand
    } else {
      # Avoid selecting duplicate comparison patterns in the D pool.
      corr<-max(vapply(chosen,function(g){
        a<-as.integer(R[,target]>R[,g])
        b<-as.integer(R[,target]>R[,cand])
        mean(a==b)
      },numeric(1)))
      if(corr<.96)chosen<-c(chosen,cand)
    }
    if(length(chosen)>=k)break
  }
  # If insufficient diversity, fill deterministically using best scores.
  if(length(chosen)<k)chosen<-c(chosen,
    head(setdiff(candidates[orderidx],chosen),k-length(chosen)))
  aa_need(length(chosen)==k&&!anyDuplicated(chosen),"anchor selection invalid")
  chosen
}

aa_node_score <- function(R,y,node,min_leaf=8L) {
  gate<-aa_call(R,node)
  counts<-tabulate(gate+1L,nbins=2L)
  if(min(counts)<min_leaf)return(Inf)
  truth<-aa_y(y)
  q<-c(aa_leaf_p(truth[gate==0]),aa_leaf_p(truth[gate==1]))
  rg_ll(truth,q[gate+1L])
}

aa_find_node <- function(R,y,pool,min_leaf=8L,
                          max_hubs=14L,k_set=c(3L,5L)) {
  targets<-aa_hubs(pool,R,max_hubs)
  records<-list()
  for(t in targets)for(k in k_set){
    a<-tryCatch(aa_anchors(R,y,pool,t,k),
                error=function(e)NULL)
    if(is.null(a))next
    node<-aa_node(t,a)
    loss<-aa_node_score(R,y,node,min_leaf=min_leaf)
    if(!is.finite(loss))next
    records[[length(records)+1L]]<-list(node=node,loss=loss)
  }
  if(length(records)==0L)return(NULL)
  scores<-vapply(records,function(x)x$loss,numeric(1))
  records[[which.min(scores)]]
}

# True top-down induction: the root is chosen on D, each eligible child
# is then optimized on its own D subset with new supervised pair ranking.
aa_train_tree <- function(X,y,pool,depth=2L) {
  R<-aa_ranks(X)
  root<-aa_find_node(R,y,pool,min_leaf=8L,max_hubs=14L)
  if(is.null(root))return(NULL)
  children<-vector("list",2L)
  gate<-aa_call(R,root$node)
  if(depth>=2L)for(b in 0:1){
    idx<-which(gate==b)
    yy<-factor(as.character(y[idx]),levels=levels(y))
    if(length(idx)<24L||min(table(yy))<5L)next
    XX<-X[idx,,drop=FALSE]
    local_pool<-tryCatch(compute_top_pairs_switchbox(
        XX,yy,top_n=600L,buffer_mult=4L),
        error=function(e)NULL)
    if(is.null(local_pool))next
    candidate<-aa_find_node(R[idx,,drop=FALSE],yy,
           local_pool,min_leaf=7L,max_hubs=10L,k_set=c(3L,5L))
    if(is.null(candidate))next
    # Modest D-only improvement must exceed fixed cost of new split.
    part<-aa_call(R[idx,,drop=FALSE],candidate$node)
    if(min(tabulate(part+1L,2L))<7L)next
    tr<-aa_y(yy)
    parentloss<-rg_ll(tr,rep(aa_leaf_p(tr),length(tr)))
    childloss<-rg_ll(tr,c(aa_leaf_p(tr[part==0]),
                          aa_leaf_p(tr[part==1]))[part+1L])
    if(parentloss-childloss>.025)
      children[[b+1L]]<-candidate$node
  }
  list(root=root$node,children=children)
}

aa_route <- function(R,tree,depth=2L) {
  g<-aa_call(R,tree$root)
  out<-g+1L
  if(depth>=2L)for(b in 0:1)if(!is.null(tree$children[[b+1L]])){
    ids<-which(g==b)
    if(length(ids)==0L)next
    local<-aa_call(R[ids,,drop=FALSE],tree$children[[b+1L]])
    out[ids]<-(b+1L)*2L+local+1L
  }
  out
}
aa_train_leaf <- function(R,y,tree,depth=2L) {
  route<-aa_route(R,tree,depth)
  yy<-aa_y(y)
  # For a missing child on a route, use parent probabilities.
  lookup<-vapply(1:6,function(i)
      aa_leaf_p(yy[route==i]),numeric(1))
  list(tree=tree,lookup=lookup,depth=depth)
}
aa_predict_leaf <- function(model,R) {
  route<-aa_route(R,model$tree,model$depth)
  as.numeric(model$lookup[route])
}

aa_features <- function(R,tree,depth=2L) {
  nodes<-list(tree$root)
  if(depth>=2L)nodes<-c(nodes,
      Filter(Negate(is.null),tree$children))
  M<-do.call(cbind,lapply(nodes,function(node) {
    z<-vapply(node$anchors,function(a)
      as.integer(R[,node$target]>R[,a]),integer(nrow(R)))
    if(is.null(dim(z)))z<-matrix(z,ncol=1L)
    z
  }))
  colnames(M)<-paste0("comparison_",seq_len(ncol(M)))
  M
}
aa_global <- function(T,y,E,tree,depth=2L) {
  tr<-aa_features(aa_ranks(T),tree,depth)
  te<-aa_features(aa_ranks(E),tree,depth)
  fit<-glmnet::glmnet(x=tr,y=aa_y(y),alpha=0,lambda=.1,
              family="binomial",standardize=FALSE)
  as.numeric(predict(fit,newx=te,type="response",s=.1))
}

aa_inner_select <- function(T,yT) {
  fold<-rg_inner(yT)
  D<-T[fold$D,,drop=FALSE];V<-T[fold$V,,drop=FALSE]
  yd<-factor(as.character(yT[fold$D]),levels=levels(yT))
  yv<-factor(as.character(yT[fold$V]),levels=levels(yT))
  pool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
  parent<-rg_parent(D,yd,V,pool,K=20L)
  tree<-aa_train_tree(D,yd,pool,depth=2L)
  fallback<-list(tree=NULL,mode="NONE",accepted=FALSE,
                 V_delta=NA_real_,V_ba_delta=NA_real_,
                 root_gene="",anchors="",children=0L,
                 reason="no_valid_anchor_root")
  if(is.null(tree))return(fallback)
  rd<-aa_ranks(D);rv<-aa_ranks(V)
  stump<-aa_train_leaf(rd,yd,tree,depth=1L)
  deep<-aa_train_leaf(rd,yd,tree,depth=2L)
  ps<-aa_predict_leaf(stump,rv)
  pd<-aa_predict_leaf(deep,rv)
  vtruth<-aa_y(yv)
  mode<-if(rg_ll(vtruth,pd)<rg_ll(vtruth,ps)-.005)
      "DEPTH2" else "STUMP"
  valprob<-if(mode=="DEPTH2")pd else ps
  delta<-rg_ll(vtruth,parent$test)-rg_ll(vtruth,valprob)
  badelta<-rg_ba(vtruth,valprob)-rg_ba(vtruth,parent$test)
  accepted<-is.finite(delta)&&delta>=.005&&badelta>=-.02
  list(tree=tree,mode=mode,accepted=accepted,
       V_delta=delta,V_ba_delta=badelta,
       root_gene=tree$root$target,
       anchors=paste(tree$root$anchors,collapse="|"),
       children=sum(!vapply(tree$children,is.null,logical(1))),
       reason="V_EVIDENCE_GATE")
}

aa_eval <- function(T,yT,E,yE) {
  sel<-aa_inner_select(T,yT)
  pool<-compute_top_pairs_switchbox(T,yT,top_n=1000L,buffer_mult=5L)
  orig<-rg_parent(T,yT,E,pool,K=20L)$test
  p_stump<-orig;p_d2<-orig;p_global<-orig
  if(!is.null(sel$tree)){
    rt<-aa_ranks(T);re<-aa_ranks(E)
    p_stump<-aa_predict_leaf(aa_train_leaf(rt,yT,sel$tree,1L),re)
    p_d2<-aa_predict_leaf(aa_train_leaf(rt,yT,sel$tree,2L),re)
    p_global<-aa_global(T,yT,E,sel$tree,depth=2L)
  }
  p_tree<-if(sel$mode=="DEPTH2")p_d2 else p_stump
  p_gated<-if(isTRUE(sel$accepted))p_tree else orig
  list(probs=list(ORIGINAL_RLR=orig,
                  ANCHOR_STUMP=p_stump,
                  ANCHOR_TREE_D2=p_d2,
                  SAME_COMPARISONS_RIDGE=p_global,
                  V_SELECTED_ANCHOR_TREE=p_tree,
                  RLR_GATED_ANCHOR_TREE=p_gated),
       selection=sel,truth=aa_y(yE))
}
