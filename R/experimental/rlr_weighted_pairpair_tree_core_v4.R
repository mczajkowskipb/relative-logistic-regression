# RLR-WT v4 — autonomous rank tree with logistic-coefficient-weighted
# TWO-pair ordinal predicates at EACH node. Original unpublished RLR20
# provides its disjoint pair panel and *coefficients*, not an inference
# dependency. A test costs precisely TWO original ternary pair comparisons.
# Max depth 2 => at most FOUR pair comparisons along any patient path.
#
# Unlike TSPDT single-pair nodes, weighted node uses RLR's learned beta.
# IMPORTANT: previously published TSPDT/GTSPDT/weighted-kTSP already
# encompass multipair tree tests; novelty not presumed on architecture
# alone: comparisons vs unit-sign and same-weight label-only necessary.
source("R/experimental/rlr_distilled_rank_tree_core_v1.R")
wt_need<-rg_need

wt_candidates <- function(parent,max_pairs=8L,mode=c("BETA","SIGN")) {
 mode<-match.arg(mode)
 allcoef<-as.numeric(stats::coef(parent$fit,s="lambda.min"))[-1L]
 wt_need(length(allcoef)==nrow(parent$selected) &&
           length(allcoef)>=2,"original RLR coefficient alignment")
 pick<-head(order(-abs(allcoef),seq_along(allcoef)),max_pairs)
 pp<-parent$selected[pick,c("gene_left","gene_right",
                            "ordered_pair"),drop=FALSE]
 beta<-allcoef[pick]
 # Normalizing by max abs coefficient preserves every zero-threshold
 # weighted test, makes diagnostics comparable across folds/datasets.
 if(mode=="SIGN")beta<-sign(beta)
 if(mode=="BETA"&&max(abs(beta))>0)beta<-beta/max(abs(beta))
 cmb<-utils::combn(seq_along(pick),2L)
 cand<-data.frame(
  gene1=as.integer(cmb[1L,]),gene2=as.integer(cmb[2L,]),
  weight1=as.numeric(beta[cmb[1L,]]),
  weight2=as.numeric(beta[cmb[2L,]]),
  stringsAsFactors=FALSE)
 rownames(pp)<-NULL
 list(pairs=pp,candidates=cand,mode=mode,
      originally_nonzero=sum(abs(allcoef)>1e-9),
      RLR_nominal=nrow(parent$selected))
}
wt_term <- function(gap,m=10L){
 ifelse(gap>=m,1L,ifelse(gap<= -m,-1L,0L))
}
wt_split_matrix <- function(X,selection){
 Z<-encode_pair_matrix(X,selection$pairs,m=10L)
 C<-selection$candidates
 ans<-vapply(seq_len(nrow(C)),function(j){
  as.integer(C$weight1[j]*Z[,C$gene1[j]] +
             C$weight2[j]*Z[,C$gene2[j]]>0)+1L
 },integer(nrow(X)))
 if(is.null(dim(ans)))ans<-matrix(ans,ncol=1L)
 colnames(ans)<-paste0("w",seq_len(ncol(ans)))
 rownames(ans)<-rownames(X)
 ans
}
wt_fit <- function(D,y,original_parent,p_teacher=NULL,
                   coeff_mode="BETA",alpha=.25,leaf_alpha=0,
                   max_depth=2L,max_internal=3L,
                   min_leaf=4L,penalty=.003) {
 wt_need(is.factor(y)&&length(levels(y))==2L&&length(y)==nrow(D),
          "invalid original train classes/shape")
 selection<-wt_candidates(original_parent,max_pairs=8L,
                          mode=coeff_mode)
 A<-wt_split_matrix(D,selection)
 n<-nrow(A);yy<-as.numeric(y==levels(y)[2L])
 if(is.null(p_teacher))p_teacher<-rep(.5,n)
 wt_need(length(p_teacher)==n&&all(is.finite(p_teacher))&&
    all(p_teacher>0&p_teacher<1),"teacher OOF not aligned")
 z<-(1-alpha)*yy+alpha*p_teacher
 zleaf<-(1-leaf_alpha)*yy+leaf_alpha*p_teacher
 prior<-(sum(zleaf)+2)/(n+4)
 splitprior<-(sum(z)+2)/(n+4)
 prob<-function(x,idx,prior){
    (sum(x[idx])+4*prior)/(length(idx)+4)
 }
 nodes<-list(list(id=1L,idx=seq_len(n),depth=0L,
   p=prob(zleaf,seq_len(n),prior),split=NA_integer_,
   children=integer()))
 leaves<-1L
 for(step in seq_len(max_internal)){
  scored<-list();k<-0L
  for(id in leaves){
   nd<-nodes[[id]]
   if(nd$depth>=max_depth||length(nd$idx)<2L*min_leaf)next
   old<-dt_entropy(z[nd$idx],prob(z,nd$idx,splitprior))
   for(feature in seq_len(ncol(A))){
    ids0<-nd$idx[A[nd$idx,feature]==1L]
    ids1<-nd$idx[A[nd$idx,feature]==2L]
    if(min(length(ids0),length(ids1))<min_leaf)next
    loss<-dt_entropy(z[ids0],prob(z,ids0,splitprior))+
          dt_entropy(z[ids1],prob(z,ids1,splitprior))
    gain<-(old-loss)/n-penalty
    if(gain<=1e-12)next
    k<-k+1L
    scored[[k]]<-list(parent=id,split=feature,
                      child=list(ids0,ids1),gain=gain)
   }
  }
  if(k==0L)break
  chosen<-scored[[which.max(vapply(scored,
                                  function(z)z$gain,numeric(1)))]]
  nd<-nodes[[chosen$parent]]
  idx<-seq.int(length(nodes)+1L,length(nodes)+2L)
  for(b in 1:2){
    ids<-chosen$child[[b]]
    nodes[[idx[b]]]<-list(id=idx[b],idx=ids,depth=nd$depth+1L,
      p=prob(zleaf,ids,prior),split=NA_integer_,children=integer())
  }
  nodes[[chosen$parent]]$split<-chosen$split
  nodes[[chosen$parent]]$children<-idx
  leaves<-c(setdiff(leaves,chosen$parent),idx)
 }
 selected_nodes<-vapply(nodes,function(z)
     if(is.na(z$split))NA_integer_ else z$split,integer(1))
 used<-sort(unique(selected_nodes[!is.na(selected_nodes)]))
 pairix<-if(length(used))unique(c(
   selection$candidates$gene1[used],
   selection$candidates$gene2[used])) else integer()
 structure(list(nodes=nodes,selector=selection,
   used_tests=used,used_pairs=pairix,
   n_internal=length(nodes)-length(leaves),n_leaves=length(leaves),
   n_unique_pairs=length(pairix),
   n_unique_genes=if(length(pairix))length(unique(c(
      selection$pairs$gene_left[pairix],
      selection$pairs$gene_right[pairix])))else 0L,
   max_depth=max_depth,alpha=alpha,leaf_alpha=leaf_alpha,
   min_leaf=min_leaf,penalty=penalty,training_n=n,
   max_compared_per_patient=2L*max_depth),
   class="rlr_weighted_pairpair_tree")
}
wt_predict <- function(model,E) {
 ranks<-rg_ranks(E)
 P<-model$selector$pairs;C<-model$selector$candidates
 wt_need(all(P$gene_left %in% colnames(ranks))&&
         all(P$gene_right %in% colnames(ranks)),
         "model's trained original gene pair probe absent")
 n<-nrow(E);probs<-numeric(n);visits<-integer(n)
 leaves<-integer(n);paths<-character(n)
 for(i in seq_len(n)){
  id<-1L;route<-character()
  repeat{
   node<-model$nodes[[id]]
   if(is.na(node$split)){
     probs[i]<-node$p;leaves[i]<-id;break
   }
   c<-C[node$split,,drop=FALSE]
   p1<-c$gene1[1L];p2<-c$gene2[1L]
   d1<-ranks[i,P$gene_right[p1]]-ranks[i,P$gene_left[p1]]
   d2<-ranks[i,P$gene_right[p2]]-ranks[i,P$gene_left[p2]]
   r1<-wt_term(d1);r2<-wt_term(d2)
   score<-c$weight1[1L]*r1+c$weight2[1L]*r2
   branch<-if(score>0)2L else 1L
   route<-c(route,sprintf("[%s %.3f; %s %.3f] state=%d",
      P$ordered_pair[p1],c$weight1[1L],
      P$ordered_pair[p2],c$weight2[1L],branch))
   visits[i]<-visits[i]+2L
   id<-node$children[[branch]]
  }
  paths[i]<-paste(route,collapse=" / ")
 }
 wt_need(all(probs>0&probs<1)&&all(visits<=4L),
          "invalid probability/path complexity")
 list(prob=probs,label=as.integer(probs>=.5),
      n_evaluated=visits,leaf=leaves,path=paths)
}
