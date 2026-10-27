# SAVS v1 — Selective Anchor-Vote Stump.
# NONPARAMETRIC two-branch interpretable rank-only tree. Every leaf holds
# a short list of one-SIDED, class-directed ordinal pair rules. A rule
# contributes a +/-1 vote ONLY WHEN its train-frozen rank inequality
# fires; otherwise it ABSTAINS. No per-leaf glmnet or minimum class count.
# All features/root/directions/calibration use original outer TRAIN only.
#
# IMPORTANT: unlike the original RLR m10 ternary sign, an opposite
# ordering is NOT an active negative vote of the same feature here.
# It is an ABSTENTION; negative votes come from distinct class-0
# associated one-sided rules. Do not conflate with RLR's exact notion.
source("R/experimental/rlr_ordinal_anchor_score_tree_v1.R")

sv_need <- rg_need
sv_truth <- os_truth
sv_m <- function(p, mode=c("RLR10","RELATIVE10PCT")){
  mode<-match.arg(mode)
  if(mode=="RLR10")return(10L)
  as.integer(max(10L,ceiling(.10*(p-1L))))
}
sv_rule_cols<-c("gene_left","gene_right","fire_sign","vote_sign","strength")
sv_empty_rules<-function(){
  data.frame(gene_left=character(),gene_right=character(),
    fire_sign=integer(),vote_sign=integer(),strength=numeric(),
    stringsAsFactors=FALSE)
}
sv_active <- function(R, rules, margin){
  if(nrow(rules)==0L)
    return(matrix(FALSE,nrow=nrow(R),ncol=0L))
  sv_need(all(rules$gene_left %in% colnames(R)) &&
          all(rules$gene_right %in% colnames(R)),
          "missing trained gene/probe at prediction")
  gap<-vapply(seq_len(nrow(rules)),function(j)
       R[,rules$gene_right[j]]-R[,rules$gene_left[j]],
       numeric(nrow(R)))
  if(is.null(dim(gap)))gap<-matrix(gap,ncol=1L)
  sweep(gap,2L,rules$fire_sign,"*")>=margin
}
sv_cast <- function(R,rules,margin){
  firings<-sv_active(R,rules,margin)
  if(ncol(firings)==0L)
    return(list(vote=integer(nrow(R)),active=integer(nrow(R)),
      bucket=integer(nrow(R))))
  s<-as.integer(firings %*% as.numeric(rules$vote_sign))
  list(vote=s,active=as.integer(rowSums(firings)),
       bucket=as.integer(sign(s)))
}
# Conditioned rule selection: prefer balanced support for BOTH class votes,
# choose at most 2 positive and 2 negative one-sided tests (4 total).
# Candidate genes cannot be repeated on one leaf path. Selection works
# with a single-class branch, falling back to original D global contrast.
sv_choose_rules <- function(R,y,pool,ids,margin=10L,max_per_class=2L,
                            excluded=character(), max_pool=350L){
  yy<-sv_truth(y)
  zz<-which(seq_len(nrow(R)) %in% ids)
  source_ids<-zz
  if(length(zz)<4L || length(unique(yy[zz]))<2L)
    source_ids<-seq_len(nrow(R))
  pos<-source_ids[yy[source_ids]==1L]
  neg<-source_ids[yy[source_ids]==0L]
  if(length(pos)==0L || length(neg)==0L)return(sv_empty_rules())
  pp<-pool[seq_len(min(nrow(pool),max_pool)),,drop=FALSE]
  pp<-pp[!pp$gene_left %in% excluded &
         !pp$gene_right %in% excluded,,drop=FALSE]
  if(nrow(pp)==0L)return(sv_empty_rules())
  gap<-rg_gaps(R,pp)
  scored<-list()
  k<-0L
  for(fire in c(1L,-1L)){
    flag<-gap*fire >= margin
    cp<-colMeans(flag[pos,,drop=FALSE])
    cn<-colMeans(flag[neg,,drop=FALSE])
    contrast<-cp-cn
    prevalence<-colMeans(flag[source_ids,,drop=FALSE])
    # Hard sparse support screen; prefer rules firing in a reasonable
    # fraction of train samples rather than one noisy outlier.
    good<-which(abs(contrast)>1e-12 &
                prevalence>=.08 & prevalence<=.92)
    if(!length(good))next
    for(j in good){
      k<-k+1L
      scored[[k]]<-data.frame(
       gene_left=as.character(pp$gene_left[j]),
       gene_right=as.character(pp$gene_right[j]),
       fire_sign=as.integer(fire),
       vote_sign=as.integer(sign(contrast[j])),
       strength=abs(contrast[j])*sqrt(prevalence[j]*(1-prevalence[j])),
       stringsAsFactors=FALSE)
    }
  }
  if(length(scored)==0L)return(sv_empty_rules())
  options<-do.call(rbind,scored)
  options<-options[order(-options$strength,options$gene_left,
              options$gene_right,options$fire_sign),,drop=FALSE]
  chosen<-integer()
  used<-excluded
  counts<-setNames(c(0L,0L),c("-1","1"))
  for(i in seq_len(nrow(options))){
    current<-options[i,]
    genes<-c(current$gene_left,current$gene_right)
    tag<-as.character(current$vote_sign)
    if(counts[[tag]]>=max_per_class||any(genes %in% used))next
    chosen<-c(chosen,i)
    used<-c(used,genes)
    counts[[tag]]<-counts[[tag]]+1L
    if(length(chosen)>=2L*max_per_class)break
  }
  if(!length(chosen))return(sv_empty_rules())
  ans<-options[chosen,sv_rule_cols,drop=FALSE]
  rownames(ans)<-NULL
  ans
}
sv_bucket_fit<-function(R,y,ids,rules,margin,global_prior,strength=4){
  yy<-sv_truth(y)[ids]
  n<-length(ids)
  branch_prior<-(sum(yy)+5*global_prior)/(n+5)
  score<-sv_cast(R[ids,,drop=FALSE],rules,margin)$bucket
  probs<-vapply(-1:1,function(b){
    ii<-which(score==b)
    (sum(yy[ii])+strength*branch_prior)/(length(ii)+strength)
  },numeric(1))
  list(probs=probs,prior=branch_prior)
}
sv_fit<-function(D,y,pool=NULL,root=NULL,fixed_rules=NULL,
                  mode="RLR10",max_targets=8L){
  sv_need(nrow(D)==length(y) && length(levels(y))==2L,
          "invalid binary training labels")
  R<-rg_ranks(D)
  margin<-sv_m(ncol(D),mode)
  if(is.null(pool))
    pool<-compute_top_pairs_switchbox(D,y,top_n=1000L,buffer_mult=5L)
  if(is.null(root)){
    propose<-os_choose_node(R,y,pool,
                 leaf_min=max(3L,as.integer(floor(.08*nrow(D)))),
                 each_class_min=0L,max_targets=max_targets)
    root<-if(is.null(propose))NULL else propose$node
  }
  if(is.null(root)) {
    branch<-rep(0L,nrow(D))
  } else {
    branch<-os_gate(R,root)
  }
  yy<-sv_truth(y)
  prior<-(sum(yy)+2)/(length(yy)+4)
  rules<-vector("list",2L)
  fits<-vector("list",2L)
  for(b in 0:1){
    ids<-which(branch==b)
    excluded<-if(is.null(root))character() else
      c(root$target,root$anchors)
    frozen<-if(is.null(fixed_rules))NULL else fixed_rules[[b+1L]]
    rules[[b+1L]]<-if(is.null(fixed_rules))
      sv_choose_rules(R,y,pool,ids,margin,excluded=excluded)
    else frozen
    fits[[b+1L]]<-sv_bucket_fit(R,y,ids,rules[[b+1L]],
                  margin,global_prior=prior)
  }
  list(root=root,rule_panels=rules,fits=fits,margin=margin,
       prior=prior,training_branches=c(sum(branch==0L),
                                      sum(branch==1L)),
       n_rules=vapply(rules,nrow,integer(1)))
}
sv_predict<-function(model,E) {
 R<-rg_ranks(E)
 branch<-if(is.null(model$root))rep(0L,nrow(E)) else
       os_gate(R,model$root)
 p<-rep(model$prior,nrow(E))
 activity<-integer(nrow(E));balance<-integer(nrow(E))
 panel_size<-integer(nrow(E))
 for(b in 0:1){
   ids<-which(branch==b)
   if(!length(ids))next
   rules<-model$rule_panels[[b+1L]]
   votes<-sv_cast(R[ids,,drop=FALSE],rules,model$margin)
   p[ids]<-model$fits[[b+1L]]$probs[votes$bucket+2L]
   activity[ids]<-votes$active
   balance[ids]<-votes$vote
   panel_size[ids]<-nrow(rules)
 }
 list(p=p,branch=branch,n_active=activity,
      n_considered=panel_size,
      n_root=if(is.null(model$root))0L else length(model$root$anchors),
      sum_votes=balance,
      abstain=as.integer(activity==0L),
      active_fraction=ifelse(panel_size>0L,activity/panel_size,0))
}
# Global comparator: exactly same nominated conditional gene-pair tests,
# each with its original direction and class vote. Deduplicate identical
# tests but do not add new genes or select from external labels.
sv_samepair_global<-function(D,y,model){
  all<-do.call(rbind,model$rule_panels)
  if(is.null(all)||nrow(all)==0L)all<-sv_empty_rules()
  if(nrow(all)){
    sig<-paste(all$gene_left,all$gene_right,
               all$fire_sign,all$vote_sign,sep="|")
    all<-all[!duplicated(sig),sv_rule_cols,drop=FALSE]
    rownames(all)<-NULL
  }
  p<-(sum(sv_truth(y))+2)/(length(y)+4)
  base<-list(root=NULL,rule_panels=list(all,sv_empty_rules()),
         margin=model$margin,prior=p,
         training_branches=c(nrow(D),0L),n_rules=c(nrow(all),0L))
  R<-rg_ranks(D)
  base$fits<-list(sv_bucket_fit(R,y,seq_len(nrow(D)),all,
                        model$margin,p),
                  list(probs=rep(p,3L),prior=p))
  base
}
sv_root_only<-function(D,y,model){
  p<-(sum(sv_truth(y))+2)/(length(y)+4)
  R<-rg_ranks(D)
  branch<-if(is.null(model$root))rep(0L,nrow(D)) else
         os_gate(R,model$root)
  panels<-list(sv_empty_rules(),sv_empty_rules())
  list(root=model$root,rule_panels=panels,margin=model$margin,
       prior=p,training_branches=c(sum(branch==0L),sum(branch==1L)),
       n_rules=c(0L,0L),
       fits=lapply(0:1,function(b){
          pr<-(sum(sv_truth(y)[branch==b])+5*p)/
                      (sum(branch==b)+5)
          list(probs=rep(pr,3L),prior=pr)
       }))
}
sv_global4<-function(D,y,pool,margin){
  R<-rg_ranks(D)
  rules<-sv_choose_rules(R,y,pool,seq_len(nrow(D)),
                        margin,excluded=character())
  m<-list(root=NULL,rule_panels=list(rules,sv_empty_rules()),
     margin=margin,prior=(sum(sv_truth(y))+2)/(nrow(D)+4),
     training_branches=c(nrow(D),0L),
     n_rules=c(nrow(rules),0L))
  m$fits<-list(sv_bucket_fit(R,y,seq_len(nrow(D)),rules,
                    margin,m$prior),
             list(probs=rep(m$prior,3L),prior=m$prior))
  m
}
