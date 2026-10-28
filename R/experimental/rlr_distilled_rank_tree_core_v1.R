# RLR-DT v1: an autonomous 2-level compact rank-relation tree DISTILLED
# from the out-of-fold probabilities of the ORIGINAL unpublished RLR20.
# At prediction only comparisons of within-patient gene RANKS are used.
# The teacher RLR and its large pair panel are NOT run at inference time.
# The same rank split candidates and node budget are used by label-only
# TSPDT-like controls. This is exploratory, not claimed a novel TSPDT.
source("R/experimental/rlr_rank_gate_core_v1.R")

dt_need <- rg_need
dt_pair_names <- function(pairs){
 dt_need(all(c("gene_left","gene_right") %in% colnames(pairs)),
         "missing original gene pair identifiers")
 if("ordered_pair" %in% colnames(pairs)) as.character(pairs$ordered_pair)
 else paste(pairs$gene_left,">",pairs$gene_right)
}
# Deterministic FALLBACK if rare-class OOF folds break glmnet's
# internal binomial machinery. Same original RLR gene-pair ternary m10
# representation, but a distinct strongly ridge-regularized logistic fit.
# This is NOT literal original RLR20 training and is reported separately.
dt_tiny_rank_ridge_teacher <- function(Tr,yr,Va,pool,K=20L,lambda=2) {
 pairs<-apply_disjointness(pool,K=K,d=2L)
 xtr<-encode_pair_matrix(Tr,pairs,m=10L)
 xva<-encode_pair_matrix(Va,pairs,m=10L)
 z<-cbind(1,xtr)
 y<-as.numeric(yr==levels(yr)[2L])
 pen<-c(0,rep(1,ncol(xtr)))
 beta<-numeric(ncol(z))
 beta[1L]<-qlogis((sum(y)+1)/(length(y)+2))
 for(it in seq_len(60L)){
  eta<-drop(z%*%beta)
  p<-plogis(pmax(-25,pmin(25,eta)))
  weight<-pmax(1e-5,p*(1-p))
  hessian<-crossprod(z,z*weight)/nrow(z)+
          diag(lambda*pen+1e-8,nrow=ncol(z))
  grad<-drop(crossprod(z,p-y))/nrow(z)+lambda*pen*beta
  delta<-tryCatch(drop(solve(hessian,grad)),error=function(e)NULL)
  if(is.null(delta) || any(!is.finite(delta)))break
  # bounded Newton steps stabilize samples with one rare-class instance
  scale<-min(1,2/max(abs(delta)))
  beta<-beta-scale*delta
  if(max(abs(scale*delta))<1e-7)break
 }
 as.numeric(plogis(pmax(-25,pmin(25,
                      cbind(1,xva)%*%beta))))
}

dt_teacher_oof <- function(D,y,K=20L,folds=3L,seed=119L) {
 dt_need(nrow(D)==length(y) && min(table(y))>=folds,
         "insufficient class support for 3fold teacher OOF")
 parts<-make_repeated_stratified_folds(y,folds,1L,seed)
 pp<-rep(NA_real_,nrow(D))
 n_fallback<-0L
 for(sp in parts){
  tr<-sp$train_idx
  va<-sp$test_idx
  yr<-factor(as.character(y[tr]),levels=levels(y))
  pool<-compute_top_pairs_switchbox(D[tr,,drop=FALSE],yr,
                                    top_n=1000L,buffer_mult=5L)
  # cv.glmnet requires at least 4 folds. The frozen original helper
  # can choose nfolds=3 in tiny OOF subsets such as Pomeroy, although
  # the full D teacher uses the original unmodified implementation.
  # In that narrow case copy the SAME switchBox/d2/m10/glmnet pipeline
  # with the smallest permitted CV-fold count of 4 (train data only).
  teacher_result<-tryCatch({
    if(safe_nfolds(yr,5L)>=4L){
      rg_parent(D[tr,,drop=FALSE],yr,D[va,,drop=FALSE],
                pool,K=K)$test
    } else {
      pairs<-apply_disjointness(pool,K=K,d=2L)
      xtr<-encode_pair_matrix(D[tr,,drop=FALSE],pairs,m=10L)
      xva<-encode_pair_matrix(D[va,,drop=FALSE],pairs,m=10L)
      set.seed(42L)
      cv<-glmnet::cv.glmnet(x=xtr,
         y=as.numeric(yr==levels(yr)[2L]),family="binomial",
         alpha=.5,nfolds=4L,type.measure="deviance")
      as.numeric(predict(cv,newx=xva,s="lambda.min",
                         type="response"))
    }
  },error=function(e)NULL)
  if(is.null(teacher_result) ||
     length(teacher_result)!=length(va) ||
     any(!is.finite(teacher_result))){
    # The frozen original RLR model on the ENTIRE D is still evaluated
    # normally; ONLY the tiny OOF fold teacher gets this marked fallback.
    teacher_result<-dt_tiny_rank_ridge_teacher(
       D[tr,,drop=FALSE],yr,D[va,,drop=FALSE],pool,K=K)
    n_fallback<-n_fallback+1L
    message("RLRTREE_SMALL_COHORT_RIDGE_TEACHER_SURROGATE=TRUE")
  }
  pp[va]<-teacher_result
 }
 dt_need(all(is.finite(pp))&&all(pp>0&pp<1),
         "not truly out-of-fold teacher predictions")
 attr(pp,"ridge_fallback_folds")<-n_fallback
 pp
}
dt_parent_candidates <- function(parent,max_candidates=20L) {
 co<-as.numeric(stats::coef(parent$fit,s="lambda.min"))[-1L]
 dt_need(length(co)==nrow(parent$selected),
         "RLR original pair/coefficient mismatch")
 # Coefficient-ranked original RLR20 panel; full selected if many zero
 # coefficients so no special labels or E-induced pair search is needed.
 ix<-order(-abs(co),seq_along(co))
 tab<-parent$selected[head(ix,max_candidates),
                      c("gene_left","gene_right","ordered_pair"),
                      drop=FALSE]
 rownames(tab)<-NULL
 tab
}
dt_test_matrix <- function(X,pairs,mode=c("BINARY","TERNARY")) {
 mode<-match.arg(mode)
 gap<-rg_gaps(rg_ranks(X),pairs)
 if(mode=="BINARY") {
  z<-ifelse(gap>0,2L,1L)
 } else {
  z<-ifelse(gap>=10,3L,ifelse(gap<= -10,1L,2L))
 }
 dim(z)<-dim(gap)
 storage.mode(z)<-"integer"
 dimnames(z)<-list(rownames(X),dt_pair_names(pairs))
 z
}
dt_entropy <- function(targ,p){
 p<-pmax(1e-8,pmin(1-1e-8,p))
 sum(-targ*log(p)-(1-targ)*log1p(-p))
}
dt_leaf_probability <- function(targ,idx,prior,strength=4) {
 (sum(targ[idx])+strength*prior)/(length(idx)+strength)
}
dt_train <- function(D,y,pairs,p_teacher=NULL,
                     alpha=0.55,leaf_alpha=alpha,mode="BINARY",
                     max_depth=2L,max_internal=3L,
                     min_leaf=4L,penalty=.003,
                     prior_strength=4) {
 dt_need(length(y)==nrow(D)&&is.factor(y)&&length(levels(y))==2L,
         "invalid original binary labels")
 dt_need(nrow(pairs)>=1L&&is.finite(alpha)&&alpha>=0&&alpha<=1 &&
         max_depth>=1 && max_internal>=0 && min_leaf>=2,
         "invalid tree complexity or distillation loss")
 yy<-as.numeric(y==levels(y)[2L])
 if(is.null(p_teacher))p_teacher<-rep(.5,length(yy))
 dt_need(length(p_teacher)==length(yy)&&
         all(is.finite(p_teacher)) &&
         all(p_teacher>0&p_teacher<1),
         "OOF teacher probabilities invalid")
 # The proper cross-entropy soft target: class truth and teacher OOF mix.
 dt_need(is.finite(leaf_alpha)&&leaf_alpha>=0&&leaf_alpha<=1,
         "invalid RLR-guided tree leaf target alpha")
 target<-(1-alpha)*yy+alpha*p_teacher
 leaf_target<-(1-leaf_alpha)*yy+leaf_alpha*p_teacher
 prior<-(sum(leaf_target)+2)/(length(leaf_target)+4)
 split_prior<-(sum(target)+2)/(length(target)+4)
 mat<-dt_test_matrix(D,pairs,mode)
 n<-nrow(mat)
 nbranch<-if(mode=="BINARY")2L else 3L
 node<-list(list(id=1L,idx=seq_len(n),depth=0L,
     p=dt_leaf_probability(leaf_target,seq_len(n),prior,prior_strength),
     split=NA_integer_,children=integer()))
 leaf_ids<-1L
 n_internal<-0L
 repeat {
  if(n_internal>=max_internal)break
  proposals<-list()
  k<-0L
  for(id in leaf_ids) {
   nd<-node[[id]]
   if(nd$depth>=max_depth || length(nd$idx)<2L*min_leaf)next
   old_p<-dt_leaf_probability(target,nd$idx,split_prior,prior_strength)
   old_cost<-dt_entropy(target[nd$idx],old_p)
   for(feature in seq_len(ncol(mat))) {
     child<-lapply(seq_len(nbranch),function(b)
        nd$idx[mat[nd$idx,feature]==b])
     if(any(lengths(child)<min_leaf))next
     child_cost<-sum(vapply(child,function(ix){
       pchild<-dt_leaf_probability(target,ix,split_prior,prior_strength)
       dt_entropy(target[ix],pchild)
     },numeric(1)))
     gain<-(old_cost-child_cost)/n - penalty
     if(gain<=1e-12)next
     k<-k+1L
     proposals[[k]]<-list(node=id,feature=feature,child=child,
                          gain=gain,raw_gain=(old_cost-child_cost)/n)
   }
  }
  if(k==0L)break
  # Prefer largest cost-complexity-adjusted gain; stable low index tie.
  ix<-which.max(vapply(proposals,function(z)z$gain,numeric(1)))
  sel<-proposals[[ix]]
  nd<-node[[sel$node]]
  newids<-seq.int(length(node)+1L,length(node)+nbranch)
  for(b in seq_len(nbranch)) {
    ids<-sel$child[[b]]
    node[[newids[b]]]<-list(id=newids[b],idx=ids,
       depth=nd$depth+1L,
       p=dt_leaf_probability(leaf_target,ids,prior,prior_strength),
       split=NA_integer_,children=integer())
  }
  node[[sel$node]]$split<-sel$feature
  node[[sel$node]]$children<-newids
  leaf_ids<-c(setdiff(leaf_ids,sel$node),newids)
  n_internal<-n_internal+1L
 }
 observed<-sort(unique(vapply(node,function(nd)
   if(is.na(nd$split))NA_integer_ else nd$split,integer(1))),
   na.last=NA)
 structure(list(nodes=node,pairs=pairs,
    mode=mode,depth=max_depth,internal=n_internal,
    n_leaves=length(leaf_ids),used_pairs=observed,
    n_unique_genes=if(length(observed))
      length(unique(c(pairs$gene_left[observed],
                      pairs$gene_right[observed]))) else 0L,
    global_prior=prior,alpha=alpha,leaf_alpha=leaf_alpha,penalty=penalty,
    min_leaf=min_leaf,n_training=nrow(D)),
    class="rlr_distilled_rank_tree")
}
dt_predict <- function(fit,E){
 ranked<-rg_ranks(E)
 dt_need(all(fit$pairs$gene_left %in% colnames(ranked)) &&
         all(fit$pairs$gene_right %in% colnames(ranked)),
         "missing trained pair gene at inference")
 n<-nrow(E);p<-numeric(n);visits<-integer(n)
 paths<-character(n);ends<-integer(n)
 for(i in seq_len(n)){
  node_id<-1L;route<-character()
  repeat{
   cur<-fit$nodes[[node_id]]
   if(is.na(cur$split)){
     p[i]<-cur$p;ends[i]<-node_id;break
   }
   feat<-cur$split
   gap<-ranked[i,fit$pairs$gene_right[feat]]-
        ranked[i,fit$pairs$gene_left[feat]]
   state<-if(fit$mode=="BINARY") {
     if(gap>0)2L else 1L
   } else {
     if(gap>=10)3L else if(gap<=-10)1L else 2L
   }
   route<-c(route,paste0(dt_pair_names(fit$pairs)[feat],
                             ":",state))
   visits[i]<-visits[i]+1L
   node_id<-cur$children[[state]]
  }
  paths[i]<-paste(route,collapse=" / ")
 }
 dt_need(all(p>0&p<1),"invalid tree leaf probability")
 list(prob=p,label=as.integer(p>=.5),depth=visits,
      leaf=ends,path=paths)
}
dt_fidelity <- function(teacher_prob,student_prob){
 dt_need(length(teacher_prob)==length(student_prob) &&
           all(is.finite(teacher_prob))&&all(is.finite(student_prob)),
           "teacher/student samples mismatch")
 list(soft_logloss=dt_entropy(teacher_prob,student_prob)/length(teacher_prob),
      teacher_brier=mean((teacher_prob-student_prob)^2),
      agreement=mean(as.integer(teacher_prob>=.5)==
                     as.integer(student_prob>=.5)))
}
