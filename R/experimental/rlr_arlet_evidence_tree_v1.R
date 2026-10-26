# ARLET v1 — probabilistic three-way relation-evidence cascade.
# A TRUE conditional tree: two confident branches, one unresolved branch
# evaluated by a second, different 3/5/7-relation logistic multitest.
# Original RLR ternary rank evidence m=10; no full logistic fitting IN leaves.
# Conditional feature nomination uses D only. EA searches over gate relations,
# NOT predictive labels E. No new cohorts, no manuscript/source mutation.
arlet_clip<-function(p)pmin(1-1e-7,pmax(1e-7,p))
arlet_need<-function(x,msg)if(!isTRUE(x))stop("ARLET_STOP: ",msg,call.=FALSE)
arlet_ll<-function(y,p) {
  p<-arlet_clip(p);mean(-y*log(p)-(1-y)*log1p(-p))
}
arlet_scorer<-function(y,p) {
  predicted<-as.integer(p>=.5)
  r0<-if(sum(y==0))mean(predicted[y==0]==0) else 0
  r1<-if(sum(y==1))mean(predicted[y==1]==1) else 0
  f0<-if(sum(predicted==0)+sum(y==0)>0)
    2*sum(predicted==0 & y==0)/(sum(predicted==0)+sum(y==0)) else 0
  f1<-if(sum(predicted==1)+sum(y==1)>0)
    2*sum(predicted==1 & y==1)/(sum(predicted==1)+sum(y==1)) else 0
  list(logloss=arlet_ll(y,p),brier=mean((p-y)^2),
       ba=(r0+r1)/2,macro_f1=(f0+f1)/2)
}
arlet_pair_pool<-function(D,y,seed=42L) {
  # Original signed switchBox top 1000 is INCLUDED as-is.
  signed<-compute_top_pairs_switchbox(D,y,top_n=1000L,
                                        buffer_mult=5L)
  g<-colnames(D)
  # Deterministic seeded random gene pairs from the ENTIRE input feature
  # universe, not the global RLR's disjoint shortlist.
  set.seed(seed)
  q<-min(450L,max(250L,length(g)*2L))
  one<-sample.int(length(g),q,replace=TRUE)
  two<-sample.int(length(g),q,replace=TRUE)
  keep<-one!=two
  extra<-data.frame(gene_left=g[one[keep]],
              gene_right=g[two[keep]],stringsAsFactors=FALSE)
  extra$ordered_pair<-paste0(extra$gene_left," > ",extra$gene_right)
  signed<-signed[,c("gene_left","gene_right","ordered_pair"),
                  drop=FALSE]
  merged<-rbind(signed,extra)
  merged<-merged[!duplicated(merged$ordered_pair),,drop=FALSE]
  rownames(merged)<-NULL
  arlet_need(nrow(merged)>50,"candidate pool too narrow")
  merged
}
arlet_encode_both<-function(D,E,candidates) {
  # Rank EACH patient across all assayed genes, exactly as original RLR.
  ranks_d<-t(apply(D,1,rank,ties.method="average"))
  ranks_e<-t(apply(E,1,rank,ties.method="average"))
  colnames(ranks_d)<-colnames(D)
  colnames(ranks_e)<-colnames(E)
  enc<-function(rk) {
    result<-vapply(seq_len(nrow(candidates)),function(j) {
      dx<-rk[,candidates$gene_right[j]]-rk[,candidates$gene_left[j]]
      ifelse(dx>=10,1,ifelse(dx<=(-10),-1,0))
    },numeric(nrow(rk)))
    if(!is.matrix(result))result<-matrix(result,ncol=nrow(candidates))
    colnames(result)<-candidates$ordered_pair
    rownames(result)<-rownames(rk)
    result
  }
  list(train=enc(ranks_d),test=enc(ranks_e))
}
arlet_select<-function(proposal,pool,size,exclude=integer()) {
  # Strict gene disjointness per multitest, no need to force disjoint genes
  # BETWEEN nodes: reusing biological genes on different branches can help.
  proposal<-unique(as.integer(proposal))
  proposal<-c(proposal,setdiff(seq_len(nrow(pool)),proposal))
  selected<-integer();used<-character()
  for(j in proposal) {
    if(j %in% exclude || j<1L || j>nrow(pool))next
    a<-pool$gene_left[[j]];b<-pool$gene_right[[j]]
    if(a==b || a %in% used || b %in% used)next
    selected<-c(selected,j);used<-c(used,a,b)
    if(length(selected)>=size)break
  }
  selected
}
arlet_lrfit<-function(X,y) {
  arlet_need(is.matrix(X)&&nrow(X)==length(y)&&ncol(X)>0,
             "small logistic design invalid")
  fallback<-mean(y)
  if(length(unique(y))<2L || !any(apply(X,2,sd)>0))
    return(list(type="intercept",p=fallback))
  fit<-tryCatch(suppressWarnings(glmnet::glmnet(
          x=X,y=y,family="binomial",alpha=.5,
          lambda=.07,standardize=FALSE,
          maxit=50000L)),error=function(e)NULL)
  if(is.null(fit))return(list(type="intercept",p=fallback))
  list(type="glmnet",fit=fit)
}
arlet_lrpredict<-function(fit,X) {
  if(fit$type=="intercept")return(rep(fit$p,nrow(X)))
  as.numeric(predict(fit$fit,newx=X,s=.07,type="response"))
}
arlet_fit<-function(R,y,pool,proposal,s,confidence,
                      second_size=7L) {
  top<-arlet_select(proposal,pool,s)
  arlet_need(length(top)>=s,"not enough diverse root relations")
  model1<-arlet_lrfit(R[,top,drop=FALSE],y)
  p_root<-arlet_lrpredict(model1,R[,top,drop=FALSE])
  root_active<-rowSums(R[,top,drop=FALSE]!=0)
  defer<-root_active<2L | (p_root>1-confidence & p_root<confidence)
  # Select SECOND STAGE on difficult training cases, but its coefficients
  # are fitted on ALL D observations to avoid extremely small leaves.
  subset<-which(defer)
  if(length(subset)<12 || length(unique(y[subset]))<2L)subset<-seq_along(y)
  cond<-vapply(seq_len(ncol(R)),function(j){
     r<-R[subset,j]
     a<-y[subset]==1
     if(!any(a)||!any(!a))return(0)
     abs(mean(r[a])-mean(r[!a]))
  },numeric(1))
  order_cond<-order(-cond,seq_along(cond))
  second<-arlet_select(order_cond,pool,second_size,exclude=top)
  arlet_need(length(second)>=3L,"second stage lacks 3 disjoint features")
  model2<-arlet_lrfit(R[,second,drop=FALSE],y)
  p2<-arlet_lrpredict(model2,R[,second,drop=FALSE])
  pred<-ifelse(defer,p2,p_root)
  list(gate=list(top=top,s=s,confidence=confidence,
                 second=second,
                 model1=model1,model2=model2),
       D_deferred=sum(defer),p_train=pred)
}
arlet_predict<-function(fit,R) {
  m<-fit$gate
  rt<-R[,m$top,drop=FALSE]
  root<-arlet_lrpredict(m$model1,rt)
  a<-rowSums(rt!=0)
  defer<-a<2L | (root>1-m$confidence & root<m$confidence)
  p2<-arlet_lrpredict(m$model2,R[,m$second,drop=FALSE])
  p<-ifelse(defer,p2,root)
  # Tested relations root always evaluated; second only if deferred.
  used<-a+ifelse(defer,rowSums(R[,m$second,drop=FALSE]!=0),0)
  list(prob=p,deferred=defer,active_evidence=used,
       evaluated=length(m$top)+defer*length(m$second),
       route=ifelse(defer,"UNRESOLVED",ifelse(root<.5,
                                     "CONFIDENT_CLASS0","CONFIDENT_CLASS1")))
}
arlet_search<-function(RD,yD,RV,yV,pool,seed=42L) {
  score<-function(proposal,s,conf) {
    f<-arlet_fit(RD,yD,pool,proposal,s,conf)
    p<-arlet_predict(f,RV)
    # Complexity measured by actual active evidence on each patient, and
    # number of evaluated pair tests. Penalty fixed before looking at E.
    objective<-arlet_ll(yV,p$prob)+
               .0008*mean(p$active_evidence)+
               .0004*mean(p$evaluated)
    list(fitness=objective,raw_ll=arlet_ll(yV,p$prob),
         s=s,conf=conf,proposal=proposal,
         val_deferral=mean(p$deferred))
  }
  # Deterministic root, based on training-only signed-score ranking.
  det<-score(1:5,5L,.70)
  set.seed(seed+23L)
  best<-det
  # REAL EA: 3 generations x 8 offspring, elitist winner retention,
  # mutation of selected relation genes or confidence/split size.
  # Bounded trial budget; V is nested inside outer training.
  pop<-list(det)
  for(i in 1:7) {
    proposal<-sample(seq_len(nrow(pool)),sample(c(3L,5L),1L))
    pop[[length(pop)+1L]]<-score(proposal,length(proposal),
                                  sample(c(.65,.75,.85),1L))
  }
  for(gen in 1:3){
    pop<-pop[order(vapply(pop,function(z)z$fitness,numeric(1)))]
    if(pop[[1L]]$fitness<best$fitness)best<-pop[[1L]]
    children<-list(pop[[1L]],pop[[2L]])
    while(length(children)<8L){
      parent<-pop[[sample(1:4,1L)]]
      s<-if(runif(1)<.25)sample(c(3L,5L),1L) else parent$s
      conf<-if(runif(1)<.4)sample(c(.65,.75,.85),1L) else parent$conf
      proposal<-parent$proposal
      if(length(proposal)<s)proposal<-c(proposal,
                        sample(seq_len(nrow(pool)),s-length(proposal)))
      proposal<-proposal[seq_len(s)]
      proposal[sample(seq_len(s),1L)]<-sample(seq_len(nrow(pool)),1L)
      children[[length(children)+1L]]<-score(proposal,s,conf)
    }
    pop<-children
  }
  pop<-pop[order(vapply(pop,function(z)z$fitness,numeric(1)))]
  if(pop[[1L]]$fitness<best$fitness)best<-pop[[1L]]
  # Winner is selected on training-internal V only; no outer E.
  list(deterministic=det,evolutionary=best,
       ea_val_gain=det$fitness-best$fitness,
       distinct_ea=best$fitness < det$fitness)
}
# Small, fair complexity-control proxy for the CUDA TSPDT, NOT an original
# TSPDT reproduction. Uses a single pair at each split (m10).
arlet_tspdt_proxy<-function(R,y,RT,depth=2L,min_leaf=7L) {
  recurse<-function(ix,ix_test,d) {
    p<-(sum(y[ix])+1)/(length(ix)+2)
    if(d>=depth || length(ix)<2L*min_leaf ||
       length(unique(y[ix]))<2L || !length(ix_test))
      return(rep(p,length(ix_test)))
    best_score<-Inf;best_j<-0L
    cands<-seq_len(min(ncol(R),150L))
    for(j in cands){
      z<-R[ix,j]>0
      if(sum(z)<min_leaf || sum(!z)<min_leaf)next
      impurity<-function(t)if(length(t))2*mean(t)*(1-mean(t)) else 1
      v<-sum(z)*impurity(y[ix][z])+sum(!z)*impurity(y[ix][!z])
      if(v<best_score){best_score<-v;best_j<-j}
    }
    if(best_j==0L)return(rep(p,length(ix_test)))
    ztr<-R[ix,best_j]>0;zte<-RT[ix_test,best_j]>0
    ans<-numeric(length(ix_test))
    ans[zte]<-recurse(ix[which(ztr)],ix_test[which(zte)],d+1L)
    ans[!zte]<-recurse(ix[which(!ztr)],ix_test[which(!zte)],d+1L)
    ans
  }
  recurse(seq_along(y),seq_len(nrow(RT)),0L)
}
