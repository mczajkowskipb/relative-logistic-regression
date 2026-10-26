# T-GRLR V1: a *single* global elastic-net RLR with branch-gated relational
# interaction terms, driven by a 3-way (-1/0/+1) biological gene-pair root.
# A joint glmnet fit (not independent tiny-leaf experts) learns sparse,
# branch-specific corrections. D-only supervised candidate selection;
# V-only model selection; E is untouched for comparisons.
# This is an exploratory hierarchical logistic interaction method; it must
# NOT be advertised as a conventional terminal-classification tree.
tgr_need<-function(ok,msg)if(!isTRUE(ok))stop("TGR_STOP: ",msg,call.=FALSE)
tgr_ll<-function(y,p){
  p<-pmin(1-1e-7,pmax(1e-7,p))
  mean(-y*log(p)-(1-y)*log1p(-p))
}
tgr_ba<-function(y,p){
  v<-as.integer(p>=.5)
  mean(c(if(any(y==0L))mean(v[y==0L]==0L)else 0,
         if(any(y==1L))mean(v[y==1L]==1L)else 0))
}
tgr_f1<-function(y,p){
  z<-as.integer(p>=.5)
  mean(vapply(c(0L,1L),function(a){
    tp<-sum(z==a & y==a)
    if((sum(z==a)+sum(y==a))==0L)return(0)
    2*tp/(sum(z==a)+sum(y==a))
  },numeric(1)))
}
tgr_enfit<-function(R,y,seed=42L){
  tgr_need(is.matrix(R)&&ncol(R)>=2L&&nrow(R)==length(y),
           "elastic-net design invalid")
  set.seed(seed)
  cv<-tryCatch(suppressWarnings(glmnet::cv.glmnet(
           x=R,y=as.integer(y),family="binomial",alpha=.5,
           nfolds=min(5L,max(3L,min(sum(y==0),sum(y==1)))),
           type.measure="deviance")),error=function(e)NULL)
  if(is.null(cv))return(list(type="mean",p=mean(y)))
  list(type="cv",cv=cv)
}
tgr_predict<-function(fit,R){
  if(fit$type=="mean")return(rep(fit$p,nrow(R)))
  as.numeric(predict(fit$cv,newx=R,s="lambda.min",type="response"))
}
tgr_allpairs<-function(X,y,seed=42L){
  # Expand candidate set beyond global K10/50, but keep EXACT m10 encoding.
  top<-compute_top_pairs_switchbox(X,factor(y),top_n=1000L,
                                  buffer_mult=5L)
  p<-top[,c("gene_left","gene_right","ordered_pair"),drop=FALSE]
  set.seed(seed)
  features<-colnames(X);n<-length(features)
  m<-min(550L,max(275L,2L*n))
  a<-sample.int(n,m,replace=TRUE)
  b<-sample.int(n,m,replace=TRUE)
  keep<-a!=b
  ex<-data.frame(gene_left=features[a[keep]],
                 gene_right=features[b[keep]],stringsAsFactors=FALSE)
  ex$ordered_pair<-paste0(ex$gene_left," > ",ex$gene_right)
  cands<-rbind(p,ex)
  cands<-cands[!duplicated(cands$ordered_pair),,drop=FALSE]
  rownames(cands)<-NULL
  tgr_need(nrow(cands)>=500L,"candidate pool too small")
  list(pool=cands,source_top=top)
}
tgr_design<-function(X,main,root,branch,candidate){
  A<-encode_pair_matrix(X,main,m=10L)
  B<-encode_pair_matrix(X,candidate,m=10L)
  base<-matrix(B[,root],ncol=1L)
  colnames(base)<-"ROOT_PAIR"
  mats<-list(A,base)
  for(s in c(-1L,0L,1L)){
    jj<-branch[[as.character(s)]]
    if(!length(jj))next
    Z<-B[,jj,drop=FALSE]*(B[,root]==s)
    colnames(Z)<-paste0("COND_",s,"_",seq_along(jj))
    mats[[length(mats)+1L]]<-Z
  }
  out<-do.call(cbind,mats)
  # Preserve exact train/test column alignment; no test-set feature ranking.
  colnames(out)<-make.unique(colnames(out))
  out
}
tgr_branch_rels<-function(B,y,candidate,root,branch_min=8L,rel_k=3L){
  out<-list("-1"=integer(),"0"=integer(),"1"=integer())
  for(s in c(-1L,0L,1L)){
    ix<-which(B[,root]==s)
    if(length(ix)<branch_min || min(table(factor(y[ix],levels=c(0,1))))<3L)
      next
    z<-B[ix,,drop=FALSE]
    aa<-y[ix]==1L
    contrast<-abs(colMeans(z[aa,,drop=FALSE])-
                  colMeans(z[!aa,,drop=FALSE]))
    contrast[!is.finite(contrast)]<-0
    contrast[[root]]<-0
    ord<-order(-contrast,seq_along(contrast))
    selected<-integer();genes<-c(candidate$gene_left[[root]],
                                 candidate$gene_right[[root]])
    for(j in ord){
      if(contrast[[j]]<=.03)break
      g1<-candidate$gene_left[[j]];g2<-candidate$gene_right[[j]]
      if(g1 %in% genes || g2 %in% genes || g1==g2)next
      selected<-c(selected,j);genes<-c(genes,g1,g2)
      if(length(selected)==rel_k)break
    }
    out[[as.character(s)]]<-selected
  }
  out
}
tgr_root_candidates<-function(B,y,candidate,main,max_roots=9L){
  # Root gate must be ternary. Require evidence of inactive/active states,
  # reasonable support, and not duplicate global RLR main features.
  sz<-nrow(B);rr<-integer()
  used<-main$ordered_pair
  for(j in seq_len(min(180L,ncol(B)))){
    v<-B[,j]
    counts<-table(factor(v,levels=c(-1,0,1)))
    if(sum(counts>=7L)<2L || sum(v==0)<4L || sum(v!=0)<9L)next
    if(candidate$ordered_pair[[j]] %in% used)next
    rr<-c(rr,j)
    if(length(rr)>=max_roots)break
  }
  rr
}
tgr_propose<-function(D,y,V,v_y,K=10L,seed=42L){
  aa<-tgr_allpairs(D,factor(y),seed)
  candidates<-aa$pool
  main<-apply_disjointness(aa$source_top,K=K,d=2L)
  tgr_need(nrow(main)>=2L,"global main shortlist too narrow")
  BD<-encode_pair_matrix(D,candidates,m=10L)
  roots<-tgr_root_candidates(BD,y,candidates,main,max_roots=9L)
  tgr_need(length(roots)>0L,"no valid ternary gate")
  M_D<-encode_pair_matrix(D,main,m=10L)
  M_V<-encode_pair_matrix(V,main,m=10L)
  nullfit<-tgr_enfit(M_D,y,seed)
  p0<-tgr_predict(nullfit,M_V)
  null_ll<-tgr_ll(v_y,p0)
  best<-list(objective=null_ll,tree=FALSE,
             validation_ll=null_ll,
             improvement=0,baseline_ll=null_ll,
             n_root_candidates=length(roots))
  for(root in roots){
    branches<-tgr_branch_rels(BD,y,candidates,root,
                          branch_min=8L,rel_k=3L)
    n_feature<-sum(lengths(branches))
    if(n_feature<2L)next
    A<-tgr_design(D,main,root,branches,candidates)
    A_V<-tgr_design(V,main,root,branches,candidates)
    tgr_need(identical(colnames(A),colnames(A_V)),
             "train V tree columns disagree")
    fit<-tgr_enfit(A,y,seed)
    p<-tgr_predict(fit,A_V)
    ll<-tgr_ll(v_y,p)
    objective<-ll+.001*n_feature+.004 # penalize hierarchical complexity
    if(objective<best$objective){
      best<-list(objective=objective,validation_ll=ll,
        tree=TRUE,root=root,branches=branches,n_feature=n_feature,
        improvement=null_ll-ll,baseline_ll=null_ll,
        n_root_candidates=length(roots))
    }
  }
  if(!isTRUE(best$tree) ||
     best$improvement<.008){
    best$tree<-FALSE
    best$validated<-FALSE
  }else best$validated<-TRUE
  best$main<-main
  best$candidates<-candidates
  best
}
tgr_refit<-function(T,y,E,proposal,seed=42L){
  if(!isTRUE(proposal$tree)){
    # A fallback uses original-nominated global pairs but is not used
    # as official original RLR baseline: that comes from original function.
    pool<-compute_top_pairs_switchbox(T,factor(y),top_n=1000L,
                                             buffer_mult=5L)
    main<-apply_disjointness(pool,K=10L,d=2L)
    A<-encode_pair_matrix(T,main,m=10L)
    B<-encode_pair_matrix(E,main,m=10L)
    fit<-tgr_enfit(A,y,seed)
    return(list(prob=tgr_predict(fit,B),tree_used=FALSE,
      n_extra=0L,n_main=ncol(A),mean_local_evaluated=ncol(A),
      mean_local_active=mean(rowSums(B!=0)),n_global_relations=ncol(A)))
  }
  A<-tgr_design(T,proposal$main,proposal$root,
                         proposal$branches,proposal$candidates)
  B<-tgr_design(E,proposal$main,proposal$root,
                         proposal$branches,proposal$candidates)
  tgr_need(identical(colnames(A),colnames(B)),
              "refit train/test tree design mismatch")
  fit<-tgr_enfit(A,y,seed)
  p<-tgr_predict(fit,B)
  root_val<-encode_pair_matrix(E,
            proposal$candidates[proposal$root,,drop=FALSE],m=10L)[,1]
  extra_each<-vapply(root_val,function(s)
    length(proposal$branches[[as.character(s)]]),integer(1))
  main_eval<-nrow(proposal$main)
  # Complexity is based on path-specific AND global base evidence.
  path_R<-B
  mask<-matrix(FALSE,nrow=nrow(E),ncol=ncol(B))
  mask[,seq_len(main_eval+1L)]<-TRUE
  offset<-main_eval+1L
  for(s in c(-1L,0L,1L)){
    n<-length(proposal$branches[[as.character(s)]])
    if(n){
      cols<-offset+seq_len(n)
      mask[root_val==s,cols]<-TRUE
      offset<-offset+n
    }
  }
  local_activity<-rowSums((B!=0)&mask)
  list(prob=p,tree_used=TRUE,n_extra=sum(lengths(proposal$branches)),
       n_main=main_eval,mean_local_evaluated=main_eval+1+mean(extra_each),
       mean_local_active=mean(local_activity),
       n_global_relations=main_eval+1+sum(lengths(proposal$branches)))
}
