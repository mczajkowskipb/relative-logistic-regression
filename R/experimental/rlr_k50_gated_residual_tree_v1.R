# K50-RGTree V1 — branch-gated RELATIONAL OFFSET corrections to original RLR.
# Full original RLR K50 probabilities enter as a LOCKED logistic offset.
# A single (not per leaf) elastic-net penalized correction fit uses
# three-way {-1,0,+1} root-dependent candidate relations.
# Every structure is nominated only from D; accepted by heldout inner V.
# If no improvement is verified, exact original RLR is returned.
rgt_need<-function(x,msg)if(!isTRUE(x))stop("RGT_STOP: ",msg,call.=FALSE)
rgt_clip<-function(p)pmin(1-1e-8,pmax(1e-8,p))
rgt_ll<-function(y,p) {
  p<-rgt_clip(p)
  mean(-y*log(p)-(1-y)*log1p(-p))
}
rgt_ba<-function(y,p) {
  z<-as.integer(p>=.5)
  0.5*(mean(z[y==0]==0)+mean(z[y==1]==1))
}
rgt_f1<-function(y,p) {
  z<-as.integer(p>=.5)
  mean(vapply(c(0L,1L),function(s) {
    den<-sum(z==s)+sum(y==s)
    if(!den)return(0)
    2*sum(z==s & y==s)/den
  },numeric(1)))
}
rgt_pseudologit<-function(p)qlogis(rgt_clip(p))
rgt_candidates<-function(D,dy,seed) {
  top<-compute_top_pairs_switchbox(D,dy,top_n=1000L,buffer_mult=5L)
  pairs<-top
  set.seed(seed)
  genes<-colnames(D);n<-length(genes)
  ns<-min(350L,max(180L,2L*n))
  a<-sample.int(n,ns,replace=TRUE)
  b<-sample.int(n,ns,replace=TRUE)
  x<-data.frame(gene_left=genes[a[a!=b]],
                gene_right=genes[b[a!=b]],stringsAsFactors=FALSE)
  x$ordered_pair<-paste0(x$gene_left," > ",x$gene_right)
  # Separate from original selected pairs by unique oriented name.
  unique_pairs<-rbind(pairs[,c("gene_left","gene_right","ordered_pair")],x)
  unique_pairs<-unique_pairs[!duplicated(unique_pairs$ordered_pair),,drop=FALSE]
  rownames(unique_pairs)<-NULL
  rgt_need(nrow(unique_pairs)>=300L,"insufficient candidate pairs")
  list(candidates=unique_pairs,top=top)
}
rgt_root_proposals<-function(B,pairs,parent_pairs,nroots=8L) {
  forbidden<-parent_pairs$ordered_pair
  stats<-vapply(seq_len(ncol(B)),function(j){
    counts<-as.integer(table(factor(B[,j],levels=c(-1,0,1))))
    if(sum(counts>=5L)<2L || pairs$ordered_pair[[j]] %in% forbidden)
      return(-Inf)
    x<-counts/sum(counts)
    -sum(x[x>0]*log(x[x>0])) +.005*counts[[2L]]
  },numeric(1))
  ranked<-which(is.finite(stats))
  if(!length(ranked))return(integer())
  head(ranked[order(-stats[ranked],ranked)],nroots)
}
rgt_branch_pairs<-function(B,dy,pbase,pairs,root,rparent,
                            minbranch=10L,each=3L) {
  # Root train-only selection. Select the evidence most associated with
  # residual y-p_RLR within each routed subgroup (not target alone).
  branches<-list("-1"=integer(),"0"=integer(),"1"=integer())
  for(s in c(-1L,0L,1L)) {
    ix<-which(B[,root]==s)
    if(length(ix)<minbranch || sum(dy[ix]==0)<3L ||
        sum(dy[ix]==1)<3L)next
    resid<-as.numeric(dy[ix]-pbase[ix])
    if(sd(resid)<1e-10)next
    mat<-B[ix,,drop=FALSE]
    contrast<-as.numeric(crossprod(resid,mat))/length(ix)
    contrast<-abs(contrast)
    contrast[!is.finite(contrast)]<-0
    contrast[[root]]<-0
    ord<-order(-contrast,seq_along(contrast))
    chosen<-integer()
    used<-unique(c(pairs$gene_left[[root]],pairs$gene_right[[root]]))
    for(j in ord) {
      if(contrast[[j]]<.02)break
      if(pairs$ordered_pair[[j]] %in% rparent$ordered_pair)next
      g<-c(pairs$gene_left[[j]],pairs$gene_right[[j]])
      if(g[[1L]]==g[[2L]] || any(g %in% used))next
      chosen<-c(chosen,j);used<-c(used,g)
      if(length(chosen)>=each)break
    }
    if(length(chosen)==each)
      branches[[as.character(s)]]<-chosen
  }
  branches
}
rgt_design<-function(X,pairs,root,branches) {
  needed<-sort(unique(c(root,unlist(branches,use.names=FALSE))))
  # Rank full gene universe once, as original RLR, select candidate pairs.
  R<-encode_pair_matrix(X,pairs[needed,,drop=FALSE],m=10L)
  colnames(R)<-as.character(needed)
  g<-R[,as.character(root)]
  # Three dedicated branch-indicator interaction blocks. Each feature is
  # r_{u,v}(x) * I(root(x)=s); intercept stays fixed at 0 correction.
  out<-list()
  for(s in c(-1L,0L,1L)) {
    jj<-branches[[as.character(s)]]
    if(length(jj)) {
      x<-R[,as.character(jj),drop=FALSE]*
                         as.integer(g==s)
      colnames(x)<-paste0("STATE_",s,"_PAIR_",jj)
      out[[length(out)+1L]]<-x
    }
  }
  if(!length(out))return(NULL)
  z<-do.call(cbind,out)
  if(!is.matrix(z))z<-matrix(z,ncol=1L)
  rownames(z)<-rownames(X)
  z
}
rgt_correct_fit<-function(Z,y,parent_p,lambda=.05) {
  if(is.null(Z) || ncol(Z)<2L ||
     !any(apply(Z,2,sd)>1e-12))return(NULL)
  fit<-tryCatch(suppressWarnings(glmnet::glmnet(
        x=Z,y=as.numeric(y),family="binomial",
        alpha=.5,lambda=lambda,offset=rgt_pseudologit(parent_p),
        intercept=FALSE,standardize=FALSE,maxit=50000L)),
      error=function(e)NULL)
  fit
}
rgt_correct_pred<-function(fit,Z,parent_p,lambda=.05) {
  if(is.null(fit))return(parent_p)
  p<-as.numeric(predict(fit,newx=Z,s=lambda,type="response",
         newoffset=rgt_pseudologit(parent_p)))
  rgt_need(length(p)==length(parent_p)&&
     all(is.finite(p))&&all(p>=0 & p<=1),
     "bad parent-offset correction prediction")
  p
}
rgt_parent_fit<-function(D,yD,V,yV,pool,folder) {
  fit_rlr_train_test(D,yD,V,yV,pairs_pool=pool,K=50L,
       m=10L,d=2L,out_dir=folder,alpha=.5,
       lambda_choice="lambda.min",seed=42L)
  shortlist<-apply_disjointness(pool,K=50L,d=2L)
  cv<-readRDS(file.path(folder,"cvfit.rds"))
  RD<-encode_pair_matrix(D,shortlist,m=10L)
  RV<-encode_pair_matrix(V,shortlist,m=10L)
  pD<-as.numeric(predict(cv,newx=RD,s="lambda.min",type="response"))
  pV<-as.numeric(predict(cv,newx=RV,s="lambda.min",type="response"))
  list(pD=pD,pV=pV,shortlist=shortlist,cv=cv)
}
rgt_select<-function(D,yD,V,yV,seed=42L,path) {
  # Candidate orientation & parent fitting: D ONLY.
  nom<-rgt_candidates(D,yD,seed)
  parent<-rgt_parent_fit(D,yD,V,yV,nom$top,path)
  BD<-encode_pair_matrix(D,nom$candidates,m=10L)
  roots<-rgt_root_proposals(BD,nom$candidates,
                            parent$shortlist,nroots=8L)
  null_ll<-rgt_ll(as.integer(yV==levels(yV)[2L]),parent$pV)
  null_ba<-rgt_ba(as.integer(yV==levels(yV)[2L]),parent$pV)
  best<-list(root=NA_integer_,branches=NULL,
             gain=0,BAchange=0,validation_ll=null_ll,
             n_nonzero_corrections=0L)
  diagnostic<-list(branch_designs_valid=0L,correction_fits_valid=0L,
    corrections_with_nonzero_coefficients=0L,
    candidates_with_positive_validation_gain=0L,
    max_unpenalized_validation_gain=NA_real_)
  for(root in roots) {
    branches<-rgt_branch_pairs(BD,
      as.integer(yD==levels(yD)[2L]),parent$pD,
      nom$candidates,root,parent$shortlist,
      minbranch=10L,each=3L)
    count<-sum(lengths(branches))
    if(count<3L)next
    diagnostic$branch_designs_valid<-
      diagnostic$branch_designs_valid+1L
    A<-rgt_design(D,nom$candidates,root,branches)
    B<-rgt_design(V,nom$candidates,root,branches)
    rgt_need(identical(colnames(A),colnames(B)),
             "D/V branch pair design mismatch")
    fit<-rgt_correct_fit(A,as.integer(yD==levels(yD)[2L]),
                          parent$pD,lambda=.05)
    if(is.null(fit))next
    diagnostic$correction_fits_valid<-
      diagnostic$correction_fits_valid+1L
    pr<-rgt_correct_pred(fit,B,parent$pV,lambda=.05)
    gain<-null_ll-rgt_ll(as.integer(yV==levels(yV)[2L]),pr)
    ba_change<-rgt_ba(as.integer(yV==levels(yV)[2L]),pr)-null_ba
    bb<-as.numeric(stats::coef(fit,s=.05))[-1L]
    nnz<-sum(abs(bb)>1e-9)
    if(nnz>0L)diagnostic$corrections_with_nonzero_coefficients<-
       diagnostic$corrections_with_nonzero_coefficients+1L
    if(gain>0)diagnostic$candidates_with_positive_validation_gain<-
       diagnostic$candidates_with_positive_validation_gain+1L
    diagnostic$max_unpenalized_validation_gain<-
       if(is.na(diagnostic$max_unpenalized_validation_gain))gain else
         max(diagnostic$max_unpenalized_validation_gain,gain)
    # Fixed predeclared validation model complexity cost.
    objective<-gain-.0015*nnz
    if(nnz>=1L && ba_change>=(-.02) &&
       objective>(best$gain-.0015*best$n_nonzero_corrections)){
      best<-list(root=root,branches=branches,gain=gain,
            BAchange=ba_change,validation_ll=null_ll-gain,
            n_nonzero_corrections=nnz)
    }
  }
  # Require a meaningful absolute V logloss improvement; do not
  # use the outer E to pick or optimize any model parameters.
  best$accept<-is.finite(best$root) &&
     best$gain>=.01 && best$BAchange>=(-.02)
  best$baseline_ll<-null_ll
  best$root_nominated<-length(roots)
  best$diagnostic<-diagnostic
  best$candidates<-nom$candidates
  best
}
rgt_refit<-function(T,yT,E,original_pool,selected,baseline,
                   original_cv) {
  if(!isTRUE(selected$accept))
    return(list(p=baseline,used_tree=FALSE,
         n_correction_nonzero=0L,global_extra=0L,
         path_extra_mean=0,root_active_rate=NA_real_,
         correction_delta_max=0))
  root<-selected$root
  bs<-selected$branches
  A<-rgt_design(T,selected$candidates,root,bs)
  B<-rgt_design(E,selected$candidates,root,bs)
  rgt_need(identical(colnames(A),colnames(B)),
           "train/test interaction columns mismatch")
  shortlist<-apply_disjointness(original_pool,K=50L,d=2L)
  RT<-encode_pair_matrix(T,shortlist,m=10L)
  pT<-as.numeric(predict(original_cv,newx=RT,
                         s="lambda.min",type="response"))
  fit<-rgt_correct_fit(A,
     as.integer(yT==levels(yT)[2L]),pT,lambda=.05)
  if(is.null(fit))return(list(p=baseline,used_tree=FALSE,
    n_correction_nonzero=0L,global_extra=0L,
    path_extra_mean=0,root_active_rate=NA_real_,
    correction_delta_max=0))
  pred<-rgt_correct_pred(fit,B,baseline,lambda=.05)
  bb<-as.numeric(stats::coef(fit,s=.05))[-1L]
  rootR<-encode_pair_matrix(E,selected$candidates[root,,drop=FALSE],m=10L)
  rx<-as.integer(rootR[,1L])
  path_count<-vapply(rx,function(s)length(bs[[as.character(s)]]),integer(1))
  list(p=pred,used_tree=TRUE,
     n_correction_nonzero=sum(abs(bb)>1e-9),
     global_extra=sum(lengths(bs)),
     path_extra_mean=mean(path_count)+1,
     root_active_rate=mean(rx!=0),
     correction_delta_max=max(abs(pred-baseline)))
}
