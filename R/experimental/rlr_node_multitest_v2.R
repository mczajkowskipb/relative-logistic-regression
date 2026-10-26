# CRIT publication gate: interpretably sparse RLR multivariate split tests (depth ONE).
# Frozen original RLR is EXACT rlmt4_parent / rlmt4_strict_oof.
# RLR is in NODE, not independent RLR in leaf; 2 terminal leaves each
# receive one strongly shrunk offset correction to global original RLR.
# V2 exploratory: allow near-pure terminal classes; keep min 12 patients per leaf.
# D-only node selection, one V confirmation, E untouched.
node_cfg<-function()list(seed=42L,alpha=.5,
  min_node=12L,min_each_class=0L,threshold_q=c(.20,.35,.50,.65,.80),
  ridge=.08,cost=.002,min_gain=.002)
node_need<-function(x,why)if(!isTRUE(x))
  stop("NODE_STOP: ",why,call.=FALSE)
node_ll<-function(y,p){
  p<-pmin(1-1e-8,pmax(1e-8,p))
  -y*log(p)-(1-y)*log1p(-p)
}
node_logit<-function(p)stats::qlogis(pmin(1-1e-8,pmax(1e-8,p)))
node_sigmoid<-function(eta)stats::plogis(eta)
node_y01<-function(y)as.numeric(y==levels(y)[2L])
node_leaf_offsets<-function(y,p,route,ridge=.08){
  node_need(length(y)==length(p)&&length(route)==length(y)&&
              all(is.finite(p)),"leaf coefficient lengths/nonfinite")
  M<-cbind(left=as.numeric(!route),right=as.numeric(route))
  offset<-node_logit(p)
  objective<-function(b) {
    eta<-as.numeric(offset+M%*%b)
    mean(pmax(eta,0)+log1p(exp(-abs(eta)))-y*eta)+
      ridge*sum(b*b)/2
  }
  gradient<-function(b) {
    eta<-as.numeric(offset+M%*%b)
    as.numeric(crossprod(M,node_sigmoid(eta)-y)/length(y)+ridge*b)
  }
  fit<-stats::optim(c(0,0),objective,gradient,method="BFGS",
                    control=list(maxit=300L,reltol=1e-10))
  node_need(fit$convergence==0L&&all(is.finite(fit$par)),
            "leaf ridge numerical convergence")
  fit$par
}
node_apply_offsets<-function(p,route,b){
  node_need(length(p)==length(route)&&length(b)==2L,"route/apply")
  as.numeric(node_sigmoid(node_logit(p)+ifelse(route,b[[2L]],b[[1L]])))
}
node_route<-function(X,model){
  if(model$kind=="SINGLE") {
    a<-model$pair$gene_left[[1L]]
    b<-model$pair$gene_right[[1L]]
    val<-as.numeric(X[[a]]>X[[b]])
  }else{
    R<-encode_pair_matrix(X,model$pairs,m=10L)
    node_need(ncol(R)==length(model$coef),"multi-node pair count")
    val<-as.numeric(model$intercept+R%*%model$coef)
  }
  node_need(all(is.finite(val)),"route score nonfinite")
  val>model$threshold
}
node_supported<-function(y,route,cfg){
  if(sum(route)<cfg$min_node||sum(!route)<cfg$min_node)return(FALSE)
  tab1<-table(factor(y[route],levels=levels(y)))
  tab0<-table(factor(y[!route],levels=levels(y)))
  all(tab1>=cfg$min_each_class)&&all(tab0>=cfg$min_each_class)
}
node_pick_threshold<-function(X,y,p_d_oof,model,cfg){
  score<-if(model$kind=="SINGLE"){
    as.numeric(X[[model$pair$gene_left[[1L]]]]>
                 X[[model$pair$gene_right[[1L]]]])
  }else{
    R<-encode_pair_matrix(X,model$pairs,m=10L)
    as.numeric(model$intercept+R%*%model$coef)
  }
  options<-if(model$kind=="SINGLE")0.5 else
    unique(as.numeric(stats::quantile(score,cfg$threshold_q,names=FALSE,
                                      type=7)))
  record<-list()
  for(t in options){
    route<-score>t
    if(!node_supported(y,route,cfg))next
    b<-node_leaf_offsets(node_y01(y),p_d_oof,route,cfg$ridge)
    candidate<-node_apply_offsets(p_d_oof,route,b)
    gain<-mean(node_ll(node_y01(y),p_d_oof)-
               node_ll(node_y01(y),candidate))
    record[[length(record)+1L]]<-list(
      threshold=t,D_gain=gain,D_left=sum(!route),D_right=sum(route),
      D_b=b)
  }
  if(!length(record))return(NULL)
  best<-which.max(vapply(record,function(z)z$D_gain,numeric(1)))
  record[[best]]
}
node_candidates<-function(D,y,p_d,parentD,kind,k,cfg){
  if(kind=="SINGLE"){
    # Broaden the possible root relation pool: its top marginal members
    # alone can form class-pure gates with too few members in both classes.
    pool<-crit7_pool(D,parentD$pairs,cfg=crit7_cfg())
    pool$gene_left<-pool$gene_a
    pool$gene_right<-pool$gene_b
    node_need(nrow(pool)>=1L,"empty D single-root comparison pool")
    val<-list()
    for(i in seq_len(nrow(pool))) {
      m<-list(kind="SINGLE",pair=pool[i,,drop=FALSE],
              threshold=.5)
      t<-node_pick_threshold(D,y,p_d,m,cfg)
      if(is.null(t))next
      val[[length(val)+1L]]<-list(model=m,choice=t,
                                  n_pool=nrow(pool))
    }
    return(val)
  }
  # RLR node: top k pairs among selected original switchBox shortlist.
  # Each small multivariate node is itself fitted as a sparse logistic
  # model on D. The global parent is not replaced or retrained.
  pairset<-parentD$pairs
  node_need(nrow(pairset)>=k,"insufficient original RLR relations")
  allR<-encode_pair_matrix(D,pairset,m=10L)
  fit_beta<-as.numeric(stats::coef(parentD$fit,
                s=parentD$config$lambda_choice))[-1L]
  node_need(length(fit_beta)==nrow(pairset),"coefficient / pair alignment")
  selection<-order(-abs(fit_beta),seq_along(fit_beta))[seq_len(k)]
  R<-allR[,selection,drop=FALSE]
  # Fixed regularization class = elastic-net; lambda chosen within D.
  set.seed(cfg$seed+1000L+k)
  small<-glmnet::cv.glmnet(x=R,y=node_y01(y),family="binomial",
      alpha=cfg$alpha,nfolds=safe_nfolds(y,3L),
      type.measure="deviance",standardize=FALSE)
  co<-as.numeric(stats::coef(small,s="lambda.min"))
  m<-list(kind="MULTI",pairs=pairset[selection,,drop=FALSE],
          intercept=co[[1L]],coef=co[-1L],threshold=NA_real_,
          k=k,alpha=cfg$alpha,lambda=small$lambda.min,
          fitted_nonzero=sum(abs(co[-1L])>1e-9))
  candidate<-node_pick_threshold(D,y,p_d,m,cfg)
  if(is.null(candidate))return(list())
  m$threshold<-candidate$threshold
  list(list(model=m,choice=candidate,n_pool=nrow(pairset)))
}
node_run<-function(T,y,E,pcfg=rlmt4_config(),cfg=node_cfg()) {
  node_need(is.factor(y)&&length(levels(y))==2L&&
       nrow(T)==length(y)&&identical(names(T),names(E)),
       "outer data / classes")
  full_parent<-rlmt4_parent(T,y,pcfg)
  outer_parent_E<-rlmt4_parent_predict(full_parent,E)$prob
  full_oof<-rlmt4_strict_oof(T,y,pcfg)
  node_need(!is.null(full_oof),"full training strict OOF missing")
  fold<-make_repeated_stratified_folds(
    y,n_folds=3L,n_repeats=1L,seed=cfg$seed+731L)[[1L]]
  d<-fold$train_idx;v<-fold$test_idx
  node_need(length(intersect(d,v))==0 &&
       identical(sort(c(d,v)),seq_len(nrow(T))),"D/V patient partition overlap")
  D<-T[d,,drop=FALSE];V<-T[v,,drop=FALSE]
  parent_d<-rlmt4_parent(D,y[d],pcfg)
  po_d<-rlmt4_strict_oof(D,y[d],pcfg)
  node_need(!is.null(po_d),"inner D strict original-RLR OOF missing")
  po_v<-rlmt4_parent_predict(parent_d,V)$prob
  schemes<-list(SINGLE_REL=c("SINGLE",1L),
                RLR_NODE_2=c("MULTI",2L),
                RLR_NODE_5=c("MULTI",5L))
  outcomes<-list(ORIGINAL_RLR=list(
      prob=outer_parent_E,accepted=NA,reason="unchanged_original_RLR",
      v_mean=NA_real_,v_lower=NA_real_,node_relations=NA_integer_,
      nonzero_node_coeff=NA_integer_,gate="",D_left=NA_integer_,
      D_right=NA_integer_))
  for(method in names(schemes)){
    spec<-schemes[[method]]
    nodes<-node_candidates(D,y[d],po_d$prob,parent_d,
                            kind=spec[[1L]],k=as.integer(spec[[2L]]),cfg)
    if(!length(nodes)){
      outcomes[[method]]<-list(
        prob=outer_parent_E,accepted=FALSE,reason="no_supported_D_node",
        v_mean=NA_real_,v_lower=NA_real_,
        node_relations=as.integer(spec[[2L]]),
        nonzero_node_coeff=NA_integer_,gate="",D_left=NA_integer_,
        D_right=NA_integer_)
      next
    }
    # D chooses one candidate node only. Single-rel pool may contain
    # several hundred; no V labels involved in its ranking.
    best<-which.max(vapply(nodes,function(e)e$choice$D_gain,numeric(1)))
    selected<-nodes[[best]]
    m<-selected$model
    routed_D<-node_route(D,m)
    routed_V<-node_route(V,m)
    bv<-node_leaf_offsets(node_y01(y[d]),po_d$prob,routed_D,
                          cfg$ridge)
    p_v_candidate<-node_apply_offsets(po_v,routed_V,bv)
    v_gain<-node_ll(node_y01(y[v]),po_v)-
            node_ll(node_y01(y[v]),p_v_candidate)
    lower<-mean(v_gain)-stats::sd(v_gain)/sqrt(length(v_gain))-cfg$cost
    accepted<-is.finite(lower)&&lower>cfg$min_gain
    pred_E<-outer_parent_E
    if(accepted){
      # Gate remains trained ONLY on D. Refit 2 simple leaf offsets on
      # ALL outer training T, using strict OOF original-RLR predictions.
      # No E data/labels in any fitted coefficient or threshold.
      bT<-node_leaf_offsets(node_y01(y),full_oof$prob,node_route(T,m),
                            cfg$ridge)
      pred_E<-node_apply_offsets(outer_parent_E,node_route(E,m),bT)
    }
    if(!accepted)
      node_need(identical(pred_E,outer_parent_E),
                "exact unchanged RLR fallback violated")
    node_need(all(is.finite(pred_E))&&all(pred_E>=0&pred_E<=1),
              "invalid E probability")
    outcomes[[method]]<-list(prob=pred_E,accepted=accepted,
      reason=if(accepted)"V_CONFIRMED" else "V_REJECTED",
      v_mean=mean(v_gain),v_lower=lower,
      node_relations=if(m$kind=="SINGLE")1L else nrow(m$pairs),
      nonzero_node_coeff=if(m$kind=="SINGLE")1L else m$fitted_nonzero,
      gate=if(m$kind=="SINGLE"){
        paste(m$pair$gene_left[[1L]],m$pair$gene_right[[1L]],sep=">")
      }else{
        paste(sprintf("%s>%s",m$pairs$gene_left,m$pairs$gene_right),
              collapse=";")
      },
      D_left=sum(!routed_D),D_right=sum(routed_D))
  }
  list(preds=outcomes,meta=list(
    D_ids=rownames(D),V_ids=rownames(V),E_ids=rownames(E),
    source_train_ids=rownames(T),original_RLR_pair_count=nrow(full_parent$pairs),
    parent_identity=full_parent$model_identity))
}
