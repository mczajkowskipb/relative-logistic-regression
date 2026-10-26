# Isolated RLR + sparse two-leaf model pilot (2 or 5 local rank relations).
# Unmodified original RLR: rlmt4_parent() / rlmt4_strict_oof().
# D selects gate/features; V once confirms each prespecified method; heldout E
# is NEVER used to choose methods. Diagnostic on EXPOSED CompCancer sources only.
sparse_cfg <- function() list(
  seed=42L,k=c(2L,5L),min_node=12L,min_class_node=3L,
  lambda_intercept=.08,lambda_slope=.25,
  confirmation_cost=.002,confirmation_min_gain=.002
)
sparse_need <- function(ok,msg) if(!isTRUE(ok))
   stop("SPARSE_PILOT_STOP: ",msg,call.=FALSE)
sparse_prob <- function(p) pmax(1e-8,pmin(1-1e-8,p))
sparse_logit <- function(p) stats::qlogis(sparse_prob(p))
sparse_ll <- function(y,prob) {
  a<-sparse_prob(prob)
  -(y*log(a)+(1-y)*log1p(-a))
}
sparse_fit <- function(y,p_offset,Z,n_intercepts,cfg=sparse_cfg()) {
  sparse_need(nrow(Z)==length(y)&&length(p_offset)==length(y),
              "design sample dimensions")
  sparse_need(ncol(Z)>=1L && all(is.finite(Z)) && 
              all(is.finite(p_offset)) && all(is.finite(y)),
              "nonfinite fit")
  off<-sparse_logit(p_offset)
  penalty<-c(rep(cfg$lambda_intercept,n_intercepts),
             rep(cfg$lambda_slope,ncol(Z)-n_intercepts))
  obj<-function(b) {
    eta<-as.numeric(off+Z%*%b)
    mean(pmax(eta,0)+log1p(exp(-abs(eta)))-y*eta)+
      sum(penalty*b*b)/2
  }
  grad<-function(b) {
    eta<-as.numeric(off+Z%*%b)
    as.numeric(crossprod(Z,plogis(eta)-y)/length(y)+penalty*b)
  }
  o<-stats::optim(rep(0,ncol(Z)),obj,grad,method="BFGS",
                  control=list(maxit=350L,reltol=1e-10))
  sparse_need(o$convergence==0L&&all(is.finite(o$par)),
              "ridge leaf/global numerical convergence")
  o$par
}
sparse_predict <- function(prob,Z,b) {
  sparse_need(ncol(Z)==length(b)&&nrow(Z)==length(prob),
              "prediction dimensions")
  as.numeric(stats::plogis(sparse_logit(prob)+
                            as.numeric(Z%*%b)))
}
sparse_y <- function(y) as.numeric(y==levels(y)[2L])
sparse_design <- function(X, gatepair, features_a, features_b,
                          kind) {
  n<-nrow(X)
  rel<-function(j) ifelse(X[[j$gene_left[[1L]]]]>X[[j$gene_right[[1L]]]],1,-1)
  getmat<-function(pairs) {
    if(!nrow(pairs))return(matrix(numeric(0),nrow=n,ncol=0))
    encode_pair_matrix(X,pairs,m=10L)
  }
  if(kind=="GLOBAL") {
    A<-getmat(features_a)
    return(cbind(intercept=rep(1,n),A))
  }
  g<-rel(gatepair)
  a<-as.numeric(g<0)
  b<-1-a
  A<-getmat(features_a)
  B<-getmat(features_b)
  cbind(left_intercept=a,right_intercept=b,
        sweep(A,1L,a,"*"),sweep(B,1L,b,"*"))
}
sparse_select <- function(D,y,poo,pd,cfg=sparse_cfg()) {
  all<-pd$pairs
  sparse_need(nrow(all)>=6L, "original D parent too few gene pairs")
  R<-encode_pair_matrix(D,all,m=10L)
  gates<-vapply(seq_len(nrow(all)),function(i){
    g<-D[[all$gene_left[[i]]]]>D[[all$gene_right[[i]]]]
    ok<-sum(g)>=cfg$min_node&&sum(!g)>=cfg$min_node&&
       all(table(factor(y[g],levels=levels(y)))>=cfg$min_class_node) &&
       all(table(factor(y[!g],levels=levels(y)))>=cfg$min_class_node)
    if(!ok)return(-Inf)
    r<-sparse_y(y)-poo
    abs(mean(r[g])-mean(r[!g]))*sqrt(sum(g)*sum(!g)/length(g))
  },numeric(1))
  sparse_need(any(is.finite(gates)),"no balanced rank gate in D")
  gate_index<-which.max(gates)
  gate<-all[gate_index,,drop=FALSE]
  g<-D[[gate$gene_left]]>D[[gate$gene_right]]
  r<-sparse_y(y)-poo
  candidates<-seq_len(nrow(all))
  candidates<-setdiff(candidates,gate_index)
  score<-function(indices) {
    vapply(candidates,function(j){
      z<-R[indices,j]
      if(stats::sd(z)==0)return(-Inf)
      abs(sum((r[indices]-mean(r[indices]))*(z-mean(z))))/
          sqrt(sum((z-mean(z))^2)+1e-8)
    },numeric(1))
  }
  a<-candidates[order(-score(which(!g)),candidates)]
  b<-candidates[order(-score(which(g)),candidates)]
  allsc<-vapply(candidates,function(j){
    z<-R[,j]
    if(stats::sd(z)==0)return(-Inf)
    abs(sum((r-mean(r))*(z-mean(z))))/
        sqrt(sum((z-mean(z))^2)+1e-8)
  },numeric(1))
  global<-candidates[order(-allsc,candidates)]
  sparse_need(length(a)>=5&&length(b)>=5&&length(global)>=5,
              "not enough leaf feature candidates")
  list(gate=gate,gate_index=gate_index,
       left=all[a,,drop=FALSE],right=all[b,,drop=FALSE],
       global=all[global,,drop=FALSE],
       leaf_sizes=c(left=sum(!g),right=sum(g)),
       original_pair_count=nrow(all))
}
sparse_spec <- function(selection,k,variant) {
  if(variant=="GLOBAL") {
    return(list(kind="GLOBAL",k=k,independent=FALSE,
          gate=NULL,fa=selection$global[seq_len(k),,drop=FALSE],
          fb=selection$global[FALSE,,drop=FALSE]))
  }
  list(kind="TREE",k=k,independent=(variant=="INDEPENDENT"),
       gate=selection$gate,
       fa=selection$left[seq_len(k),,drop=FALSE],
       fb=selection$right[seq_len(k),,drop=FALSE])
}
sparse_one <- function(T,y,E,methods=c("GLOBAL_2","GLOBAL_5",
                                       "SHARED_2","SHARED_5",
                                       "INDEPENDENT_2","INDEPENDENT_5"),
                       parent_cfg=rlmt4_config(),
                       cfg=sparse_cfg()) {
  sparse_need(identical(names(T),names(E))&&
              length(levels(y))==2L && nrow(T)==length(y),
              "outer X/y schema")
  original<-rlmt4_parent(T,y,parent_cfg)
  pb<-rlmt4_parent_predict(original,E)$prob
  foo<-rlmt4_strict_oof(T,y,parent_cfg)
  sparse_need(!is.null(foo),"outer original RLR strict OOF failed")
  folds<-make_repeated_stratified_folds(y,n_folds=3L,
                                       n_repeats=1L,
                                       seed=cfg$seed+731L)
  d<-folds[[1L]]$train_idx
  v<-folds[[1L]]$test_idx
  sparse_need(!length(intersect(d,v))&&
              identical(sort(c(d,v)),seq_len(nrow(T))),
              "D/V overlap or coverage")
  D<-T[d,,drop=FALSE]
  V<-T[v,,drop=FALSE]
  pd<-rlmt4_parent(D,y[d],parent_cfg)
  dood<-rlmt4_strict_oof(D,y[d],parent_cfg)
  sparse_need(!is.null(dood),"D-only original RLR strict OOF failed")
  pv<-rlmt4_parent_predict(pd,V)$prob
  sel<-sparse_select(D,y[d],dood$prob,pd,cfg)
  output<-list(ORIGINAL_RLR=list(prob=pb,accepted=NA,
                                v_mean=NA_real_,v_lower=NA_real_,
                                reason="unchanged_exact_parent"))
  outmeta<-list(root_pair=paste(sel$gate$gene_left,sel$gate$gene_right,
                                sep=">"),
                original_active_pair_budget=parent_cfg$K,
                leaf_sizes=sel$leaf_sizes,
                train_ids=rownames(T),E_ids=rownames(E),
                D_ids=rownames(D),V_ids=rownames(V),
                candidate_original_pairs=sel$original_pair_count)
  for(method in methods) {
    part<-strsplit(method,"_",fixed=TRUE)[[1L]]
    variant<-part[[1L]]
    k<-as.integer(part[[2L]])
    sparse_need(k %in% cfg$k &&
                variant %in% c("GLOBAL","SHARED","INDEPENDENT"),
                "unknown arm")
    spec<-sparse_spec(sel,k,variant)
    Zd<-sparse_design(D,spec$gate,spec$fa,spec$fb,spec$kind)
    Zv<-sparse_design(V,spec$gate,spec$fa,spec$fb,spec$kind)
    Zt<-sparse_design(T,spec$gate,spec$fa,spec$fb,spec$kind)
    Ze<-sparse_design(E,spec$gate,spec$fa,spec$fb,spec$kind)
    n_i<-if(variant=="GLOBAL")1L else 2L
    po_d<-if(spec$independent)rep(.5,length(d)) else dood$prob
    po_v<-if(spec$independent)rep(.5,length(v)) else pv
    po_t<-if(spec$independent)rep(.5,nrow(T)) else foo$prob
    po_e<-if(spec$independent)rep(.5,nrow(E)) else pb
    trained<-tryCatch({
       fit_d<-sparse_fit(sparse_y(y[d]),po_d,Zd,n_i,cfg)
       pv_candidate<-sparse_predict(po_v,Zv,fit_d)
       gains<-sparse_ll(sparse_y(y[v]),pv)-
                 sparse_ll(sparse_y(y[v]),pv_candidate)
       lower<-mean(gains)-stats::sd(gains)/sqrt(length(gains))-
                 cfg$confirmation_cost
       pass<-is.finite(lower) && lower>cfg$confirmation_min_gain
       prob<-pb
       if(pass) {
         fit_all<-sparse_fit(sparse_y(y),po_t,Zt,n_i,cfg)
         prob<-sparse_predict(po_e,Ze,fit_all)
       }
       list(prob=prob,accepted=pass,v_mean=mean(gains),
            v_lower=lower,
            reason=if(pass)"confirmed_on_V" else "rejected_on_V",
            fit_dims=c(ncol(Zd),n_i,k),
            model_pairs=list(gate=spec$gate,fa=spec$fa,fb=spec$fb))
    },error=function(e) {
       list(prob=pb,accepted=FALSE,v_mean=NA_real_,
            v_lower=NA_real_,
            reason=paste0("fitting_failure:",conditionMessage(e)),
            fit_dims=c(ncol(Zd),n_i,k),
            model_pairs=list(gate=spec$gate,fa=spec$fa,fb=spec$fb))
    })
    sparse_need(length(trained$prob)==nrow(E) &&
                all(is.finite(trained$prob)) &&
                all(trained$prob>=0&trained$prob<=1),
                "invalid heldout prediction vector")
    if(!isTRUE(trained$accepted))
      sparse_need(identical(pb,trained$prob),
                  "baseline fallback must be bit-exact")
    output[[method]]<-trained
  }
  list(preds=output,meta=outmeta)
}
