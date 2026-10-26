# Rank-invariant RLR-gated stump V1. Development experiment only.
# Source original RLR separately (R/new/07_world_protocol_utils.R).
# No raw abundance, ratios, cross-sample absolute cutoffs, or E-based selection.

rg_need <- function(ok, message) if (!isTRUE(ok)) stop("RANK_GATE_STOP: ", message, call.=FALSE)

rg_ranks <- function(X, schema=NULL) {
  rg_need(is.data.frame(X) && nrow(X)>0L && ncol(X)>1L,
          "numeric nonempty data.frame required")
  rg_need(all(vapply(X, function(z) is.numeric(z)&&all(is.finite(z)),logical(1))),
          "invalid/nonfinite gene values")
  rg_need(!anyDuplicated(names(X)), "duplicate gene identifiers")
  if (!is.null(schema)) rg_need(identical(names(X),schema),"schema/order mismatch")
  r <- t(vapply(seq_len(nrow(X)), function(i)
    rank(unlist(X[i,,drop=FALSE],use.names=FALSE), ties.method="average"),
    numeric(ncol(X))))
  dimnames(r) <- list(rownames(X), names(X))
  r
}

rg_pairs <- function(X, pool, max_pool=650L, max_extra=700L) {
  base <- pool[seq_len(min(nrow(pool),as.integer(max_pool))),
               c("gene_left","gene_right"),drop=FALSE]
  # An explicit D-only non-TSP candidate expansion: high-MAD genes.
  # MAD is used for feature SELECTION on training samples, never as prediction input.
  mad <- vapply(X, stats::mad, numeric(1))
  top <- names(sort(mad,decreasing=TRUE))[seq_len(min(48L,ncol(X)))]
  extra <- if (length(top)>1L) {
    cmb <- utils::combn(top,2L)
    data.frame(gene_left=cmb[1,],gene_right=cmb[2,],
               stringsAsFactors=FALSE)
  } else base[0,,drop=FALSE]
  if (nrow(extra)>max_extra) {
    set.seed(42L)
    extra <- extra[sort(sample.int(nrow(extra),max_extra)),,drop=FALSE]
  }
  pairs <- rbind(base,extra)
  key <- vapply(seq_len(nrow(pairs)), function(i)
    paste(sort(c(pairs$gene_left[i],pairs$gene_right[i])),collapse="|"),
    character(1))
  pairs <- pairs[!duplicated(key),,drop=FALSE]
  pairs$ordered_pair <- paste0(pairs$gene_left," > ",pairs$gene_right)
  rownames(pairs)<-NULL
  rg_need(nrow(pairs)>2L,"candidate pool too small")
  pairs
}

rg_gaps <- function(ranks,pairs) {
  rg_need(all(pairs$gene_left %in% colnames(ranks)) &&
            all(pairs$gene_right %in% colnames(ranks)),
          "missing gate/partner gene")
  z <- vapply(seq_len(nrow(pairs)),function(i)
    as.numeric(ranks[,pairs$gene_right[i]]-ranks[,pairs$gene_left[i]]),
    numeric(nrow(ranks)))
  if (is.null(dim(z))) z <- matrix(z,ncol=1L)
  dimnames(z) <- list(rownames(ranks),pairs$ordered_pair)
  z
}

rg_gate <- function(gap, mode, n_genes) {
  rg_need(mode %in% c("ORDER","WIDE"), "unknown fixed gate type")
  if (mode=="ORDER") as.numeric(gap>0) else
    as.numeric(gap >= max(10,ceiling(.15*(n_genes-1))))
}
rg_ternary <- function(gap) ifelse(gap>=10,1,ifelse(gap<=-10,-1,0))
rg_expit <- function(z) plogis(pmin(30,pmax(-30,z)))
rg_ll <- function(y,p) {
  q<-pmin(1-1e-15,pmax(1e-15,p))
  mean(-y*log(q)-(1-y)*log1p(-q))
}
rg_ba <- function(y,p) {
  predicted<-as.integer(p>=.5)
  mean(c(mean(predicted[y==0L]==0L),mean(predicted[y==1L]==1L)))
}

rg_parent <- function(D,y,E,pool,K=20L) {
  rg_need(identical(names(D),names(E)),"parent schema mismatch")
  rg_need(nrow(D)==length(y) && length(unique(as.character(y)))==2L,
          "parent training labels invalid")
  pairs <- apply_disjointness(pool,K=K,d=2L)
  rg_need(nrow(pairs)>1L,"original RLR shortlist too short")
  train <- encode_pair_matrix(D,pairs,m=10L)
  test <- encode_pair_matrix(E,pairs,m=10L)
  set.seed(42L)
  cv <- glmnet::cv.glmnet(x=train,y=as.numeric(y==levels(y)[2L]),
         family="binomial",alpha=.5,nfolds=safe_nfolds(y,5L),
         type.measure="deviance")
  list(train=as.numeric(predict(cv,newx=train,s="lambda.min",type="response")),
       test=as.numeric(predict(cv,newx=test,s="lambda.min",type="response")),
       fit=cv, selected=pairs)
}

# Search for conditional residual evidence. All label use occurs on D, never V/E.
# Scores are screening approximations, NOT unbiased prospective estimates.
rg_candidates <- function(X,y,prob,pool,max_root=48L,max_return=8L) {
  rg_need(nrow(X)==length(y)&&length(prob)==length(y),"D alignment")
  rk <- rg_ranks(X)
  pairs<-rg_pairs(X,pool)
  delta<-rg_gaps(rk,pairs)
  Z<-rg_ternary(delta)
  rg_need(all(Z %in% c(-1,0,1)),"invalid ternary matrix")
  # Roots sampled across the ranked pool; not only first few global winners.
  primary<-unique(pmax(1L,as.integer(round(seq(1,min(650L,nrow(pairs)),
                                                   length.out=max_root)))))
  y01<-as.integer(y==levels(y)[2L])
  e<-y01-prob
  w<-pmax(1e-5,prob*(1-prob))
  out<-list(); j<-0L
  for(k in primary) for(mode in c("ORDER","WIDE")) {
    g<-rg_gate(delta[,k],mode,ncol(X))
    a<-which(g==0); b<-which(g==1)
    if (min(length(a),length(b))<12L) next
    if (min(table(factor(y01[a],levels=0:1)),
            table(factor(y01[b],levels=0:1)))<3L) next
    interaction <- Z * (g-mean(g))
    numer <- as.numeric(crossprod(e,interaction))^2
    denom <- as.numeric(crossprod(w,interaction^2))+1
    scores<-numer/denom
    scores[k]<- -Inf
    if (!any(is.finite(scores))) next
    partner<-which.max(scores)
    j<-j+1L
    out[[j]]<-data.frame(
      root=k,partner=partner,gate=mode,score=scores[partner],
      root_n_left=length(a),root_n_right=length(b),
      stringsAsFactors=FALSE)
  }
  rg_need(j>0L, "no balanced rank root; use exact RLR fallback")
  cand<-do.call(rbind,out)
  cand<-cand[order(-cand$score,cand$root,cand$partner,cand$gate),,drop=FALSE]
  cand<-head(cand,max_return)
  rownames(cand)<-NULL
  list(candidates=cand,pairs=pairs)
}

rg_design <- function(X,pairs,root,partner,gate) {
  rk<-rg_ranks(X)
  d<-rg_gaps(rk,pairs[c(root,partner),,drop=FALSE])
  g<-rg_gate(d[,1],gate,ncol(X))
  z<-rg_ternary(d[,2])
  cbind(root=g,partner=z,interaction=g*z)
}

# Ridge logistic correction with FROZEN parent logit offset.
# Equal strength penalty and identical predictors in the global control,
# except the contextual interaction, which is penalized double.
rg_correct <- function(Ztrain,y,base_train,Ztest,base_test,tree=TRUE) {
  cols<-if(tree) 1:3 else 1:2
  x<-Ztrain[,cols,drop=FALSE]
  ex<-Ztest[,cols,drop=FALSE]
  rg_need(all(is.finite(x)) && all(is.finite(ex)), "invalid design")
  f<-glmnet::glmnet(x=x,y=as.integer(y==levels(y)[2L]),
         family="binomial",alpha=0,lambda=.10,standardize=FALSE,
         penalty.factor=if(tree)c(1,1,2) else c(1,1),
         intercept=FALSE,offset=qlogis(pmin(1-1e-10,pmax(1e-10,base_train))))
  as.numeric(predict(f,newx=ex,s=.10,type="response",
       newoffset=qlogis(pmin(1-1e-10,pmax(1e-10,base_test)))))
}

# Split outer T into inner D/V with class support, without touching outer E.
rg_inner <- function(y) {
  set.seed(420L)
  vi<-unlist(lapply(split(seq_along(y),y), function(ix)
    sample(ix,max(3L,as.integer(floor(.25*length(ix)))))),
    use.names=FALSE)
  list(D=setdiff(seq_along(y),vi),V=sort(vi))
}

# D-selected candidates are V-tested against BOTH base RLR and same-pair
# global non-tree correction. A success must add contextual (not merely main)
# benefit, and cannot reduce balanced accuracy by more than 0.02.
rg_select <- function(T,yT) {
  v<-rg_inner(yT)
  D<-T[v$D,,drop=FALSE]; V<-T[v$V,,drop=FALSE]
  yD<-factor(as.character(yT[v$D]),levels=levels(yT))
  yV<-factor(as.character(yT[v$V]),levels=levels(yT))
  poolD<-compute_top_pairs_switchbox(D,yD,top_n=1000L,buffer_mult=5L)
  parent<-rg_parent(D,yD,V,poolD,K=20L)
  yV01<-as.integer(yV==levels(yV)[2L])
  fallback<-list(accepted=FALSE,root=NA_integer_,partner=NA_integer_,
       gate="NONE",pair_table=NULL,delta_vs_rlr=NA_real_,
       delta_vs_global=NA_real_,candidate_count=0L,reason="no_valid_gate")
  search<-tryCatch(rg_candidates(D,yD,parent$train,poolD),
                   error=function(e) {message(conditionMessage(e));NULL})
  if(is.null(search))return(fallback)
  best<-NULL
  for(i in seq_len(nrow(search$candidates))) {
    opt<-search$candidates[i,,drop=FALSE]
    xd<-rg_design(D,search$pairs,opt$root,opt$partner,opt$gate)
    xv<-rg_design(V,search$pairs,opt$root,opt$partner,opt$gate)
    pt<-tryCatch(rg_correct(xd,yD,parent$train,xv,parent$test,tree=TRUE),
                 error=function(e)NULL)
    pg<-tryCatch(rg_correct(xd,yD,parent$train,xv,parent$test,tree=FALSE),
                 error=function(e)NULL)
    if(is.null(pt)||is.null(pg))next
    dt<-rg_ll(yV01,parent$test)-rg_ll(yV01,pt)
    dg<-rg_ll(yV01,pg)-rg_ll(yV01,pt)
    ba<-rg_ba(yV01,pt)-rg_ba(yV01,parent$test)
    candidate<-list(accepted=dt>=.005 && dg>=.002 && ba>=-.02,
          root=opt$root,partner=opt$partner,gate=opt$gate,
          pair_table=search$pairs,delta_vs_rlr=dt,
          delta_vs_global=dg,candidate_count=nrow(search$candidates),
          reason="V_GATE",ba_delta=ba)
    if (is.null(best)||dt>best$delta_vs_rlr) best<-candidate
  }
  if(is.null(best))return(fallback)
  best$accepted<-best$delta_vs_rlr>=.005 &&
                 best$delta_vs_global>=.002 &&
                 best$ba_delta>=-.02
  best
}

rg_evaluate <- function(T,yT,E,yE) {
  select<-rg_select(T,yT)
  poolT<-compute_top_pairs_switchbox(T,yT,top_n=1000L,buffer_mult=5L)
  original<-rg_parent(T,yT,E,poolT,K=20L)
  yt<-as.integer(yE==levels(yE)[2L])
  pg<-original$test; pt<-original$test
  if(!is.null(select$pair_table)){
    xd<-rg_design(T,select$pair_table,select$root,select$partner,select$gate)
    xe<-rg_design(E,select$pair_table,select$root,select$partner,select$gate)
    pg<-tryCatch(rg_correct(xd,yT,original$train,xe,original$test,FALSE),
                 error=function(e){message(conditionMessage(e));original$test})
    pt<-tryCatch(rg_correct(xd,yT,original$train,xe,original$test,TRUE),
                 error=function(e){message(conditionMessage(e));original$test})
  }
  gated<-if(isTRUE(select$accepted))pt else original$test
  list(probs=list(ORIGINAL_RLR=original$test,
                  GLOBAL_MAIN_CONTROL=pg,
                  RANK_STUMP_UNGATED=pt,
                  RANK_STUMP_V_GATED=gated),
       selection=select,truth=yt,selected_original=original$selected)
}
