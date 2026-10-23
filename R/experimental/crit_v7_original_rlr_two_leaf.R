# CRIT-v7 — real, interpretable rank-split MODEL TREE on EXACT ORIGINAL RLR.
# New experimental method; not a proof of superiority on biomedical data.
# Source original parent utilities then rlmt_v4_rank_repair.R before this file.
#
# Each accepted decision stump is genuinely conditional:
#   logit p = logit p_ORIGINAL_RLR + b_false * (1-z) + b_true*z
#             + gamma_false*(1-z)*r_false + gamma_true*z*r_true,
# with z=I(g_a>g_b), and potentially DIFFERENT selected gene pairs
# r_false and r_true in its two small leaf experts. Each leaf has ONE
# relational predictor. This is a 2-leaf model tree, not a probability
# stump and not a globally additive rank-margin correction.
#
# Full supervised route/leaf selection takes place strictly within
# OUTER TRAIN. Separate inner D/V controls split acceptance.
# Parent D-OOF and full-outer-OOF RLR wrappers use independently
# refit original switchBox pair selection, not in-sample RLR labels.
# Exact heldout ORIGINAL_RLR parity is asserted externally.
#
# WARNING: v2-v6 19-cohort exposures are DEVELOPMENT ONLY.
# This 1-split tree can be expressed as logistic interaction features.
# A separate non-tree interaction control MUST be reported.

crit7_modes <- c("ORIGINAL_RLR","GLOBAL_CAL",
                 "GLOBAL_INTERACTION","CRIT_TREE")
crit7_cfg <- function()list(
  seed=42L,
  proposal_n=256L,
  unsupervised_top_genes=12L,
  shortlist=6L,
  min_leaf_n=10L,
  min_leaf_per_class=2L,
  min_leaf_relation_side=3L,
  lambda_intercept=.10,
  lambda_slope=.20,
  max_abs_correction=2.5,
  confirmation_SE_mult=1.0,
  confirmation_cost=.002,
  confirmation_min_gain=.002)

crit7_need <- function(ok,msg) {
  if(!isTRUE(ok))stop("CRIT7_STOP: ",msg,call.=FALSE)
  invisible(TRUE)
}
crit7_logit <- function(p)
  stats::qlogis(pmax(1e-8,pmin(1-1e-8,p)))

# Candidate proposal does NOT require marginal label association:
#  - original RLR pair genes are included, but no restriction to K10;
#  - an unsupervised top-gene rank-variability pool can capture XOR pairs;
#  - deterministic globally sampled pairs supplement it.
# All unsupervised statistics and random choices come from TRAIN D only.
crit7_pool <- function(X,original_pairs,cfg=crit7_cfg()) {
  nm <- names(X)
  p <- length(nm)
  rankmat <- rlmt4_rank_matrix(X)
  variability <- apply(rankmat,2L,stats::var)
  top <- order(-variability,seq_along(variability))[
    seq_len(min(p,cfg$unsupervised_top_genes))]
  combos <- utils::combn(top,2L)
  from_unsup <- data.frame(
    gene_a=nm[combos[1L,]],gene_b=nm[combos[2L,]],
    stringsAsFactors=FALSE)
  from_random <- rlmt_dc_random_pairs(nm,cfg$proposal_n,
                                       seed=cfg$seed+617L)
  from_parent <- data.frame(
    gene_a=as.character(original_pairs$gene_left),
    gene_b=as.character(original_pairs$gene_right),
    stringsAsFactors=FALSE)
  pool <- rbind(from_unsup,from_parent,
                from_random[,c("gene_a","gene_b"),drop=FALSE])
  # One undirected gene pair per input. Signs are still interpreted
  # unambiguously as "gene_a > gene_b" for every tree rule.
  ia <- match(pool$gene_a,nm)
  ib <- match(pool$gene_b,nm)
  left <- pmin(ia,ib)
  right <- pmax(ia,ib)
  idx <- !duplicated(paste(left,right,sep=":")) & left!=right
  pool <- data.frame(gene_a=nm[left[idx]],
                     gene_b=nm[right[idx]],
                     stringsAsFactors=FALSE)
  pool$pool_index <- seq_len(nrow(pool))
  pool
}
crit7_relations <- function(X,pool) {
  g <- vapply(seq_len(nrow(pool)),function(j)
    ifelse(X[[pool$gene_a[[j]]]]>X[[pool$gene_b[[j]]]],1,-1),
    numeric(nrow(X)))
  if(is.null(dim(g)))g <- matrix(g,nrow=nrow(X))
  storage.mode(g) <- "double"
  colnames(g) <- sprintf("pair%03d",seq_len(ncol(g)))
  rownames(g) <- rownames(X)
  g
}
crit7_pair_disjoint <- function(pool,a,b) {
  length(intersect(c(pool$gene_a[[a]],pool$gene_b[[a]]),
                   c(pool$gene_a[[b]],pool$gene_b[[b]])))==0L
}
crit7_leaf_feasible <- function(z,y,cfg) {
  if(sum(z)<cfg$min_leaf_n ||
     sum(!z)<cfg$min_leaf_n)return(FALSE)
  both <- function(s) {
    counts <- table(factor(y[s],levels=levels(y)))
    all(counts>=cfg$min_leaf_per_class)
  }
  both(z) && both(!z)
}
crit7_design <- function(gate,rel_false,rel_true) {
  crit7_need(length(gate)==length(rel_false) &&
               length(gate)==length(rel_true),
               "tree design length mismatch")
  left <- as.numeric(gate<0)
  right <- as.numeric(gate>0)
  cbind(leaf_false=left,leaf_true=right,
        false_rank_pair=left*rel_false,true_rank_pair=right*rel_true)
}
crit7_interaction_design <- function(gate,rel) {
  cbind(intercept=rep(1,length(gate)),gate=gate,
        model_rel=rel,gate_times_model_rel=gate*rel)
}
crit7_fit_offset <- function(y,p_oof,design,cfg=crit7_cfg()) {
  yy <- rlmt_y01(y)
  x <- as.matrix(design)
  crit7_need(nrow(x)==length(yy) && ncol(x)>0L &&
              all(is.finite(x)) &&
              length(p_oof)==length(yy) &&
              all(is.finite(p_oof)) && all(p_oof>=0 & p_oof<=1),
              "invalid RLR offset/fit data")
  penal <- if(ncol(x)==1L)cfg$lambda_intercept else
    c(rep(cfg$lambda_intercept,2L),
      rep(cfg$lambda_slope,ncol(x)-2L))
  off <- crit7_logit(p_oof)
  objective <- function(theta) {
    eta <- as.numeric(off+x%*%theta)
    mean(pmax(eta,0)+log1p(exp(-abs(eta)))-yy*eta)+
      .5*sum(penal*theta^2)
  }
  grad <- function(theta) {
    eta <- as.numeric(off+x%*%theta)
    as.numeric(crossprod(x,stats::plogis(eta)-yy)/length(yy)+
                 penal*theta)
  }
  opt <- stats::optim(rep(0,ncol(x)),objective,grad,
    method="BFGS",
    control=list(maxit=350L,reltol=1e-12))
  crit7_need(opt$convergence==0L &&
                all(is.finite(opt$par)) &&
                max(abs(grad(opt$par)))<1e-5,
                "ridge offset coefficient fit failure")
  as.numeric(opt$par)
}
crit7_correction_prob <- function(prob,design,theta,cfg=crit7_cfg()) {
  d <- as.matrix(design)
  crit7_need(ncol(d)==length(theta)&&
               nrow(d)==length(prob),
               "predictor dimensions invalid")
  delta <- as.numeric(d%*%theta)
  cap <- cfg$max_abs_correction
  as.numeric(stats::plogis(crit7_logit(prob)+
               pmin(cap,pmax(-cap,delta))))
}
crit7_logloss <- function(y,p)
  rlmt_logloss(rlmt_y01(y),p)
crit7_fold_gain <- function(y,baseline,candidate) {
  yy <- rlmt_y01(y)
  b <- pmax(1e-15,pmin(1-1e-15,baseline))
  c <- pmax(1e-15,pmin(1-1e-15,candidate))
  -(yy*log(b)+(1-yy)*log1p(-b))+
     (yy*log(c)+(1-yy)*log1p(-c))
}
crit7_leaf_pair <- function(g,yres,pool,excluded,z,cfg) {
  ix <- which(z)
  if(length(ix)<cfg$min_leaf_n)return(NULL)
  keep <- which(vapply(seq_len(ncol(g)),function(j)
    crit7_pair_disjoint(pool,excluded,j) &&
    min(sum(g[ix,j]>0),sum(g[ix,j]<0))>=
      cfg$min_leaf_relation_side,logical(1)))
  if(!length(keep))return(NULL)
  centered_res <- yres[ix]-mean(yres[ix])
  values <- vapply(keep,function(j) {
    r <- g[ix,j]
    rc <- r-mean(r)
    abs(sum(centered_res*rc))/sqrt(sum(rc^2)+1e-8)
  },numeric(1))
  keep[[which.max(values)]]
}
crit7_nominate <- function(g,y,p_oo,pool,cfg) {
  n <- nrow(g)
  k <- ncol(g)
  crit7_need(length(y)==n&&length(p_oo)==n&&
                k==nrow(pool),"selection inputs invalid")
  residual <- rlmt_y01(y)-p_oo
  # Interaction score, exactly label-dependent ONLY in discovery D.
  cross <- abs(crossprod(g,g*as.numeric(residual)))
  accepted_gate <- vapply(seq_len(k),function(j)
    crit7_leaf_feasible(g[,j]>0,y,cfg),logical(1))
  candidate <- list()
  chosen_gate <- integer()
  # To avoid O(k^2) R loops, take strongest entries first; each
  # gate appears at most once and must have a disjoint relation.
  order_ix <- order(as.vector(cross),decreasing=TRUE)
  for(id in order_ix) {
    i <- ((id-1L)%%k)+1L
    j <- ((id-1L)%/%k)+1L
    if(!accepted_gate[[i]] || !crit7_pair_disjoint(pool,i,j) ||
       i %in% chosen_gate)next
    z <- g[,i]>0
    lf <- crit7_leaf_pair(g,residual,pool,i,!z,cfg)
    lt <- crit7_leaf_pair(g,residual,pool,i,z,cfg)
    if(is.null(lf)||is.null(lt))next
    chosen_gate <- c(chosen_gate,i)
    candidate[[length(candidate)+1L]] <-
      list(gate=i,rel_false=lf,rel_true=lt,
           interaction_rel=j,screen_score=cross[i,j])
    if(length(candidate)>=cfg$shortlist)break
  }
  candidate
}

crit7_select <- function(Xd,yd,pood,Xdval,yval,pval,parentD,
                         cfg=crit7_cfg()) {
  pool <- crit7_pool(Xd,parentD$pairs,cfg)
  gD <- crit7_relations(Xd,pool)
  gV <- crit7_relations(Xdval,pool)
  candidates <- crit7_nominate(gD,yd,pood,pool,cfg)
  best <- NULL
  best_lcb <- -Inf
  all <- list()
  for(candidate in candidates) {
    i <- candidate$gate
    fd <- crit7_design(gD[,i],gD[,candidate$rel_false],
                          gD[,candidate$rel_true])
    fv <- crit7_design(gV[,i],gV[,candidate$rel_false],
                          gV[,candidate$rel_true])
    theta <- crit7_fit_offset(yd,pood,fd,cfg)
    pp <- crit7_correction_prob(pval,fv,theta,cfg)
    gains <- crit7_fold_gain(yval,pval,pp)
    mg <- mean(gains)
    se <- if(length(gains)>1L)stats::sd(gains)/
      sqrt(length(gains)) else Inf
    conservative <- mg-cfg$confirmation_SE_mult*se-
      cfg$confirmation_cost
    all[[length(all)+1L]]<-data.frame(
      gate_index=i,rel_false_index=candidate$rel_false,
      rel_true_index=candidate$rel_true,
      non_tree_interaction_index=candidate$interaction_rel,
      discovery_score=candidate$screen_score,
      confirmation_mean_gain=mg,confirmation_se=se,
      confirmation_lower=conservative,
      stringsAsFactors=FALSE)
    if(is.finite(conservative)&&conservative>best_lcb) {
      best_lcb<-conservative
      best<-candidate
    }
  }
  accepted<-!is.null(best) &&
    best_lcb>cfg$confirmation_min_gain
  list(accepted=accepted,selected=best,
       candidate_pool=pool,
       confirmation=if(length(all))do.call(rbind,all) else
         data.frame(),best_confirmed_lower=best_lcb,
       evaluated_candidate_count=length(candidates),
       inner_holdout_n=length(yval))
}

crit7_fit <- function(X,y,cfg=crit7_cfg(),
                      parent_cfg=rlmt4_config()) {
  original <- rlmt4_parent(X,y,parent_cfg)
  result <- list(original=original,cfg=cfg,
      original_schema=names(X),selected=NULL,
      gate=NULL,pair_false=NULL,pair_true=NULL,
      confirm=NULL,tree_coef=NULL,interaction_coef=NULL,
      interaction_pair=NULL,cal_coef=NULL,
      model_identity="CRIT7_ORIGINAL_RLR_CONDITIONAL_RANK_TREE")
  # Full train OOF is the only basis of final correction coefficients.
  full_oof <- rlmt4_strict_oof(X,y,parent_cfg)
  if(!is.null(full_oof))
    result$cal_coef<-crit7_fit_offset(y,full_oof$prob,
                            matrix(1,nrow(X),1),cfg)
  # For honest D/V we insist on enough class support to refit
  # ORIGINAL RLR and its own strict 3fold OOF within D.
  parts <- make_repeated_stratified_folds(
        y,n_folds=3L,n_repeats=1L,seed=cfg$seed+731L)
  d<-parts[[1L]]$train_idx
  v<-parts[[1L]]$test_idx
  result$partition <- list(train=rownames(X)[d],
                           confirmation=rownames(X)[v])
  if(min(table(y[d]))<12L||min(table(y[v]))<2L ||
     is.null(full_oof)){
    result$reason<-"STRICT_SUPPORT_OR_OOF_INSUFFICIENT"
    return(result)
  }
  parent_d<-rlmt4_parent(X[d,,drop=FALSE],y[d],parent_cfg)
  dOOF<-rlmt4_strict_oof(X[d,,drop=FALSE],y[d],parent_cfg)
  if(is.null(dOOF)) {
    result$reason<-"DISCOVERY_RLR_OOF_UNAVAILABLE"
    return(result)
  }
  pval<-rlmt4_parent_predict(parent_d,X[v,,drop=FALSE])$prob
  selection<-crit7_select(X[d,,drop=FALSE],y[d],dOOF$prob,
        X[v,,drop=FALSE],y[v],pval,parent_d,cfg)
  result$confirm<-selection
  if(is.null(selection$selected)){
    result$reason<-"NO_FEASIBLE_RANK_INTERACTION_GATE"
    return(result)
  }
  pair<-selection$candidate_pool
  choose<-selection$selected
  result$gate<-pair[choose$gate,c("gene_a","gene_b"),drop=FALSE]
  result$pair_false<-pair[choose$rel_false,c("gene_a","gene_b"),drop=FALSE]
  result$pair_true<-pair[choose$rel_true,c("gene_a","gene_b"),drop=FALSE]
  result$interaction_pair<-pair[choose$interaction_rel,
                              c("gene_a","gene_b"),drop=FALSE]
  # The GLOBAL_INTERACTION control is fit regardless of whether
  # the actual tree's D/V gain satisfies split acceptance.
  gate_full<-crit7_relations(X,result$gate)[,1L]
  rint<-crit7_relations(X,result$interaction_pair)[,1L]
  global_features<-crit7_interaction_design(gate_full,rint)
  result$interaction_coef<-crit7_fit_offset(
              y,full_oof$prob,global_features,cfg)
  if(!selection$accepted) {
    result$reason<-"VALIDATION_DID_NOT_ACCEPT_TREE"
    return(result)
  }
  rf<-crit7_relations(X,result$pair_false)[,1L]
  rt<-crit7_relations(X,result$pair_true)[,1L]
  tree_features<-crit7_design(gate_full,rf,rt)
  result$tree_coef<-crit7_fit_offset(y,full_oof$prob,
                                   tree_features,cfg)
  result$selected<-TRUE
  result$reason<-"CONFIRMED_GENUINE_TWO_LEAF_TREE"
  result
}
crit7_predict <- function(model,X) {
  root <- rlmt4_parent_predict(model$original,X)$prob
  ans <- list(ORIGINAL_RLR=root)
  ans$GLOBAL_CAL <- if(is.null(model$cal_coef))root else
    crit7_correction_prob(root,matrix(1,nrow(X),1),
                          model$cal_coef,model$cfg)
  if(is.null(model$interaction_coef)) {
    ans$GLOBAL_INTERACTION<-root
  } else {
    gate<-crit7_relations(X,model$gate)[,1L]
    rr<-crit7_relations(X,model$interaction_pair)[,1L]
    ans$GLOBAL_INTERACTION <- crit7_correction_prob(root,
       crit7_interaction_design(gate,rr),
       model$interaction_coef,model$cfg)
  }
  actual <- !is.null(model$tree_coef)
  if(!actual) {
    ans$CRIT_TREE<-root
    leaf_id<-rep(0L,nrow(X))
  } else {
    gate<-crit7_relations(X,model$gate)[,1L]
    rf<-crit7_relations(X,model$pair_false)[,1L]
    rt<-crit7_relations(X,model$pair_true)[,1L]
    ans$CRIT_TREE<-crit7_correction_prob(root,
                 crit7_design(gate,rf,rt),
                 model$tree_coef,model$cfg)
    leaf_id<-ifelse(gate>0,2L,1L)
  }
  crit7_need(identical(names(ans),crit7_modes)&&
       all(vapply(ans,function(p)
        length(p)==nrow(X) && all(is.finite(p)) &&
        all(p>=0&p<=1),logical(1))),
        "CRIT7 output probabilities invalid")
  list(probs=ans,leaf_id=leaf_id,split_accepted=actual)
}
