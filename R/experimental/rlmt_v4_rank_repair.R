# RLMT-v4 — interpretable rank-repair models on EXACT ORIGINAL RLR.
# This is a NEW method family, NOT proof of out-of-sample superiority.
# Source after 07_world_protocol_utils.R and 28_rlmt_utils.R,
# 40_rlmt_dc_utils.R (for deterministic random candidate proposals),
# 42_rlmt_dc_run.R (for utilities).
#
# Three distinct RELATIONAL residual features plus exact nested no-rule
# calibrator. All supervised steps use OUTER TRAIN ONLY. No test labels.
#
# CAL:       logit(p_original) + beta
# CONF:      CAL + gamma * 4*p*(1-p)*(I(a>b)-training_frequency)
# ACTIVITY:  CAL + gamma * (1-mean active original RLR pairs) *
#                                  (I(a>b)-training_frequency)
# SOFT_GAP:  CAL + gamma * (tanh((rank(a)-rank(b))/10)-train_mean)
#
# p_original is the literal original paper's RLR fitting procedure,
# not the earlier stratifed-fold RLR_ALIGNED_V2 baseline.
# One parent and one strict 3-fold out-of-fold parent fitted per fold.
# Pair selection uses OOF (y-p) correlation on the training set,
# predeclared 4096 random pairs plus original 10 selected pairs.
# Beta and gamma penalized logistic MLE with fixed λ(β)=0.02,
# λ(γ)=0.10, gamma shrunk to zero; no small-sample leaf cv.glmnet.

rlmt4_modes <- c("ORIGINAL_RLR","GLOBAL_CAL","CONF_ROUTE",
                 "ACTIVITY_ROUTE","SOFT_GAP_ROUTE")
rlmt4_config <- function()list(seed=42L,
  pool=1000L,K=10L,m=10L,d=2L,alpha=0.5,
  lambda_choice="lambda.min",proposal_n=4096L,
  min_side=8L,rank_scale=10,lambda_beta=.02,lambda_gamma=.10)
rlmt4_require <- function(x,msg) {
  if(!isTRUE(x))stop("RLMT4_STOP: ",msg,call.=FALSE)
  invisible(TRUE)
}
rlmt4_parent <- function(X,y,cfg=rlmt4_config()) {
  rlmt4_require(is.data.frame(X) && is.factor(y) &&
    length(levels(y))==2L && length(y)==nrow(X) &&
    all(vapply(X,function(a)is.numeric(a)&&all(is.finite(a)),logical(1))),
    "parent X/y schema")
  pool <- compute_top_pairs_switchbox(
    X_train=X,y_train=y,top_n=cfg$pool,buffer_mult=5L,
    out_dir=NULL,save_all=FALSE)
  pairs <- apply_disjointness(pool,K=cfg$K,d=cfg$d)
  rlmt4_require(nrow(pairs)>0L,"no original RLR pairs")
  R <- encode_pair_matrix(X,pairs,m=cfg$m)
  set.seed(cfg$seed)
  cvfit <- glmnet::cv.glmnet(
    x=R,y=rlmt_y01(y),family="binomial",
    alpha=cfg$alpha,nfolds=safe_nfolds(y,5),
    type.measure="deviance")
  structure(list(fit=cvfit,pairs=pairs,
    gene_names=names(X),seed=cfg$seed,
    model_identity="FROZEN_ORIGINAL_RLR_K10_M10_D2",
    config=cfg),class="rlmt4_original_parent")
}
rlmt4_parent_predict <- function(parent,X) {
  rlmt4_require(identical(names(X),parent$gene_names),
                "parent predictor gene schema mismatch")
  design <- encode_pair_matrix(X,parent$pairs,
                               m=parent$config$m)
  prob <- as.numeric(stats::predict(parent$fit,newx=design,
                 s=parent$config$lambda_choice,type="response"))
  rlmt4_require(length(prob)==nrow(X) && all(is.finite(prob)) &&
                all(prob>=0 & prob<=1),
                "original RLR prediction invalid")
  activity <- 1-rowMeans(design!=0)
  list(prob=prob,inactivity=activity)
}
rlmt4_strict_oof <- function(X,y,cfg=rlmt4_config()) {
  folds <- make_repeated_stratified_folds(
    y,n_folds=3L,n_repeats=1L,seed=cfg$seed+101L)
  need <- min(vapply(folds,function(f)
    min(table(y[f$train_idx])),numeric(1)))
  if(need<5L)return(NULL)
  pred <- rep(NA_real_,nrow(X))
  inactivity <- pred
  count <- integer(nrow(X))
  for(f in folds) {
    train <- f$train_idx
    val <- f$test_idx
    rlmt4_require(length(intersect(train,val))==0L,
                  "train/valid OOF overlap")
    parent <- rlmt4_parent(X[train,,drop=FALSE],y[train],cfg)
    p <- rlmt4_parent_predict(parent,X[val,,drop=FALSE])
    pred[val] <- p$prob
    inactivity[val] <- p$inactivity
    count[val] <- count[val]+1L
  }
  rlmt4_require(all(count==1L) &&
                 all(is.finite(pred)) &&
                 all(is.finite(inactivity)),
                 "honest OOF coverage incomplete")
  list(prob=pred,inactivity=inactivity,
       partition=lapply(folds,function(f)
         list(train=rownames(X)[f$train_idx],
              heldout=rownames(X)[f$test_idx])))
}
rlmt4_logit <- function(p) stats::qlogis(pmax(1e-8,pmin(1-1e-8,p)))
rlmt4_fit_correction <- function(y,p_oof,feature=NULL,cfg=rlmt4_config()) {
  yy <- rlmt_y01(y)
  n <- length(yy)
  rlmt4_require(length(p_oof)==n && all(is.finite(p_oof)) &&
                all(p_oof>=0&p_oof<=1),"bad OOF correction input")
  use_route <- !is.null(feature)
  if(use_route)
    rlmt4_require(length(feature)==n &&
                  all(is.finite(feature)),
                  "bad residual-feature length/values")
  x <- if(use_route) cbind(1,as.numeric(feature)) else
                       matrix(1,nrow=n,ncol=1L)
  penalty <- if(use_route)
               c(cfg$lambda_beta,cfg$lambda_gamma) else
               cfg$lambda_beta
  off <- rlmt4_logit(p_oof)
  objective <- function(theta) {
    eta <- as.numeric(off+x%*%theta)
    # stable logistic negative log likelihood: softplus(eta)-y*eta
    mean(pmax(eta,0)+log1p(exp(-abs(eta)))-yy*eta)+
        sum(penalty*theta^2)/2
  }
  gradient <- function(theta) {
    eta <- as.numeric(off+x%*%theta)
    as.numeric(crossprod(x,stats::plogis(eta)-yy)/n+
               penalty*theta)
  }
  opt <- stats::optim(rep(0,ncol(x)),objective,gradient,
                method="BFGS",control=list(maxit=300L,reltol=1e-12))
  rlmt4_require(opt$convergence==0L &&
    all(is.finite(opt$par)) &&
    max(abs(gradient(opt$par)))<1e-6,
    "penalized correction fit did not converge")
  list(beta=unname(opt$par[[1L]]),
       gamma=if(use_route)unname(opt$par[[2L]]) else 0,
       penalty_beta=cfg$lambda_beta,
       penalty_gamma=cfg$lambda_gamma)
}
rlmt4_rank_matrix <- function(X) {
  A <- t(apply(X,1L,rank,ties.method="average"))
  colnames(A) <- names(X)
  rownames(A) <- rownames(X)
  A
}
rlmt4_route_feature <- function(X,p,inactivity,
                                 gene_a,gene_b,mode,
                                 frequency=NA_real_,gap_mean=NA_real_,
                                 rankmat=NULL,cfg=rlmt4_config()) {
  rlmt4_require(gene_a %in% names(X) && gene_b %in% names(X) &&
                gene_a!=gene_b,"route pair genes not present")
  if(mode=="SOFT_GAP_ROUTE") {
    if(is.null(rankmat))rankmat <- rlmt4_rank_matrix(X)
    feature <- tanh((rankmat[,gene_a]-rankmat[,gene_b])/
                     cfg$rank_scale)
    return(if(is.finite(gap_mean))feature-gap_mean else feature)
  }
  route <- as.numeric(X[[gene_a]]>X[[gene_b]])
  if(!is.finite(frequency))return(route)
  gate <- switch(mode,
      CONF_ROUTE=4*p*(1-p),
      ACTIVITY_ROUTE=inactivity,
      stop("RLMT4_STOP: unsupported route feature"))
  rlmt4_require(all(is.finite(gate)) &&
                all(gate>=0 & gate<=1),
                "route reliability gate invalid")
  as.numeric(gate*(route-frequency))
}
rlmt4_candidates <- function(X,parent,cfg=rlmt4_config()) {
  random <- rlmt_dc_random_pairs(names(X),
                  cfg$proposal_n,seed=cfg$seed+201L)
  extra <- data.frame(
    pool_index=seq_len(nrow(parent$pairs))+nrow(random),
    gene_a=as.character(parent$pairs$gene_left),
    gene_b=as.character(parent$pairs$gene_right))
  all <- rbind(random,extra)
  keys <- paste(all$gene_a,all$gene_b,sep="\t")
  all[!duplicated(keys),,drop=FALSE]
}
rlmt4_nominate <- function(X,y,oof,parent,mode,
                           cfg=rlmt4_config()) {
  candidates <- rlmt4_candidates(X,parent,cfg)
  yy <- rlmt_y01(y)
  residual <- yy-oof$prob
  ranks <- if(mode=="SOFT_GAP_ROUTE")
             rlmt4_rank_matrix(X) else NULL
  best <- NULL
  highest <- -Inf
  for(k in seq_len(nrow(candidates))) {
    a <- candidates$gene_a[[k]]
    b <- candidates$gene_b[[k]]
    z <- as.numeric(X[[a]]>X[[b]])
    n1 <- sum(z)
    n0 <- length(z)-n1
    if(n1<cfg$min_side || n0<cfg$min_side)next
    pi <- mean(z)
    gap <- if(mode=="SOFT_GAP_ROUTE")
      mean(tanh((ranks[,a]-ranks[,b])/cfg$rank_scale)) else NA_real_
    feature <- rlmt4_route_feature(X,oof$prob,oof$inactivity,
               a,b,mode,frequency=pi,gap_mean=gap,
               rankmat=ranks,cfg=cfg)
    v <- sum(feature^2)
    if(v<1e-9)next
    score <- abs(sum(residual*feature))/sqrt(v)
    if(score>highest+1e-12) {
      highest <- score
      best <- list(gene_a=a,gene_b=b,
        pi=pi,gap_mean=gap,train_true=n1,train_false=n0,
        residual_score=score,pool_index=candidates$pool_index[[k]])
    }
  }
  list(pair=best,candidate_n=nrow(candidates))
}
rlmt4_fit <- function(X,y,cfg=rlmt4_config()) {
  base <- rlmt4_parent(X,y,cfg)
  oof <- rlmt4_strict_oof(X,y,cfg)
  if(is.null(oof))return(list(
    parent=base,cfg=cfg,fit=NULL,
    routes=list(),model_identity="ORIGINAL_RLR_ONLY_INSUFFICIENT_OOF",
    reason="strict_original_RLR_OOF_requires_class_5_per_inner_train"))
  cal <- rlmt4_fit_correction(y,oof$prob,cfg=cfg)
  routes <- list()
  for(method in c("CONF_ROUTE","ACTIVITY_ROUTE","SOFT_GAP_ROUTE")) {
    nominated <- rlmt4_nominate(X,y,oof,base,method,cfg)
    if(is.null(nominated$pair)) {
      routes[[method]] <- list(pair=NULL,coefficients=cal,
                         candidate_n=nominated$candidate_n)
      next
    }
    pair <- nominated$pair
    feature <- rlmt4_route_feature(X,oof$prob,oof$inactivity,
       pair$gene_a,pair$gene_b,method,
       frequency=pair$pi,gap_mean=pair$gap_mean,cfg=cfg)
    fit <- rlmt4_fit_correction(y,oof$prob,feature,cfg)
    routes[[method]] <- list(pair=pair,coefficients=fit,
                            candidate_n=nominated$candidate_n)
  }
  list(parent=base,cfg=cfg,fit=cal,routes=routes,
       model_identity="ORIGINAL_RLR_CONFIDENCE_GATED_RANK_REPAIR",
       oof_training_sample_ids=rownames(X),
       oof_prob_training=oof$prob,
       oof_partition=oof$partition)
}
rlmt4_predict <- function(model,X) {
  p <- rlmt4_parent_predict(model$parent,X)
  probs <- list(ORIGINAL_RLR=p$prob)
  if(is.null(model$fit)) {
    for(k in rlmt4_modes[-1L])probs[[k]] <- p$prob
  } else {
    cal <- model$fit
    probs$GLOBAL_CAL <- stats::plogis(
      rlmt4_logit(p$prob)+cal$beta)
    for(mode in c("CONF_ROUTE","ACTIVITY_ROUTE","SOFT_GAP_ROUTE")) {
      route <- model$routes[[mode]]
      if(is.null(route$pair)) {
        probs[[mode]] <- probs$GLOBAL_CAL
        next
      }
      pair <- route$pair
      feature <- rlmt4_route_feature(X,p$prob,p$inactivity,
        pair$gene_a,pair$gene_b,mode,
        frequency=pair$pi,gap_mean=pair$gap_mean,cfg=model$cfg)
      probs[[mode]] <- stats::plogis(
        rlmt4_logit(p$prob)+route$coefficients$beta+
                    route$coefficients$gamma*feature)
    }
  }
  rlmt4_require(identical(names(probs),rlmt4_modes) &&
                 all(vapply(probs,function(x)
                   length(x)==nrow(X) && all(is.finite(x)) &&
                   all(x>=0 & x<=1),logical(1))),
                 "invalid or incomplete v4 predictions")
  probs
}
