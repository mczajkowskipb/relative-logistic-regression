# RLMT-v6 (development): capture CONTINUOUS within-state rank margin
# information discarded by ternary rank coding of ORIGINAL paper RLR.
#
# No new genes selected beyond ORIGINAL_RLR <=10 signed TSP pairs.
# One original RLR parent; original full training fit unchanged.
# Its exact frozen archive probabilities are gated by the real pilot.
#
# FROZEN modes:
# ORIGINAL_RLR           literal original published model
# GLOBAL_CAL             penalized shared intercept only
# ONEPAIR_MARGIN         single selected among <=10 original pairs,
#                        with original-ternary-state-orthogonalized rank
#                        distance, 1 extra γ and ridge shrinkage
# CONSENSUS_MARGIN       predetermined weighted within-bin margin summary
#                        using fixed ORIGINAL RLR signed coefficients
# ALLPAIR_MARGIN         <=10 within-bin rank margin residuals, ridge only
#
# Important methodological caveat: original <=10 pairs are selected using
# all outer-training labels. OOF scores come from *independently trained*
# parents, but candidate features selected using all OUTER TRAIN labels.
# This is not feature-selection-fully-crossfitted. NO outer TEST leakage.
# The external heldout fold evaluation is necessary to assess optimism.
# v4/v5 19 cohorts already exposed; pilot and any 19×5 = development.

rlmt6_modes <- c("ORIGINAL_RLR","GLOBAL_CAL",
                 "ONEPAIR_MARGIN","CONSENSUS_MARGIN","ALLPAIR_MARGIN")
rlmt6_config <- function()list(
  lambda_beta=.02,lambda_gamma=.15,
  min_rms=.02,clip=3,original_m=10L,
  correction_cap=2.0)

rlmt6_require <- function(x,msg) {
  if(!isTRUE(x))stop("RLMT6_STOP: ",msg,call.=FALSE)
  invisible(TRUE)
}
rlmt6_logit <- function(p) stats::qlogis(
  pmax(1e-8,pmin(1-1e-8,p)))

rlmt6_rank_within_original_state_fit <- function(X,parent,cfg=rlmt6_config()) {
  pairs <- parent$pairs
  rlmt6_require(nrow(pairs)>0L&&nrow(pairs)<=10L &&
     all(c("gene_left","gene_right")%in%names(pairs)),
     "original selected K10 pair list invalid")
  rlmt6_require(cfg$original_m==parent$config$m,
     "original published ternary pair threshold differs")
  R <- encode_pair_matrix(X,pairs,m=cfg$original_m)
  ranks <- rlmt4_rank_matrix(X)
  k <- nrow(pairs)
  p <- ncol(X)
  raw <- matrix(0,nrow(X),k)
  means <- matrix(0,nrow=3,ncol=k,
                  dimnames=list(c("-1","0","1"),NULL))
  scale <- numeric(k)
  for(j in seq_len(k)) {
    a <- as.character(pairs$gene_left[[j]])
    b <- as.character(pairs$gene_right[[j]])
    raw[,j] <- (ranks[,b]-ranks[,a])/max(1,p-1L)
    cat <- R[,j]
    for(t in c(-1L,0L,1L)) {
      present <- which(cat==t)
      means[as.character(t),j] <-
        if(length(present))mean(raw[present,j]) else 0
    }
    residue <- raw[,j]-means[as.character(cat),j]
    scale[j] <- max(cfg$min_rms,sqrt(mean(residue^2)))
  }
  spec <- list(original_pair_list=pairs,
    training_genes=names(X),feature_n=p,
    original_m=cfg$original_m,conditional_means=means,
    scales=scale,clip=cfg$clip,
    interpretation="within_original_ternary_category_rank_distance")
  class(spec) <- "rlmt6_margin_transformer"
  M <- rlmt6_rank_within_original_state_apply(X,spec)
  rlmt6_require(ncol(M)==k&&all(is.finite(M)),
                 "training rank margins are invalid")
  list(transformer=spec,train=M)
}

rlmt6_rank_within_original_state_apply <- function(X,spec) {
  rlmt6_require(inherits(spec,"rlmt6_margin_transformer") &&
      identical(names(X),spec$training_genes),
      "rank distance transformer schema mismatch")
  R <- encode_pair_matrix(X,spec$original_pair_list,m=spec$original_m)
  ranks <- rlmt4_rank_matrix(X)
  p <- spec$feature_n
  k <- nrow(spec$original_pair_list)
  out <- matrix(0,nrow(X),k)
  for(j in seq_len(k)) {
    a <- as.character(spec$original_pair_list$gene_left[[j]])
    b <- as.character(spec$original_pair_list$gene_right[[j]])
    x <- (ranks[,b]-ranks[,a])/max(1,p-1L)
    v <- (x-spec$conditional_means[as.character(R[,j]),j])/
             spec$scales[[j]]
    out[,j] <- pmax(-spec$clip,pmin(spec$clip,v))
  }
  colnames(out) <- as.character(spec$original_pair_list$ordered_pair)
  rownames(out) <- rownames(X)
  rlmt6_require(all(is.finite(out)),"rank distances not finite")
  out
}

rlmt6_fit_offset <- function(y,p_oof,margins=NULL,
                             cfg=rlmt6_config()) {
  yy <- rlmt_y01(y)
  n <- length(y)
  rlmt6_require(length(p_oof)==n&&
    all(is.finite(p_oof))&&all(p_oof>=0&p_oof<=1),
    "crossfit parent probability invalid")
  extra <- if(is.null(margins))matrix(0,n,0) else
    as.matrix(margins)
  rlmt6_require(is.numeric(extra)&&nrow(extra)==n &&
    all(is.finite(extra)),
    "rank-margin columns invalid")
  X <- cbind(1,extra)
  penalties <- c(cfg$lambda_beta,
                 rep(cfg$lambda_gamma,ncol(extra)))
  offset <- rlmt6_logit(p_oof)
  objective <- function(coef) {
    eta <- as.numeric(offset+X%*%coef)
    mean(pmax(eta,0)+log1p(exp(-abs(eta)))-yy*eta)+
      sum(penalties*coef^2)/2
  }
  gradient <- function(coef) {
    eta <- as.numeric(offset+X%*%coef)
    as.numeric(crossprod(X,stats::plogis(eta)-yy)/n+
               penalties*coef)
  }
  opt <- stats::optim(rep(0,ncol(X)),objective,gradient,
     method="BFGS",
     control=list(maxit=600L,reltol=1e-13))
  rlmt6_require(opt$convergence==0 &&
      all(is.finite(opt$par))&&
      max(abs(gradient(opt$par)))<5e-6,
      "penalized fixed feature correction optimizer failed")
  list(beta=opt$par[[1L]],gamma=unname(opt$par[-1L]),
       lambda_beta=cfg$lambda_beta,
       lambda_gamma=cfg$lambda_gamma)
}

rlmt6_consensus_weights <- function(parent) {
  beta <- as.numeric(stats::coef(parent$fit,
                   s=parent$config$lambda_choice))[-1L]
  k <- nrow(parent$pairs)
  rlmt6_require(length(beta)==k&&all(is.finite(beta)),
     "original RLR coefficient order or length differs")
  if(sum(abs(beta))<=1e-10)return(rep(0,k))
  beta/sum(abs(beta))
}

rlmt6_fit <- function(X,y,parent_cfg=rlmt4_config(),
                      cfg=rlmt6_config()) {
  parent <- rlmt4_parent(X,y,parent_cfg)
  transform <- rlmt6_rank_within_original_state_fit(X,parent,cfg)
  oo <- rlmt4_strict_oof(X,y,parent_cfg)
  if(is.null(oo)) {
    return(list(parent=parent,transform=transform$transformer,
      fixed_margin_weights=rlmt6_consensus_weights(parent),
      root_only=TRUE,coefficients=NULL,selected_original_pair=NA_integer_,
      reason="INSUFFICIENT_STRICT_3FOLD_PARENT_OOF",
      model_identity="RLMT6_ORIGINAL_RLR_FALLBACK"))
  }
  M <- transform$train
  yy <- rlmt_y01(y)
  r <- yy-oo$prob
  # Deterministic within-original-pair nomination; only <=10 choices.
  score <- vapply(seq_len(ncol(M)),function(j)
    abs(sum(r*M[,j]))/sqrt(sum(M[,j]^2)+1e-9),numeric(1))
  best <- which.max(score)
  weights <- rlmt6_consensus_weights(parent)
  comp <- as.numeric(M%*%weights)
  fit <- list(
    GLOBAL_CAL=rlmt6_fit_offset(y,oo$prob,cfg=cfg),
    ONEPAIR_MARGIN=rlmt6_fit_offset(y,oo$prob,
                         M[,best,drop=FALSE],cfg),
    CONSENSUS_MARGIN=rlmt6_fit_offset(y,oo$prob,
                         matrix(comp,ncol=1L),cfg),
    ALLPAIR_MARGIN=rlmt6_fit_offset(y,oo$prob,M,cfg))
  list(parent=parent,transform=transform$transformer,
       fixed_margin_weights=weights,
       selected_original_pair=as.integer(best),
       selected_original_gene_a=as.character(
         parent$pairs$gene_left[[best]]),
       selected_original_gene_b=as.character(
         parent$pairs$gene_right[[best]]),
       feature_screen_scores=score,root_only=FALSE,
       coefficients=fit,
       oof_training_parent_probs=oo$prob,
       oof_training_sample_ids=rownames(X),
       oof_partition=oo$partition,
       model_identity="RLMT6_ORIGINAL_PAIR_WITHIN_CATEGORY_MARGIN")
}

rlmt6_predict <- function(model,X) {
  parent <- rlmt4_parent_predict(model$parent,X)
  out <- list(ORIGINAL_RLR=parent$prob)
  if(isTRUE(model$root_only)) {
    for(method in rlmt6_modes[-1L])out[[method]]<-parent$prob
    return(out)
  }
  M <- rlmt6_rank_within_original_state_apply(X,model$transform)
  extra <- list(
    GLOBAL_CAL=matrix(0,nrow(M),0),
    ONEPAIR_MARGIN=M[,model$selected_original_pair,drop=FALSE],
    CONSENSUS_MARGIN=matrix(
      as.numeric(M%*%model$fixed_margin_weights),ncol=1),
    ALLPAIR_MARGIN=M)
  off <- rlmt6_logit(parent$prob)
  for(method in names(extra)) {
    coeff <- model$coefficients[[method]]
    correction <- rep(coeff$beta,nrow(M))
    if(length(coeff$gamma))
      correction <- correction+as.numeric(extra[[method]]%*%coeff$gamma)
    # cap the full model's deviation from parent in log-odds;
    # the same cap is used in every method.
    correction <- pmax(-2,pmin(2,correction))
    out[[method]]<-as.numeric(stats::plogis(off+correction))
  }
  rlmt6_require(identical(names(out),rlmt6_modes) &&
    all(vapply(out,function(v)length(v)==nrow(X) &&
       all(is.finite(v)) && all(v>=0&v<=1),logical(1))),
    "v6 prediction inventory or values invalid")
  out
}
