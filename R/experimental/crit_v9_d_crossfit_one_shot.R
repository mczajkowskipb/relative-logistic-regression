# CRIT-v9: 6 D-nominated candidates -> D-only coefficient-crossfit
# ranking -> ONE V-confirmed candidate, 5 method-specific V gates.
# A scientific development candidate, not publication-ready inference.
# V is never used to choose a candidate, although V is used for each
# of the five predeclared arm acceptance decisions.
# D score is not an unbiased CV estimate because feature nomination
# originally used all D labels. Its only function is D-only ranking.
# CRIT-v8 experimental: ONE D-only nomination, separate ONE-SHOT V gates.
# EXPOSED_RETROSPECTIVE_DEV_ONLY, no inference or external confirmation.
# Source original utilities, rlmt_v4_rank_repair.R and CRIT-v7 first.
crit9_modes <- c("ORIGINAL_RLR","GLOBAL_CAL","GATE_ONLY",
                 "ADDITIVE_PAIRS","GLOBAL_INTERACTION","CONDITIONAL_TREE")
crit9_cfg <- function() {
  z <- crit7_cfg()
  z$shortlist <- 6L
  z
}
crit9_need <- function(ok,msg) {
  if(!isTRUE(ok))stop("CRIT9_STOP: ",msg,call.=FALSE)
  invisible(TRUE)
}
crit9_designs <- function(X,model) {
  n <- nrow(X)
  ans <- list(GLOBAL_CAL=matrix(1,nrow=n,ncol=1L,
                                 dimnames=list(NULL,"intercept")))
  if(is.null(model$gate))return(ans)
  g <- crit7_relations(X,model$gate)[,1L]
  f <- crit7_relations(X,model$pair_false)[,1L]
  t <- crit7_relations(X,model$pair_true)[,1L]
  h <- crit7_relations(X,model$interaction_pair)[,1L]
  ds <- crit7_design(g,f,t)
  ans$GATE_ONLY <- ds[,1:2,drop=FALSE]
  ans$ADDITIVE_PAIRS <- cbind(
       ds[,1:2,drop=FALSE],pair_false=f,pair_true=t)
  ans$GLOBAL_INTERACTION <- crit7_interaction_design(g,h)
  ans$CONDITIONAL_TREE <- ds
  ans
}

crit9_rank_D_candidates <- function(D,y,pood,parentD,cfg) {
  pool <- crit7_pool(D,parentD$pairs,cfg)
  g <- crit7_relations(D,pool)
  candidate <- crit7_nominate(g,y,pood,pool,cfg)
  if(!length(candidate))
    return(list(pool=pool,candidates=candidate,
                ranking=data.frame(),selected_rank=NULL))
  folds <- make_repeated_stratified_folds(
       y,n_folds=3L,n_repeats=1L,seed=cfg$seed+1907L)
  need <- function(ok,msg)crit9_need(ok,msg)
  all <- list()
  for(j in seq_along(candidate)) {
    x <- candidate[[j]]
    feat <- crit7_design(g[,x$gate],g[,x$rel_false],
                          g[,x$rel_true])
    gains <- rep(NA_real_,length(y))
    for(f in folds) {
      tr <- f$train_idx
      te <- f$test_idx
      need(!length(intersect(tr,te)),
           "within-D ranking fold patient overlap")
      theta <- crit7_fit_offset(y[tr],pood[tr],
                               feat[tr,,drop=FALSE],cfg)
      out <- crit7_correction_prob(
                 pood[te],feat[te,,drop=FALSE],theta,cfg)
      gains[te] <- crit7_fold_gain(y[te],pood[te],out)
    }
    need(all(is.finite(gains)),
         "D ranking failed full patient coverage")
    avg <- mean(gains)
    se <- stats::sd(gains)/sqrt(length(gains))
    all[[j]]<-data.frame(D_rank=j,
       discovery_screen_score=as.numeric(x$screen_score),
       D_cv_mean_gain=avg,D_cv_SE=se,
       D_cv_lower=avg-se-cfg$confirmation_cost,
       stringsAsFactors=FALSE)
  }
  ranking<-do.call(rbind,all)
  # Tie break by fixed discovery rank; all inputs come from D only.
  selected <- which.max(ranking$D_cv_lower)
  need(length(selected)==1L && is.finite(ranking$D_cv_lower[[selected]]),
       "no finite D-ranked candidate")
  list(pool=pool,candidates=candidate,
       ranking=ranking,selected_rank=as.integer(selected))
}

crit9_fit <- function(X,y,cfg=crit9_cfg(),parent_cfg=rlmt4_config()) {
  crit9_need(is.data.frame(X)&&is.factor(y)&&
         nrow(X)==length(y)&&length(levels(y))==2L,
         "input schema")
  crit9_need(cfg$shortlist==6L,
         "v9 discovery shortlist must be exactly six; only selected one enters V")
  root <- rlmt4_parent(X,y,parent_cfg)
  model <- list(original=root,cfg=cfg,gate=NULL,
        pair_false=NULL,pair_true=NULL,interaction_pair=NULL,
        accepted=setNames(rep(FALSE,length(crit9_modes)-1L),crit9_modes[-1L]),
        coefficients=setNames(vector("list",length(crit9_modes)-1L),crit9_modes[-1L]),
        validation=data.frame(),
        nomination_n=0L,selected_D_rank=NULL,ranking=data.frame(),partition=NULL,
        reason="START",schema=names(X),
        model_identity="CRIT9_D_CROSSFIT_SINGLE_V_FIVE_ARM_CONTROL")
  full_oof <- rlmt4_strict_oof(X,y,parent_cfg)
  if(is.null(full_oof)) {
    model$reason <- "STRICT_FULL_TRAIN_OOF_UNAVAILABLE"
    return(model)
  }
  parts <- make_repeated_stratified_folds(
            y,n_folds=3L,n_repeats=1L,seed=cfg$seed+731L)[[1L]]
  d <- parts$train_idx
  v <- parts$test_idx
  crit9_need(!length(intersect(d,v)) &&
      identical(sort(c(d,v)),seq_len(nrow(X))),
      "D/V overlap or missing training patients")
  model$partition <- list(discovery=rownames(X)[d],
        confirmation=rownames(X)[v])
  if(min(table(y[d]))<12L || min(table(y[v]))<2L) {
    model$reason <- "INSUFFICIENT_DISCOVERY_OR_CONFIRMATION_SUPPORT"
    return(model)
  }
  D <- X[d,,drop=FALSE]
  V <- X[v,,drop=FALSE]
  pd <- rlmt4_parent(D,y[d],parent_cfg)
  dood <- rlmt4_strict_oof(D,y[d],parent_cfg)
  if(is.null(dood)) {
    model$reason <- "STRICT_DISCOVERY_RLR_OOF_UNAVAILABLE"
    return(model)
  }
  pv <- rlmt4_parent_predict(pd,V)$prob
  # DISCOVERY D ONLY: screen up to six, choose one by within-D
  # cross-fitted ridge gain. Nothing from V is used to choose features.
  selection <- crit9_rank_D_candidates(D,y[d],dood$prob,pd,cfg)
  pool <- selection$pool
  nominal <- selection$candidates
  model$ranking <- selection$ranking
  model$nomination_n <- length(nominal)
  model$selected_D_rank <- selection$selected_rank
  if(!is.null(selection$selected_rank)) {
    c <- nominal[[selection$selected_rank]]
    model$gate <- pool[c$gate,c("gene_a","gene_b"),drop=FALSE]
    model$pair_false <- pool[c$rel_false,c("gene_a","gene_b"),drop=FALSE]
    model$pair_true <- pool[c$rel_true,c("gene_a","gene_b"),drop=FALSE]
    model$interaction_pair <- pool[c$interaction_rel,c("gene_a","gene_b"),drop=FALSE]
  }
  designsD <- crit9_designs(D,model)
  designsV <- crit9_designs(V,model)
  designsT <- crit9_designs(X,model)
  model$validation <- data.frame(
    method=crit9_modes[-1L],available=FALSE,accepted=FALSE,
    validation_mean=NA_real_,validation_se=NA_real_,
    validation_lower=NA_real_,stringsAsFactors=FALSE)
  for(i in seq_along(model$validation$method)) {
    mode <- model$validation$method[[i]]
    if(is.null(designsD[[mode]]))next
    fitD <- crit7_fit_offset(y[d],dood$prob,designsD[[mode]],cfg)
    val <- crit7_correction_prob(pv,designsV[[mode]],fitD,cfg)
    gains <- crit7_fold_gain(y[v],pv,val)
    mg <- mean(gains)
    se <- if(length(gains)>1L)stats::sd(gains)/sqrt(length(gains)) else Inf
    lower <- mg - cfg$confirmation_SE_mult*se -
                       cfg$confirmation_cost
    # Each control separately passes/fails exactly one predeclared
    # test on V. There is no winner picked across the five arms on V.
    yes <- is.finite(lower) && lower>cfg$confirmation_min_gain
    model$validation$available[[i]] <- TRUE
    model$validation$accepted[[i]] <- yes
    model$validation$validation_mean[[i]] <- mg
    model$validation$validation_se[[i]] <- se
    model$validation$validation_lower[[i]] <- lower
    if(yes) {
      model$coefficients[[mode]] <-
          crit7_fit_offset(y,full_oof$prob,designsT[[mode]],cfg)
      model$accepted[[mode]] <- TRUE
    }
  }
  model$reason <- if(is.null(model$gate)) "NO_SINGLE_D_CANDIDATE" else
                    "D_CROSSFIT_ONE_NOMINATION_ALL_V_DECISIONS_RECORDED"
  model
}
crit9_predict <- function(model,X) {
  crit9_need(identical(names(X),model$schema),
            "predictor schema mismatch")
  root <- rlmt4_parent_predict(model$original,X)$prob
  ds <- crit9_designs(X,model)
  probs <- list(ORIGINAL_RLR=root)
  for(mode in crit9_modes[-1L]) {
    beta <- model$coefficients[[mode]]
    if(is.null(beta)) {
      probs[[mode]] <- root
    } else {
      crit9_need(isTRUE(model$accepted[[mode]]) &&
                 !is.null(ds[[mode]]),
                 "unaccepted fitted model or missing features")
      probs[[mode]] <- crit7_correction_prob(
             root,ds[[mode]],beta,model$cfg)
    }
  }
  crit9_need(identical(names(probs),crit9_modes) &&
    all(vapply(probs,function(p)
      length(p)==nrow(X) && all(is.finite(p)) &&
      all(p>=0&p<=1),logical(1))),"invalid prediction vectors")
  leaf_id <- if(is.null(model$gate))rep(0L,nrow(X)) else
    ifelse(crit7_relations(X,model$gate)[,1L]>0,2L,1L)
  list(probs=probs,leaf_id=leaf_id,
       accepted=model$accepted,validation=model$validation)
}
