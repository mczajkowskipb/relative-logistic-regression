# CRIT-v8 experimental: ONE D-only nomination, separate ONE-SHOT V gates.
# EXPOSED_RETROSPECTIVE_DEV_ONLY, no inference or external confirmation.
# Source original utilities, rlmt_v4_rank_repair.R and CRIT-v7 first.
crit8_modes <- c("ORIGINAL_RLR","GLOBAL_CAL","GATE_ONLY",
                 "ADDITIVE_PAIRS","GLOBAL_INTERACTION","CONDITIONAL_TREE")
crit8_cfg <- function() {
  z <- crit7_cfg()
  z$shortlist <- 1L
  z
}
crit8_need <- function(ok,msg) {
  if(!isTRUE(ok))stop("CRIT8_STOP: ",msg,call.=FALSE)
  invisible(TRUE)
}
crit8_designs <- function(X,model) {
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
crit8_fit <- function(X,y,cfg=crit8_cfg(),parent_cfg=rlmt4_config()) {
  crit8_need(is.data.frame(X)&&is.factor(y)&&
         nrow(X)==length(y)&&length(levels(y))==2L,
         "input schema")
  crit8_need(cfg$shortlist==1L,
         "protocol prohibits multiple V-evaluated nominations")
  root <- rlmt4_parent(X,y,parent_cfg)
  model <- list(original=root,cfg=cfg,gate=NULL,
        pair_false=NULL,pair_true=NULL,interaction_pair=NULL,
        accepted=setNames(rep(FALSE,length(crit8_modes)-1L),crit8_modes[-1L]),
        coefficients=setNames(vector("list",length(crit8_modes)-1L),crit8_modes[-1L]),
        validation=data.frame(),
        nomination_n=0L,partition=NULL,
        reason="START",schema=names(X),
        model_identity="CRIT8_SINGLE_NOMINATION_SIX_ARM_CONTROL")
  full_oof <- rlmt4_strict_oof(X,y,parent_cfg)
  if(is.null(full_oof)) {
    model$reason <- "STRICT_FULL_TRAIN_OOF_UNAVAILABLE"
    return(model)
  }
  parts <- make_repeated_stratified_folds(
            y,n_folds=3L,n_repeats=1L,seed=cfg$seed+731L)[[1L]]
  d <- parts$train_idx
  v <- parts$test_idx
  crit8_need(!length(intersect(d,v)) &&
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
  # One candidate completely determined by discovery D. No V-driven
  # shortlist ranking, no V-driven choice of predictor pair or structure.
  pool <- crit7_pool(D,pd$pairs,cfg)
  gD <- crit7_relations(D,pool)
  nominal <- crit7_nominate(gD,y[d],dood$prob,pool,cfg)
  crit8_need(length(nominal)<=1L,
       "multiple validation candidates are forbidden")
  model$nomination_n <- length(nominal)
  if(length(nominal)==1L) {
    c <- nominal[[1L]]
    model$gate <- pool[c$gate,c("gene_a","gene_b"),drop=FALSE]
    model$pair_false <- pool[c$rel_false,c("gene_a","gene_b"),drop=FALSE]
    model$pair_true <- pool[c$rel_true,c("gene_a","gene_b"),drop=FALSE]
    model$interaction_pair <- pool[c$interaction_rel,c("gene_a","gene_b"),drop=FALSE]
  }
  designsD <- crit8_designs(D,model)
  designsV <- crit8_designs(V,model)
  designsT <- crit8_designs(X,model)
  model$validation <- data.frame(
    method=crit8_modes[-1L],available=FALSE,accepted=FALSE,
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
                    "SINGLE_D_NOMINATION_ALL_V_DECISIONS_RECORDED"
  model
}
crit8_predict <- function(model,X) {
  crit8_need(identical(names(X),model$schema),
            "predictor schema mismatch")
  root <- rlmt4_parent_predict(model$original,X)$prob
  ds <- crit8_designs(X,model)
  probs <- list(ORIGINAL_RLR=root)
  for(mode in crit8_modes[-1L]) {
    beta <- model$coefficients[[mode]]
    if(is.null(beta)) {
      probs[[mode]] <- root
    } else {
      crit8_need(isTRUE(model$accepted[[mode]]) &&
                 !is.null(ds[[mode]]),
                 "unaccepted fitted model or missing features")
      probs[[mode]] <- crit7_correction_prob(
             root,ds[[mode]],beta,model$cfg)
    }
  }
  crit8_need(identical(names(probs),crit8_modes) &&
    all(vapply(probs,function(p)
      length(p)==nrow(X) && all(is.finite(p)) &&
      all(p>=0&p<=1),logical(1))),"invalid prediction vectors")
  leaf_id <- if(is.null(model$gate))rep(0L,nrow(X)) else
    ifelse(crit7_relations(X,model$gate)[,1L]>0,2L,1L)
  list(probs=probs,leaf_id=leaf_id,
       accepted=model$accepted,validation=model$validation)
}
