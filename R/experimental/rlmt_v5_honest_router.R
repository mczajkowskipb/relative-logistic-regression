# RLMT-v5 — HONEST TRAIN-ONLY ROUTER TO ORIGINAL RLR.
# DISCOVERY stage is based on *already exposed* v4 development failures,
# so any positive result is a retrospective result, NOT independent proof.
#
# Instead of forcing an extra relational correction for every dataset,
# fit three full ORIGINAL_RLR-anchored alternatives on each inner training
# fold, evaluate their predictions on DISJOINT inner heldout patients,
# and select a low-complexity mixture with the exact original RLR.
# Full parent is fit on all outer-training observations afterward.
#
# Candidate families are the FROZEN v4 algorithms, not newly tuned;
# the meta-router and blend weights are the new methodology.
# All data-dependent decisions are made within OUTER TRAIN ONLY.
#
# A single model has at most 10 original RLR rank pairs PLUS one named
# correction pair, β/γ correction coefficients and one shrinkage weight.
# Weight=0 guarantees EXACT original RLR predictions.

rlmt5_candidates <- c("GLOBAL_CAL","CONF_ROUTE","ACTIVITY_ROUTE")
rlmt5_weights <- c(0.25,0.50,0.75,1.00)

rlmt5_config <- function() list(
  seed=42L,
  inner_folds=3L,
  inner_seed=1403L,
  selected_candidates=rlmt5_candidates,
  blend_grid=rlmt5_weights,
  se_penalty=0.50,
  selection_cost=0.00025,
  min_improvement=0.00025)

rlmt5_loss <- function(y,prob) {
  rlmt4_require(all(is.finite(prob))&&
    length(prob)==length(y)&&all(prob>=0&prob<=1),
    "v5 invalid probability for paired loss")
  rlmt_logloss(rlmt_y01(y),prob)
}

rlmt5_decide <- function(y,fold_id,predictions,cfg=rlmt5_config()) {
  candidates <- cfg$selected_candidates
  modes <- c("ORIGINAL_RLR",candidates)
  rlmt4_require(identical(names(predictions),modes),
                "v5 CV table order or contents diverged")
  n <- length(y)
  folds <- sort(unique(fold_id))
  rlmt4_require(length(folds)==cfg$inner_folds &&
    all(table(fold_id)>0) &&
    all(vapply(predictions,function(p)
      length(p)==n&&all(is.finite(p))&&all(p>=0&p<=1),
      logical(1))),"bad heldout evaluation grid")
  baseline <- predictions$ORIGINAL_RLR
  # Compute paired *fold mean* logloss improvements. The score
  # subtracts an SE penalty across truly heldout validation folds
  # and a frozen "try-a-more-complex-model" charge. Weighted blends
  # are strictly inside [0,1], including exact w=0 default.
  report <- list()
  k <- 0L
  for(method in candidates)for(w in cfg$blend_grid) {
    pr <- (1-w)*baseline+w*predictions[[method]]
    gain <- vapply(folds,function(f) {
      ix <- which(fold_id==f)
      rlmt5_loss(y[ix],baseline[ix])-rlmt5_loss(y[ix],pr[ix])
    },numeric(1))
    fold_average <- mean(gain)
    se <- stats::sd(gain)/sqrt(length(gain))
    conservative <- fold_average-cfg$se_penalty*se-
                        cfg$selection_cost
    k <- k+1L
    report[[k]] <- data.frame(
      method=method,weight=w,mean_gain=fold_average,
      gain_se=se,conservative_gain=conservative,
      inner_fold_gains=paste(sprintf("%+.9f",gain),collapse=";"),
      stringsAsFactors=FALSE)
  }
  tab <- do.call(rbind,report)
  best <- which.max(tab$conservative_gain)
  accept <- is.finite(tab$conservative_gain[[best]]) &&
    tab$conservative_gain[[best]]>cfg$min_improvement
  selected <- if(accept)as.character(tab$method[[best]]) else
                "ORIGINAL_RLR"
  w <- if(accept)as.numeric(tab$weight[[best]]) else 0
  list(model_identity="TRAIN_ONLY_3FOLD_ONE_SE_SHRUNK_ROUTER",
       selected=selected,weight=w,
       evidence=tab,
       best_inner_gain=tab$mean_gain[[best]],
       best_inner_gain_se=tab$gain_se[[best]],
       best_penalized_gain=tab$conservative_gain[[best]],
       accepted=accept,config=cfg)
}

rlmt5_fit <- function(X,y,config=rlmt5_config(),
                      parent_cfg=rlmt4_config()) {
  rlmt4_require(is.factor(y)&&
    nrow(X)==length(y)&&length(levels(y))==2L,
    "v5 train labels or input invalid")
  # When the minority class is too small for nested parent
  # crossfitting, fall back to original RLR rather than invent
  # a leaky replacement classifier.
  inner <- make_repeated_stratified_folds(
    y,config$inner_folds,1L,seed=config$inner_seed)
  counts <- vapply(inner,function(f)
    min(table(y[f$train_idx])),numeric(1))
  fit <- NULL
  if(all(counts>=8L)) {
    n <- nrow(X)
    score <- setNames(lapply(
       c("ORIGINAL_RLR",config$selected_candidates),
       function(.)rep(NA_real_,n)),
       c("ORIGINAL_RLR",config$selected_candidates))
    fold_id <- integer(n)
    full_cover <- integer(n)
    for(i in seq_along(inner)) {
      f <- inner[[i]]
      a <- f$train_idx
      v <- f$test_idx
      rlmt4_require(length(intersect(a,v))==0L,
                    "v5 internal discovery-validation leakage")
      local <- rlmt4_fit(X[a,,drop=FALSE],y[a],parent_cfg)
      holdout <- rlmt4_predict(local,X[v,,drop=FALSE])
      for(method in names(score))
        score[[method]][v] <- holdout[[method]]
      fold_id[v] <- as.integer(i)
      full_cover[v] <- full_cover[v]+1L
    }
    rlmt4_require(all(full_cover==1L) &&
      all(vapply(score,function(z)all(is.finite(z)),logical(1))),
      "v5 strict CV prediction coverage broken")
    fit <- rlmt5_decide(y,fold_id,score,config)
    fit$cross_fitted_sample_ids <- rownames(X)
    fit$inner_fold_ids <- fold_id
    fit$inner_original_proba <- score$ORIGINAL_RLR
    fit$inner_selected_proba <-
      score[[fit$selected]]
    fit$inner_partitions <- lapply(inner,function(f)
      list(train=rownames(X)[f$train_idx],
           validation=rownames(X)[f$test_idx]))
  } else {
    fit <- list(model_identity="STRICT_INNER_SUPPORT_ORIGINAL_ONLY",
        selected="ORIGINAL_RLR",weight=0,accepted=FALSE,
        reason="MIN_CLASS_LT_8_IN_INNER_TRAIN",
        config=config)
  }
  # Refit once on FULL OUTER TRAIN. This parent is exactly the RLR
  # paper baseline, *independently checked against archived predictions*
  # by the external runner. No heldout test label anywhere here.
  full <- rlmt4_fit(X,y,parent_cfg)
  if(fit$selected!="ORIGINAL_RLR")
    rlmt4_require(fit$selected %in% names(rlmt4_predict(
      full,X[seq_len(1L),,drop=FALSE])),
      "v5 selected method not present in full fitted model")
  structure(list(router=fit,full_bundle=full,
      model_identity="RLMT5_HONEST_CV_GATED_RANK_REPAIR",
      originally_fitted_training_names=rownames(X)),
      class="rlmt5_model")
}

rlmt5_predict <- function(model,X) {
  rlmt4_require(inherits(model,"rlmt5_model"),
                "v5 wrong model object type")
  all <- rlmt4_predict(model$full_bundle,X)
  selection <- model$router$selected
  weight <- model$router$weight
  rlmt4_require(selection %in% names(all) &&
                is.finite(weight) && weight>=0 && weight<=1,
                "v5 incomplete saved router")
  original <- all$ORIGINAL_RLR
  correction <- all[[selection]]
  predict <- as.numeric((1-weight)*original+weight*correction)
  rlmt4_require(all(is.finite(predict)) &&
     all(predict>=0&predict<=1),
     "v5 invalid model probabilities")
  if(weight==0) rlmt4_require(identical(predict,original),
    "v5 zero router must be exact ORIGINAL RLR")
  list(prob=predict,original_prob=original,
       model=selection,weight=weight,
       active_pair=if(selection %in% names(model$full_bundle$routes))
         model$full_bundle$routes[[selection]]$pair else NULL,
       candidates=all)
}
