# RLMT-v2 RESEARCH PROTOTYPE: two separate hypotheses H1/H2.
# Uses existing switchBox/glmnet. NEVER publishes/changes original RLR.
# Source after R/new/07_world_protocol_utils.R, 28_rlmt_utils.R,
# 40_rlmt_dc_utils.R and 42_rlmt_dc_run.R.
#
# H1 = SAME aligned parent/leaf learner procedure for D and T;
#      original DC nomination + D/V 1SE independent confirmation.
# H2 = SAME aligned parent + leaf-specific fixed-ridge OFFSET correction;
#      NO internal offset-lambda CV, parent offsets generated OOF on S.
# Both retain frozen candidate count 4096 and min child (16,6/class).
# This intentionally creates new parent identity "RLR_ALIGNED", NOT RLR.

rlmt_v2_conf <- function() {
  z <- rlmt_dc_defaults(list())
  if (z$leaf$K!=10L || z$leaf$m!=10L ||
      z$leaf$d!=2L || z$leaf$pool!=1000L ||
      z$leaf$alpha!=0.5 || z$leaf$min_child_n!=16L ||
      z$leaf$min_child_class!=6L || z$leaf$seed!=42L ||
      z$dc$proposal_n!=4096L || z$dc$confirmation_min_n!=20L)
    stop("V2_STOP: frozen v2/source constants changed")
  z
}

rlmt_v2_check_xy <- function(X,y) {
  if (!is.data.frame(X) || !is.factor(y) ||
      length(levels(y))!=2L || anyNA(y) ||
      nrow(X)!=length(y) || ncol(X)<10L ||
      !all(vapply(X,is.numeric,logical(1))) ||
      !all(vapply(X,function(v)all(is.finite(v)),logical(1))) ||
      anyDuplicated(colnames(X)) || anyDuplicated(rownames(X)) ||
      is.null(rownames(X)) || any(table(y)<3L))
    stop("V2_STOP: invalid binary numeric expression cohort")
  invisible(TRUE)
}

# A(S): EXACTLY ONE train-only fitted procedure, irrespective of S=D/T.
rlmt_v2_fit_parent <- function(X,y,conf=rlmt_v2_conf()$leaf) {
  rlmt_v2_check_xy(X,y)
  if (any(table(y)<5L)) stop("V2_STOP: parent inner CV needs >=5/class")
  pool <- compute_top_pairs_switchbox(
    X,y,top_n=conf$pool,buffer_mult=5L,
    out_dir=NULL,save_all=FALSE)
  pairs <- apply_disjointness(pool,K=conf$K,d=conf$d)
  if (nrow(pairs)<1L) stop("V2_STOP: aligned parent has no features")
  encoded <- encode_pair_matrix(X,pairs,m=conf$m)
  foldid <- rlmt_stratified_foldid(y,nfolds=5L,seed=conf$seed)
  fit <- glmnet::cv.glmnet(
    x=encoded,y=rlmt_y01(y),family="binomial",
    alpha=conf$alpha,nfolds=5L,foldid=foldid,
    type.measure="deviance")
  list(identity="RLR_ALIGNED_V2",fit=fit,pairs=pairs,pool=pool,
       foldid=foldid,lambda=conf$lambda_choice,
       gene_names=names(X),class_levels=levels(y),
       fitting_procedure="rlmt_v2_fit_parent",
       selected_lambda=as.numeric(fit[[conf$lambda_choice]]))
}
rlmt_v2_predict_parent <- function(mod,X,conf=rlmt_v2_conf()$leaf) {
  if (!identical(names(X),mod$gene_names))
    stop("V2_STOP: aligned parent feature names/order mismatch")
  p <- as.numeric(stats::predict(mod$fit,
    newx=encode_pair_matrix(X,mod$pairs,m=conf$m),
    s=mod$lambda,type="response"))
  if (length(p)!=nrow(X) || anyNA(p) ||
      any(!is.finite(p)) || any(p<0|p>1))
    stop("V2_STOP: invalid aligned parent predictions")
  p
}

# Training-only OOF probability p_{-fold}(X_i) for every S observation.
# Never feed in-sample parent probabilities into H2 offset fitting.
rlmt_v2_parent_oof <- function(X,y,conf=rlmt_v2_conf()$leaf) {
  rlmt_v2_check_xy(X,y)
  if (any(table(y)<9L)) stop("V2_STOP: >=9/class required for honest OOF")
  folds <- make_repeated_stratified_folds(
    y,n_folds=3L,n_repeats=1L,seed=conf$seed+101L)
  p <- rep(NA_real_,nrow(X))
  count <- integer(nrow(X))
  for (f in folds) {
    tr <- f$train_idx
    va <- f$test_idx
    if (length(intersect(tr,va))>0L ||
        any(table(y[tr])<5L)) stop("V2_STOP: invalid OOF fold support")
    fit <- rlmt_v2_fit_parent(X[tr,,drop=FALSE],y[tr],conf)
    p[va] <- rlmt_v2_predict_parent(fit,X[va,,drop=FALSE],conf)
    count[va] <- count[va]+1L
  }
  if (any(count!=1L) || anyNA(p))
    stop("V2_STOP: OOF predictions do not cover S exactly once")
  list(prob=p,folds=folds)
}

rlmt_v2_route <- function(X,candidate) {
  if (!is.data.frame(candidate) || nrow(candidate)!=1L ||
      !all(c("gene_a","gene_b") %in% names(candidate)))
    stop("V2_STOP: invalid nominated pair")
  a <- as.character(candidate$gene_a[[1L]])
  b <- as.character(candidate$gene_b[[1L]])
  if (!a %in% names(X) || !b %in% names(X) || a==b)
    stop("V2_STOP: nominated features not in cohort")
  z <- X[[a]]>X[[b]]
  if(anyNA(z)) stop("V2_STOP: NA route")
  z
}

# H1 leaf: identical local reranking/pool/stratified glmnet on S=D or T.
rlmt_v2_fit_h1_leaf <- function(X,y,pool,conf) {
  rlmt_fit_leaf(X,y,pool,conf,rerank=TRUE,seed=conf$seed+702L)
}

rlmt_v2_logit <- function(p,epsilon=1e-8) {
  if (!is.numeric(p) || anyNA(p) || any(p<0|p>1))
    stop("V2_STOP: invalid logistic offset probability")
  stats::qlogis(pmax(epsilon,pmin(1-epsilon,p)))
}

# H2: fixed ridge penalty lambda=1 (not optimized on prior patient
# cohorts or within-sample leaf CV). No intercept and no unpenalized
# coefficients: parent baseline is a fixed OOF logit offset.
rlmt_v2_fit_h2_leaf <- function(X,y,pool,parent_oof_prob,conf) {
  if (any(table(y)<6L) || nrow(X)<16L ||
      length(parent_oof_prob)!=nrow(X))
    stop("V2_STOP: offset leaf class support/OOF mismatch")
  pairs <- apply_disjointness(rerank_pool_in_node(X,y,pool),
                              K=conf$K,d=conf$d)
  if (nrow(pairs)==0L) stop("V2_STOP: empty offset leaf shortlist")
  R <- encode_pair_matrix(X,pairs,m=conf$m)
  offset <- rlmt_v2_logit(parent_oof_prob)
  fit <- glmnet::glmnet(
    x=R,y=rlmt_y01(y),family="binomial",offset=offset,
    alpha=0,lambda=1,intercept=FALSE,standardize=FALSE,
    penalty.factor=rep(1,ncol(R)),control=list(thresh=1e-8))
  list(identity="H2_FIXED_RIDGE_OFFSET_V2",fit=fit,pairs=pairs,
       lambda=1,alpha=0,intercept=FALSE,standardize=FALSE,
       penalty_factor=rep(1,ncol(R)),parent_oof=parent_oof_prob,
       offset_epsilon=1e-8)
}
rlmt_v2_predict_h2_leaf <- function(mod,X,parent_prob,conf) {
  if (length(parent_prob)!=nrow(X))
    stop("V2_STOP: heldout parent/leaf length mismatch")
  R <- encode_pair_matrix(X,mod$pairs,m=conf$m)
  p <- as.numeric(stats::predict(mod$fit,newx=R,s=mod$lambda,
     type="response",newoffset=rlmt_v2_logit(parent_prob)))
  if (anyNA(p) || any(!is.finite(p)) ||
      any(p<0|p>1)) stop("V2_STOP: invalid offset-corrected probability")
  p
}

# For a pre-nominated split, perform the SAME learning procedure on D or T.
rlmt_v2_fit_tree <- function(X,y,candidate,mode,conf) {
  if (!mode %in% c("H1","H2")) stop("V2_STOP: unknown method")
  z <- rlmt_v2_route(X,candidate)
  if (!rlmt_children_ok(z,y,conf))
    stop("V2_STOP: candidate violates both-classes training support")
  parent <- rlmt_v2_fit_parent(X,y,conf)
  ids <- list(which(z),which(!z))
  parent_oof <- if(mode=="H2") rlmt_v2_parent_oof(X,y,conf) else NULL
  # Reuse the SAME full S-only pool that fitted the aligned root.
  # Keep this pool in the parent artifact so D/T procedure can be audited.
  fullpool <- parent$pool
  leaves <- lapply(ids,function(ix) {
    if(mode=="H1") {
      rlmt_v2_fit_h1_leaf(X[ix,,drop=FALSE],y[ix],fullpool,conf)
    } else {
      rlmt_v2_fit_h2_leaf(X[ix,,drop=FALSE],y[ix],fullpool,
                          parent_oof$prob[ix],conf)
    }
  })
  list(parent=parent,leaves=leaves,mode=mode,
       nominated_pair=candidate[,c("gene_a","gene_b"),drop=FALSE],
       oof=if(is.null(parent_oof)) NULL else list(
         sample_id=rownames(X),prob=parent_oof$prob,
         fold_ids=lapply(parent_oof$folds,function(f)
           list(train=rownames(X)[f$train_idx],
                validation=rownames(X)[f$test_idx]))),
       training_route_count=as.integer(table(factor(z,levels=c(FALSE,TRUE)))),
       train_class_levels=levels(y))
}
rlmt_v2_predict_tree <- function(model,X,conf=rlmt_v2_conf()$leaf) {
  p0 <- rlmt_v2_predict_parent(model$parent,X,conf)
  z <- rlmt_v2_route(X,model$nominated_pair)
  p <- numeric(nrow(X))
  for(i in 1:2) {
    idx <- which(if(i==1L) z else !z)
    if(!length(idx))next
    p[idx] <- if(model$mode=="H1") {
      rlmt_predict_leaf(model$leaves[[i]],X[idx,,drop=FALSE],conf)$prob
    } else {
      rlmt_v2_predict_h2_leaf(model$leaves[[i]],
                              X[idx,,drop=FALSE],p0[idx],conf)
    }
  }
  if(anyNA(p) || any(!is.finite(p)) || any(p<0|p>1))
    stop("V2_STOP: invalid final tree probabilities")
  list(prob=p,parent=p0,route=ifelse(z,2L,3L))
}

# Independent training-only 1SE confirmation; no test labels used
# for any fitting, split decision or resulting probabilities.
rlmt_v2_run <- function(X_train,y_train,X_test,y_test,
                        out_dir,mode=c("H1","H2")) {
  mode <- match.arg(mode)
  opt <- rlmt_v2_conf()
  conf <- opt$leaf
  dc <- opt$dc
  rlmt_v2_check_xy(X_train,y_train)
  if(!is.data.frame(X_test) ||
     !identical(names(X_test),names(X_train)) ||
     nrow(X_test)!=length(y_test) ||
     any(!vapply(X_test,function(x)is.numeric(x)&&all(is.finite(x)),
                logical(1))) ||
     !is.factor(y_test) ||
     !identical(levels(y_test),levels(y_train)) ||
     anyNA(y_test) || anyDuplicated(rownames(X_test)))
    stop("V2_STOP: invalid heldout schema")
  if(file.exists(out_dir) || dir.exists(out_dir))
    stop("V2_STOP: refusing to overwrite any prior run")
  ensure_dir(out_dir)
  decision <- list(method=paste0("RLMT_V2_",mode),
       parent_identity="RLR_ALIGNED_V2",accepted=FALSE,
       reason="insufficient_training_support",
       confirmation=NULL,nominated_pair=NULL)
  if(nrow(X_train)>=3L*conf$min_child_n &&
     all(table(y_train)>=3L*conf$min_child_class)) {
    phase <- rlmt_dc_partition(y_train,seed=conf$seed)
    D <- phase$discovery
    V <- phase$confirmation
    Xd <- X_train[D,,drop=FALSE]
    yd <- y_train[D]
    Xv <- X_train[V,,drop=FALSE]
    write.csv(data.frame(sample_id=rownames(Xd)),
      file.path(out_dir,"discovery_ids.csv"),row.names=FALSE)
    write.csv(data.frame(sample_id=rownames(Xv)),
      file.path(out_dir,"confirmation_ids.csv"),row.names=FALSE)
    if(nrow(Xd)>=2L*conf$min_child_n &&
       all(table(yd)>=2L*conf$min_child_class) &&
       length(V)>=dc$confirmation_min_n) {
      # Common H1/H2 nomination: ORIGINAL DC discovery residual contrast.
      parent_oof <- rlmt_v2_parent_oof(Xd,yd,conf)
      proposals <- rlmt_dc_random_pairs(
        names(Xd),dc$proposal_n,seed=conf$seed+201L)
      chosen <- rlmt_dc_propose(Xd,yd,parent_oof$prob,
                                proposals,conf)
      write.csv(chosen$audit,file.path(out_dir,"discovery_audit.csv"),
        row.names=FALSE)
      decision$feasible_candidates <- chosen$feasible_n
      if(!is.null(chosen$candidate)) {
        choice <- chosen$candidate
        decision$nominated_pair <- list(
          gene_a=choice$gene_a[[1L]],
          gene_b=choice$gene_b[[1L]],
          score=choice$residual_screen_score[[1L]])
        decision$reason <- "candidate_evaluated"
        # One frozen D candidate: fit on D with EXACT SAME helper as on T.
        # V labels are not passed until rlmt_dc_confirm.
        D_model <- tryCatch(rlmt_v2_fit_tree(Xd,yd,choice,mode,conf),
                            error=function(e)e)
        if(inherits(D_model,"error")) {
          decision$reason <- paste0("candidate_D_fit_failed:",
                                   conditionMessage(D_model))
        } else {
          pV <- rlmt_v2_predict_tree(D_model,Xv,conf)
          confit <- rlmt_dc_confirm(
            y_train[V],pV$parent,pV$prob,
            min_n=dc$confirmation_min_n)
          decision$confirmation <- confit
          decision$accepted <- isTRUE(confit$accepted)
          decision$reason <- if(decision$accepted)
              "independent_train_only_confirmation_accepted" else
              "independent_train_only_confirmation_rejected"
        }
      } else decision$reason <- "no_feasible_candidate"
    } else decision$reason <- "insufficient_discovery_support"
  }
  if(decision$accepted) {
    selected <- data.frame(
      gene_a=decision$nominated_pair$gene_a,
      gene_b=decision$nominated_pair$gene_b)
    model <- rlmt_v2_fit_tree(X_train,y_train,selected,mode,conf)
    final <- rlmt_v2_predict_tree(model,X_test,conf)
    saveRDS(model,file.path(out_dir,"complete_fitted_tree.rds"))
    if(!is.null(model$oof))
      saveRDS(model$oof,file.path(out_dir,"parent_oof_training.rds"))
  } else {
    parent <- rlmt_v2_fit_parent(X_train,y_train,conf)
    final <- list(prob=rlmt_v2_predict_parent(parent,X_test,conf),
                  route=rep(1L,nrow(X_test)))
    final$parent <- final$prob
    saveRDS(parent,file.path(out_dir,"complete_aligned_parent.rds"))
  }
  if(!decision$accepted && !identical(final$parent,final$prob))
    stop("V2_STOP: rejected split not bit-exact aligned RLR parent")
  write.csv(data.frame(sample_id=rownames(X_test),
    y_true=as.character(y_test),
    y_pred=ifelse(final$prob>=.5,levels(y_test)[2L],levels(y_test)[1L]),
    prob=final$prob,parent_prob=final$parent,
    leaf_id=final$route),
    file.path(out_dir,"predictions_test.csv"),row.names=FALSE)
  rlmt_dc_write_json(decision,file.path(out_dir,"decision.json"))
  metric <- list(method=paste0("RLMT_V2_",mode),
     parent="RLR_ALIGNED_V2",
     accepted=decision$accepted,
     logloss=rlmt_logloss(rlmt_y01(y_test),final$prob),
     brier=mean((rlmt_y01(y_test)-final$prob)^2),
     parent_logloss=rlmt_logloss(rlmt_y01(y_test),final$parent))
  rlmt_dc_write_json(metric,file.path(out_dir,"metrics_test.json"))
  invisible(list(prob=final$prob,parent=final$parent,
       route=final$route,decision=decision,metrics=metric))
}
