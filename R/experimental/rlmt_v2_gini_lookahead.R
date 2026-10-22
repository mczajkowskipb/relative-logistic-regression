# RLMT v2: two *real* RLR-leaf tree selectors, separately from H2 offset.
# TSPDT_GINI: binary Gini reduction, 4096 deterministic rel. proposals.
# RLMT_LA: top20 residual-screen proposals, real inner-OOF RLR-leaf
#           refits on discovery only, then independent V confirmation.
#
# Both use exactly rlmt_v2_fit_tree(..., mode="H1") on D and T.
# Superset of, but NEVER modifies, the frozen RLR and v2 H1/H2 code.
# Source after rlmt_v2_aligned_offset.R.

rlmt_la_gini <- function(y) {
  if (!length(y)) return(0)
  p <- as.numeric(table(y))/length(y)
  1-sum(p^2)
}

rlmt_la_candidates <- function(X,y,conf,dc,method) {
  if(!method %in% c("TSPDT_GINI","RLMT_LA"))
    stop("RLMT_LA_STOP: unsupported candidate selection mode")
  pairs <- rlmt_dc_random_pairs(
    names(X),n_pairs=dc$proposal_n,seed=conf$seed+201L)
  if(method=="TSPDT_GINI") {
    audit <- vector("list",nrow(pairs))
    root <- rlmt_la_gini(y)
    for(i in seq_len(nrow(pairs))) {
      z <- rlmt_v2_route(X,pairs[i,,drop=FALSE])
      admissible <- rlmt_children_ok(z,y,conf)
      benefit <- if(admissible) root-
        (sum(z)*rlmt_la_gini(y[z])+
         sum(!z)*rlmt_la_gini(y[!z]))/length(y) else NA_real_
      audit[[i]] <- data.frame(
        pool_index=pairs$pool_index[[i]],
        gene_a=pairs$gene_a[[i]],gene_b=pairs$gene_b[[i]],
        feasible=admissible,impurity_gain=benefit)
    }
    all <- do.call(rbind,audit)
    keep <- which(all$feasible & is.finite(all$impurity_gain))
    if(!length(keep))return(list(candidate=NULL,audit=all,
                                 screening="Gini impurity reduction"))
    selected <- keep[order(-all$impurity_gain[keep],
                             all$pool_index[keep])][1L]
    return(list(candidate=all[selected,,drop=FALSE],audit=all,
                screening="Gini impurity reduction"))
  }
  parent <- rlmt_v2_parent_oof(X,y,conf)
  ranked <- rlmt_dc_propose(X,y,parent$prob,pairs,conf)
  admissible <- which(ranked$audit$feasible &
                      is.finite(ranked$audit$residual_screen_score))
  if(!length(admissible))
    return(list(candidate=NULL,audit=ranked$audit,
                screening="Residual rank shortlist then full leaf OOF"))
  budget <- min(length(admissible),as.integer(conf$screen_top_B))
  top <- admissible[order(
    -ranked$audit$residual_screen_score[admissible],
    ranked$audit$pool_index[admissible])][seq_len(budget)]
  shortlisted <- ranked$audit[top,,drop=FALSE]
  if(nrow(shortlisted)>20L)
    stop("RLMT_LA_STOP: look-ahead budget > 20 candidates")

  # Cache inner-training supervised candidate pools once per OOF fold.
  # Importantly, never rank using D-wide labels inside heldout folds.
  cache <- lapply(parent$folds,function(f) {
    tr <- f$train_idx
    list(train=tr,test=f$test_idx,
         pool=compute_top_pairs_switchbox(
           X[tr,,drop=FALSE],y[tr],top_n=conf$pool,
           buffer_mult=5L,out_dir=NULL,save_all=FALSE))
  })
  scored <- vector("list",nrow(shortlisted))
  y01 <- rlmt_y01(y)
  for(k in seq_len(nrow(shortlisted))) {
    candidate <- shortlisted[k,,drop=FALSE]
    predictions <- rep(NA_real_,nrow(X))
    complete <- TRUE
    reason <- ""
    for(j in seq_along(cache)) {
      part <- cache[[j]]
      tr <- part$train
      va <- part$test
      route <- rlmt_v2_route(X[tr,,drop=FALSE],candidate)
      if(!rlmt_children_ok(route,y[tr],conf)) {
        complete <- FALSE
        reason <- "inner_train_child_support"
        break
      }
      left <- tr[route]
      right <- tr[!route]
      result <- tryCatch({
        m1 <- rlmt_v2_fit_h1_leaf(
          X[left,,drop=FALSE],y[left],part$pool,conf)
        m2 <- rlmt_v2_fit_h1_leaf(
          X[right,,drop=FALSE],y[right],part$pool,conf)
        z <- rlmt_v2_route(X[va,,drop=FALSE],candidate)
        prob <- rep(NA_real_,length(va))
        if(any(z)) prob[z] <- rlmt_predict_leaf(
          m1,X[va[z],,drop=FALSE],conf)$prob
        if(any(!z)) prob[!z] <- rlmt_predict_leaf(
          m2,X[va[!z],,drop=FALSE],conf)$prob
        prob
      },error=function(e)e)
      if(inherits(result,"error") || anyNA(result) ||
         any(!is.finite(result)) || any(result<0|result>1)) {
        complete <- FALSE
        reason <- if(inherits(result,"error"))
          conditionMessage(result) else "invalid_oof_predictions"
        break
      }
      predictions[va] <- result
    }
    if(complete && anyNA(predictions))
      stop("RLMT_LA_STOP: incomplete OOF prediction grid")
    gain <- if(complete) {
      p0 <- pmin(1-1e-15,pmax(1e-15,parent$prob))
      p1 <- pmin(1-1e-15,pmax(1e-15,predictions))
      per <- -(y01*log(p0)+(1-y01)*log1p(-p0))+
             (y01*log(p1)+(1-y01)*log1p(-p1))
      c(mean(per),stats::sd(per)/sqrt(length(per)))
    } else c(NA_real_,NA_real_)
    scored[[k]] <- data.frame(
      pool_index=candidate$pool_index[[1L]],
      gene_a=candidate$gene_a[[1L]],
      gene_b=candidate$gene_b[[1L]],
      residual_screen_score=candidate$residual_screen_score[[1L]],
      full_r_lr_leaf_oof=complete,
      mean_gain=gain[[1L]],
      se_gain=gain[[2L]],
      lower_gain=gain[[1L]]-gain[[2L]],
      failure_reason=reason,stringsAsFactors=FALSE)
  }
  audit <- do.call(rbind,scored)
  eligible <- which(audit$full_r_lr_leaf_oof &
                    is.finite(audit$lower_gain) &
                    audit$mean_gain>0 & audit$lower_gain>0)
  selected <- if(length(eligible)) {
    pick <- eligible[order(-audit$lower_gain[eligible],
                           audit$pool_index[eligible])][1L]
    audit[pick,,drop=FALSE]
  } else NULL
  list(candidate=selected,audit=audit,
       screened_pairs=nrow(pairs),lookahead_fitted=nrow(audit),
       screening="strict discovery OOF RLR leaf lookahead")
}

rlmt_la_run <- function(X_train,y_train,X_test,y_test,out_dir,
                         method=c("RLMT_LA","TSPDT_GINI")) {
  method <- match.arg(method)
  opt <- rlmt_v2_conf()
  conf <- opt$leaf
  dc <- opt$dc
  rlmt_v2_check_xy(X_train,y_train)
  if(!is.data.frame(X_test) ||
     !identical(names(X_train),names(X_test)) ||
     nrow(X_test)!=length(y_test) ||
     !is.factor(y_test) ||
     !identical(levels(y_test),levels(y_train)) ||
     anyNA(y_test) ||
     !all(vapply(X_test,function(x)
       is.numeric(x)&&all(is.finite(x)),logical(1))) ||
     is.null(rownames(X_test)) ||
     anyDuplicated(rownames(X_test)))
    stop("RLMT_LA_STOP: invalid test input schema")
  if(file.exists(out_dir) || dir.exists(out_dir))
    stop("RLMT_LA_STOP: output root exists")
  ensure_dir(out_dir)
  dec <- list(method=method,parent="RLR_ALIGNED_V2",
              status="NO_SPLIT",reason="insufficient_training_support",
              candidate=NULL,confirmation=NULL,
              candidate_budget=dc$proposal_n,
              lookahead_budget=conf$screen_top_B)
  if(nrow(X_train)>=3L*conf$min_child_n &&
     all(table(y_train)>=3L*conf$min_child_class)) {
    phase <- rlmt_dc_partition(y_train,seed=conf$seed)
    D <- phase$discovery
    V <- phase$confirmation
    write.csv(data.frame(sample_id=rownames(X_train)[D]),
              file.path(out_dir,"discovery_ids.csv"),row.names=FALSE)
    write.csv(data.frame(sample_id=rownames(X_train)[V]),
              file.path(out_dir,"confirmation_ids.csv"),row.names=FALSE)
    Xd <- X_train[D,,drop=FALSE]
    yd <- y_train[D]
    Xv <- X_train[V,,drop=FALSE]
    if(nrow(Xd)>=2L*conf$min_child_n &&
       all(table(yd)>=2L*conf$min_child_class) &&
       length(V)>=dc$confirmation_min_n) {
      chosen <- rlmt_la_candidates(Xd,yd,conf,dc,method)
      write.csv(chosen$audit,file.path(out_dir,"candidate_audit.csv"),
                row.names=FALSE)
      dec$discovery_screen <- chosen$screening
      if(!is.null(chosen$candidate)) {
        cand <- chosen$candidate
        dec$candidate <- list(gene_a=as.character(cand$gene_a[[1L]]),
                              gene_b=as.character(cand$gene_b[[1L]]))
        fitted <- tryCatch(rlmt_v2_fit_tree(
          Xd,yd,cand,"H1",conf),error=function(e)e)
        if(inherits(fitted,"error")) {
          dec$reason <- paste0("discovery_RLR_leaf_fit_failed:",
                               conditionMessage(fitted))
        } else {
          val <- rlmt_v2_predict_tree(fitted,Xv,conf)
          confirm <- rlmt_dc_confirm(y_train[V],val$parent,val$prob,
                                     min_n=dc$confirmation_min_n)
          dec$confirmation <- confirm
          if(isTRUE(confirm$accepted)) {
            dec$status <- "ACCEPTED"
            dec$reason <- "independent_training_confirmation"
          } else dec$reason <- "confirmation_rejected"
        }
      } else dec$reason <- "no_feasible_discovery_candidate"
    } else dec$reason <- "discovery_or_confirmation_support_insufficient"
  }
  if(identical(dec$status,"ACCEPTED")) {
    pair <- data.frame(gene_a=dec$candidate$gene_a,
                       gene_b=dec$candidate$gene_b)
    model <- rlmt_v2_fit_tree(X_train,y_train,pair,"H1",conf)
    final <- rlmt_v2_predict_tree(model,X_test,conf)
    saveRDS(model,file.path(out_dir,"complete_tree.rds"))
  } else {
    model <- rlmt_v2_fit_parent(X_train,y_train,conf)
    pp <- rlmt_v2_predict_parent(model,X_test,conf)
    final <- list(prob=pp,parent=pp,
                  route=rep(1L,nrow(X_test)))
    saveRDS(model,file.path(out_dir,"complete_parent.rds"))
  }
  if(identical(dec$status,"NO_SPLIT") &&
     !identical(final$parent,final$prob))
    stop("RLMT_LA_STOP: fallback is not exact aligned parent")
  write.csv(data.frame(sample_id=rownames(X_test),
    y_true=as.character(y_test),
    y_pred=ifelse(final$prob>=0.5,levels(y_test)[2L],
                   levels(y_test)[1L]),
    prob=final$prob,parent_prob=final$parent,
    leaf_id=final$route),
    file.path(out_dir,"predictions_test.csv"),row.names=FALSE)
  rlmt_dc_write_json(dec,file.path(out_dir,"decision.json"))
  rlmt_dc_write_json(list(method=method,
    status=dec$status,
    logloss=rlmt_logloss(rlmt_y01(y_test),final$prob),
    aligned_parent_logloss=rlmt_logloss(
      rlmt_y01(y_test),final$parent),
    brier=mean((rlmt_y01(y_test)-final$prob)^2)),
    file.path(out_dir,"metrics_test.json"))
  invisible(list(prob=final$prob,parent=final$parent,
                 route=final$route,decision=dec))
}
