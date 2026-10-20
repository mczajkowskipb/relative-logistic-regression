# RLMT-DC-v0: one-split full training workflow (NEW method).
# Must be sourced after 07_world_protocol_utils.R, 28_rlmt_utils.R
# and 40_rlmt_dc_utils.R. No changes to published RLR / RLMT v1.

rlmt_dc_write_json <- function(value, path) {
  jsonlite::write_json(value, path, auto_unbox=TRUE, pretty=TRUE)
}

rlmt_dc_record <- function(path, state) {
  ensure_dir(file.path(path,"dc"))
  rlmt_dc_write_json(state,file.path(path,"dc","decision.json"))
}

rlmt_dc_root <- function(X_train,y_train,X_test,y_test,global_pool,
                         conf,out_dir,reason,extra_metrics) {
  method <- "RLMT_DC_v0_K10_m10_d2_pool1000_L2"
  rlmt_root_output(X_train,y_train,X_test,y_test,global_pool,
                   conf,out_dir,method,extra_metrics,reason)
  file <- file.path(out_dir,"metrics_test.json")
  metrics <- jsonlite::read_json(file,simplifyVector=TRUE)
  metrics$method_family <- "RLMT-DC-v0"
  metrics$dc_stop_reason <- reason
  rlmt_dc_write_json(metrics,file)
  invisible(list(accepted=FALSE,reason=reason))
}

rlmt_dc_write_two_leaf_result <- function(
    X_train,y_train,X_test,y_test,global_pool,conf,
    out_dir,proposal,confirmation,extra_metrics) {
  a <- as.character(proposal$gene_a[[1L]])
  b <- as.character(proposal$gene_b[[1L]])
  z_train <- X_train[[a]] > X_train[[b]]
  if (!rlmt_children_ok(z_train,y_train,conf))
    stop("Full-training candidate does not meet unchanged class support")
  z_test <- X_test[[a]] > X_test[[b]]
  if (anyNA(z_test)) stop("NA in final test routing")

  tr <- list(which(z_train),which(!z_train))
  te <- list(which(z_test),which(!z_test))
  # Both complete outer-training leaf refits must succeed. No alternative
  # candidate or test-dependent model fallback is attempted here.
  leaves <- lapply(seq_along(tr),function(i) {
    rlmt_fit_leaf(X_train[tr[[i]],,drop=FALSE],
                  y_train[tr[[i]]],global_pool,conf,
                  seed=conf$seed+800L+i)
  })
  prob <- rep(NA_real_,nrow(X_test))
  local_active_n <- rep(NA_real_,nrow(X_test))
  for (i in seq_along(te)) {
    ix <- te[[i]]
    rlmt_save_leaf(leaves[[i]],file.path(out_dir,"leaf_models",
                                         paste0("node_",i+1L)))
    if (length(ix)) {
      pred <- rlmt_predict_leaf(leaves[[i]],
                                X_test[ix,,drop=FALSE],conf)
      prob[ix] <- pred$prob
      local_active_n[ix] <- pred$local_active_n
    }
  }
  if (anyNA(prob) || any(!is.finite(prob)) ||
      any(prob < 0 | prob > 1)) stop("Invalid final probability")
  y01 <- rlmt_y01(y_test)
  pred <- factor(ifelse(prob>=0.5,levels(y_test)[2L],
                        levels(y_test)[1L]),levels=levels(y_test))
  tp <- sum(pred==levels(y_test)[2L] & y_test==levels(y_test)[2L])
  fp <- sum(pred==levels(y_test)[2L] & y_test==levels(y_test)[1L])
  fn <- sum(pred==levels(y_test)[1L] & y_test==levels(y_test)[2L])
  precision <- if (tp+fp==0L) 0 else tp/(tp+fp)
  recall <- if (tp+fn==0L) 0 else tp/(tp+fn)
  f1 <- if (precision+recall==0) 0 else
    2*precision*recall/(precision+recall)
  active <- unique(unlist(lapply(leaves,function(m)
    m$active$ordered_pair),use.names=FALSE))
  genes <- unique(unlist(lapply(leaves,function(m)
    c(m$active$gene_left,m$active$gene_right)),
    use.names=FALSE))
  model_name <- "RLMT_DC_v0_K10_m10_d2_pool1000_L2"
  met <- c(list(
    method_family="RLMT-DC-v0",method=model_name,
    requested_K=conf$K,m=conf$m,d=conf$d,pool=conf$pool,
    alpha=conf$alpha,lambda_choice=conf$lambda_choice,
    max_leaves=2L,realized_leaves=2L,
    auc=as.numeric(pROC::auc(y01,prob)),f1=f1,
    logloss=rlmt_logloss(y01,prob),
    brier=calc_brier(prob,y01),ece=calc_ece(prob,y01),
    global_active_pairs=length(active),unique_genes=length(genes),
    local_active_mean=mean(local_active_n),
    local_active_median=median(local_active_n),
    local_active_iqr=IQR(local_active_n),
    dc_confirmation_gain=confirmation$gain,
    dc_confirmation_se_gain=confirmation$se_gain,
    dc_confirmation_lower_gain=confirmation$lower_gain
  ),extra_metrics)
  write.csv(data.frame(
    sample_id=rownames(X_test),y_true=as.character(y_test),
    y_pred=as.character(pred),score=prob,prob=prob),
    file.path(out_dir,"predictions_test.csv"),row.names=FALSE)
  write.csv(data.frame(sample_id=rownames(X_test),
                       local_active_n=local_active_n),
    file.path(out_dir,"local_activity_test.csv"),row.names=FALSE)
  rlmt_dc_write_json(met,file.path(out_dir,"metrics_test.json"))

  node <- data.frame(
    node_id=c(1L,2L,3L),parent_id=c(NA_integer_,1L,1L),
    depth=c(0L,1L,1L),is_leaf=c(FALSE,TRUE,TRUE),
    n=c(length(y_train),length(tr[[1L]]),length(tr[[2L]])),
    class_0=c(sum(y_train==levels(y_train)[1L]),
              sum(y_train[tr[[1L]]]==levels(y_train)[1L]),
              sum(y_train[tr[[2L]]]==levels(y_train)[1L])),
    class_1=c(sum(y_train==levels(y_train)[2L]),
              sum(y_train[tr[[1L]]]==levels(y_train)[2L]),
              sum(y_train[tr[[2L]]]==levels(y_train)[2L])),
    active_relations=c(NA_integer_,
                       nrow(leaves[[1L]]$active),
                       nrow(leaves[[2L]]$active)),
    stringsAsFactors=FALSE)
  write.csv(node,file.path(out_dir,"nodes.csv"),row.names=FALSE)
  write.csv(data.frame(sample_id=rownames(X_test),
                       leaf_id=ifelse(z_test,2L,3L)),
            file.path(out_dir,"test_routes.csv"),row.names=FALSE)
  rlmt_dc_write_json(list(
    method=model_name,size=2L,max_leaves=2L,
    root_id=1L,accepted_split=TRUE,
    split=list(
      gene_a=a,gene_b=b,
      discovery_score=as.numeric(proposal$residual_screen_score),
      confirmation_mean_gain=confirmation$gain,
      confirmation_se_gain=confirmation$se_gain,
      confirmation_lower_gain=confirmation$lower_gain,
      left_child=2L,right_child=3L),
    nodes=list(list(node_id=1L,is_leaf=FALSE),
               list(node_id=2L,is_leaf=TRUE),
               list(node_id=3L,is_leaf=TRUE))
  ),file.path(out_dir,"tree.json"))
  invisible(list(accepted=TRUE,reason="confirmed"))
}

rlmt_dc_run <- function(X_train,y_train,X_test,y_test,
                        outer_training_pool,config=list(),out_dir,
                        extra_metrics=list()) {
  if (dir.exists(out_dir)) stop("RLMT-DC refuses to overwrite model output")
  if (!is.data.frame(X_train) || !is.data.frame(X_test) ||
      !identical(names(X_train),names(X_test)) ||
      !is.factor(y_train) || !is.factor(y_test) ||
      length(levels(y_train))!=2L ||
      !identical(levels(y_train),levels(y_test)) ||
      nrow(X_train)!=length(y_train) ||
      nrow(X_test)!=length(y_test) ||
      anyNA(X_train) || anyNA(X_test) ||
      any(!vapply(X_train,is.numeric,logical(1))) ||
      any(!vapply(X_test,is.numeric,logical(1))) ||
      !identical(rownames(X_train),make.unique(rownames(X_train))) ||
      !identical(rownames(X_test),make.unique(rownames(X_test)))) {
    stop("Invalid full outer-training/test schema")
  }
  if (!all(vapply(X_train,function(x)all(is.finite(x)),logical(1))) ||
      !all(vapply(X_test,function(x)all(is.finite(x)),logical(1))))
    stop("RLMT-DC requires finite numeric expression values")
  opt <- rlmt_dc_defaults(config)
  conf <- opt$leaf
  dc <- opt$dc
  if (nrow(outer_training_pool) < conf$K) {
    stop("Outer-training RLR pool is too small")
  }

  ensure_dir(out_dir)
  ensure_dir(file.path(out_dir,"dc"))
  fallback <- function(reason,context=list()) {
    rlmt_dc_record(out_dir,c(list(
      method="RLMT-DC-v0",status="NO_SPLIT",
      reason=reason,
      discovery_n=length(phase$discovery %||% integer()),
      confirmation_n=length(phase$confirmation %||% integer())
    ),context))
    rlmt_dc_root(X_train,y_train,X_test,y_test,outer_training_pool,
                 conf,out_dir,reason,extra_metrics)
  }

  phase <- NULL
  if (nrow(X_train) < 3L*conf$min_child_n ||
      any(table(y_train)<3L*conf$min_child_class)) {
    return(fallback("insufficient_outer_training_support"))
  }
  phase <- rlmt_dc_partition(y_train,seed=conf$seed)
  D <- phase$discovery
  V <- phase$confirmation
  write.csv(data.frame(sample_id=rownames(X_train)[D]),
            file.path(out_dir,"dc","discovery_ids.csv"),row.names=FALSE)
  write.csv(data.frame(sample_id=rownames(X_train)[V]),
            file.path(out_dir,"dc","confirmation_ids.csv"),row.names=FALSE)
  yD <- y_train[D]
  X_D <- X_train[D,,drop=FALSE]
  X_V <- X_train[V,,drop=FALSE]
  if (nrow(X_D) < 2L*conf$min_child_n ||
      any(table(yD)<2L*conf$min_child_class) ||
      length(V)<dc$confirmation_min_n) {
    return(fallback("insufficient_discovery_or_confirmation_support"))
  }

  # DISCOVERY POOL: recomputed on D only (not inherited from T).
  disc_pool <- compute_top_pairs_switchbox(
    X_train=X_D,y_train=yD,top_n=conf$pool,
    buffer_mult=5L,out_dir=file.path(out_dir,"dc","discovery_pool"),
    save_all=FALSE)
  if (nrow(disc_pool)<conf$K)
    stop("Training-only discovery RLR pool was unexpectedly too small")
  discovery_inner <- make_repeated_stratified_folds(
    yD,dc$discovery_oof_folds,1L,seed=conf$seed+101L)
  parent_oof <- rlmt_parent_oof(X_D,yD,disc_pool,conf,discovery_inner)
  proposals <- rlmt_dc_random_pairs(
    names(X_D),n_pairs=dc$proposal_n,seed=conf$seed+201L)
  choice <- rlmt_dc_propose(
    X_D,yD,parent_oof,proposals,conf)
  write.csv(choice$audit,file.path(out_dir,"dc","discovery_candidates.csv"),
            row.names=FALSE)
  if (is.null(choice$candidate))
    return(fallback("no_discovery_support_feasible_split",
                    list(proposed_n=nrow(proposals),
                         feasible_n=choice$feasible_n)))

  # The one nominated split is frozen before any V labels are read.
  proposal <- choice$candidate
  write.csv(proposal,file.path(out_dir,"dc","nominated_split.csv"),
            row.names=FALSE)
  fit <- tryCatch(
    rlmt_dc_parent_and_child_on_validation(
      X_D,yD,X_V,proposal,disc_pool,conf),
    error=function(e)e)
  if (inherits(fit,"error")) {
    return(fallback("discovery_leaf_fit_failed",
                    list(message=conditionMessage(fit),
                         proposed_n=nrow(proposals),
                         feasible_n=choice$feasible_n,
                         gene_a=proposal$gene_a,
                         gene_b=proposal$gene_b)))
  }

  # Read V labels ONCE to decide acceptance. Do not choose second candidate.
  confirmation <- rlmt_dc_confirm(
    y_train[V],fit$parent_prob,fit$split_prob,
    min_n=dc$confirmation_min_n)
  rlmt_dc_record(out_dir,list(
    method="RLMT-DC-v0",
    status=if (confirmation$accepted) "ACCEPTED" else "REJECTED",
    discovery_n=length(D),confirmation_n=length(V),
    proposed_n=nrow(proposals),feasible_n=choice$feasible_n,
    nominated_gene_a=proposal$gene_a,
    nominated_gene_b=proposal$gene_b,
    nomination_score=proposal$residual_screen_score,
    confirmation=confirmation))

  if (!confirmation$accepted) {
    return(rlmt_dc_root(X_train,y_train,X_test,y_test,
                        outer_training_pool,conf,out_dir,
                        "confirmation_gain_not_above_1SE",extra_metrics))
  }
  # The selected split and acceptance are now fixed.
  rlmt_dc_write_two_leaf_result(
    X_train,y_train,X_test,y_test,outer_training_pool,conf,
    out_dir,proposal,confirmation,extra_metrics)
}
