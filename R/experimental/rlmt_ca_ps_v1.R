# RLMT-CA-PS-v1 — EXPERIMENTAL class-aware, parent-shrunk 2-leaf RLR tree.
# Scientific identity: NEW method. NOT published RLR, RLMT-DC-v0 or TSPDT-2011.
# Source after 07_world_protocol_utils.R, 28_rlmt_utils.R,
# 40_rlmt_dc_utils.R, and 42_rlmt_dc_run.R.
#
# No outer-test labels enter nomination, split screening, mixing weights
# or the internal independent confirmation gate.
#
# Defined once before any biomedical fit; NOT tuned on breast2 outcomes:
# - proposal: unchanged deterministic 4096 pairs from frozen RLMT-DC;
# - feasibility: each child >=16 samples AND >=6 from EACH class;
# - score: discovery residual-contrast times sqrt(min balancedness index);
# - child balance index B=4*n0*n1/(n0+n1)^2, 0<=B<=1;
# - shrink weight g=4*n0*n1/[n*(n+2*K)] (K=10 from RLR method).
# - probability p = (1-g)*p_parent + g*p_child, using node g.
# - acceptance: independent confirmation gain-SE > 0 on MIXED predictions.
# - final root RLR is fitted by the unmodified ORIGINAL RLR function.
#
# Balancedness is a soft support reliability penalty. Deliberately NOT
# mandatory 50:50 class mixing; biologically real class enrichment is allowed.

rlmt_ca_ps_config <- function(config=list()) {
  if (length(config)!=0L)
    stop("RLMT_CA_PS_V1: fixed research contract; no tuning overrides")
  opt <- rlmt_dc_defaults(list())
  if (opt$leaf$K != 10L || opt$leaf$min_child_n != 16L ||
      opt$leaf$min_child_class != 6L ||
      opt$dc$proposal_n != 4096L ||
      opt$dc$confirmation_min_n != 20L)
    stop("RLMT_CA_PS_V1: frozen source configuration changed")
  opt
}

rlmt_ca_ps_node_balance <- function(n0,n1) {
  if (length(n0)!=1L || length(n1)!=1L ||
      !is.finite(n0) || !is.finite(n1) ||
      n0<0 || n1<0 ||
      n0!=as.integer(n0) || n1!=as.integer(n1))
    stop("RLMT_CA_PS_V1: invalid nonnegative class counts")
  if (n0+n1==0) return(0)
  as.numeric(4*n0*n1/(n0+n1)^2)
}

rlmt_ca_ps_node_weight <- function(n0,n1,K=10L) {
  if (length(K)!=1L || !is.finite(K) ||
      K!=as.integer(K) || K<1L)
    stop("RLMT_CA_PS_V1: invalid shrinkage complexity K")
  n <- n0+n1
  b <- rlmt_ca_ps_node_balance(n0,n1)
  as.numeric(b*n/(n+2*K))
}

rlmt_ca_ps_probs <- function(parent_prob,child_prob,z_test,
                             w_true,w_false) {
  if (!is.numeric(parent_prob) || !is.numeric(child_prob) ||
      length(parent_prob)!=length(child_prob) ||
      !is.logical(z_test) || length(z_test)!=length(parent_prob) ||
      anyNA(z_test) || anyNA(parent_prob) || anyNA(child_prob) ||
      any(!is.finite(parent_prob)) || any(!is.finite(child_prob)) ||
      any(parent_prob<0|parent_prob>1) ||
      any(child_prob<0|child_prob>1) ||
      length(w_true)!=1L || length(w_false)!=1L ||
      !is.finite(w_true) || !is.finite(w_false) ||
      w_true<0 || w_true>1 || w_false<0 || w_false>1)
    stop("RLMT_CA_PS_V1: invalid independent shrink inputs")
  w <- ifelse(z_test,w_true,w_false)
  out <- parent_prob+w*(child_prob-parent_prob)
  if (anyNA(out) || any(!is.finite(out)) || any(out<0|out>1))
    stop("RLMT_CA_PS_V1: invalid shrunk probabilities")
  as.numeric(out)
}

rlmt_ca_ps_propose <- function(X,y,parent_oof,candidate_pairs,conf) {
  if (!is.data.frame(X) || !is.factor(y) ||
      length(levels(y))!=2L || nrow(X)!=length(y) ||
      !is.numeric(parent_oof) || length(parent_oof)!=nrow(X) ||
      anyNA(parent_oof) || any(!is.finite(parent_oof)) ||
      any(parent_oof<0|parent_oof>1) ||
      !all(c("pool_index","gene_a","gene_b") %in% names(candidate_pairs)))
    stop("RLMT_CA_PS_V1: invalid discovery input")
  residual <- rlmt_y01(y)-parent_oof
  rows <- vector("list",nrow(candidate_pairs))
  for (i in seq_len(nrow(candidate_pairs))) {
    a <- as.character(candidate_pairs$gene_a[i])
    b <- as.character(candidate_pairs$gene_b[i])
    if (!a %in% names(X) || !b %in% names(X) || a==b)
      stop("RLMT_CA_PS_V1: bad candidate gene")
    z <- X[[a]]>X[[b]]
    if (anyNA(z)) stop("RLMT_CA_PS_V1: missing split route")
    cleft <- as.integer(table(factor(y[z],levels=levels(y))))
    cright <- as.integer(table(factor(y[!z],levels=levels(y))))
    nl <- sum(cleft)
    nr <- sum(cright)
    feasible <- (nl>=conf$min_child_n &&
                 nr>=conf$min_child_n &&
                 all(cleft>=conf$min_child_class) &&
                 all(cright>=conf$min_child_class))
    bal_l <- rlmt_ca_ps_node_balance(cleft[1L],cleft[2L])
    bal_r <- rlmt_ca_ps_node_balance(cright[1L],cright[2L])
    contrast <- if (feasible) {
      abs(mean(residual[z])-mean(residual[!z]))*
        sqrt(nl*nr/length(y))
    } else NA_real_
    penalty <- sqrt(min(bal_l,bal_r))
    adjusted <- if (feasible) contrast*penalty else NA_real_
    rows[[i]] <- data.frame(
      pool_index=as.integer(candidate_pairs$pool_index[i]),
      gene_a=a,gene_b=b,feasible=feasible,
      n_left=nl,n_right=nr,
      left_class0=cleft[1L],left_class1=cleft[2L],
      right_class0=cright[1L],right_class1=cright[2L],
      class_balance_left=bal_l,class_balance_right=bal_r,
      unadjusted_residual_score=contrast,
      balance_penalty=penalty,
      residual_screen_score=adjusted,
      stringsAsFactors=FALSE)
  }
  if (!length(rows)) return(list(candidate=NULL,audit=data.frame(),
                                feasible_n=0L))
  tab <- do.call(rbind,rows)
  eligible <- which(tab$feasible & is.finite(tab$residual_screen_score))
  if (!length(eligible)) return(list(candidate=NULL,audit=tab,
                                     feasible_n=0L))
  rank <- eligible[order(-tab$residual_screen_score[eligible],
                          tab$pool_index[eligible])]
  list(candidate=tab[rank[1L],,drop=FALSE],
       audit=tab,feasible_n=length(eligible))
}

rlmt_ca_ps_support <- function(z,y,K) {
  if (!is.logical(z) || length(z)!=length(y) ||
      anyNA(z) || !is.factor(y) ||
      length(levels(y))!=2L || anyNA(y))
    stop("RLMT_CA_PS_V1: invalid split support")
  tab_true <- as.integer(table(factor(y[z],levels=levels(y))))
  tab_false <- as.integer(table(factor(y[!z],levels=levels(y))))
  list(counts_true=tab_true,counts_false=tab_false,
       weight_true=rlmt_ca_ps_node_weight(tab_true[1L],tab_true[2L],K),
       weight_false=rlmt_ca_ps_node_weight(tab_false[1L],tab_false[2L],K))
}

rlmt_ca_ps_check_data <- function(Xtr,ytr,Xte,yte,pool) {
  if (!is.data.frame(Xtr) || !is.data.frame(Xte) ||
      nrow(Xtr)!=length(ytr) || nrow(Xte)!=length(yte) ||
      !is.factor(ytr) || !is.factor(yte) ||
      !identical(levels(ytr),levels(yte)) ||
      length(levels(ytr))!=2L || anyNA(ytr) || anyNA(yte) ||
      !identical(names(Xtr),names(Xte)) ||
      anyDuplicated(names(Xtr)) ||
      any(!vapply(Xtr,is.numeric,logical(1))) ||
      any(!vapply(Xte,is.numeric,logical(1))) ||
      any(!vapply(Xtr,function(x)all(is.finite(x)),logical(1))) ||
      any(!vapply(Xte,function(x)all(is.finite(x)),logical(1))) ||
      !is.data.frame(pool) ||
      !all(c("gene_left","gene_right","ordered_pair") %in% names(pool)))
    stop("RLMT_CA_PS_V1: training/test schema mismatch")
  invisible(TRUE)
}

rlmt_ca_ps_metrics <- function(y,prob) {
  y01 <- rlmt_y01(y)
  label <- ifelse(prob>=0.5,levels(y)[2L],levels(y)[1L])
  pred <- factor(label,levels=levels(y))
  metrics <- vapply(levels(y),function(cls) {
    tp <- sum(pred==cls & y==cls)
    fp <- sum(pred==cls & y!=cls)
    fn <- sum(pred!=cls & y==cls)
    c(recall=if(sum(y==cls)>0L) tp/sum(y==cls) else NA_real_,
      f1=if(2*tp+fp+fn>0L) 2*tp/(2*tp+fp+fn) else NA_real_)
  },numeric(2))
  list(class=pred,accuracy=mean(pred==y),
       balanced_accuracy=if(anyNA(metrics[1L,])) NA_real_ else mean(metrics[1L,]),
       macro_f1=if(anyNA(metrics[2L,])) NA_real_ else mean(metrics[2L,]),
       logloss=rlmt_logloss(y01,prob),
       brier=mean((y01-prob)^2))
}

rlmt_ca_ps_v1_run <- function(X_train,y_train,X_test,y_test,
                              outer_training_pool,out_dir,
                              config=list(),extra_metrics=list()) {
  opt <- rlmt_ca_ps_config(config)
  conf <- opt$leaf
  dc <- opt$dc
  rlmt_ca_ps_check_data(X_train,y_train,X_test,y_test,
                        outer_training_pool)
  if (dir.exists(out_dir))
    stop("RLMT_CA_PS_V1: refusing to overwrite any existing model root")
  if (nrow(outer_training_pool)<conf$K)
    stop("RLMT_CA_PS_V1: outer-training pair pool too small")
  ensure_dir(out_dir)
  ensure_dir(file.path(out_dir,"ca_ps"))
  phase <- NULL
  nominated <- NULL
  confirmed <- NULL
  why <- "insufficient_outer_support"
  if (nrow(X_train)>=3L*conf$min_child_n &&
      all(table(y_train)>=3L*conf$min_child_class)) {
    phase <- rlmt_dc_partition(y_train,seed=conf$seed)
    D <- phase$discovery
    V <- phase$confirmation
    write.csv(data.frame(sample_id=rownames(X_train)[D]),
              file.path(out_dir,"ca_ps","discovery_ids.csv"),
              row.names=FALSE)
    write.csv(data.frame(sample_id=rownames(X_train)[V]),
              file.path(out_dir,"ca_ps","confirmation_ids.csv"),
              row.names=FALSE)
    Xd <- X_train[D,,drop=FALSE]
    yd <- y_train[D]
    Xv <- X_train[V,,drop=FALSE]
    if (nrow(Xd)>=2L*conf$min_child_n &&
        all(table(yd)>=2L*conf$min_child_class) &&
        length(V)>=dc$confirmation_min_n) {
      disc_pool <- compute_top_pairs_switchbox(
        Xd,yd,top_n=conf$pool,buffer_mult=5L,
        out_dir=file.path(out_dir,"ca_ps","discovery_pool"),
        save_all=FALSE)
      discovery_cv <- make_repeated_stratified_folds(
        yd,dc$discovery_oof_folds,1L,seed=conf$seed+101L)
      parent_oof <- rlmt_dc_parent_oof_strict(
        Xd,yd,conf,discovery_cv)
      proposals <- rlmt_dc_random_pairs(
        names(Xd),n_pairs=dc$proposal_n,seed=conf$seed+201L)
      candidate <- rlmt_ca_ps_propose(
        Xd,yd,parent_oof,proposals,conf)
      write.csv(candidate$audit,file.path(out_dir,"ca_ps",
                                         "discovery_candidates.csv"),
                row.names=FALSE)
      nominated <- candidate$candidate
      why <- if(is.null(nominated)) "no_class_support_feasible_candidate"
             else "nominated"
      if (!is.null(nominated)) {
        write.csv(nominated,file.path(out_dir,"ca_ps",
                                     "nominated_split.csv"),
                  row.names=FALSE)
        # Nomination is now frozen; confirmation labels not yet accessed.
        predv <- tryCatch(
          rlmt_dc_parent_and_child_on_validation(
            Xd,yd,Xv,nominated,disc_pool,conf),
          error=function(e)e)
        if (inherits(predv,"error")) {
          why <- paste0("discovery_leaf_fit_failed:",
                        conditionMessage(predv))
        } else {
          zd <- Xd[[nominated$gene_a[[1L]]]] >
                Xd[[nominated$gene_b[[1L]]]]
          zv <- Xv[[nominated$gene_a[[1L]]]] >
                Xv[[nominated$gene_b[[1L]]]]
          w <- rlmt_ca_ps_support(zd,yd,conf$K)
          mixed <- rlmt_ca_ps_probs(
            predv$parent_prob,predv$split_prob,zv,
            w$weight_true,w$weight_false)
          confirmed <- rlmt_dc_confirm(y_train[V],
            predv$parent_prob,mixed,
            min_n=dc$confirmation_min_n)
          confirmed$weight_true_D <- w$weight_true
          confirmed$weight_false_D <- w$weight_false
          why <- if(isTRUE(confirmed$accepted))
            "confirmed_mixed_gain_above_1se" else
            "confirmation_mixed_gain_not_above_1se"
        }
      }
    } else {
      why <- "insufficient_discovery_or_confirmation_support"
    }
  }

  # Exact original published RLR parent, no alternative approximation.
  # Saves a complete unmodified reference fit under a DISTINCT subdirectory.
  parent_dir <- file.path(out_dir,"rlr_reference")
  fit_rlr_train_test(
    X_train,y_train,X_test,y_test,outer_training_pool,
    K=conf$K,m=conf$m,d=conf$d,out_dir=parent_dir,
    alpha=conf$alpha,lambda_choice=conf$lambda_choice,
    seed=conf$seed,extra_metrics=extra_metrics)
  refpred <- read.csv(file.path(parent_dir,"predictions_test.csv"),
                      stringsAsFactors=FALSE,
                      colClasses=c(sample_id="character",y_true="character",
                                   y_pred="character"))
  if (!identical(as.character(refpred$sample_id),rownames(X_test)) ||
      !identical(as.character(refpred$y_true),as.character(y_test)) ||
      anyNA(refpred$prob) ||
      any(!is.finite(refpred$prob)))
    stop("RLMT_CA_PS_V1: original RLR reference schema invalid")
  final_p <- as.numeric(refpred$prob)
  accepted <- !is.null(confirmed) && isTRUE(confirmed$accepted)
  support <- NULL
  route <- rep(1L,nrow(X_test))
  if (accepted) {
    a <- as.character(nominated$gene_a[[1L]])
    b <- as.character(nominated$gene_b[[1L]])
    z_train <- X_train[[a]]>X_train[[b]]
    z_test <- X_test[[a]]>X_test[[b]]
    if (!rlmt_children_ok(z_train,y_train,conf))
      stop("RLMT_CA_PS_V1: final node class support invariant")
    support <- rlmt_ca_ps_support(z_train,y_train,conf$K)
    train_ix <- list(which(z_train),which(!z_train))
    leaf <- lapply(seq_along(train_ix),function(i) {
      ix <- train_ix[[i]]
      rlmt_fit_leaf(X_train[ix,,drop=FALSE],y_train[ix],
                    outer_training_pool,conf,rerank=TRUE,
                    seed=conf$seed+800L+i)
    })
    leaf_p <- final_p
    for (i in 1:2) {
      test_ix <- which(if (i==1L) z_test else !z_test)
      rlmt_save_leaf(leaf[[i]],file.path(out_dir,"leaf_models",
                                      paste0("node_",i+1L)))
      if (length(test_ix)) {
        leaf_p[test_ix] <- rlmt_predict_leaf(
          leaf[[i]],X_test[test_ix,,drop=FALSE],conf)$prob
      }
    }
    final_p <- rlmt_ca_ps_probs(
      final_p,leaf_p,z_test,
      support$weight_true,support$weight_false)
    route <- ifelse(z_test,2L,3L)
  }
  result <- rlmt_ca_ps_metrics(y_test,final_p)
  write.csv(data.frame(sample_id=rownames(X_test),
                       y_true=as.character(y_test),
                       y_pred=as.character(result$class),
                       prob=final_p,score=final_p,leaf_id=route),
            file.path(out_dir,"predictions_test.csv"),row.names=FALSE)
  write.csv(data.frame(sample_id=rownames(X_test),leaf_id=route),
            file.path(out_dir,"test_routes.csv"),row.names=FALSE)
  decision <- list(
    method="RLMT_CA_PS_V1_K10_m10_d2_pool1000_L2",
    identity="new_class_aware_parent_shrunk_not_original_RLMT_DC",
    status=if(accepted) "ACCEPTED" else "NO_SPLIT",
    reason=why,
    discovery_n=if(is.null(phase)) 0L else length(phase$discovery),
    confirmation_n=if(is.null(phase)) 0L else length(phase$confirmation),
    nominated_pair=if(is.null(nominated)) NULL
      else list(gene_a=nominated$gene_a[[1L]],
                gene_b=nominated$gene_b[[1L]],
                adjusted_score=nominated$residual_screen_score[[1L]],
                balance_left=nominated$class_balance_left[[1L]],
                balance_right=nominated$class_balance_right[[1L]]),
    confirmation=confirmed,
    final_leaf_support=support,
    realized_leaves=if(accepted) 2L else 1L,
    decision_probability_threshold=0.5)
  rlmt_dc_write_json(decision,file.path(out_dir,"ca_ps","decision.json"))
  metrics <- c(list(
    method_family="RLMT-CA-PS-v1",
    method="RLMT_CA_PS_V1_K10_m10_d2_pool1000_L2",
    realized_leaves=decision$realized_leaves,
    accepted_split=accepted,
    balanced_accuracy=result$balanced_accuracy,
    macro_f1=result$macro_f1,
    accuracy=result$accuracy,
    logloss=result$logloss,brier=result$brier,
    probability_rule="parent_plus_support_shrunk_leaf_difference",
    score_rule="discovery_residual_contrast_times_child_balance_factor"),
    extra_metrics)
  rlmt_dc_write_json(metrics,file.path(out_dir,"metrics_test.json"))
  invisible(list(status="OK",accepted_split=accepted,
                 realized_leaves=decision$realized_leaves,
                 metrics=metrics,decision=decision,
                 prob=final_p,class=result$class))
}
