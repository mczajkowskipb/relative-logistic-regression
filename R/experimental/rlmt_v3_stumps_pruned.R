# RLMT v3 — deliberately SIMPLE and genuinely splittable rank-relation
# probability trees. ISOLATED RESEARCH, NOT the 2011 full TSPDT replica.
# Run after 07_world_protocol_utils, 28_rlmt_utils, 40_rlmt_dc_utils,
# 42_rlmt_dc_run and rlmt_v2_aligned_offset.
#
# Three predeclared methods:
#  STUMP_TSP_BETA: always fit 1 relational stump when size-eligible,
#                   Beta-smoothed leaf probabilities, no leaf glmnet.
#  STUMP_RIDGE_OFFSET: same TSP stump, 2 ridge-regularized coefficients
#                   on branch indicators with aligned parent OOF offset.
#  TREE_D2_PRUNED_BETA: top-down TSP-scored tree, max depth 2,
#                   independent D/V chooses depth 1 vs depth 2;
#                   leaves smoothed toward parent nodes.
#
# Source pair nomination is supervised strictly inside training, never
# using outer-test labels. No split acceptance by 1SE vs aligned parent.
# This is an explicit *change of model family*, not posthoc relaxation
# of the earlier 19x5 study; full externally independent validation needed.

rlmt_v3_config <- function() {
  list(seed=42L,proposal_n=4096L,min_leaf_n=8L,
       beta_strength=16,logistic_ridge=1,max_depth=2L)
}

rlmt_v3_check <- function(ok,context) {
  if(!isTRUE(ok))stop("RLMT_V3_STOP: ",context,call.=FALSE)
}

rlmt_v3_beta <- function(y01,prior,strength=16) {
  rlmt_v3_check(length(y01)>0L && all(y01 %in% c(0,1)) &&
                 length(prior)==1L && is.finite(prior) &&
                 prior>0 && prior<1 && strength>0,
                 "invalid Beta-binomial leaf prior/data")
  (sum(y01)+strength*prior)/(length(y01)+strength)
}

# TSP *training-only* label-contrast, not Gini and not original
# 2011 TSPDT with k-TSP internal tests and C4.5 pessimistic pruning.
rlmt_v3_choose_pair <- function(X,y,params=rlmt_v3_config(),
                                node_seed=params$seed) {
  if(length(unique(y))!=2L ||
     nrow(X)<2L*params$min_leaf_n) return(NULL)
  yy <- rlmt_y01(y)
  positive <- which(yy==1L)
  negative <- which(yy==0L)
  if(!length(positive) || !length(negative))return(NULL)
  candidates <- rlmt_dc_random_pairs(
    names(X),params$proposal_n,seed=node_seed)
  best <- NULL
  winner <- -Inf
  tol <- 1e-12
  for(i in seq_len(nrow(candidates))) {
    a <- candidates$gene_a[[i]]
    b <- candidates$gene_b[[i]]
    z <- X[[a]]>X[[b]]
    n1 <- sum(z)
    n0 <- length(z)-n1
    if(n1<params$min_leaf_n || n0<params$min_leaf_n) next
    tsp <- abs(mean(z[positive])-mean(z[negative]))
    # Pool order breaks ties deterministically.
    if(tsp>winner+tol) {
      winner <- tsp
      best <- candidates[i,c("gene_a","gene_b"),drop=FALSE]
      best$tsp_score <- tsp
      best$pool_index <- candidates$pool_index[[i]]
    }
  }
  best
}

rlmt_v3_train_node <- function(X,y,max_depth=1L,
                               params=rlmt_v3_config(),
                               level=0L,
                               prior=NULL) {
  yy <- rlmt_y01(y)
  p0 <- if(is.null(prior)) (sum(yy)+1)/(length(yy)+2) else prior
  p <- rlmt_v3_beta(yy,p0,params$beta_strength)
  result <- list(prob=p,n=length(yy),positive=sum(yy),
                 depth=level,pair=NULL,children=NULL)
  if(level>=max_depth)return(result)
  chosen <- rlmt_v3_choose_pair(
    X,y,params,node_seed=params$seed+level*701L+17L)
  if(is.null(chosen))return(result)
  z <- rlmt_v2_route(X,chosen)
  if(sum(z)<params$min_leaf_n ||
     sum(!z)<params$min_leaf_n)
    stop("RLMT_V3_STOP: violating minimum leaf sample count")
  result$pair <- chosen
  result$children <- list(
    true=rlmt_v3_train_node(
      X[z,,drop=FALSE],y[z],max_depth,params,
      level=level+1L,prior=p),
    false=rlmt_v3_train_node(
      X[!z,,drop=FALSE],y[!z],max_depth,params,
      level=level+1L,prior=p))
  result
}

rlmt_v3_predict_node <- function(node,X,max_depth=Inf) {
  n <- nrow(X)
  if(is.null(node$pair) || node$depth>=max_depth)
    return(rep(node$prob,n))
  route <- rlmt_v2_route(X,node$pair)
  ans <- numeric(n)
  if(any(route))
    ans[route] <- rlmt_v3_predict_node(
      node$children$true,X[route,,drop=FALSE],max_depth)
  if(any(!route))
    ans[!route] <- rlmt_v3_predict_node(
      node$children$false,X[!route,,drop=FALSE],max_depth)
  ans
}

rlmt_v3_depth_count <- function(node) {
  if(is.null(node$pair))return(0L)
  as.integer(1L+
    rlmt_v3_depth_count(node$children$true)+
    rlmt_v3_depth_count(node$children$false))
}

rlmt_v3_postpruned_tree <- function(X,y,params=rlmt_v3_config()) {
  phase <- rlmt_dc_partition(y,seed=params$seed)
  D <- phase$discovery
  V <- phase$confirmation
  rlmt_v3_check(length(D)>0L && length(V)>0L &&
                length(intersect(D,V))==0L,
                "D/V patient groups overlap")
  initial <- rlmt_v3_train_node(
    X[D,,drop=FALSE],y[D],2L,params)
  v1 <- rlmt_v3_predict_node(
    initial,X[V,,drop=FALSE],1L)
  v2 <- rlmt_v3_predict_node(
    initial,X[V,,drop=FALSE],2L)
  y01 <- rlmt_y01(y[V])
  perloss <- function(p) {
    q <- pmax(1e-15,pmin(1-1e-15,p))
    -(y01*log(q)+(1-y01)*log1p(-q))
  }
  gain <- perloss(v1)-perloss(v2)
  av <- mean(gain)
  se <- if(length(gain)>1L)
    stats::sd(gain)/sqrt(length(gain)) else Inf
  if(!is.finite(se))se <- Inf
  initial_deeper <- rlmt_v3_depth_count(initial)>1L
  keep_depth <- if(initial_deeper &&
                   is.finite(av) && av>se && av>0) 2L else 1L
  # After depth choice, refit only on OUTER TRAIN (never outer test).
  fitted <- rlmt_v3_train_node(X,y,keep_depth,params)
  list(tree=fitted,chosen_depth=keep_depth,
       initial_deeper=initial_deeper,
       validation_depth2_minus_depth1_logloss=-av,
       validation_depth2_gain_se=se,
       discovery_ids=rownames(X)[D],
       confirmation_ids=rownames(X)[V])
}

rlmt_v3_fit <- function(X,y,method,params=rlmt_v3_config()) {
  rlmt_v2_check_xy(X,y)
  rlmt_v3_check(
    method %in% c("STUMP_TSP_BETA",
                  "STUMP_RIDGE_OFFSET",
                  "TREE_D2_PRUNED_BETA"),
    "unknown new-tree variant")
  global <- rlmt_v2_fit_parent(X,y)
  if(method=="STUMP_TSP_BETA") {
    tree <- rlmt_v3_train_node(X,y,1L,params)
    return(list(method=method,params=params,parent=global,
                tree=tree,model_identity="SMOOTHED_TSP_STUMP"))
  }
  if(method=="TREE_D2_PRUNED_BETA") {
    pruned <- rlmt_v3_postpruned_tree(X,y,params)
    return(list(method=method,params=params,parent=global,
                pruned=pruned,tree=pruned$tree,
                model_identity="VALIDATION_PRUNED_TSP_TREE"))
  }
  # Ridge logistic uses ONLY two branch indicators, NOT 10 selected
  # high-dimensional RLR features in each small leaf.
  chosen <- rlmt_v3_choose_pair(X,y,params,params$seed+17L)
  if(is.null(chosen))return(list(
    method=method,params=params,parent=global,
    chosen=NULL,fit=NULL,model_identity="OFFSET_STUMP_NO_ROUTE"))
  oof <- rlmt_v2_parent_oof(X,y)$prob
  z <- rlmt_v2_route(X,chosen)
  R <- cbind(route_false=as.numeric(!z),
             route_true=as.numeric(z))
  fit <- glmnet::glmnet(
    x=R,y=rlmt_y01(y),family="binomial",
    offset=rlmt_v2_logit(oof),
    lambda=params$logistic_ridge,alpha=0,
    intercept=FALSE,standardize=FALSE,
    penalty.factor=c(1,1))
  list(method=method,params=params,parent=global,
       chosen=chosen,fit=fit,model_identity="RLR_PLUS_TWO_RIDGE_LOGIT_OFFSETS",
       training_parent_oof=oof)
}

rlmt_v3_predict <- function(model,X) {
  params <- model$params
  p0 <- rlmt_v2_predict_parent(model$parent,X)
  if(model$method=="STUMP_RIDGE_OFFSET") {
    if(is.null(model$chosen))return(list(
      prob=p0,parent=p0,leaf_id=rep(1L,nrow(X))))
    route <- rlmt_v2_route(X,model$chosen)
    R <- cbind(route_false=as.numeric(!route),
               route_true=as.numeric(route))
    p <- as.numeric(stats::predict(
      model$fit,newx=R,s=params$logistic_ridge,
      newoffset=rlmt_v2_logit(p0),type="response"))
    leafid <- ifelse(route,2L,3L)
  } else {
    p <- rlmt_v3_predict_node(model$tree,X)
    if(is.null(model$tree$pair)) {
      leafid <- rep(1L,nrow(X))
    } else {
      route <- rlmt_v2_route(X,model$tree$pair)
      leafid <- ifelse(route,2L,3L)
    }
  }
  rlmt_v3_check(length(p)==nrow(X) &&
                 all(is.finite(p)) &&
                 all(p>=0 & p<=1),
                 "invalid tree probability")
  list(prob=p,parent=p0,leaf_id=leafid)
}

rlmt_v3_run <- function(X_train,y_train,X_test,y_test,
                        out_dir,method) {
  rlmt_v2_check_xy(X_train,y_train)
  rlmt_v3_check(is.data.frame(X_test) &&
       identical(names(X_train),names(X_test)) &&
       nrow(X_test)==length(y_test) &&
       !anyNA(X_test) &&
       all(vapply(X_test,function(z)is.numeric(z)&&
                  all(is.finite(z)),logical(1))) &&
       is.factor(y_test) &&
       identical(levels(y_test),levels(y_train)) &&
       !anyNA(y_test) &&
       !is.null(rownames(X_test)) &&
       !anyDuplicated(rownames(X_test)) &&
       !length(intersect(rownames(X_train),rownames(X_test))),
       "invalid heldout data schema/patient separation")
  rlmt_v3_check(!file.exists(out_dir) &&
                 !dir.exists(out_dir),
                 "no-overwrite protection")
  model <- rlmt_v3_fit(X_train,y_train,method)
  result <- rlmt_v3_predict(model,X_test)
  ensure_dir(out_dir)
  saveRDS(model,file.path(out_dir,"complete_model.rds"))
  write.csv(data.frame(
     sample_id=rownames(X_test),y_true=as.character(y_test),
     y_pred=ifelse(result$prob>=.5,levels(y_test)[2L],
                   levels(y_test)[1L]),
     prob=result$prob,parent_prob=result$parent,
     leaf_id=result$leaf_id),
     file.path(out_dir,"predictions_test.csv"),row.names=FALSE)
  m <- list(method=method,
     root_identity="RLR_ALIGNED_V2",
     n_internal_nodes=if(!is.null(model$tree))
        rlmt_v3_depth_count(model$tree) else
        as.integer(!is.null(model$chosen)),
     chosen_depth=if(!is.null(model$pruned))
        model$pruned$chosen_depth else 1L,
     logloss=rlmt_logloss(rlmt_y01(y_test),result$prob),
     brier=mean((rlmt_y01(y_test)-result$prob)^2),
     parent_logloss=rlmt_logloss(rlmt_y01(y_test),result$parent),
     scientific_stage="EXPOSED_RETROSPECTIVE_DEV_ONLY")
  rlmt_dc_write_json(m,file.path(out_dir,"metrics_test.json"))
  invisible(c(result,list(model=model,metrics=m)))
}
