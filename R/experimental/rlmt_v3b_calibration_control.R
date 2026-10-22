# RLMT-v3B: calibration-controlled compact relational stumps.
# This is a NEW exploratory model FAMILY following the already EXPOSED
# 19×5 and v3 2×3 results. Comparisons are retrospective diagnostics,
# not an independent proof of superiority.
#
# For each OUTER TRAIN fold:
#  (1) fit ONE identical RLR_ALIGNED_V2 parent and strict parent OOF;
#  (2) fit a single uniform-offset ridge intercept as CALIBRATION CONTROL;
#  (3) fit TSP top-contrast one-rule 2-coefficient ridge;
#  (4) fit Gini top-gain one-rule 2-coefficient ridge;
#  (5) fit OOF-residual-contrast one-rule 2-coefficient ridge.
# Every candidate route is x_i > x_j and is chosen without outer test.
# No child cv.glmnet; no 1SE gate forbidding stumps.
#
# Sourced after R/experimental/rlmt_v3_stumps_pruned.R.
rlmt_v3b_modes <- c("GLOBAL_RIDGE_CAL",
                    "STUMP_TSP_RIDGE",
                    "STUMP_GINI_RIDGE",
                    "STUMP_RESID_RIDGE")

rlmt_v3b_global_fit <- function(y,parent_oof,lambda=1) {
  # glmnet uses 1/n negative binomial log-likelihood + lambda/2*beta^2
  # (alpha=0); this single-parameter Newton fit avoids invalid one-column
  # glmnet designs and uses the SAME lambda=1 penalty as the route models.
  offset <- rlmt_v2_logit(parent_oof)
  yy <- rlmt_y01(y)
  rlmt_v3b_check <- rlmt_v3_check
  rlmt_v3b_check(length(offset)==length(yy) &&
                  all(is.finite(offset)) &&
                  lambda>0,"global ridge offset training input")
  beta <- 0
  for(step in seq_len(80L)) {
    p <- stats::plogis(offset+beta)
    grad <- mean(p-yy)+lambda*beta
    hess <- mean(p*(1-p))+lambda
    delta <- grad/hess
    beta <- beta-delta
    if(abs(delta)<1e-13) break
  }
  rlmt_v3b_check(is.finite(beta) &&
                  abs(mean(stats::plogis(offset+beta)-yy)+
                      lambda*beta)<1e-10,
                  "global offset penalized logistic optimizer convergence")
  list(model_identity="GLOBAL_ONE_COEFFICIENT_RIDGE_CONTROL",
       offset_intercept=beta,lambda=lambda)
}

rlmt_v3b_global_predict <- function(cal,parent_probs) {
  p <- stats::plogis(rlmt_v2_logit(parent_probs)+cal$offset_intercept)
  rlmt_v3_check(all(is.finite(p))&&all(p>=0&p<=1),
                "invalid global ridge probabilities")
  p
}

rlmt_v3b_proposals <- function(X,y,parent_oof,
                                params=rlmt_v3_config()) {
  rlmt_v3_check(length(parent_oof)==nrow(X) &&
                 all(is.finite(parent_oof)) &&
                 all(parent_oof>=0 & parent_oof<=1),
                 "parent OOF probability screening mismatch")
  yy <- rlmt_y01(y)
  res <- yy-parent_oof
  p0 <- mean(yy)
  root_gini <- 2*p0*(1-p0)
  index <- rlmt_dc_random_pairs(
    names(X),params$proposal_n,seed=params$seed+17L)
  best <- list(TSP=NULL,GINI=NULL,RESID=NULL)
  high <- c(TSP=-Inf,GINI=-Inf,RESID=-Inf)
  admissible <- 0L
  class0 <- which(yy==0L)
  class1 <- which(yy==1L)
  for(i in seq_len(nrow(index))) {
    a <- index$gene_a[[i]]
    b <- index$gene_b[[i]]
    z <- X[[a]]>X[[b]]
    n1 <- sum(z); n0 <- length(z)-n1
    if(n1<params$min_leaf_n || n0<params$min_leaf_n) next
    admissible <- admissible+1L
    p_true <- mean(yy[z])
    p_false <- mean(yy[!z])
    gini <- root_gini-
      (n1*2*p_true*(1-p_true)+n0*2*p_false*(1-p_false))/
         length(yy)
    tsp <- abs(mean(z[class1])-mean(z[class0]))
    residual <- abs(mean(res[z])-mean(res[!z]))*
      sqrt(n1*n0/length(yy))
    scores <- c(TSP=tsp,GINI=gini,RESID=residual)
    for(type in names(scores)) {
      if(scores[[type]]>high[[type]]+1e-12) {
        high[[type]] <- scores[[type]]
        best[[type]] <- data.frame(
           gene_a=a,gene_b=b,pool_index=index$pool_index[[i]],
           score=scores[[type]],n_true=n1,n_false=n0,
           train_true_pos=sum(yy[z]),
           train_false_pos=sum(yy[!z]),
           stringsAsFactors=FALSE)
      }
    }
  }
  list(best=best,feasible_pairs=admissible,
       scanned_pairs=nrow(index))
}

rlmt_v3b_two_offset_fit <- function(X,y,parent_oof,candidate,
                                     lambda=1) {
  if(is.null(candidate))return(NULL)
  route <- rlmt_v2_route(X,candidate)
  design <- cbind(route_false=as.numeric(!route),
                  route_true=as.numeric(route))
  fit <- glmnet::glmnet(
     x=design,y=rlmt_y01(y),family="binomial",
     offset=rlmt_v2_logit(parent_oof),
     lambda=lambda,alpha=0,intercept=FALSE,
     standardize=FALSE,penalty.factor=c(1,1))
  list(model_identity="TWO_ROUTE_COEFFICIENT_LOGIT_RIDGE",
       fit=fit,candidate=candidate,lambda=lambda)
}

rlmt_v3b_two_offset_predict <- function(mod,X,parent) {
  if(is.null(mod))return(parent)
  z <- rlmt_v2_route(X,mod$candidate)
  design <- cbind(route_false=as.numeric(!z),
                  route_true=as.numeric(z))
  p <- as.numeric(stats::predict(
    mod$fit,newx=design,s=mod$lambda,
    type="response",newoffset=rlmt_v2_logit(parent)))
  rlmt_v3_check(all(is.finite(p))&&all(p>=0&p<=1),
                "invalid route ridge predictions")
  p
}

rlmt_v3b_fit_bundle <- function(X,y,params=rlmt_v3_config()) {
  rlmt_v2_check_xy(X,y)
  # Strict RLR parent OOF requires >=9 per class. In extremely small
  # cohorts the only defensible PREDECLARED outcome is the exact
  # one-leaf parent for ALL four methods, never a leaked/pseudo OOF
  # vector and never dropping the dataset from the independent audit.
  if(any(table(y)<9L)) {
    parent <- rlmt_v2_fit_parent(X,y)
    return(list(parent=parent,training_parent_oof=NULL,
      proposals=list(best=list(TSP=NULL,GINI=NULL,RESID=NULL),
        feasible_pairs=0L,scanned_pairs=0L),
      calibrator=list(
        model_identity="NO_CALIBRATION_LOW_CLASS_SUPPORT",
        offset_intercept=0,lambda=params$logistic_ridge),
      TSP=NULL,GINI=NULL,RESID=NULL,params=params,
      fallback_reason="MIN_CLASS_LT_9_STRICT_OOF_NOT_FEASIBLE"))
  }
  # Reuse *exact* v3 STUMP_RIDGE_OFFSET method where available,
  # rather than silently replace the initially promising algorithm.
  base <- rlmt_v3_fit(X,y,"STUMP_RIDGE_OFFSET",params)
  # On rare datasets no route meets min_leaf_n. The original v3
  # stump then correctly falls back to its aligned parent without
  # computing OOF. Still fit the calibration control fairly.
  if(is.null(base$training_parent_oof)) {
    rlmt_v3_check(is.null(base$chosen) && is.null(base$fit),
                  "missing OOF despite a fitted route")
    base$training_parent_oof <- rlmt_v2_parent_oof(X,y)$prob
  }
  rlmt_v3_check(!is.null(base$parent) &&
                 length(base$training_parent_oof)==nrow(X),
                 "v3 TSP stump/parent-OOF prerequisite")
  proposed <- rlmt_v3b_proposals(
     X,y,base$training_parent_oof,params)
  nominated <- proposed$best$TSP
  if(!is.null(base$chosen) && !is.null(nominated)) {
    rlmt_v3_check(
     identical(as.character(base$chosen$gene_a[[1L]]),
               as.character(nominated$gene_a[[1L]])) &&
     identical(as.character(base$chosen$gene_b[[1L]]),
               as.character(nominated$gene_b[[1L]])),
     "TSP branch differs from frozen first v3 method")
  }
  # TSP model is exactly the v3 fitted glmnet+route:
  tsp <- if(is.null(base$fit)) NULL else list(
       fit=base$fit,candidate=base$chosen,
       lambda=params$logistic_ridge,
       model_identity="EXACT_V3_TSP_RIDGE_STUMP")
  list(parent=base$parent,training_parent_oof=base$training_parent_oof,
       proposals=proposed,
       calibrator=rlmt_v3b_global_fit(
            y,base$training_parent_oof,params$logistic_ridge),
       TSP=tsp,
       GINI=rlmt_v3b_two_offset_fit(
            X,y,base$training_parent_oof,
            proposed$best$GINI,params$logistic_ridge),
       RESID=rlmt_v3b_two_offset_fit(
            X,y,base$training_parent_oof,
            proposed$best$RESID,params$logistic_ridge),
       params=params)
}

rlmt_v3b_predict_bundle <- function(fit,X) {
  p0 <- rlmt_v2_predict_parent(fit$parent,X)
  predicted <- list(
    GLOBAL_RIDGE_CAL=rlmt_v3b_global_predict(fit$calibrator,p0),
    STUMP_TSP_RIDGE=rlmt_v3b_two_offset_predict(fit$TSP,X,p0),
    STUMP_GINI_RIDGE=rlmt_v3b_two_offset_predict(fit$GINI,X,p0),
    STUMP_RESID_RIDGE=rlmt_v3b_two_offset_predict(fit$RESID,X,p0))
  rlmt_v3_check(all(vapply(predicted,function(p)
           length(p)==nrow(X)&&all(is.finite(p))&&
           all(p>=0&p<=1),logical(1))),
           "bundle incomplete predictions")
  list(parent=p0,prob=predicted)
}
