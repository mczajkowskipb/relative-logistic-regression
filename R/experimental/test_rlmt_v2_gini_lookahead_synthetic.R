# RLR leaf RDS tests for new (a) TSPDT_GINI and (b) RLMT_LA modes.
# Entirely synthetic, with explicitly DISCLOSED separate test-only
# nomination + independent confirmation Boolean forced paths.
for(pkg in c("switchBox","glmnet","jsonlite","pROC")) {
  if(!requireNamespace(pkg,quietly=TRUE))
    stop("RLMT_LA_SYNTH_STOP: missing package ",pkg)
}
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v2_aligned_offset.R")
source("R/experimental/rlmt_v2_gini_lookahead.R")
check <- function(x,msg) {
  if(!isTRUE(x))stop("RLMT_LA_SYNTH_STOP: ",msg,call.=FALSE)
  message("PASS ",msg)
}
fake <- function(n,seed,prefix) {
  set.seed(seed)
  p <- 42L
  m <- matrix(rnorm(n*p),nrow=n,ncol=p)
  colnames(m) <- sprintf("rel_gene_%02d",seq_len(p))
  z <- m[,1L]>m[,2L]
  eta <- .5*(m[,3L]-m[,4L])+
       ifelse(z,1,-1)*1.7*(m[,5L]-m[,6L])
  y <- factor(ifelse(runif(n)<plogis(eta),"P","N"),
              levels=c("N","P"))
  X <- as.data.frame(m,check.names=FALSE)
  rownames(X) <- sprintf("%s_%04d",prefix,seq_len(n))
  list(X=X,y=y)
}
a <- fake(180L,27361L,"train")
b <- fake(48L,27362L,"test")
conf <- rlmt_v2_conf()$leaf
dc <- rlmt_v2_conf()$dc
check(all(table(a$y)>=65L),"synthetic train class support")
check(rlmt_la_gini(factor(c("N","P")))==.5,
      "binary Gini exact")
check(rlmt_la_gini(factor(c("N","N","N"),
       levels=c("N","P")))==0,
      "pure-node Gini exactly zero")
phase <- rlmt_dc_partition(a$y,seed=conf$seed)
Xd <- a$X[phase$discovery,,drop=FALSE]
yd <- a$y[phase$discovery]
choice0 <- rlmt_la_candidates(Xd,yd,conf,dc,"TSPDT_GINI")
check(nrow(choice0$audit)==as.integer(min(dc$proposal_n,
      ncol(Xd)*(ncol(Xd)-1L)/2)),
      "deterministic Gini proposal inventory")
check(!is.null(choice0$candidate),
      "Gini viable training-only relation")
choice1 <- rlmt_la_candidates(Xd,yd,conf,dc,"RLMT_LA")
check(nrow(choice1$audit)>0 &&
      nrow(choice1$audit)<=conf$screen_top_B &&
      all(choice1$audit$pool_index %in%
          choice0$audit$pool_index),
      "lookahead at most top-20 actual RLR-child fits")
if(!is.null(choice1$candidate)) {
  check(choice1$candidate$full_r_lr_leaf_oof &&
      choice1$candidate$mean_gain>0 &&
      choice1$candidate$lower_gain>0,
      "lookahead nominates only validated discovery OOF gain")
}
for(method in c("TSPDT_GINI","RLMT_LA")) {
  path <- tempfile(pattern=paste0("rlmt_v2_",method,"_natural_"))
  fit <- rlmt_la_run(a$X,a$y,b$X,b$y,path,method)
  saved <- read.csv(file.path(path,"predictions_test.csv"),
                    stringsAsFactors=FALSE)
  check(length(fit$prob)==nrow(b$X) &&
        isTRUE(all.equal(as.numeric(saved$prob),fit$prob,
                         tolerance=1e-12)) &&
        isTRUE(all.equal(as.numeric(saved$parent_prob),fit$parent,
                         tolerance=1e-12)),
        paste0(method,"_natural_RLR_leaf_execution"))
  if(identical(fit$decision$status,"ACCEPTED")) {
    check(file.exists(file.path(path,"complete_tree.rds")),
          paste0(method,"_natural_accepted_complete_RDS"))
  } else check(identical(fit$parent,fit$prob) &&
        file.exists(file.path(path,"complete_parent.rds")),
        paste0(method,"_natural_exact_parent_fallback"))
  unlink(path,recursive=TRUE,force=TRUE)
}
# TEST ONLY: force a feasible route chosen by the training data and
# allow acceptance after actual D/V confirmation calculation.
# NO scientific gain/conclusion is inferred from this forced path.
real_select <- rlmt_la_candidates
real_confirm <- rlmt_dc_confirm
force_counts <- 0L
rlmt_la_candidates <- function(X,y,conf,dc,method) {
  candidates <- real_select(X,y,conf,dc,"TSPDT_GINI")
  candidates$candidate <- candidates$audit[
    candidates$audit$feasible, ,drop=FALSE][1L,,drop=FALSE]
  check(nrow(candidates$candidate)==1L &&
        isTRUE(candidates$candidate$feasible),
        "forced_test_candidate_comes_from_feasible_train_D")
  candidates
}
rlmt_dc_confirm <- function(y_confirm,parent_prob,split_prob,min_n=20L) {
  result <- real_confirm(y_confirm,parent_prob,split_prob,min_n)
  check(is.finite(result$gain) && result$n==length(y_confirm),
        "real_confirmation_loss_measured_before_test_force")
  force_counts <<- force_counts + 1L
  result$accepted <- TRUE
  result
}
for(method in c("TSPDT_GINI","RLMT_LA")) {
  path <- tempfile(pattern=paste0("rlmt_v2_",method,"_forced_"))
  out <- rlmt_la_run(a$X,a$y,b$X,b$y,path,method)
  check(identical(out$decision$status,"ACCEPTED"),
        paste0(method,"_two_leaf_real_training"))
  model <- readRDS(file.path(path,"complete_tree.rds"))
  replay <- rlmt_v2_predict_tree(model,b$X,conf)
  check(isTRUE(all.equal(out$prob,replay$prob,tolerance=0)) &&
        isTRUE(all.equal(out$parent,replay$parent,tolerance=0)) &&
        all(replay$route %in% c(2L,3L)),
        paste0(method,"_complete_parent_and_both_RLR_leaf_RDS_replay"))
  # Model probability predicts using only X, never y_test.
  reversed <- factor(rev(as.character(b$y)),levels=levels(b$y))
  check(!identical(reversed,b$y) &&
        identical(rlmt_v2_predict_tree(model,b$X,conf)$prob,out$prob),
        paste0(method,"_test_outcome_ignored_by_saved_predictor"))
  unlink(path,recursive=TRUE,force=TRUE)
}
check(force_counts==2L,"two_forced_real_independent_confirmation_calls")
message("RLMT_V2_GINI_LA_SYNTHETIC_REAL_GATE=PASS")
message("FORCED_NOMINATION_AND_CONFIRMATION_TEST_ONLY=TRUE")
message("NO_BIOMEDICAL_DATA=TRUE")
