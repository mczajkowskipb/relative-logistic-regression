# All four v3B models are executed for real on synthetic train/test.
# No mocks, no forced acceptances, no test-data dependent selector.
for(pkg in c("switchBox","glmnet","jsonlite","pROC"))
  if(!requireNamespace(pkg,quietly=TRUE))
    stop("V3B_SYNTH_STOP: required package absent ",pkg)
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v2_aligned_offset.R")
source("R/experimental/rlmt_v3_stumps_pruned.R")
source("R/experimental/rlmt_v3b_calibration_control.R")
check <- function(ok,msg) {
 if(!isTRUE(ok))stop("V3B_SYNTH_STOP: ",msg)
 message("PASS ",msg)
}
fake <- function(n,seed,prefix) {
  set.seed(seed)
  m <- matrix(rnorm(n*40L),nrow=n,ncol=40L)
  colnames(m) <- sprintf("gene_%02d",seq_len(40L))
  route <- m[,1L]>m[,2L]
  response <- stats::plogis(
     ifelse(route,1.3,-1.3)+.7*(m[,3L]-m[,4L]))
  y <- factor(ifelse(runif(n)<response,"P","N"),
              levels=c("N","P"))
  X <- as.data.frame(m,check.names=FALSE)
  rownames(X) <- sprintf("%s_%04d",prefix,seq_len(n))
  list(X=X,y=y)
}
train <- fake(180L,3031L,"train")
test <- fake(60L,3032L,"test")
check(all(table(train$y)>=60L),"both class support")
model <- rlmt_v3b_fit_bundle(train$X,train$y)
p <- rlmt_v3b_predict_bundle(model,test$X)
check(identical(names(p$prob),rlmt_v3b_modes),
      "predeclared four methods returned")
check(all(vapply(p$prob,function(z)
    length(z)==nrow(test$X) && all(z>0&z<1),logical(1))),
    "four methods have valid probabilities")
check(!is.null(model$TSP) &&
      !is.null(model$GINI) &&
      !is.null(model$RESID),
    "real TSP/Gini/residual splits built")
check(identical(model$TSP$candidate$gene_a[[1L]],
                model$proposals$best$TSP$gene_a[[1L]]) &&
      identical(model$TSP$candidate$gene_b[[1L]],
                model$proposals$best$TSP$gene_b[[1L]]),
      "exact original v3 TSP route is preserved")
check(model$proposals$feasible_pairs>0L &&
      model$proposals$scanned_pairs==
         min(4096L,ncol(train$X)*(ncol(train$X)-1L)/2),
      "train-only 4096 relation candidate inventory")
one <- rlmt_v3b_global_fit(train$y,
    model$training_parent_oof,lambda=1)
check(identical(one$offset_intercept,
                model$calibrator$offset_intercept),
      "one-coefficient ridge calibration deterministic")
check(all(abs(rlmt_v3b_global_predict(one,p$parent)-
              p$prob$GLOBAL_RIDGE_CAL)<1e-14),
      "calibrated p test replay")
old <- rlmt_v3_predict(
   list(method="STUMP_RIDGE_OFFSET",
        params=model$params,
        parent=model$parent,
        fit=model$TSP$fit,
        chosen=model$TSP$candidate),test$X)
check(isTRUE(all.equal(
       p$prob$STUMP_TSP_RIDGE,old$prob,
       tolerance=0)),
      "STUMP_TSP_RIDGE is EXACT v3 method prediction")
f <- tempfile(fileext=".rds")
saveRDS(model,f)
again <- rlmt_v3b_predict_bundle(readRDS(f),test$X)
unlink(f)
check(identical(p,again),"all four predictions bit-exact after RDS replay")
check(length(intersect(rownames(train$X),rownames(test$X)))==0,
      "outer train/test sample independence")
message("RLMT_V3B_SYNTHETIC_REAL_FITS=PASS")
message("NO_BIOMEDICAL_DATA=TRUE")
