# RLMT-v5: actual REAL synthetic original RLR/switchBox/glmnet training.
for(pkg in c("switchBox","glmnet","jsonlite","pROC"))
  if(!requireNamespace(pkg,quietly=TRUE))
    stop("RLMT5_SYNTH_STOP: missing ",pkg)
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v4_rank_repair.R")
source("R/experimental/rlmt_v5_honest_router.R")
check <- function(x,msg){
 if(!isTRUE(x))stop("RLMT5_SYNTH_STOP: ",msg,call.=FALSE)
 message("PASS ",msg)
}
sim <- function(n,seed,mark) {
 set.seed(seed)
 p <- 42L
 m <- matrix(stats::rnorm(n*p),ncol=p)
 colnames(m) <- sprintf("gene_%02d",seq_len(p))
 prob <- stats::plogis(1.1*(m[,1L]>m[,2L])-
   0.7*(m[,3L]>m[,4L])+0.4*m[,5L]-0.5*m[,6L])
 y <- factor(ifelse(stats::runif(n)<prob,"P","N"),
   levels=c("N","P"))
 X <- as.data.frame(m,check.names=FALSE)
 rownames(X) <- sprintf("%s_%03d",mark,seq_len(n))
 list(X=X,y=y)
}
tr <- sim(150L,7201L,"train")
te <- sim(50L,7202L,"test")
check(min(table(tr$y))>=40L,"real train adequate class support")
model <- rlmt5_fit(tr$X,tr$y)
pred <- rlmt5_predict(model,te$X)
check(inherits(model,"rlmt5_model") &&
      length(pred$prob)==nrow(te$X) &&
      all(pred$prob>0 & pred$prob<1),
      "real v5 model and heldout probabilities")
check(identical(model$full_bundle$parent$model_identity,
          "FROZEN_ORIGINAL_RLR_K10_M10_D2"),
      "root is manuscript-original RLR, not aligned surrogate")
check(model$router$selected %in%
      c("ORIGINAL_RLR",rlmt5_candidates) &&
      model$router$weight>=0 &&
      model$router$weight<=1,
      "honest v5 calibrated router has legal mode/weight")
check(length(model$router$inner_partitions)==3L &&
      all(vapply(model$router$inner_partitions,function(f)
       length(intersect(f$train,f$validation))==0L,logical(1))),
       "inner validation samples disjoint from fitted inner models")
check(nrow(model$router$evidence)==12L,
      "fixed 3 mode x 4 weight CV evidence inventory")
check(all(is.finite(model$router$evidence$conservative_gain)),
      "finite honest validation evidence retained")
check(all(abs(pred$prob-
      ((1-pred$weight)*pred$original_prob+
       pred$weight*pred$candidates[[pred$model]]))<1e-14),
      "saved v5 prediction exactly equals declared mixture")

# Explicit negative case: all candidate correctors are worse
# than unchanged original RLR -> the router must REJECT ALL.
yy <- factor(rep(c("N","P"),each=36L),levels=c("N","P"))
root <- ifelse(yy=="P",.90,.10)
wrong <- 1-root
model_score <- list(ORIGINAL_RLR=root,
 GLOBAL_CAL=wrong,CONF_ROUTE=wrong,ACTIVITY_ROUTE=wrong)
reject <- rlmt5_decide(yy,rep(1:3,24L),model_score)
check(reject$selected=="ORIGINAL_RLR" &&
      reject$weight==0 && !reject$accepted,
      "exact RLR fallback for uniformly harmful corrections")

# Explicit positive case: internally heldout corrector is
# demonstrably better than a weak but valid original RLR.
base <- ifelse(yy=="P",.56,.44)
good <- ifelse(yy=="P",.80,.20)
score <- list(ORIGINAL_RLR=base,
 GLOBAL_CAL=base,CONF_ROUTE=good,ACTIVITY_ROUTE=base)
accept <- rlmt5_decide(yy,rep(1:3,24L),score)
check(accept$selected=="CONF_ROUTE" &&
      accept$weight>0 && accept$accepted,
      "genuine conservative train-only gate accepts strong gain")
f <- tempfile(fileext=".rds")
saveRDS(model,f)
pred_saved <- rlmt5_predict(readRDS(f),te$X)
unlink(f)
check(identical(pred$prob,pred_saved$prob) &&
      identical(pred$model,pred_saved$model) &&
      identical(pred$weight,pred_saved$weight),
      "full model bit-exact RDS reload")
message("RLMT5_SYNTH_REAL_GUARDED_CV=PASS")
message("BIOMEDICAL_DATA_USED=FALSE")
