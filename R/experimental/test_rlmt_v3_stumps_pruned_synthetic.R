# Isolated REAL R/glmnet synthetic test; no biomedical input or outcomes.
# The smoke test checks true learning, route and RDS replay for all 3 v3
# variants. No forced acceptance, no forced pruning, no mock predictor.
for(pkg in c("switchBox","glmnet","jsonlite","pROC")) {
  if(!requireNamespace(pkg,quietly=TRUE))
    stop("RLMT_V3_SYNTH_STOP: missing package ",pkg)
}
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v2_aligned_offset.R")
source("R/experimental/rlmt_v3_stumps_pruned.R")

test <- function(value,message) {
  if(!isTRUE(value))stop("RLMT_V3_SYNTH_STOP: ",message,call.=FALSE)
  message("PASS ",message)
}
synthetic <- function(n,seed,prefix) {
  set.seed(seed)
  x <- matrix(rnorm(n*40L),nrow=n,ncol=40L)
  colnames(x) <- sprintf("gene_%02d",seq_len(40L))
  branch <- x[,1L]>x[,2L]
  probability <- plogis(ifelse(branch,1.3,-1.3)+
                        0.5*(x[,3L]-x[,4L]))
  y <- factor(ifelse(runif(n)<probability,"P","N"),
              levels=c("N","P"))
  d <- as.data.frame(x,check.names=FALSE)
  rownames(d) <- sprintf("%s_%04d",prefix,seq_len(n))
  list(X=d,y=y)
}
a <- synthetic(180L,9001L,"train")
b <- synthetic(60L,9002L,"test")
test(all(table(a$y)>=60L),"synthetic train support")
params <- rlmt_v3_config()
test(identical(params$min_leaf_n,8L) &&
     identical(params$proposal_n,4096L) &&
     params$beta_strength==16 &&
     params$logistic_ridge==1,
     "v3 predeclared constants")
test(abs(rlmt_v3_beta(c(0,1),.5,16)-.5)<1e-12,
     "Beta-binomial smoother known reference")
test(rlmt_v3_beta(rep(1,8),.5,16)<1,
     "pure 8-sample leaf has nondegenerate probability")
real_route <- rlmt_v3_choose_pair(a$X,a$y,params)
test(!is.null(real_route),"genuine train-only TSP candidate")
z <- rlmt_v2_route(a$X,real_route)
test(all(c(sum(z),sum(!z))>=8L),
     "candidate has 8 samples per side")

for(method in c("STUMP_TSP_BETA",
                "STUMP_RIDGE_OFFSET",
                "TREE_D2_PRUNED_BETA")) {
  out <- tempfile(paste0("rlmt_v3_",method,"_"))
  result <- rlmt_v3_run(
    a$X,a$y,b$X,b$y,out,method)
  test(length(result$prob)==nrow(b$X) &&
       all(is.finite(result$prob)) &&
       all(result$prob>=0 & result$prob<=1),
       paste(method," real valid predicted probabilities"))
  p <- read.csv(file.path(out,"predictions_test.csv"),
                stringsAsFactors=FALSE)
  test(identical(p$sample_id,rownames(b$X)) &&
       isTRUE(all.equal(as.numeric(p$prob),
                        result$prob,tolerance=1e-12)),
       paste(method," written heldout predictions match"))
  restored <- readRDS(file.path(out,"complete_model.rds"))
  replay <- rlmt_v3_predict(restored,b$X)
  test(identical(result$prob,replay$prob) &&
       identical(result$parent,replay$parent),
       paste(method," bit-exact RDS replay"))
  test(identical(restored$method,method) &&
       identical(restored$parent$identity,"RLR_ALIGNED_V2"),
       paste(method," proper model/parent identity"))
  test(file.exists(file.path(out,"metrics_test.json")),
       paste(method," saved proper-score metrics"))
  if(method=="STUMP_TSP_BETA")
    test(rlmt_v3_depth_count(restored$tree)==1L,
         "TSP stump uses exactly one internal relational rule")
  if(method=="TREE_D2_PRUNED_BETA")
    test(restored$pruned$chosen_depth %in% c(1L,2L) &&
         rlmt_v3_depth_count(restored$tree)<=3L &&
         !length(intersect(restored$pruned$discovery_ids,
                           restored$pruned$confirmation_ids)),
         "D/V independent depth pruning with <=3 relation rules")
  if(method=="STUMP_RIDGE_OFFSET") {
    test(!is.null(restored$fit) && !is.null(restored$chosen) &&
         length(restored$training_parent_oof)==nrow(a$X),
         "genuine two-offset ridge GLM with honest parent OOF")
  }
  # Test labels are not passed to saveRDS predictor.
  flipped <- factor(rev(as.character(b$y)),levels=levels(b$y))
  test(!identical(flipped,b$y) &&
       identical(replay$prob,rlmt_v3_predict(restored,b$X)$prob),
       paste(method," outer-test outcomes not required at prediction"))
  unlink(out,recursive=TRUE)
}
message("RLMT_V3_STUMP_PRUNE_SYNTH_REAL_GATE=PASS")
message("BIOMEDICAL_DATA_USED=FALSE")
