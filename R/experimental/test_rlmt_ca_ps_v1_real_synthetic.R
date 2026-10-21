# Real ALGORITHM integration smoke, synthetic data only.
# Needs the EXISTING project dependencies installed: glmnet, pROC,
# jsonlite, SWAP and all original RLR dependencies.
# Unlike unit/mocked test, runs the actual original RLR parent,
# discovery SWAP feature selection, glmnet leaf CV and actual confirmation.
#
# Does not use real omics/ARFF datasets or inspect any prior held-out results.
for (pkg in c("glmnet","pROC","jsonlite","SWAP")) {
  if (!requireNamespace(pkg,quietly=TRUE))
    stop("RLMT_CA_PS_REAL_SMOKE_STOP: missing existing dependency ",pkg)
}
suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
  library(SWAP)
})
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_ca_ps_v1.R")
set.seed(8361L)
p <- 42L
columns <- sprintf("syn_gene_%03d",seq_len(p))
make_data <- function(n,prefix) {
  X <- matrix(rnorm(n*p),nrow=n,ncol=p)
  colnames(X) <- columns
  route <- X[,1L]>X[,2L]
  eta <- -.25+.85*(X[,3L]-X[,4L])+
         ifelse(route,1,-1)*.65*(X[,5L]-X[,6L])
  prob <- plogis(eta)
  label <- factor(ifelse(runif(n)<prob,"P","N"),
                  levels=c("N","P"))
  X <- as.data.frame(X,check.names=FALSE)
  rownames(X) <- sprintf("%s_%03d",prefix,seq_len(n))
  list(X=X,y=label)
}
training <- make_data(180L,"train")
testing <- make_data(64L,"holdout")
stopifnot(all(table(training$y)>=45L),
          all(table(testing$y)>0L),
          identical(colnames(training$X),colnames(testing$X)))
global_pool <- compute_top_pairs_switchbox(
  training$X,training$y,top_n=1000L,buffer_mult=5L,
  out_dir=NULL,save_all=FALSE)
stopifnot(nrow(global_pool)>=10L)
out <- tempfile(pattern="rlmt_ca_ps_real_synthetic_")
result <- rlmt_ca_ps_v1_run(
  training$X,training$y,testing$X,testing$y,
  outer_training_pool=global_pool,
  out_dir=out,
  config=list(),
  extra_metrics=list(dataset="synthetic_only_not_omics"))
stopifnot(result$status=="OK",
          result$realized_leaves %in% c(1L,2L),
          length(result$prob)==nrow(testing$X),
          all(result$prob>=0 & result$prob<=1),
          identical(as.character(result$class),
                    ifelse(result$prob>=.5,"P","N")),
          file.exists(file.path(out,"metrics_test.json")),
          file.exists(file.path(out,"rlr_reference",
                               "predictions_test.csv")),
          file.exists(file.path(out,"prediction_components_test.csv")))
components <- read.csv(
  file.path(out,"prediction_components_test.csv"),
  stringsAsFactors=FALSE)
stopifnot(isTRUE(all.equal(
  components$mixed_prob,
  components$parent_prob+
    components$shrink_weight*
    (components$child_prob-components$parent_prob),
  tolerance=1e-12)))
if (!result$accepted_split) {
  stopifnot(identical(result$prob,as.numeric(components$parent_prob)),
            all(components$shrink_weight==0),
            result$realized_leaves==1L)
} else {
  stopifnot(all(components$shrink_weight>0 &
                components$shrink_weight<1),
            file.exists(file.path(out,"leaf_models",
                                  "node_2","model_complete.rds")),
            file.exists(file.path(out,"leaf_models",
                                  "node_3","model_complete.rds")))
}
metrics <- jsonlite::read_json(file.path(out,"metrics_test.json"),
                               simplifyVector=TRUE)
stopifnot(is.finite(metrics$logloss),
          is.finite(metrics$brier),
          metrics$method_family=="RLMT-CA-PS-v1")
cat("RLMT_CA_PS_REAL_SYNTHETIC=PASS\n")
cat("ACCEPTED_SPLIT=",result$accepted_split,"\n",sep="")
cat("REALIZED_LEAVES=",result$realized_leaves,"\n",sep="")
cat("HELDOUT_SYNTH_LOGLOSS=",format(metrics$logloss,digits=8),"\n",sep="")
cat("NO_BIOMEDICAL_DATA=TRUE\n")
unlink(out,recursive=TRUE,force=TRUE)
