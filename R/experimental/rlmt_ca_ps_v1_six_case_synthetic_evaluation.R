# Fixed, non-adaptive, REAL RLR/CA-PS-v1 comparison on SYNTHETIC DATA ONLY.
#
# 3 fixed data-generating regimes x 2 fixed independent train/test seeds.
# No tuning, no early stopping, no cherry-picked success. All six cases
# must pass. Outer holdout labels are generated only AFTER training data;
# they never enter split selection, model fitting or shrink weights.
# The old 19 publication cohorts are NEVER loaded by this script.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L)
  stop("SYNTH_SUITE_STOP: provide a distinct new output directory")
OUT <- normalizePath(args[[1L]],mustWork=FALSE)
if (dir.exists(OUT) || file.exists(OUT))
  stop("SYNTH_SUITE_STOP: refusing overwrite of synthetic suite root")
for (pkg in c("glmnet","pROC","jsonlite","switchBox","foreign")) {
  if (!requireNamespace(pkg,quietly=TRUE))
    stop("SYNTH_SUITE_STOP: missing dependency ",pkg)
}
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_ca_ps_v1.R")

SCENARIOS <- list(
  balanced_heterogeneous=list(intercept=0,heterogeneity=1.0,
                              group_intercept=.5),
  minority_heterogeneous=list(intercept=-.95,heterogeneity=1.0,
                              group_intercept=.5),
  balanced_null_tree=list(intercept=0,heterogeneity=0,
                          group_intercept=0)
)
SEEDS <- c(1139L,2267L)
P <- 42L
NTRAIN <- 240L
NTEST <- 96L

make_cohort <- function(n,seed,scenario,id) {
  set.seed(seed)
  A <- matrix(rnorm(n*P),nrow=n,ncol=P)
  colnames(A) <- sprintf("syn_%03d",seq_len(P))
  z <- A[,1L]>A[,2L]
  shared <- .7*(A[,3L]-A[,4L])
  contrast <- ifelse(z,1,-1)*(A[,5L]-A[,6L])
  theta <- scenario$intercept + shared +
    scenario$heterogeneity*contrast +
    ifelse(z,scenario$group_intercept,-scenario$group_intercept)
  yy <- factor(ifelse(runif(n)<plogis(theta),"P","N"),
               levels=c("N","P"))
  XX <- as.data.frame(A,check.names=FALSE)
  rownames(XX) <- sprintf("%s_%03d",id,seq_len(n))
  list(X=XX,y=yy)
}

metric <- function(y,p) {
  yy <- as.numeric(y==levels(y)[2L])
  pred <- ifelse(p>=.5,levels(y)[2L],levels(y)[1L])
  bal <- mean(vapply(levels(y),function(cls)
    mean(pred[y==cls]==cls),numeric(1)))
  c(accuracy=mean(pred==as.character(y)),
    balanced_accuracy=bal,
    logloss=rlmt_logloss(yy,p),
    brier=mean((p-yy)^2))
}

dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
jsonlite::write_json(list(
  method="RLMT_CA_PS_V1_K10_m10_d2_pool1000_L2",
  scientific_status="SYNTHETIC_PROTOTYPE_NOT_EXTERNAL_VALIDATION",
  dataset_type="INDEPENDENT_SIMULATED_GAUSSIAN_GENES",
  scenarios=SCENARIOS,seeds=SEEDS,n_train=NTRAIN,
  n_test=NTEST,n_features=P,
  no_human_omics_or_patient_data=TRUE,
  no_adaptive_scenario_selection=TRUE,
  primary_metric="heldout_unweighted_logloss",
  secondary_metrics=c("accuracy","balanced_accuracy","brier"),
  comparator="unmodified original RLR parent from the same T-only run"
),file.path(OUT,"suite_protocol.json"),pretty=TRUE,auto_unbox=TRUE,digits=12)
manifest <- vector("list",length(SCENARIOS)*length(SEEDS))
j <- 0L
for (name in names(SCENARIOS)) for (seed in SEEDS) {
  j <- j+1L
  scenario <- SCENARIOS[[name]]
  tr <- make_cohort(NTRAIN,seed,scenario,"train")
  te <- make_cohort(NTEST,seed+100000L,scenario,"test")
  if (any(table(tr$y)<18L) || any(table(te$y)<8L))
    stop("SYNTH_SUITE_STOP: finite fixed simulation has too few classes: ",
         name,":",seed,"; do not silently resample")
  pool <- compute_top_pairs_switchbox(
    tr$X,tr$y,top_n=1000L,buffer_mult=5L,
    out_dir=NULL,save_all=FALSE)
  if (nrow(pool)<10L)
    stop("SYNTH_SUITE_STOP: too few training-only relations")
  model_dir <- file.path(OUT,"models",name,paste0("seed_",seed))
  fit <- rlmt_ca_ps_v1_run(
    tr$X,tr$y,te$X,te$y,outer_training_pool=pool,
    out_dir=model_dir,extra_metrics=list(
      scenario=name,synthetic_seed=seed,
      dataset_type="synthetic_gaussian"))
  original <- read.csv(file.path(model_dir,"rlr_reference",
                                 "predictions_test.csv"),
                       stringsAsFactors=FALSE)
  p1 <- as.numeric(fit$prob)
  p0 <- as.numeric(original$prob)
  if (!identical(as.character(original$sample_id),rownames(te$X)) ||
      !identical(as.character(original$y_true),as.character(te$y)) ||
      length(p1)!=NTEST || length(p0)!=NTEST ||
      anyNA(p0) || anyNA(p1) ||
      any(!is.finite(p0)) || any(!is.finite(p1)) ||
      any(p0<0 | p0>1) || any(p1<0 | p1>1))
    stop("SYNTH_SUITE_STOP: paired heldout mismatch ",name,":",seed)
  m0 <- metric(te$y,p0)
  m1 <- metric(te$y,p1)
  components <- read.csv(file.path(model_dir,
                                   "prediction_components_test.csv"),
                         stringsAsFactors=FALSE)
  if (nrow(components)!=NTEST ||
      !isTRUE(all.equal(as.numeric(components$parent_prob),p0,
                        tolerance=1e-12)) ||
      !isTRUE(all.equal(as.numeric(components$mixed_prob),p1,
                        tolerance=1e-12)) ||
      !isTRUE(all.equal(as.numeric(components$mixed_prob),
        as.numeric(components$parent_prob)+
          as.numeric(components$shrink_weight)*
          (as.numeric(components$child_prob)-
           as.numeric(components$parent_prob)),
        tolerance=1e-12)))
    stop("SYNTH_SUITE_STOP: serialized mixture mismatch ",name,":",seed)
  if (!fit$accepted_split &&
      (!identical(p1,p0) || any(components$shrink_weight!=0)))
    stop("SYNTH_SUITE_STOP: unsplit case lost bit exact RLR parity")
  d <- m1-m0
  manifest[[j]] <- data.frame(
    scenario=name,seed=seed,
    n_train=NTRAIN,n_test=NTEST,
    positive_train=sum(tr$y=="P"),
    positive_test=sum(te$y=="P"),
    accepted=fit$accepted_split,
    realized_leaves=fit$realized_leaves,
    rlr_accuracy=unname(m0["accuracy"]),
    ca_accuracy=unname(m1["accuracy"]),
    delta_accuracy=unname(d["accuracy"]),
    rlr_balanced_accuracy=unname(m0["balanced_accuracy"]),
    ca_balanced_accuracy=unname(m1["balanced_accuracy"]),
    delta_balanced_accuracy=unname(d["balanced_accuracy"]),
    rlr_logloss=unname(m0["logloss"]),
    ca_logloss=unname(m1["logloss"]),
    delta_logloss=unname(d["logloss"]),
    rlr_brier=unname(m0["brier"]),
    ca_brier=unname(m1["brier"]),
    delta_brier=unname(d["brier"]),
    stringsAsFactors=FALSE)
  write.csv(do.call(rbind,manifest[seq_len(j)]),
    file.path(OUT,"incremental_manifest.csv"),row.names=FALSE)
  cat(sprintf(
    "SYNTH_CASE=%s seed=%d accepted=%s dLogloss=%+.6f dBrier=%+.6f\n",
    name,seed,fit$accepted_split,d[["logloss"]],d[["brier"]]))
}
all <- do.call(rbind,manifest)
stopifnot(nrow(all)==6L,
          length(unique(paste(all$scenario,all$seed)))==6L,
          all(is.finite(all$delta_logloss)),
          all(is.finite(all$delta_brier)))
write.csv(all,file.path(OUT,"paired_six_scenarios.csv"),row.names=FALSE)
agg <- do.call(rbind,lapply(split(all,all$scenario),function(x)
  data.frame(scenario=x$scenario[1L],n=nrow(x),
             accepted=sum(x$accepted),
             mean_delta_logloss=mean(x$delta_logloss),
             mean_delta_brier=mean(x$delta_brier),
             mean_delta_balanced_accuracy=mean(x$delta_balanced_accuracy),
             stringsAsFactors=FALSE)))
write.csv(agg,file.path(OUT,"scenario_means.csv"),row.names=FALSE)
cat("=== ALL FIXED SYNTHETIC SCENARIOS ===\n")
print(agg,row.names=FALSE,digits=7)
cat("RLMT_CA_PS_SIX_CASE_SYNTHETIC_EVALUATION=PASS\n")
cat("SYNTH_ACCEPTED=",sum(all$accepted),"/6\n",sep="")
cat("SYNTH_MEAN_DELTA_LOGLOSS=",
    sprintf("%+.8f",mean(all$delta_logloss)),"\n",sep="")
cat("NO_BIOMEDICAL_DATA=TRUE\n")
