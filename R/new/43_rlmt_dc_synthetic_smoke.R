# Synthetic-only full RLMT-DC-v0 fit/predict/persist smoke.
# Does not read CompCancer, Singh, Gordon, publication results, or external cohorts.
# Run: Rscript --vanilla R/new/43_rlmt_dc_synthetic_smoke.R
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")

verify <- function(ok,reason) if(!isTRUE(ok))
  stop("RLMT_DC_SMOKE_FAIL: ",reason,call.=FALSE)

set.seed(739L)
n <- 210L
p <- 48L
X <- as.data.frame(matrix(rnorm(n*p),nrow=n,ncol=p))
colnames(X) <- sprintf("synthetic_gene_%03d",seq_len(p))
rownames(X) <- sprintf("synthetic_sample_%03d",seq_len(n))
region <- rep(c(FALSE,TRUE),length.out=n)
region <- sample(region,n,replace=FALSE)
X[[1L]] <- ifelse(region,2,-2)+rnorm(n,sd=0.10)
X[[2L]] <- rnorm(n,sd=0.10)
# Distinct class rules in two regions: deterministic seed, not a claim of
# superiority. The response of each region depends on different pairs.
signal <- ifelse(region,
                 2.5*(X[[3L]]-X[[4L]]),
                 2.5*(X[[5L]]-X[[6L]]))
y <- factor(ifelse(stats::runif(n)<plogis(signal),"PR","N"),
            levels=c("N","PR"))
verify(all(table(y)>35L),"synthetic binary support")
fold <- make_repeated_stratified_folds(y,5L,1L,seed=42L)[[1L]]
Xtr <- X[fold$train_idx,,drop=FALSE]
ytr <- y[fold$train_idx]
Xte <- X[fold$test_idx,,drop=FALSE]
yte <- y[fold$test_idx]
verify(!length(intersect(rownames(Xtr),rownames(Xte))),
       "outer split contamination")

base <- tempfile("rlmt_dc_synthetic_")
dir.create(base)
pool <- compute_top_pairs_switchbox(Xtr,ytr,top_n=1000L,
                                     buffer_mult=5L,
                                     out_dir=file.path(base,"outer_pool"),
                                     save_all=FALSE)
verify(nrow(pool)>=10L,"synthetic outer RLR pool")

# Published RLR comparison, run on the EXACT same synthetic outer split.
rlr_dir <- file.path(base,"reference_rlr")
dir.create(rlr_dir)
fit_rlr_train_test(
  X_train=Xtr,y_train=ytr,X_test=Xte,y_test=yte,
  pairs_pool=pool,K=10L,m=10L,d=2L,out_dir=rlr_dir,
  alpha=0.5,lambda_choice="lambda.min",seed=42L)

dc_dir <- file.path(base,"rlmt_dc")
result <- rlmt_dc_run(Xtr,ytr,Xte,yte,pool,
                      config=list(),out_dir=dc_dir,
                      extra_metrics=list(dataset_id="SYNTHETIC_ONLY",
                                         split_id=fold$split_id))
verify(file.exists(file.path(dc_dir,"dc","decision.json")),
       "missing pretest selection decision artifact")
dec <- jsonlite::read_json(file.path(dc_dir,"dc","decision.json"),
                           simplifyVector=TRUE)
tree <- jsonlite::read_json(file.path(dc_dir,"tree.json"),
                            simplifyVector=TRUE)
m <- jsonlite::read_json(file.path(dc_dir,"metrics_test.json"),
                         simplifyVector=TRUE)
verify(isTRUE(tree$accepted_split)==isTRUE(result$accepted) &&
       tree$size %in% c(1L,2L),"tree/result mismatch")
verify(identical(m$method_family,"RLMT-DC-v0") &&
       m$realized_leaves==tree$size,
       "wrong or missing separate method identifier")

readp <- function(d) read.csv(file.path(d,"predictions_test.csv"),
                              stringsAsFactors=FALSE,
                              colClasses=c(sample_id="character",
                                           y_true="character",
                                           y_pred="character"))
pred <- readp(dc_dir)
ref <- readp(rlr_dir)

# Hard invariant: max_leaves=1 retains EXACT unchanged published RLR,
# regardless of whether the new two-leaf method nominates a split.
l1_dir <- file.path(base,"reference_rlmt_L1")
rlmt_root_output(Xtr,ytr,Xte,yte,pool,
                 rlmt_defaults(list()),l1_dir,
                 "RLMT_K10_m10_d2_pool1000_L1",
                 list(dataset_id="SYNTHETIC_ONLY"),
                 "max_leaves=1")
p_l1 <- readp(l1_dir)
c_ref <- read.csv(file.path(rlr_dir,"coefficients.csv"),
                  stringsAsFactors=FALSE)
c_l1 <- read.csv(file.path(l1_dir,"coefficients.csv"),
                 stringsAsFactors=FALSE)
verify(identical(ref$sample_id,p_l1$sample_id) &&
       identical(ref$y_pred,p_l1$y_pred) &&
       isTRUE(all.equal(ref$prob,p_l1$prob,tolerance=1e-12)) &&
       identical(c_ref,c_l1),
       "frozen RLR / max_leaves=1 parity")
cat("RLMT_DC_SYNTHETIC_L1_PARITY=PASS\n")

verify(identical(pred$sample_id,rownames(Xte)) &&
       identical(pred$y_true,as.character(yte)) &&
       all(is.finite(pred$prob)) && all(pred$prob>=0 & pred$prob<=1),
       "wrong synthetic test identities or probabilities")
y01 <- as.numeric(yte=="PR")
calc <- rlmt_logloss(y01,pred$prob)
verify(abs(calc-m$logloss)<5.1e-5,"logged synthetic test logloss")

if (!result$accepted) {
  verify(identical(ref$y_pred,pred$y_pred) &&
         isTRUE(all.equal(ref$prob,pred$prob,tolerance=1e-12)),
         "no-split result must reproduce frozen RLR")
  verify(identical(as.integer(tree$size),1L),
         "rejected split must be a single leaf")
  cat("RLMT_DC_SYNTHETIC_ROOT_PARITY=PASS\n")
} else {
  verify(tree$size==2L &&
         isTRUE(dec$confirmation$accepted) &&
         is.finite(dec$confirmation$lower_gain) &&
         dec$confirmation$lower_gain>0,
         "accepted split must pass single heldout 1SE criterion")
  route <- read.csv(file.path(dc_dir,"test_routes.csv"),
                    stringsAsFactors=FALSE)
  a <- as.character(tree$split$gene_a)
  b <- as.character(tree$split$gene_b)
  verify(identical(as.integer(route$leaf_id),
                   as.integer(ifelse(Xte[[a]]>Xte[[b]],2L,3L))),
         "test routing mismatch")
  for (id in c(2L,3L)) {
    idx <- which(route$leaf_id==id)
    leafdir <- file.path(dc_dir,"leaf_models",paste0("node_",id))
    co <- read.csv(file.path(leafdir,"coefficients.csv"),
                   stringsAsFactors=FALSE)
    if (!length(idx)) next
    parts <- strsplit(co$ordered_pair," > ",fixed=TRUE)
    verify(all(lengths(parts)==2L),"unparseable saved leaf pair")
    pairs <- data.frame(
      gene_left=vapply(parts,function(x)x[1L],character(1)),
      gene_right=vapply(parts,function(x)x[2L],character(1)),
      ordered_pair=co$ordered_pair,stringsAsFactors=FALSE)
    R <- encode_pair_matrix(Xte[idx,,drop=FALSE],pairs,m=10L)
    fit <- readRDS(file.path(leafdir,"cvfit.rds"))
    replay <- as.numeric(predict(fit,newx=R,
                                 s="lambda.min",type="response"))
    verify(isTRUE(all.equal(replay,pred$prob[idx],
                            tolerance=1e-12)),
           paste("serialized leaf prediction replay",id))
  }
  cat("RLMT_DC_SYNTHETIC_ACCEPTED_REPLAY=PASS\n")
}
cat("RLMT_DC_SYNTHETIC_PIPELINE=PASS\n")
cat("ACCEPTED_SPLIT=",isTRUE(result$accepted),"\n",sep="")
cat("SYNTHETIC_WORKDIR=",base,"\n",sep="")
