# Actual end-to-end synthetic gate, NO biomedical inputs.
# Also proves parent matches the *literal published RLR wrapper*,
# not the distinct aligned-v2 parent used in earlier experiments.
for(pkg in c("switchBox","glmnet","jsonlite","pROC")) {
  if(!requireNamespace(pkg,quietly=TRUE))
    stop("RLMT4_SYNTH_STOP: missing ",pkg)
}
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v4_rank_repair.R")
test <- function(ok,msg) {
 if(!isTRUE(ok))stop("RLMT4_SYNTH_STOP: ",msg,call.=FALSE)
 message("PASS ",msg)
}
fake <- function(n,seed,name) {
 set.seed(seed)
 p <- 44L
 mat <- matrix(rnorm(n*p),ncol=p)
 colnames(mat) <- sprintf("gene_%02d",seq_len(p))
 zz <- mat[,1L]>mat[,2L]
 eta <- .6*(mat[,3L]-mat[,4L]) +
   ifelse(zz,1,-1)*1.3*(mat[,5L]-mat[,6L])
 y <- factor(ifelse(runif(n)<stats::plogis(eta),"P","N"),
             levels=c("N","P"))
 X <- as.data.frame(mat,check.names=FALSE)
 rownames(X) <- sprintf("%s_%04d",name,seq_len(n))
 list(X=X,y=y)
}
a <- fake(160L,701L,"train")
b <- fake(54L,702L,"test")
test(min(table(a$y))>=55L,"synthetic class support")
fit <- rlmt4_fit(a$X,a$y)
p <- rlmt4_predict(fit,b$X)
test(identical(names(p),rlmt4_modes) &&
     all(vapply(p,function(q)length(q)==nrow(b$X) &&
       all(q>0&q<1),logical(1))),
     "all five real model predictions valid")
test(identical(fit$parent$model_identity,
    "FROZEN_ORIGINAL_RLR_K10_M10_D2"),
    "parent original RLR identity")
test(length(fit$oof_prob_training)==nrow(a$X) &&
     length(fit$oof_partition)==3L,
     "strict original-parent 3-fold OOF computed")
test(all(vapply(fit$oof_partition,function(f)
      length(intersect(f$train,f$heldout))==0L,logical(1))),
      "no OOF fold patient overlap")
test(all(vapply(fit$routes,function(r)
      !is.null(r$pair) &&
      r$pair$train_true>=8L &&
      r$pair$train_false>=8L,logical(1))),
      "all three genuine training-only relation gates selected")
test(all(vapply(fit$routes,function(r)
      is.finite(r$coefficients$gamma),logical(1))),
      "penalized logistic gamma fitted")

# Compare exact reproducibility of original supervised pool shortlist,
# cv.glmnet fitting and original wrapper probabilities.
pool <- compute_top_pairs_switchbox(
  a$X,a$y,top_n=1000L,buffer_mult=5L,
  out_dir=NULL,save_all=FALSE)
out <- tempfile("rlmt4_original_parity_")
fit_rlr_train_test(a$X,a$y,b$X,b$y,
  pairs_pool=pool,K=10L,m=10L,d=2L,
  out_dir=out,alpha=.5,
  lambda_choice="lambda.min",seed=42L)
original <- read.csv(file.path(out,"predictions_test.csv"),
            stringsAsFactors=FALSE)
test(identical(as.character(original$sample_id),
               rownames(b$X)) &&
     isTRUE(all.equal(as.numeric(original$prob),
                       p$ORIGINAL_RLR,tolerance=1e-12)),
     "exact original manuscript RLR predictor parity")
unlink(out,recursive=TRUE)
s <- tempfile(fileext=".rds")
saveRDS(fit,s)
r <- rlmt4_predict(readRDS(s),b$X)
unlink(s)
test(identical(p,r),
     "five model predictions bit-exact RDS replay")

# Changes to absolute expression units preserving within-patient ranks
# do not change RELATIONAL predictions (including rank-distance gate).
trans <- as.data.frame(lapply(b$X,function(v) 3*v+11),
                       check.names=FALSE)
names(trans) <- names(b$X)
rownames(trans) <- rownames(b$X)
pt <- rlmt4_predict(fit,trans)
test(all(vapply(names(p),function(k)
    isTRUE(all.equal(p[[k]],pt[[k]],tolerance=1e-10)),
    logical(1))),
    "within-patient expression-order preserving affine invariance")
test(!anyNA(p$SOFT_GAP_ROUTE),
     "soft rank gap stable")
message("RLMT4_SYNTHETIC_REAL_RLR_PARITY_GATE=PASS")
message("NO_BIOMEDICAL_DATA_USED=TRUE")
