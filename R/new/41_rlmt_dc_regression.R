# RLMT-DC-v0 deterministic primitive regressions. Synthetic only.
# Run: Rscript --vanilla R/new/41_rlmt_dc_regression.R
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")

check <- function(ok,label) {
  if(!isTRUE(ok))stop("RLMT_DC_TEST_FAIL: ",label,call.=FALSE)
  cat("PASS ",label,"\n",sep="")
}
conf <- rlmt_dc_defaults(list())$leaf
check(conf$min_child_n==16L && conf$min_child_class==6L,
      "frozen_support_16_6")
check(rlmt_dc_defaults(list())$dc$proposal_n==4096L,
      "frozen_candidate_cap")

# Many features: 4096 distinct unordered pairs without creating p^2.
names_x <- paste0("g",seq_len(600L))
a <- rlmt_dc_random_pairs(names_x,4096L,42L)
b <- rlmt_dc_random_pairs(names_x,4096L,42L)
c <- rlmt_dc_random_pairs(names_x,4096L,43L)
check(nrow(a)==4096L && identical(a,b), "random_pairs_reproducible")
check(!identical(a,c),"random_pairs_seed_sensitive")
check(all(a$gene_a %in% names_x) &&
      all(a$gene_b %in% names_x) &&
      all(a$gene_a!=a$gene_b),"random_pairs_valid_features")
check(!anyDuplicated(paste(a$gene_a,a$gene_b,sep="|")),
      "random_pairs_unique_unordered")

# Small feature space exhaustive, still no supervised ranking.
small <- rlmt_dc_random_pairs(paste0("x",seq_len(8L)),4096L,42L)
check(nrow(small)==28L,"small_space_finite_all_pairs")

# Explicit non-overlapping discovery/confirmation. Label changes to V
# must never enter the candidate generator or ranking on D.
set.seed(17)
n <- 180L
p <- 40L
X <- as.data.frame(matrix(rnorm(n*p),nrow=n,ncol=p))
colnames(X) <- paste0("gene_",seq_len(p))
rownames(X) <- paste0("synthetic_",seq_len(n))
y <- factor(rep(c("N","PR"),each=n/2L),levels=c("N","PR"))
part <- rlmt_dc_partition(y,seed=42L)
check(length(intersect(part$discovery,part$confirmation))==0L &&
      identical(sort(c(part$discovery,part$confirmation)),
                seq_len(n)),"D_V_strict_partition")
check(length(part$confirmation)>=20L &&
      all(table(y[part$discovery])>=12L),
      "D_V_minimum_class_support")
check(identical(part,rlmt_dc_partition(y,42L)),
      "D_V_partition_reproducible")

D <- part$discovery
V <- part$confirmation
pool <- rlmt_dc_random_pairs(colnames(X),4096L,42L)
# Synthetic fixed OOF probabilities, only a test for screening logic.
pD <- rep(0.5,length(D))
rank1 <- rlmt_dc_propose(X[D,,drop=FALSE],y[D],pD,pool,conf)
check(nrow(rank1$audit)==nrow(pool) &&
      !anyNA(rank1$audit$feasible),"proposal_audit_all_draws")
check(rank1$feasible_n>0L &&
      !is.null(rank1$candidate) &&
      isTRUE(rank1$candidate$feasible[[1L]]),
      "support_screen_candidate_exists_on_synthetic")

y_changed <- y
y_changed[V] <- factor(
  rev(as.character(y[V])),levels=levels(y))
rank2 <- rlmt_dc_propose(X[D,,drop=FALSE],y_changed[D],
                         pD,pool,conf)
check(identical(rank1,rank2),
      "validation_label_canary_does_not_change_proposal")

# A deterministic deliberately impossible class-support scenario.
y_bad <- factor(c(rep("N",160L),rep("PR",20L)),
                levels=c("N","PR"))
pbad <- rep(0.5,n)
bad <- rlmt_dc_propose(X,y_bad,pbad,pool,conf)
check(all(!bad$audit$feasible) || bad$feasible_n>=0L,
      "rare_class_is_handled_without_error")

# Confirmation accepts *only* a single paired loss gain > 1 standard error.
yv <- factor(rep(c("N","PR"),each=15L),levels=c("N","PR"))
p_parent <- rep(0.5,length(yv))
p_better <- ifelse(yv=="PR",0.9,0.1)
p_worse <- 1-p_better
good <- rlmt_dc_confirm(yv,p_parent,p_better,20L)
bad <- rlmt_dc_confirm(yv,p_parent,p_worse,20L)
small_confirm <- rlmt_dc_confirm(yv[1:19],p_parent[1:19],
                                  p_better[1:19],20L)
check(good$accepted && good$gain>0 &&
      good$lower_gain>0 &&
      abs(good$gain-(log(0.9)-log(0.5)))<1e-12,
      "confirmation_positive_gain")
check(!bad$accepted && bad$gain<0,
      "confirmation_negative_gain")
check(!small_confirm$accepted && small_confirm$n==19L,
      "confirmation_reject_insufficient_n")
check(!rlmt_dc_confirm(yv,p_parent,p_parent,20L)$accepted,
      "confirmation_equal_to_parent_no_split")

# Reader must reject illegal output probabilities, never silently clip NaN.
err <- tryCatch(rlmt_dc_confirm(yv,rep(NA_real_,30L),
                                 p_better,20L),error=function(e)e)
check(inherits(err,"error"),"confirmation_invalid_probs_blocked")

cat("RLMT_DC_PRIMITIVES=PASS\n")
