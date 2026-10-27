source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_margin_aware_budget_core_v1.R")
set.seed(20261027L)
X<-as.data.frame(matrix(rnorm(100L*90L),nrow=100L,ncol=90L))
names(X)<-sprintf("g%03d",seq_len(ncol(X)))
rownames(X)<-sprintf("s%03d",seq_len(nrow(X)))
rk<-rg_ranks(X)
y<-factor(as.integer(rk[,"g001"]>rk[,"g002"] |
    (rk[,"g003"]>rk[,"g004"] & rk[,"g005"]>rk[,"g006"])),
    levels=0:1)
stopifnot(min(table(y))>3)
cmb<-utils::combn(names(X)[1:30],2)
pool<-data.frame(gene_left=cmb[1,],gene_right=cmb[2,],
    stringsAsFactors=FALSE)
pool$ordered_pair<-paste(pool$gene_left,">",pool$gene_right)
s<-mb_rank_margin_scores(X,y,pool,m=10L)
stopifnot(nrow(s)==nrow(pool),
    all(is.finite(s$score)),
    all(s$score>=0),
    all(s$pct_zero>=0 & s$pct_zero<=1),
    all(s$active_class0>=0&s$active_class0<=1),
    all(s$active_class1>=0&s$active_class1<=1))
ordered<-mb_reorder(X,y,pool,m=10L)
stopifnot(identical(sort(ordered$ordered_pair),sort(pool$ordered_pair)),
          all(diff(ordered$margin_aware_score)<=1e-12))
Y<-X
for(i in seq_len(nrow(X))) {
    a<-unlist(X[i,,drop=FALSE],use.names=FALSE)
    Y[i,]<-exp(a/runif(1,2,5))*runif(1,.5,2)
}
s2<-mb_rank_margin_scores(Y,y,pool,m=10L)
stopifnot(isTRUE(all.equal(s,s2,tolerance=1e-12)),
          identical(mb_reorder(X,y,pool,m=10L)$ordered_pair,
                    mb_reorder(Y,y,pool,m=10L)$ordered_pair))
stopifnot(mb_budget(0)==2L,mb_budget(16.6)==17L,
          mb_budget(25.55)==26L,mb_budget(100)==50L)
cat("MB_RLR_STRICT_PER_PATIENT_MONOTONE_INVARIANCE=PASS\n")
cat("MB_RLR_EXACT_ORIGINAL_M10_FEATURE_REORDER_SCORES=PASS\n")
cat("MB_RLR_TRAIN_ONLY_EFFECTIVE_K_BOUNDS=PASS\n")
cat("MB_RLR_UNIT_COMPLETE=PASS\n")
