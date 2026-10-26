source("R/experimental/rlr_rank_contrast_core_v1.R")
set.seed(818L)
X<-as.data.frame(matrix(rnorm(120L*60L),nrow=120L,ncol=60L))
names(X)<-sprintf("gene_%03d",seq_len(ncol(X)))
rownames(X)<-sprintf("sample_%03d",seq_len(nrow(X)))
X1<-X
for(i in seq_len(nrow(X))) {
  a<-exp(rnorm(1,sd=.8));b<-rnorm(1,sd=12)
  X1[i,]<-unlist(X[i,,drop=FALSE],use.names=FALSE)*a+b
}
P<-data.frame(gene_left=sprintf("gene_%03d",seq(1,11,by=2)),
              gene_right=sprintf("gene_%03d",seq(2,12,by=2)))
A<-rg_pairs(X,P)
B<-rg_pairs(X1,P)
stopifnot(identical(A,B))
a<-rc_design(X,A,1L,2L,3L)
b<-rc_design(X1,B,1L,2L,3L)
stopifnot(identical(unname(a),unname(b)))
stopifnot(all(a[,"root"] %in% 0:1))
stopifnot(all(a[,"partner"] %in% c(-1,0,1)))
cat("RANK_CONTRAST_WITHIN_PATIENT_MONOTONE_INVARIANCE=PASS\n")

# Binary context signal entirely determined by frozen patient-wise ranks.
rk<-rg_ranks(X)
d<-rg_gaps(rk,A)
root<-as.integer(d[,1]>d[,2])
partner<-as.integer(d[,3]>0)
y<-factor(as.integer(xor(as.logical(root),as.logical(partner))),
          levels=0:1)
S<-rc_candidates(X,y,rep(.5,nrow(X)),P,
                 max_roots=20L,max_contrasts=50L,max_return=5L)
stopifnot(nrow(S$candidates)>=1L,
          all(is.finite(S$candidates$score)),
          all(S$candidates$root_n_left>=12L),
          all(S$candidates$root_n_right>=12L))
cat("RANK_CONTRAST_SYNTHETIC_PAIR_GAP_PROPOSALS=PASS\n")
cat("RANK_CONTRAST_ALL_UNIT_TESTS=PASS\n")
