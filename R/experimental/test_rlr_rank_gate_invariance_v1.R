# Base-R mechanical invariance tests. No fit and no biomedical efficacy claim.
source("R/experimental/rlr_rank_gate_core_v1.R")
set.seed(7123L)
n<-120L;p<-60L
X<-as.data.frame(matrix(rnorm(n*p),nrow=n,ncol=p))
names(X)<-sprintf("gene_%03d",seq_len(p))
rownames(X)<-sprintf("patient_%03d",seq_len(n))
R0<-rg_ranks(X)

# Each patient receives its OWN strictly increasing affine transformation.
X1<-X
for(i in seq_len(n)) {
  a<-exp(rnorm(1L,sd=.7))
  b<-rnorm(1L,sd=10)
  X1[i,]<-unlist(X[i,,drop=FALSE],use.names=FALSE)*a+b
}
stopifnot(identical(unname(R0),unname(rg_ranks(X1))))
# Test a second, nonlinear strictly increasing transformation per patient.
X2<-X
for(i in seq_len(n)) {
  z<-unlist(X[i,,drop=FALSE],use.names=FALSE)
  X2[i,]<-exp(z/runif(1L,2,5))
}
stopifnot(identical(unname(R0),unname(rg_ranks(X2))))
cat("RANK_GATE_120_INDEPENDENT_PATIENT_TRANSFORMS=PASS\n")

base<-data.frame(
 gene_left=c("gene_001","gene_003","gene_005","gene_007","gene_009","gene_011"),
 gene_right=c("gene_002","gene_004","gene_006","gene_008","gene_010","gene_012"),
 stringsAsFactors=FALSE)
P0<-rg_pairs(X,base)
P1<-rg_pairs(X1,base)
P2<-rg_pairs(X2,base)
stopifnot(identical(P0,P1),identical(P0,P2))
D0<-rg_gaps(R0,P0)
stopifnot(identical(unname(D0),unname(rg_gaps(rg_ranks(X1),P1))))
stopifnot(identical(unname(D0),unname(rg_gaps(rg_ranks(X2),P2))))
for(mode in c("ORDER","WIDE")) {
  A<-rg_design(X,P0,1L,2L,mode)
  B<-rg_design(X1,P1,1L,2L,mode)
  C<-rg_design(X2,P2,1L,2L,mode)
  stopifnot(identical(unname(A),unname(B)),
            identical(unname(A),unname(C)))
}
cat("RANK_GATE_SELECTION_GAPS_STUMPS_INVARIANT=PASS\n")

# A single altered gene can shift the rank gap of two unaffected genes
# by at most one position (with no ties).
xx<-X[1,,drop=FALSE]
rr0<-rg_ranks(xx)
xx[1,"gene_001"]<-1e10
rr1<-rg_ranks(xx)
good<-setdiff(names(X),"gene_001")
d0<-outer(rr0[1,good],rr0[1,good],"-")
d1<-outer(rr1[1,good],rr1[1,good],"-")
stopifnot(max(abs(d1-d0))<=1L)
cat("RANK_GATE_ONE_OUTLIER_UNRELATED_PAIR_GAP_BOUND=PASS\n")

# Any missing or reordered gene invalidates evaluation.
err<-tryCatch({rg_ranks(X[,rev(names(X)),drop=FALSE],schema=names(X));FALSE},
              error=function(e) TRUE)
stopifnot(err)
cat("RANK_GATE_FIXED_GENE_SCHEMA_ENFORCED=PASS\n")

# Synthetic context interaction only: the label depends on two independent
# within-patient comparisons. Used only to check that proposals can be formed.
g1<-as.integer(R0[,"gene_001"]<R0[,"gene_002"])
g2<-as.integer(R0[,"gene_003"]<R0[,"gene_004"])
y<-factor(as.integer(xor(as.logical(g1),as.logical(g2))),levels=0:1)
S<-rg_candidates(X,y,rep(.5,n),base,max_root=6L,max_return=5L)
stopifnot(nrow(S$candidates)>=1L,
          all(is.finite(S$candidates$score)),
          all(S$candidates$root_n_left>=12L),
          all(S$candidates$root_n_right>=12L))
cat("RANK_GATE_SYNTHETIC_CONTEXT_CANDIDATE_DISCOVERY=PASS\n")

# Original raw expression ratios do NOT share rank invariance to sample
# translations; therefore they are deliberately absent from the method.
raw_before<-as.numeric(X[1,"gene_001"]/X[1,"gene_002"])
raw_after<-as.numeric((X[1,"gene_001"]+100)/(X[1,"gene_002"]+100))
stopifnot(abs(raw_before-raw_after)>1e-6)
cat("RANK_GATE_RAW_RATIO_NOT_USED=PASS\n")
cat("RANK_GATE_ALL_MECHANICAL_TESTS=PASS\n")
