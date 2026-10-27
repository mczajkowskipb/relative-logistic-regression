source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_quantile_ordinal_margin_core_v1.R")
set.seed(20261027L)
N<-100L;P<-80L
X<-as.data.frame(matrix(rnorm(N*P),nrow=N,ncol=P))
names(X)<-sprintf("g%03d",seq_len(P))
rownames(X)<-sprintf("p%03d",seq_len(N))
rk<-rg_ranks(X)
y<-factor(as.integer(rk[,"g001"]>rk[,"g002"]),
          levels=0:1)
pairs<-data.frame(gene_left=c("g001","g003","g005"),
 gene_right=c("g002","g004","g006"),stringsAsFactors=FALSE)
pairs$ordered_pair<-paste0(pairs$gene_left," > ",pairs$gene_right)
m0<-qr_threshold(X,pairs,0)
stopifnot(identical(as.numeric(m0),c(10,10,10)),
          identical(qr_encode(X,pairs,m0),encode_pair_matrix(X,pairs,m=10L)))
for(q in c(.25,.5)){
 mar<-qr_threshold(X,pairs,q)
 train<-qr_encode(X,pairs,mar)
 stopifnot(all(mar>=10),all(train %in% c(-1,0,1)))
 Y<-X
 for(i in seq_len(N)){
  a<-unlist(X[i,,drop=FALSE],use.names=FALSE)
  Y[i,]<-exp(a/runif(1L,3,5))*runif(1L,.5,1.5)
 }
 stopifnot(identical(mar,qr_threshold(Y,pairs,q)),
           identical(train,qr_encode(Y,pairs,mar)))
}
cat("QRLR_Q0_IDENTICAL_ORIGINAL_SIGNED_M10_RELATIONS=PASS\n")
cat("QRLR_PER_PATIENT_NONLINEAR_MONOTONE_INVARIANT=PASS\n")
cat("QRLR_NO_LABELS_USED_TO_ESTIMATE_PAIR_MARGIN=PASS\n")
cat("QRLR_MECHANICAL_TESTS=PASS\n")
