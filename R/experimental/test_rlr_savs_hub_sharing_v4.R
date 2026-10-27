source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_savs_hub_sharing_v4.R")
set.seed(8711L)
n<-110L; p<-70L
X<-as.data.frame(matrix(rnorm(n*p),nrow=n,ncol=p))
names(X)<-sprintf("g%02d",seq_len(p))
rownames(X)<-sprintf("P%03d",seq_len(n))
R<-rg_ranks(X)
root<-os_node("g01",c("g02","g03","g04"),2L)
lab<-as.integer((R[,"g01"]>R[,"g05"]) |
                ((R[,"g01"]>R[,"g06"]) & (R[,"g01"]>R[,"g07"])))
y<-factor(lab,levels=0:1)
stopifnot(min(table(y))>10L)
# Only child candidates reuse the SAME g01 used at root, but not the exact
# g01/g02..g04 root pair. A hub-excluding model must have ZERO leaf rules.
pool<-data.frame(gene_left=rep("g01",25),
                 gene_right=sprintf("g%02d",8:32),
                 stringsAsFactors=FALSE)
pool$ordered_pair<-paste(pool$gene_left,pool$gene_right,sep=" > ")
original<-sv_fit(X,y,pool,root=root,mode="RLR10")
shared<-sv4_fit(X,y,pool,root,mode="RLR10")
stopifnot(sum(original$n_rules)==0L,
          sum(shared$n_rules)>0L,
          sum(shared$n_path_gene_overlap)>0L)
for(b in 1:2) if(nrow(shared$rule_panels[[b]])){
  z<-shared$rule_panels[[b]]
  stopifnot(length(unique(c(z$gene_left,z$gene_right)))==
                      2L*nrow(z))
  for(j in seq_len(nrow(z)))stopifnot(
    !any((root$target==z$gene_left[j]&root$anchors==z$gene_right[j])|
         (root$target==z$gene_right[j]&root$anchors==z$gene_left[j])))
}
cat("SAVS4_ROOT_HUB_ALLOWED_BUT_IDENTICAL_PAIR_FORBIDDEN=PASS\n")
p1<-sv_predict(shared,X)
Z<-X
for(i in seq_len(n)){
 a<-unlist(X[i,,drop=FALSE],use.names=FALSE)
 Z[i,]<-exp(a/runif(1,2,5))*runif(1,.4,3)
}
p2<-sv_predict(shared,Z)
stopifnot(identical(p1$p,p2$p),
          identical(p1$n_active,p2$n_active),
          identical(p1$branch,p2$branch),
          all(p1$n_active<=p1$n_considered),
          all(p1$n_considered<=4L),
          all(p1$n_root==3L))
cat("SAVS4_PATHWISE_REUSE_PATIENT_MONOTONE_INVARIANT=PASS\n")
cat("SAVS4_ALL_MECHANICAL_TESTS=PASS\n")
