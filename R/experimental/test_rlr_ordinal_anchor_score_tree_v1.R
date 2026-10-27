# Standalone mechanical tests: ordinal, patientwise invariance, exact routes.
source("R/experimental/rlr_ordinal_anchor_score_tree_v1.R")
set.seed(919L)
n<-130L;p<-80L
X<-as.data.frame(matrix(rnorm(n*p),nrow=n,ncol=p))
names(X)<-sprintf("gene_%03d",seq_len(p))
rownames(X)<-sprintf("patient_%03d",seq_len(n))
r0<-aa_ranks(X)
X2<-X
for(i in seq_len(n)){
  v<-unlist(X[i,,drop=FALSE],use.names=FALSE)
  # independent nonlinear, strictly increasing patient transform
  X2[i,]<-exp(v/runif(1L,2,6))*runif(1L,.1,2)
}
r2<-aa_ranks(X2)
stopifnot(identical(unname(r0),unname(r2)))
root<-os_node("gene_001",
             c("gene_002","gene_003","gene_004","gene_005","gene_006"),
             2L)
child0<-os_node("gene_007",c("gene_008","gene_009","gene_010"),1L)
child1<-os_node("gene_011",c("gene_012","gene_013","gene_014"),2L)
tree<-list(root=root,children=list(child0,child1))
for(node in c(list(root,child0,child1))){
  a<-os_score(r0,node)
  b<-os_score(r2,node)
  stopifnot(identical(a,b),all(a>=0L),all(a<=length(node$anchors)),
            identical(as.integer(a>=node$threshold),os_gate(r0,node)))
}
stopifnot(identical(os_route(r0,tree,2L),os_route(r2,tree,2L)))
stopifnot(all(os_route(r0,tree,2L) %in% 3:6))
stopifnot(length(unique(os_route(r0,tree,2L)))==4L)
cat("OASDT_PATIENTWISE_MONOTONE_INVARIANCE=PASS\n")
cat("OASDT_REAL_TWO_LEVEL_PATHS=PASS\n")

# Majority is precisely median comparison with odd panel, but score
# thresholds 1/2/4/5 retain information beyond its binary majority.
major<-os_node(root$target,root$anchors,3L)
stopifnot(identical(os_gate(r0,major),
     as.integer(r0[,root$target]>
       apply(r0[,root$anchors,drop=FALSE],1L,median))))
stopifnot(length(unique(os_score(r0,root)))>=4L)
cat("OASDT_MAJORITY_MEDIAN_AND_MULTI_STATE=PASS\n")

# One changed anchor comparison changes count by at most one;
# a changed hub might flip all reference comparisons.
z0<-os_score(r0,root)
mod<-r0
mod[1,root$anchors[1L]]<-ncol(r0)+10
z1<-os_score(mod,root)
stopifnot(abs(z1[1L]-z0[1L])<=1L,
          identical(z0[-1L],z1[-1L]))
cat("OASDT_SINGLE_REFERENCE_BOUNDED=PASS\n")

# Fixed gene schema must be a hard requirement at inference.
stopped<-tryCatch({
   rg_ranks(X2[,rev(names(X2)),drop=FALSE],schema=names(X))
   FALSE
},error=function(e)TRUE)
stopifnot(stopped)
cat("OASDT_GENE_SCHEMA_LOCK=PASS\n")

# A deliberately planted 2-level route-dependent phenotype should be
# learnable by leaf probability estimation given the TRUE tree.
rte<-os_route(r0,tree,2L)
y<-factor(as.integer(rte %in% c(4L,6L)),levels=0:1)
f<-os_fit_leaf(r0,y,tree,2L)
prob<-os_predict_leaf(f,r0)
stopifnot(length(prob)==n,
          all(prob>=0 & prob<=1),
          rg_ba(os_truth(y),prob)>.9,
          rg_ll(os_truth(y),prob)<rg_ll(os_truth(y),rep(.5,n)))
cat("OASDT_PLANTED_DEPTH2_LEAF_SCORER=PASS\n")

# Non-tree comparator must remain estimable with only one root score.
one_node<-list(root=root,children=list(NULL,NULL))
glm<-os_global_fit(r0,y,one_node,depth=1L)
pp<-os_global_predict(r2,glm)
stopifnot(isTRUE(glm$padded),length(pp)==n,
          all(is.finite(pp)),all(pp>=0&pp<=1))
cat("OASDT_SINGLE_SCORE_INERT_DUMMY_GLMNET=PASS\n")
cat("OASDT_ALL_MECHANICAL_TESTS=PASS\n")
