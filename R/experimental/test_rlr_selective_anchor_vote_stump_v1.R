source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_selective_anchor_vote_stump_v1.R")
set.seed(8701L)
n<-120L;pp<-75L
X<-as.data.frame(matrix(rnorm(n*pp),nrow=n,ncol=pp))
names(X)<-sprintf("p%03d",seq_len(pp))
rownames(X)<-sprintf("sample%03d",seq_len(n))
R<-rg_ranks(X)
root<-os_node("p001",c("p002","p003","p004"),2L)
gate<-os_gate(R,root)
f1<-as.integer(R[,"p012"]>R[,"p015"])
f2<-as.integer(R[,"p016"]>R[,"p021"])
y<-factor(as.integer(ifelse(gate==0L,f1,f2)),levels=0:1)
stopifnot(min(table(y))>=10L,
          all(table(factor(gate,levels=0:1))>=15L))
pool<-data.frame(
 gene_left=sprintf("p%03d",c(12,16,23,26,29,32,35,38,41,44,47,50,53)),
 gene_right=sprintf("p%03d",c(15,21,24,27,30,33,36,39,42,45,48,51,54)),
 stringsAsFactors=FALSE)
model<-sv_fit(X,y,pool=pool,root=root,mode="RLR10")
pred<-sv_predict(model,X)
stopifnot(length(pred$p)==nrow(X),
 all(pred$p>0 & pred$p<1),
 all(pred$n_active<=pred$n_considered),
 all(pred$n_active>=0),
 all(pred$n_considered<=4L),
 all(pred$n_root==3L))
cat("SAVS_COMPACT_ROOT_AND_MAX_FOUR_BRANCH_RULES=PASS\n")
# One-sided semantics: opposite signed rank comparison ABSTAINS,
# unlike the ternary RLR where an opposite ordering is -1 and ACTIVE.
r1<-data.frame(gene_left="p012",gene_right="p015",fire_sign=1L,
 vote_sign=1L,strength=1)
fire<-sv_active(R,r1,10L)[,1]
gap<-R[,"p015"]-R[,"p012"]
stopifnot(identical(unname(fire),unname(gap>=10)),
 all(!fire[gap<=-10]))
cat("SAVS_ONE_SIDED_ABSTENTION_NOT_ORIGINAL_RLR_ZERO=PASS\n")
U<-X
for(i in seq_len(n)) {
 a<-unlist(X[i,,drop=FALSE],use.names=FALSE)
 U[i,]<-exp(a/runif(1L,2,7))*runif(1L,.25,3)
}
stopifnot(isTRUE(all.equal(unname(rg_ranks(X)),unname(rg_ranks(U)),tolerance=0)))
ppred<-sv_predict(model,U)
stopifnot(identical(pred$p,ppred$p),
 identical(pred$n_active,ppred$n_active),
 identical(pred$branch,ppred$branch))
global<-sv_samepair_global(X,y,model)
gg<-sv_predict(global,X)
gg2<-sv_predict(global,U)
stopifnot(identical(gg$p,gg2$p),
 all(gg$n_active<=gg$n_considered),
 all(gg$n_root==0L))
root_only<-sv_root_only(X,y,model)
stopifnot(all(sv_predict(root_only,X)$n_active==0L))
cat("SAVS_PATIENTWISE_STRICT_MONOTONE_INVARIANCE=PASS\n")
cat("SAVS_SAME_NOMINATED_RULES_NON_TREE_CONTROL=PASS\n")
# No trainable logistic or per-leaf multiclass minimum: even with one
# minority observation in training, the one-sided posterior remains finite.
rare<-factor(c(rep(0L,n-1L),1L),levels=0:1)
small<-sv_fit(X,rare,pool=pool,root=root,mode="RELATIVE10PCT")
stopifnot(all(is.finite(sv_predict(small,X)$p)))
cat("SAVS_SMALL_MINORITY_NO_PER_BRANCH_GLMNET=PASS\n")
cat("SAVS_ALL_MECHANICAL_TESTS=PASS\n")
