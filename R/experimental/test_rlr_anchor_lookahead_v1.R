source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_anchor_lookahead_child_rlr_v1.R")
set.seed(621L)
n<-144L;p<-90L
X<-as.data.frame(matrix(rnorm(n*p),nrow=n,ncol=p))
names(X)<-sprintf("g%03d",seq_len(p))
rownames(X)<-sprintf("patient%03d",seq_len(n))
U<-X
for(i in seq_len(n)){
  z<-unlist(X[i,,drop=FALSE],use.names=FALSE)
  U[i,]<-exp(z/runif(1L,2,5))*exp(rnorm(1L))
}
r0<-aa_ranks(X);r1<-aa_ranks(U)
stopifnot(identical(unname(r0),unname(r1)))
root<-os_node("g001",c("g002","g003","g004"),2L)
gt<-os_gate(r0,root)
z<-as.integer(r0[,"g005"]>r0[,"g012"])
y<-factor(as.integer(xor(gt,z)),levels=0:1)
stopifnot(min(table(y))>40)
p0<-rep(.5,n)
pairs<-data.frame(gene_left=c("g005","g007","g009"),
 gene_right=c("g012","g019","g026"),stringsAsFactors=FALSE)
pairs$ordered_pair<-paste(pairs$gene_left,pairs$gene_right,sep=" > ")
fixed<-list(pairs,pairs)
model<-la_fit_structure(X,y,p0,root,pairs_fixed=fixed)
a<-la_predict_tree(X,p0,model)
b<-la_predict_tree(U,p0,model)
stopifnot(all(is.finite(a)),all(a>=0 & a<=1),
          isTRUE(all.equal(a,b,tolerance=0)))
cat("LOOKAHEAD_FROZEN_CHILD_RLR_SAMPLEWISE_INVARIANCE=PASS\n")

g0<-la_global_control(X,y,p0,X,p0,model)
g1<-la_global_control(X,y,p0,U,p0,model)
stopifnot(all(is.finite(g0)),isTRUE(all.equal(g0,g1,tolerance=0)))
v0<-la_gate_control(X,y,p0,X,p0,root)
v1<-la_gate_control(X,y,p0,U,p0,root)
stopifnot(all(is.finite(v0)),isTRUE(all.equal(v0,v1,tolerance=0)))
cat("LOOKAHEAD_GLOBAL_AND_GATE_CONTROLS_INVARIANT=PASS\n")

# Candidate ranking must not use raw-expression MAD.
pool<-do.call(rbind,lapply(seq(1L,19L,by=2L),function(i)
  data.frame(gene_left=names(X)[i],
             gene_right=names(X)[(i+1L):min(i+15L,p)],
             stringsAsFactors=FALSE)))
pool$ordered_pair<-paste(pool$gene_left,pool$gene_right,sep=" > ")
zprob<-plogis(rnorm(n,sd=.15))
ranks0<-la_roots(X,y,pool,zprob,max_hubs=8L,top_budget=5L)
ranks1<-la_roots(U,y,pool,zprob,max_hubs=8L,top_budget=5L)
stopifnot(length(ranks0$candidates)==length(ranks1$candidates),
          isTRUE(all.equal(ranks0$cheap_scores,
                           ranks1$cheap_scores,tolerance=0)))
if(length(ranks0$candidates)>0L)for(i in seq_along(ranks0$candidates)){
  stopifnot(identical(ranks0$candidates[[i]]$node,
                      ranks1$candidates[[i]]$node))
}
cat("LOOKAHEAD_OOF_RESIDUAL_SCREEN_RANK_INVARIANCE=PASS\n")

# The tree requires exactly the same schema at prediction time, preserving
# rank margin m=10 (all 90 genes, not just selected pair genes).
bad<-tryCatch({rg_ranks(U[,rev(names(U)),drop=FALSE],schema=names(X));FALSE},
              error=function(e)TRUE)
stopifnot(bad)
cat("LOOKAHEAD_FIXED_FEATURE_UNIVERSE=PASS\n")
cat("LOOKAHEAD_ALL_MECHANICAL_TESTS=PASS\n")
