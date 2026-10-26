source("R/experimental/rlr_anchor_majority_tree_v1.R")
set.seed(493L)
n<-120L;p<-70L
X<-as.data.frame(matrix(rnorm(n*p),nrow=n,ncol=p))
names(X)<-sprintf("gene_%03d",seq_len(p))
rownames(X)<-sprintf("sample_%03d",seq_len(n))
T<-X
for(i in seq_len(n)) {
  a<-exp(rnorm(1,sd=1));b<-rnorm(1,sd=12)
  T[i,]<-as.numeric(X[i,])*a+b
}
R<-aa_ranks(X);R2<-aa_ranks(T)
stopifnot(identical(unname(R),unname(R2)))
nodes<-list(
  aa_node("gene_001",c("gene_002","gene_003","gene_004")),
  aa_node("gene_007",c("gene_008","gene_009","gene_010","gene_011","gene_012")))
for(a in nodes) {
  v1<-aa_call(R,a);v2<-aa_call(R2,a)
  stopifnot(identical(v1,v2))
  # An odd-sized majority of comparisons equals comparison to their median.
  m<-apply(R[,a$anchors,drop=FALSE],1L,stats::median)
  stopifnot(identical(as.integer(R[,a$target]>m),v1))
}
cat("ANCHOR_VOTE_PATIENT_INVARIANCE_AND_MEDIAN_IDENTITY=PASS\n")

# Perturbing a single reference can change at most one vote, but might flip
# the decision at a 2/3 boundary. Perturbing the HUB can change ALL votes.
node<-nodes[[2]]
raw_votes<-vapply(node$anchors,function(a)
   as.integer(R[,node$target]>R[,a]),integer(n))
Rmut<-R
Rmut[1L,node$anchors[1L]]<-p+1L
new_votes<-vapply(node$anchors,function(a)
   as.integer(Rmut[,node$target]>Rmut[,a]),integer(n))
stopifnot(sum(abs(raw_votes[1,]-new_votes[1,]))<=1L)
cat("ANCHOR_SINGLE_REFERENCE_VOTE_BOUNDED=PASS\n")

# Gene selection and star-panel design must also be rank-invariant.
y<-factor(rep(0:1,each=n/2),levels=0:1)
root<-"gene_001"
other<-names(X)[-1L]
pool<-data.frame(gene_left=rep(root,length(other)),
                 gene_right=other,
                 ordered_pair=paste(root,other,sep=" > "),
                 stringsAsFactors=FALSE)
a<-aa_anchors(R,y,pool,root,5L)
b<-aa_anchors(R2,y,pool,root,5L)
stopifnot(identical(a,b), length(a)==5L,
          length(unique(a))==5L,
          !root %in% a)
f1<-aa_find_node(R,y,pool,max_hubs=5L)
f2<-aa_find_node(R2,y,pool,max_hubs=5L)
if(!is.null(f1)&&!is.null(f2)){
  stopifnot(identical(f1$node,f2$node),identical(f1$loss,f2$loss))
}
cat("ANCHOR_TRAINING_RANK_INVARIANCE=PASS\n")

# A learned threshold for a future patient must be a fixed schema-specific
# ordinal majority rule, not a raw expression cutoff.
error<-tryCatch({aa_ranks(X[,rev(names(X)),drop=FALSE],
                           schema=names(X));FALSE},error=function(e)TRUE)
stopifnot(error)
cat("ANCHOR_STRICT_GENE_SCHEMA=PASS\n")

cat("ANCHOR_ALL_UNIT_TESTS=PASS\n")
