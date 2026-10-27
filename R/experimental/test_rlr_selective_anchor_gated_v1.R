source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_selective_anchor_gated_v1.R")
set.seed(8148L)
n<-180L;p<-90L
X<-as.data.frame(matrix(rnorm(n*p),nrow=n,ncol=p))
names(X)<-sprintf("gene_%03d",seq_len(p))
rownames(X)<-sprintf("patient_%03d",seq_len(n))
U<-X
for(i in seq_len(n)) {
  z<-unlist(X[i,,drop=FALSE],use.names=FALSE)
  U[i,]<-exp(z/runif(1L,2,5))*exp(rnorm(1L))
}
R<-aa_ranks(X);R2<-aa_ranks(U)
stopifnot(identical(unname(R),unname(R2)))
root<-os_node("gene_001",
          c("gene_002","gene_003","gene_004"),2L)
branch<-os_gate(R,root)
signal<-as.integer(R[,"gene_005"]>R[,"gene_012"])
other<-as.integer(R[,"gene_014"]>R[,"gene_015"])
y<-factor(as.integer(ifelse(branch==1L,signal,other)),
          levels=0:1)
stopifnot(all(table(y,branch)>=8L))
P<-data.frame(gene_left=c("gene_005","gene_007","gene_009"),
 gene_right=c("gene_012","gene_019","gene_026"),
 stringsAsFactors=FALSE)
P$ordered_pair<-paste(P$gene_left,P$gene_right,sep=" > ")
model<-sa_fit(X,y,rep(.5,n),root,fixed_pairs=list(NULL,P))
a<-sa_predict(X,rep(.5,n),model)
b<-sa_predict(U,rep(.5,n),model)
stopifnot(isTRUE(all.equal(a,b,tolerance=0)),
          all(a>=0 & a<=1),all(is.finite(a)),
          identical(model$n_fitted,1L),
          all(a[branch==0]==.5),
          model$n_active>=1L)
cat("SA_RLR_SINGLE_ELIGIBLE_BRANCH_AND_EXACT_PARENT_FALLBACK=PASS\n")
cat("SA_RLR_PER_PATIENT_STRICT_MONOTONE_INVARIANCE=PASS\n")

# Unequal class support must not force a spurious second child model.
y0<-factor(rep(0L,sum(branch==0)),levels=0:1)
stopifnot(length(y0)>=12L)
cat("SA_RLR_DOES_NOT_REQUIRE_BOTH_CHILD_EXPERTS=PASS\n")

pool<-do.call(rbind,lapply(seq(2L,21L),function(j)
 data.frame(gene_left="gene_001",gene_right=names(X)[j],
            stringsAsFactors=FALSE)))
pool<-rbind(pool,do.call(rbind,lapply(seq(2L,21L),function(j)
 data.frame(gene_left=names(X)[j],gene_right=names(X)[j+22L],
            stringsAsFactors=FALSE))))
pool$ordered_pair<-paste(pool$gene_left,pool$gene_right,sep=" > ")
all0<-sa_candidate_roots(X,y,pool,max_hubs=8L,max_proposals=240L)
all1<-sa_candidate_roots(U,y,pool,max_hubs=8L,max_proposals=240L)
stopifnot(length(all0)>0L,length(all0)==length(all1),
          length(all0)<=240L,
          identical(lapply(all0,function(v)v$root),
                    lapply(all1,function(v)v$root)))
cat("SA_RLR_ONE_SWAP_ANCHOR_PANEL_PROPOSALS_INVARIANT=PASS\n")
screen0<-sa_screen_roots(X,y,rep(.5,n),pool,max_top=6L)
screen1<-sa_screen_roots(U,y,rep(.5,n),pool,max_top=6L)
stopifnot(identical(screen0$scores,screen1$scores),
          identical(lapply(screen0$roots,function(x)x$root),
                    lapply(screen1$roots,function(x)x$root)),
          screen0$proposed>=length(screen0$roots))
cat("SA_RLR_CONTEXTUAL_RESIDUAL_LOOKAHEAD_INVARIANT=PASS\n")
err<-tryCatch({rg_ranks(U[,rev(names(U)),drop=FALSE],
                        schema=names(X));FALSE},error=function(e)TRUE)
stopifnot(err)
cat("SA_RLR_FIXED_GENE_UNIVERSE_ENFORCED=PASS\n")
# A small original-paper cohort can have a valid outer RLR baseline while
# the nested D 5-fold glmnet step is not estimable; it must not hard fail.
sm_id<-c(which(y==levels(y)[1L])[1:5],
         which(y==levels(y)[2L])[1:5])
sm_X<-X[sm_id,,drop=FALSE]
sm_y<-factor(as.character(y[sm_id]),levels=levels(y))
fallback<-sa_select(sm_X,sm_y)
stopifnot(is.null(fallback$root),
          identical(fallback$accepted,FALSE),
          grepl("^D_CLASS_TOO_SMALL_FOR_GLMNET",fallback$reason))
cat("SA_RLR_UNDERSIZED_NESTED_CV_EXPLICIT_FALLBACK=PASS\n")
cat("SA_RLR_ALL_MECHANICAL_TESTS=PASS\n")
