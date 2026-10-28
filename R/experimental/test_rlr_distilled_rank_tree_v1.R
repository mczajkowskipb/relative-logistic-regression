source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_distilled_rank_tree_core_v1.R")
set.seed(20261028L)
n<-96L;p<-68L
D<-as.data.frame(matrix(rnorm(n*p),nrow=n,ncol=p))
names(D)<-sprintf("gene%03d",seq_len(p))
rownames(D)<-sprintf("patient%03d",seq_len(n))
R<-rg_ranks(D)
# Nonlinear conditional rank truth: not expressible as ONE pair direction.
y<-factor(as.integer(ifelse(R[,"gene001"]>R[,"gene002"],
                        R[,"gene003"]>R[,"gene004"],
                        R[,"gene005"]>R[,"gene006"])),levels=0:1)
stopifnot(min(table(y))>=15)
pairs<-data.frame(gene_left=sprintf("gene%03d",seq(1,29,2)),
 gene_right=sprintf("gene%03d",seq(2,30,2)),
 stringsAsFactors=FALSE)
pairs$ordered_pair<-paste0(pairs$gene_left," > ",pairs$gene_right)
teacher<-pmin(.99,pmax(.01,.2+.6*as.numeric(y)))
allfit<-list(
 dt_train(D,y,pairs,teacher,.55,mode="BINARY"),
 dt_train(D,y,pairs,teacher,.55,mode="TERNARY"),
 dt_train(D,y,pairs,teacher,0,mode="BINARY"),
 dt_train(D,y,pairs,teacher,0,mode="TERNARY"))
for(model in allfit){
 res<-dt_predict(model,D)
 stopifnot(model$internal<=3L,model$n_leaves<=7L,
           model$n_unique_genes<=6L,
           all(res$depth<=2L),all(res$depth>=0L),
           all(res$prob>0 & res$prob<1),
           all(model$used_pairs<=nrow(pairs)),
           is.finite(dt_fidelity(teacher,res$prob)$soft_logloss))
 scaled<-D
 for(i in seq_len(n)){
   z<-unlist(D[i,,drop=FALSE],use.names=FALSE)
   scaled[i,]<-exp(z/runif(1,2,4))*runif(1,.4,2)
 }
 pp<-dt_predict(model,scaled)
 stopifnot(identical(res$prob,pp$prob),
           identical(res$depth,pp$depth),
           identical(res$path,pp$path),
           identical(res$leaf,pp$leaf))
}
# No teacher required at inference: fitted model contains only local
# selected ordinal node comparisons and leaf class probabilities.
stopifnot(!any(vapply(allfit,function(m)
    "teacher" %in% names(m) || "p_teacher" %in% names(m),logical(1))))
# A weak or zero gain cannot create a tree just because 3 nodes requested.
random_y<-factor(rep(c(0L,1L),length.out=n),levels=0:1)
nullfit<-dt_train(D,random_y,pairs,alpha=0,mode="BINARY",
                  penalty=20,max_internal=3L)
stopifnot(nullfit$internal==0L,
          all(dt_predict(nullfit,D)$depth==0L))
cat("RLRTREE_SYNTHETIC_MAX3_INTERNAL_MAX2_PATHS=PASS\n")
cat("RLRTREE_STRICT_PER_PATIENT_MONOTONE_TRANSFORM_INVARIANCE=PASS\n")
cat("RLRTREE_INFERENCE_DOES_NOT_RUN_OR_STORE_RLR_TEACHER=PASS\n")
cat("RLRTREE_PENALTY_CAN_PRUNE_TO_SINGLE_LEAF=PASS\n")
cat("RLRTREE_MECHANICAL_TESTS=PASS\n")
