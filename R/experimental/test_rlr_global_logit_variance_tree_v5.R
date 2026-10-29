source("R/experimental/rlr_global_logit_variance_tree_core_v5.R")
set.seed(20261029L)
n <- 90L
p <- 36L
X <- as.data.frame(matrix(rnorm(n*p),nrow=n))
names(X) <- sprintf("g%02d",seq_len(p))
rownames(X) <- sprintf("patient%03d",seq_len(n))
r <- rrv5_ranks(X)
y <- factor(as.integer((r[,"g02"]>r[,"g01"]) |
                       (r[,"g04"]>r[,"g03"])),levels=0:1)
stopifnot(min(table(y))>=10L)
pairs <- data.frame(gene_left=sprintf("g%02d",seq(1L,15L,2L)),
                    gene_right=sprintf("g%02d",seq(2L,16L,2L)))
# Synthetic proxy scores only verify mechanics; NOT real RLR OOF.
oo <- pmin(.99,pmax(.01,.25+.50*as.integer(y==levels(y)[2L])+
                         rnorm(n,sd=.05)))
fit0 <- rrv5_fit(X,y,pairs,weight=0,min_leaf=6L)
fit1 <- rrv5_fit(X,y,pairs,p_oof=oo,weight=.25,min_leaf=6L)
for (fit in list(fit0,fit1)) {
  pr <- rrv5_predict(fit,X)
  stopifnot(fit$n_internal<=3L,fit$n_unique_pairs<=3L,
            fit$n_unique_genes<=6L,
            max(pr$n_evaluated)<=2L,
            all(pr$n_evaluated %in% 0:2),
            !("teacher" %in% names(fit)),
            !("p_oof" %in% names(fit)),
            length(fit$pairs$gene_left)==fit$n_unique_pairs,
            is.finite(fit$objective),fit$objective<=1+1e-10)
  scaled <- X
  for (i in seq_len(n)) {
    v <- as.numeric(unlist(X[i,,drop=FALSE],use.names=FALSE))
    scaled[i,] <- exp(v/runif(1L,2,4))*runif(1L,.5,1.5)
  }
  pp <- rrv5_predict(fit,scaled)
  stopifnot(identical(pr$prob,pp$prob),
            identical(pr$n_evaluated,pp$n_evaluated),
            identical(pr$path,pp$path))
}
# Independent enumeration of every root and both possible secondary
# tests for the LABEL-only global loss, against fit0's reported optimum.
A <- rrv5_matrix(X,pairs)
yy <- as.integer(y==levels(y)[2L])
prior <- (sum(yy)+1)/(n+2)
ce <- function(idx) {
  pr <- (sum(yy[idx])+4*prior)/(length(idx)+4)
  -sum(yy[idx]*log(pr)+(1-yy[idx])*log1p(-pr))
}
den <- ce(seq_len(n))
cost <- function(ix) ce(ix)/den
best <- cost(seq_len(n))
candidate <- function(idx,feature) {
  if (feature==0L) return(list(loss=cost(idx),nodes=0L,checks=0))
  lo <- idx[A[idx,feature]==0L]; hi <- idx[A[idx,feature]==1L]
  if (min(length(lo),length(hi))<6L) return(NULL)
  list(loss=cost(lo)+cost(hi),nodes=1L,checks=length(idx)/n)
}
for (root in seq_len(ncol(A))) {
  l <- which(A[,root]==0L); h <- which(A[,root]==1L)
  if (min(length(l),length(h))<6L) next
  for (a in 0:ncol(A)) for (b in 0:ncol(A)) {
    if (a==root || b==root) next
    ca <- candidate(l,a);cb <- candidate(h,b)
    if (is.null(ca)||is.null(cb))next
    s <- ca$loss+cb$loss+
         .005*(1+ca$nodes+cb$nodes)+
         .005*(1+ca$checks+cb$checks)
    best <- min(best,s)
  }
}
stopifnot(abs(best-fit0$objective)<1e-9)
cat("RLR_GLVT_V5_EXHAUSTIVE_GLOBAL_OPTIMUM_BRUTE_FORCE=PASS\n")
cat("RLR_GLVT_V5_MAX_TWO_REAL_PAIRS_MONOTONE_INVARIANCE=PASS\n")
cat("RLR_GLVT_V5_TEACHER_FREE_TEST_INFERENCE=PASS\n")
cat("RLR_GLVT_V5_UNIT_TESTS=PASS\n")
