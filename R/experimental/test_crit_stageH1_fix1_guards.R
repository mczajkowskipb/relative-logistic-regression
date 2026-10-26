# Base-R preflight test, no historical datasets or patient data.
source("R/experimental/crit_stageH1_pool_guard_fix1.R")
assert_fail<-function(expr,label,expected_pattern) {
  msg<-tryCatch({force(expr);"NO_ERROR"},error=function(e)conditionMessage(e))
  if(identical(msg,"NO_ERROR") || !grepl(expected_pattern,msg,fixed=TRUE))
     stop("STAGEH1_FIX1_TEST_STOP: ",label,": ",msg,call.=FALSE)
  message("PASS ",label," -> ",msg)
}
a<-rbind(1:320,2:321)
b<-rbind(1:496,2:497)
stopifnot(identical(critH1_need_pool(a,"FIXED_RANDOM"),invisible(TRUE)),
          identical(critH1_need_pool(a,"SEED_RANDOM"),invisible(TRUE)),
          identical(critH1_need_pool(a,"AWARE_GRAPH"),invisible(TRUE)),
          identical(critH1_need_pool(b,"EXHAUSTIVE32"),invisible(TRUE)))
message("PASS 320/320 capped strategies and 496/496 exhaustive")
assert_fail(critH1_need_pool(b,"FIXED_RANDOM"),
            "over budget","pool budget exceeded")
dups<-a
dups[,320L]<-a[,1L]
assert_fail(critH1_need_pool(dups,"AWARE_GRAPH"),
            "duplicate candidate","duplicated unordered candidate")
rev<-a
rev[,1L]<-a[2:1,1L]
assert_fail(critH1_need_pool(rev,"SEED_RANDOM"),
            "noncanonical pair","noncanonical unordered")
assert_fail(critH1_need_pool(matrix(1:100,nrow=5L),"SEED_RANDOM"),
            "bad shape","matrix shape invalid")

set.seed(20261026)
n<-120L
p<-32L
g<-matrix(stats::rnorm(n*p),nrow=n)
pairs<-utils::combn(seq_len(p),2L)
sel<-sample.int(ncol(pairs),320L)
pool<-pairs[,sel,drop=FALSE]
R<-vapply(seq_len(ncol(pool)),function(j)
     ifelse(g[,pool[1L,j]]>g[,pool[2L,j]],1,-1),numeric(n))
score<-abs(crossprod(R,R*sample(c(-.5,.5),n,replace=TRUE)))
idx<-which(upper.tri(score),arr.ind=TRUE)
slow<-vapply(seq_len(nrow(idx)),function(j){
  length(intersect(pool[,idx[j,1L]],pool[,idx[j,2L]]))==0L
},logical(1))
a<-idx[,1L];b<-idx[,2L]
fast<-pool[1L,a]!=pool[1L,b] &
      pool[1L,a]!=pool[2L,b] &
      pool[2L,a]!=pool[1L,b] &
      pool[2L,a]!=pool[2L,b]
stopifnot(identical(slow,fast))
message("PASS vectorized disjoint-pair checker exactly matches original loop")
message("STAGEH1_FIX1_POOL_GUARD_AND_VECTOR_SCREEN_TESTS=PASS")
