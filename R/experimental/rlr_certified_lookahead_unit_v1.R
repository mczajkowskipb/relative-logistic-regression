# Synthetic exhaustive unit test for *mathematical* certification only.
# The scientific outcomes are evaluated EXCLUSIVELY on original paper cohorts.
source("R/experimental/rlr_certified_lookahead_tree_v1.R")
states<-expand.grid(rep(list(c(-1,0,1)),4L))
R<-as.matrix(states)
storage.mode(R)<-"numeric"
b<-c(1.3,-.8,.45,-.3)
intercept<--.6
cases<-list(
  list(type="static"),
  list(type="lookahead",root=3L,
       second=c("-1"=1L,"0"=2L,"1"=4L)),
  list(type="lookahead",root=2L,
       second=c("-1"=3L,"0"=1L,"1"=4L)))
full<-plogis(intercept+as.numeric(R%*%b))
for(policy in cases) {
 v<-cl2_eval(R,b,intercept,policy,epsilon=.05)
 stopifnot(all(abs(v$prob-full)<=.0250001),
           all(full>=v$lo-1e-8 & full<=v$hi+1e-8),
           all(v$hi-v$lo<=.0500001),
           identical(as.integer(v$prob>=.5),as.integer(full>=.5)),
           all(v$evaluated<=4L),
           all(v$active_n<=v$evaluated))
}
g<-cl2_learn(R[1:54,,drop=FALSE],R[55:81,,drop=FALSE],
             b,intercept,epsilon=.05,max_roots=4L,
             min_v_gain=.25)
stopifnot(g$chosen$type %in% c("static","lookahead"),
          is.finite(g$V_gain),
          all(vapply(cases,function(policy)
            all(cl2_eval(R,b,intercept,policy,epsilon=.05)$prob>=0),
            logical(1))))
cat("CL2_ALL_81_TERNARY_VECTORS_FULL_CLASS_AND_POSTERIOR_CERTIFIED=PASS\n")
cat("CL2_TRAIN_ONLY_OPTIMIZER_UNIT_SMOKE=PASS\n")
