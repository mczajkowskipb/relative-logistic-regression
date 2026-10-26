# Certified RLR Early-Exit Tree (C-RET), exploratory complementary baseline.
# Original RLR K50 is fitted LITERALLY with original publication function.
# This procedure changes only order of testing its learned ternary relations.
# Since r_j(x) is -1,0,+1, partial logit and remaining abs(weights) give
# PROVABLE class and posterior bounds. No tests/recalibration on heldout E.
# Not a novel decision tree classifier until prior-art & method comparison.
cret_need<-function(ok,why)if(!isTRUE(ok))stop("CRET_STOP: ",why,call.=FALSE)
cret_evaluate<-function(R,beta,intercept,order,epsilon=.05) {
  n<-nrow(R);k<-length(beta)
  cret_need(ncol(R)==k&&length(order)==k&&
             identical(sort(order),seq_len(k)),"coefficient/order mismatch")
  # Eliminating exactly-zero regression coefficients changes no predictions.
  active<-order[abs(beta[order])>1e-12]
  if(!length(active)) {
    p<-rep(plogis(intercept),n)
    return(list(prob=p,lo=p,hi=p,certified=rep(TRUE,n),
      class=as.integer(p>=.5),evaluated=rep(0L,n),
      active_n=rep(0L,n),n_global_active=0L))
  }
  p<-numeric(n);lo<-numeric(n);hi<-numeric(n)
  num<-integer(n);active_n<-integer(n);ok<-logical(n)
  # remaining contribution bounded by sum abs(beta), even for |r|<1.
  remain<-sum(abs(beta[active]))
  partial<-rep(intercept,n)
  for(t in seq_along(active)){
    j<-active[[t]]
    remain<-max(0,remain-abs(beta[[j]]))
    pending<-which(!ok)
    if(!length(pending))break
    rv<-R[pending,j]
    partial[pending]<-partial[pending]+rv*beta[[j]]
    num[pending]<-num[pending]+1L
    active_n[pending]<-active_n[pending]+as.integer(rv!=0L)
    lower<-plogis(partial[pending]-remain)
    upper<-plogis(partial[pending]+remain)
    # If certificate margin excludes 0.5 and probability interval <=epsilon,
    # the returned probability is within epsilon/2 of original for each patient.
    c0<-upper<.5
    c1<-lower>.5
    certified<-(c0|c1)&(upper-lower<=epsilon+1e-12)
    if(any(certified)){
      id<-pending[certified]
      lo[id]<-lower[certified];hi[id]<-upper[certified]
      p[id]<-(lower[certified]+upper[certified])/2
      ok[id]<-TRUE
    }
  }
  if(any(!ok)){
    id<-which(!ok)
    # All remaining relations are evaluated, so exact probability is known.
    truth<-plogis(intercept+
      as.numeric(R[id,,drop=FALSE]%*%beta))
    lo[id]<-truth;hi[id]<-truth;p[id]<-truth
    num[id]<-length(active)
    active_n[id]<-rowSums(R[id,active,drop=FALSE]!=0)
    ok[id]<-TRUE
  }
  # Strict numerical certificate; original exact logit must be contained
  # in all probability intervals including observed E patients.
  original<-plogis(intercept+as.numeric(R%*%beta))
  cret_need(all(original>=lo-1e-8&original<=hi+1e-8),
            "probability bound certificate not valid")
  cret_need(all(abs(p-original)<=epsilon/2+1e-7),
            "approximation exceeds certificate")
  cret_need(identical(as.integer(p>=.5),as.integer(original>=.5)),
            "original model exact classification not preserved")
  list(prob=p,lo=lo,hi=hi,certified=ok,
       class=as.integer(p>=.5),evaluated=num,active_n=active_n,
       n_global_active=length(active))
}
cret_optimize<-function(R,y,beta,intercept,seed=42L,epsilon=.05) {
  # Inference objective learned strictly on outer TRAIN. Since order does not
  # modify fitted RLR, it is safe to reuse T in ordering optimization.
  active<-which(abs(beta)>1e-12)
  order1<-c(active[order(-abs(beta[active]),active)],
            which(abs(beta)<=1e-12))
  if(length(active)<3L) return(list(original=order1,ea=order1,
        train_gain=0L))
  eval_score<-function(o)mean(cret_evaluate(R,beta,intercept,o,
                                      epsilon)$evaluated)
  baseline<-eval_score(order1)
  # EA changes only order of relation tests. 32 proposals, 4 generations
  # x 8 and elitist survivor, no external test labels.
  set.seed(seed+601L)
  pop<-list(order1)
  for(i in 1:7)pop[[length(pop)+1L]]<-
     c(sample(active),setdiff(seq_along(beta),active))
  best<-order1;best_score<-baseline
  for(gen in seq_len(4L)){
    scores<-vapply(pop,eval_score,numeric(1))
    ix<-order(scores)
    if(scores[ix[1L]]<best_score){
      best<-pop[[ix[1L]]];best_score<-scores[ix[1L]]
    }
    survive<-pop[ix[seq_len(2L)]]
    children<-survive
    while(length(children)<8L){
      o<-pop[[sample(ix[seq_len(4L)],1L)]]
      positions<-sample(seq_along(active),2L)
      o[positions]<-rev(o[positions])
      children[[length(children)+1L]]<-o
    }
    pop<-children
  }
  scores<-vapply(pop,eval_score,numeric(1))
  if(min(scores)<best_score){
    best<-pop[[which.min(scores)]]
    best_score<-min(scores)
  }
  list(original=order1,ea=best,
       train_gain=baseline-best_score)
}
