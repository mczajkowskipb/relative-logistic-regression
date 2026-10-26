# C-RET-L2: globally optimized depth-two conditional ternary relational tree.
# Full original RLR coefficients remain unchanged. The first pair test is
# chosen globally; observed -1/0/+1 routes to a different second test.
# Every branch subsequently uses the coefficient-ordered suffix. Optimize
# FULL patient path certification cost on the outer training-only D part.
# V gates acceptance; E never selects features, structure or hyperparameters.
# Exact original RLR predicted classes and <=0.025 probability error
# are guaranteed for EVERY input by original ternary beta remainder bounds.
cl2_need<-function(ok,msg)if(!isTRUE(ok))stop("CL2_STOP: ",msg,call.=FALSE)
cl2_order<-function(beta)order(-abs(beta),seq_along(beta))
cl2_eval<-function(R,beta,intercept,policy,epsilon=.05) {
  cl2_need(is.matrix(R)&&ncol(R)==length(beta)&&
    all(is.finite(beta))&&is.finite(intercept),"predictor dimensions")
  active<-which(abs(beta)>1e-12)
  q<-length(active)
  n<-nrow(R)
  full<-plogis(intercept+as.numeric(R%*%beta))
  if(q==0L)return(list(prob=full,lo=full,hi=full,
    evaluated=rep(0L,n),active_n=rep(0L,n),
    stopped=rep(TRUE,n),branch=rep("INTERCEPT",n),
    second=rep(NA_integer_,n),paths=rep("",n),global_n=q))
  fallback_order<-cl2_order(beta)
  fallback_order<-fallback_order[abs(beta[fallback_order])>1e-12]
  if(policy$type=="lookahead") {
    root<-as.integer(policy$root)
    second<-as.integer(policy$second)
    names(second)<-c("-1","0","1")
    cl2_need(root %in% active && length(second)==3L &&
         all(second %in% active)&&
         all(second!=root),"tree root/second malformed")
  }
  pr<-numeric(n);lo<-numeric(n);hi<-numeric(n)
  ev<-integer(n);act<-integer(n);route<-character(n)
  second_observed<-rep(NA_integer_,n)
  paths<-character(n);stopped<-logical(n)
  for(i in seq_len(n)){
    if(policy$type=="static"){
      order_i<-fallback_order
      route[[i]]<-"STATIC"
    } else{
      v<-as.integer(R[i,root])
      route[[i]]<-as.character(v)
      j<-second[[route[[i]]]]
      second_observed[[i]]<-j
      order_i<-c(root,j,setdiff(fallback_order,c(root,j)))
    }
    remain<-sum(abs(beta[active]))
    partial<-intercept
    for(t in seq_along(order_i)){
      j<-order_i[[t]]
      partial<-partial+beta[[j]]*R[i,j]
      remain<-max(0,remain-abs(beta[[j]]))
      ev[[i]]<-t
      act[[i]]<-act[[i]]+as.integer(R[i,j]!=0)
      l<-plogis(partial-remain)
      u<-plogis(partial+remain)
      # Early stop after any prefix only when full posterior interval
      # is one-sided about 0.5 AND has maximum width <=epsilon.
      valid<-((l>.5)||(u<.5))&&(u-l<=epsilon+1e-12)
      if(valid || t==q){
        lo[[i]]<-l;hi[[i]]<-u;pr[[i]]<-(l+u)/2
        stopped[[i]]<-valid && t<q
        paths[[i]]<-paste(order_i[seq_len(t)],collapse=":")
        break
      }
    }
  }
  cl2_need(all(full>=lo-1e-8&full<=hi+1e-8),
           "original full RLR probability not certified")
  cl2_need(all(abs(pr-full)<=epsilon/2+1e-7),
           "probability error exceeds certificate width")
  cl2_need(identical(as.integer(pr>=.5),as.integer(full>=.5)),
           "exact original class not preserved")
  cl2_need(all(ev<=q & act<=ev),"invalid relational path length")
  list(prob=pr,lo=lo,hi=hi,evaluated=ev,active_n=act,
       stopped=stopped,branch=route,second=second_observed,
       paths=paths,global_n=q)
}
cl2_learn<-function(RD,RV,beta,intercept,epsilon=.05,
                     max_roots=10L,min_v_gain=.25) {
  cl2_need(is.matrix(RD)&&is.matrix(RV)&&ncol(RD)==ncol(RV),
           "D/V matrix mismatch")
  active<-which(abs(beta)>1e-12)
  static<-list(type="static")
  Dstatic<-cl2_eval(RD,beta,intercept,static,epsilon)
  Vstatic<-cl2_eval(RV,beta,intercept,static,epsilon)
  no_tree<-list(type="static",root=NA_integer_,second=c(-1L,-1L,-1L))
  if(length(active)<3L)return(list(
    chosen=no_tree,best_train=no_tree,train_gain=0,
    V_gain=0,V_baseline=mean(Vstatic$evaluated),
    V_best=mean(Vstatic$evaluated),roots_considered=0L,
    n_branches_modified=0L))
  ranked<-cl2_order(beta)
  ranked<-ranked[abs(beta[ranked])>1e-12]
  # End-to-end restricted policy class: one root + three potentially
  # DIFFERENT second pair tests, then globally fixed |beta|-ordered suffix.
  # Exhaustive search over every second pair for each observed root state,
  # nested in a limited set of high-magnitude candidate root pairs.
  roots<-head(ranked,max_roots)
  best<-list(cost=mean(Dstatic$evaluated),policy=no_tree,modified=0L)
  for(root in roots){
    branch_choice<-integer(3)
    names(branch_choice)<-c("-1","0","1")
    default_second<-head(setdiff(ranked,root),1L)
    for(s in c(-1L,0L,1L)){
      idx<-which(RD[,root]==s)
      if(!length(idx)){
        branch_choice[[as.character(s)]]<-default_second
        next
      }
      # Root+second candidate pair orders evaluated against full
      # posterior certificate on training samples in this branch.
      branch_scores<-vapply(setdiff(ranked,root),function(candidate){
        proposal<-list(type="lookahead",root=root,
           second=setNames(rep(candidate,3L),c("-1","0","1")))
        value<-cl2_eval(RD[idx,,drop=FALSE],beta,intercept,
                        proposal,epsilon)
        mean(value$evaluated)
      },numeric(1))
      options<-setdiff(ranked,root)
      # Mild structural penalty discourages low-support routes picking a
      # new feature for an optimistic one-sample improvement.
      branch_scores<-branch_scores+
        .10*as.numeric(options!=default_second)/
           max(1L,length(idx))
      ix<-which.min(branch_scores)
      branch_choice[[as.character(s)]]<-options[[ix]]
    }
    tree<-list(type="lookahead",root=root,second=branch_choice)
    pred<-cl2_eval(RD,beta,intercept,tree,epsilon)
    n_modified<-sum(branch_choice!=default_second)
    # Penalize size/branching: fewer true alternatives is preferred unless
    # there are real expected RELATIONAL TEST reductions on D.
    objective<-mean(pred$evaluated)+.05*n_modified
    if(objective<best$cost-1e-9)
      best<-list(cost=objective,policy=tree,modified=n_modified)
  }
  valid<-cl2_eval(RV,beta,intercept,best$policy,epsilon)
  Vgain<-mean(Vstatic$evaluated)-mean(valid$evaluated)
  # A genuine conditional tree must select DIFFERENT second tests after
  # at least two observed ternary root outcomes. A reordered fixed sequence
  # is not counted as a tree-specific scientific improvement.
  state_diverse<-FALSE
  if(best$policy$type=="lookahead"){
    seen<-sort(unique(as.integer(RD[,best$policy$root])))
    state_diverse<-length(seen)>=2L &&
      length(unique(best$policy$second[as.character(seen)]))>=2L
  }
  accept<-best$policy$type=="lookahead" && state_diverse &&
          Vgain>=min_v_gain
  chosen<-if(accept)best$policy else no_tree
  # Note: full original RLR was trained on T (which contains V).
  # This is a *policy validation* gate only, not independent revalidation
  # of the original RLR predictor or an external-data claim.
  list(chosen=chosen,best_train=best$policy,
       train_gain=mean(Dstatic$evaluated)-best$cost,
       V_gain=Vgain,
       V_baseline=mean(Vstatic$evaluated),
       V_best=mean(valid$evaluated),
       roots_considered=length(roots),
       n_branches_modified=best$modified)
}
