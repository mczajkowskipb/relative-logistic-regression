# C-RET-Context V1 — a genuine state-dependent three-way relational test tree.
# After seeing ternary outcome {-1,0,+1} of previous pairs, choose the next
# pair using train-only conditional observed rank patterns; no refit of RLR.
# Class and probability interval guarantees are valid for ANY ordering.
cret_context_eval<-function(Rtrain,RE,beta,intercept,epsilon=.05,shrink=8) {
  stopifnot(is.matrix(Rtrain),is.matrix(RE),
            ncol(RE)==ncol(Rtrain),ncol(RE)==length(beta),
            all(is.finite(beta)),is.finite(intercept))
  n<-nrow(RE);k<-ncol(RE)
  active<-which(abs(beta)>1e-12)
  q<-length(active)
  original<-plogis(intercept+as.numeric(RE%*%beta))
  if(q==0L)return(list(prob=original,lo=original,hi=original,
       evaluated=rep(0L,n),active_n=rep(0L,n),
       paths=rep("",n),first=rep(NA_integer_,n),
       second=rep(NA_integer_,n),exact_class=TRUE))
  # Pre-compute unconditional empirical ternary probabilities (train only).
  T<-Rtrain[,active,drop=FALSE]
  prior<-rbind(colMeans(T==(-1)),colMeans(T==0),colMeans(T==1))
  p<-numeric(n);lo<-numeric(n);hi<-numeric(n)
  evaln<-integer(n);actn<-integer(n)
  path<-character(n);first<-rep(NA_integer_,n);second<-rep(NA_integer_,n)
  for(i in seq_len(n)){
    pending<-active
    bound<-sum(abs(beta[active]))
    L<-intercept
    matched<-rep(TRUE,nrow(Rtrain))
    seen<-integer()
    for(step in seq_len(q)){
      if(!length(pending))break
      scores<-numeric(length(pending))
      matchedn<-sum(matched)
      for(ii in seq_along(pending)){
        j<-pending[[ii]]
        bnew<-max(0,bound-abs(beta[[j]]))
        w<-vapply(c(-1,0,1),function(v){
          hi0<-plogis(L+beta[[j]]*v+bnew)
          lo0<-plogis(L+beta[[j]]*v-bnew)
          hi0-lo0
        },numeric(1))
        z<-match(j,active)
        # Exact conditional subtype counts shrunk to unconditional
        # frequencies, robust when n at current path is very small.
        qj<-Rtrain[matched,j]
        count<-c(sum(qj==(-1)),sum(qj==0),sum(qj==1))
        freq<-(count+shrink*prior[,z])/(matchedn+shrink)
        scores[[ii]]<-sum(freq*w)
      }
      # Always deterministic; tie is resolved by larger |beta| then index.
      ordering<-order(scores,-abs(beta[pending]),pending)
      chosen<-pending[[ordering[[1L]]]]
      if(step==1L)first[[i]]<-chosen
      if(step==2L)second[[i]]<-chosen
      r<-RE[i,chosen]
      L<-L+r*beta[[chosen]]
      bound<-max(0,bound-abs(beta[[chosen]]))
      matched<-matched & Rtrain[,chosen]==r
      seen<-c(seen,chosen)
      pending<-setdiff(pending,chosen)
      evaln[[i]]<-step
      if(r!=0)actn[[i]]<-actn[[i]]+1L
      lower<-plogis(L-bound)
      upper<-plogis(L+bound)
      good<-(lower>.5 || upper<.5) &&
            upper-lower<=epsilon+1e-12
      if(good || !length(pending)){
        lo[[i]]<-lower;hi[[i]]<-upper
        p[[i]]<-(lower+upper)/2
        break
      }
    }
    path[[i]]<-paste(seen,collapse=":")
  }
  # All symbolic certificates must contain RLR's full posterior.
  stopifnot(all(original>=lo-1e-8 & original<=hi+1e-8),
            all(abs(p-original)<=epsilon/2+1e-7),
            identical(as.integer(p>=.5),as.integer(original>=.5)),
            all(evaln<=q),all(actn<=evaln))
  list(prob=p,lo=lo,hi=hi,evaluated=evaln,active_n=actn,
       first=first,second=second,paths=path,
       exact_class=TRUE)
}
