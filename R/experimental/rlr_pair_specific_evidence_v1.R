# RLR within-sample rank comparison encodings.
# Paper baseline stays byte-identical: encode_pair_matrix(train, pairs, m=10).
# Only the ternary encoding step is modified; pairs/alpha/CV stay fixed.
# Pair-specific thresholds and scales are estimated from outer TRAIN X only.
rlr_margin_need<-function(cond,msg) if(!isTRUE(cond))
  stop("RLR_MARGIN_STOP: ",msg,call.=FALSE)

rlr_pair_gaps<-function(X,pairs) {
  rlr_margin_need(is.data.frame(X)&&nrow(X)>0L&&nrow(pairs)>0L,
                  "empty input or candidate pairs")
  rk<-t(apply(X,1L,rank,ties.method="average"))
  colnames(rk)<-names(X)
  D<-vapply(seq_len(nrow(pairs)),function(j) {
    a<-pairs$gene_left[[j]]
    b<-pairs$gene_right[[j]]
    rlr_margin_need(a %in% names(X) && b %in% names(X),
                    "gene name not found in input matrix")
    as.numeric(rk[,b]-rk[,a])
  },numeric(nrow(X)))
  if(is.null(dim(D)))D<-matrix(D,ncol=1L)
  colnames(D)<-pairs$ordered_pair
  rownames(D)<-rownames(X)
  rlr_margin_need(all(is.finite(D)),"non-finite rank gap")
  D
}

rlr_margin_fit<-function(X,pairs,mode) {
  rlr_margin_need(mode %in% c("TERNARY_M10","ADAPTIVE_Q25","GRADED_M10_Q75"),
                  "unknown frozen mode")
  delta<-rlr_pair_gaps(X,pairs)
  k<-ncol(delta)
  threshold<-rep(10,k)
  strength_scale<-rep(1,k)
  if(mode=="ADAPTIVE_Q25"){
    threshold<-vapply(seq_len(k),function(j){
      max(1,round(as.numeric(stats::quantile(abs(delta[,j]),probs=.25,
                     names=FALSE,type=7))))
    },numeric(1))
  }
  if(mode=="GRADED_M10_Q75"){
    strength_scale<-vapply(seq_len(k),function(j){
      max(1,as.numeric(stats::quantile(pmax(abs(delta[,j])-10+1,0),
                                     probs=.75,names=FALSE,type=7)))
    },numeric(1))
  }
  list(mode=mode,pairs=pairs,gene_schema=names(X),
       margin=threshold,scale=strength_scale)
}

rlr_margin_transform<-function(model,X) {
  rlr_margin_need(identical(names(X),model$gene_schema),
                  "features and order changed between train and test")
  D<-rlr_pair_gaps(X,model$pairs)
  delta<-abs(D)
  threshold<-matrix(rep(model$margin,each=nrow(D)),nrow=nrow(D))
  if(model$mode=="GRADED_M10_Q75"){
    denominator<-matrix(rep(model$scale,each=nrow(D)),nrow=nrow(D))
    magnitude<-pmax(0,pmin(1,(delta-threshold+1)/denominator))
  }else{
    magnitude<-(delta>=threshold)*1
  }
  z<-sign(D)*magnitude
  storage.mode(z)<-"double"
  rlr_margin_need(all(z<=1&z>=-1)&all(is.finite(z)),
                  "encoding exceeds -1/+1")
  if(model$mode=="TERNARY_M10"){
    exact<-encode_pair_matrix(X,model$pairs,m=10L)
    rlr_margin_need(isTRUE(all.equal(unname(z),unname(exact),
                       check.attributes=FALSE,tolerance=0)),
                    "frozen original RLR pair encoding mismatch")
  }
  colnames(z)<-model$pairs$ordered_pair
  rownames(z)<-rownames(X)
  z
}

rlr_margin_one<-function(T,y,E,pool,K,mode,seed=42L) {
  pairs<-apply_disjointness(pool,K=K,d=2L)
  rlr_margin_need(nrow(pairs)>0L,"empty gene-disjoint shortlist")
  model<-rlr_margin_fit(T,pairs,mode)
  RT<-rlr_margin_transform(model,T)
  RE<-rlr_margin_transform(model,E)
  set.seed(seed)
  cv<-glmnet::cv.glmnet(
      x=RT,y=as.numeric(y==levels(y)[2L]),family="binomial",
      alpha=.5,nfolds=safe_nfolds(y,5),type.measure="deviance")
  p<-as.numeric(stats::predict(cv,newx=RE,s="lambda.min",
                               type="response"))
  beta<-as.numeric(stats::coef(cv,s="lambda.min"))[-1L]
  active<-which(beta!=0)
  rlr_margin_need(all(is.finite(p))&&length(p)==nrow(E)&&
                    all(p>=0&p<=1),"prediction invalid")
  list(prob=p,mode=mode,K=K,real_K=nrow(pairs),
       globally_active=length(active),
       mean_local_active=mean(rowSums(RE[,active,drop=FALSE]!=0)),
       margin_mean=mean(model$margin),
       margin_min=min(model$margin),margin_max=max(model$margin),
       train_active_mean=mean(RT!=0),test_active_mean=mean(RE!=0),
       lambda_min=cv$lambda.min,
       pairs=pairs,params=model)
}
