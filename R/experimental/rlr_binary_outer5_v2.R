# Outer 5-fold (one full repeat) matched-pair controlled comparator benchmark.
# Original RLR code is SOURCED UNCHANGED. Public Schliep reconstructed source.
# The original paper-processed ARFF SHA256 parity is NOT verified.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_rank_gate_core_v1.R")
need<-function(ok,msg)if(!isTRUE(ok))stop("OUTER5_STOP: ",msg,call.=FALSE)
datasets<-c("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
  "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
  "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
ds<-Sys.getenv("RLR_DATASET_ID",unset="")
need(ds%in%datasets,"invalid dataset name")
file<-file.path("results/rlr_full11_same_pairs_v1",ds,"source_data",paste0(ds,".csv"))
out<-file.path("results/rlr_binary_outer5_v2",ds)
need(file.exists(file)&&file.exists(file.path(out,"SOURCE_SHA.txt"))&&
       !file.exists(file.path(out,"fold_metrics.csv")),"immutable run inputs")
z<-read.csv(file,check.names=FALSE,stringsAsFactors=FALSE)
need(all(c("class","sample_id")%in%names(z)),"bad sample schema")
X<-as.data.frame(z[,setdiff(names(z),c("sample_id","class")),drop=FALSE],
                 check.names=FALSE)
rownames(X)<-z$sample_id;y<-factor(z$class)
need(nlevels(y)==2L&&nrow(X)==length(y)&&!anyDuplicated(rownames(X))&&
       all(vapply(X,function(t)is.numeric(t)&&all(is.finite(t)),logical(1))),
     "corrupt source")
sp<-make_repeated_stratified_folds(y,5L,1L,42L)
need(length(sp)==5L && identical(sort(unlist(lapply(sp,function(s)s$test_idx))),
                                seq_len(nrow(X))),"outer five coverage invalid")
methods<-c("RLR20_ORIGINAL","BINARY_ENET20_PURIST_INSPIRED",
           "BINARY_MCP20_DECAF_INSPIRED")
encode_binary<-function(mat,pairs){
  v<-vapply(seq_len(nrow(pairs)),
            function(j)as.numeric(mat[[pairs$gene_left[j]]] >
                                      mat[[pairs$gene_right[j]]]),
            numeric(nrow(mat)))
  if(is.null(dim(v)))v<-matrix(v,ncol=1L)
  v
}
ll<-function(t,p){q<-pmax(1e-15,pmin(1-1e-15,p));mean(-t*log(q)-(1-t)*log1p(-q))}
f1<-function(t,h,cls){tp<-sum(t==cls&h==cls);fp<-sum(t!=cls&h==cls)
  fn<-sum(t==cls&h!=cls);if(2*tp+fp+fn==0L)0 else 2*tp/(2*tp+fp+fn)}
auc<-function(t,p){n1<-sum(t==1L);n0<-sum(t==0L)
  need(n1>0L&&n0>0L,"test fold lacks a class")
  (sum(rank(p,ties.method="average")[t==1L])-n1*(n1+1)/2)/(n1*n0)}
ece<-function(t,p){bins<-pmin(10L,as.integer(floor(p*10))+1L)
  sum(vapply(seq_len(10L),function(b){ix<-which(bins==b)
    if(!length(ix))0 else length(ix)*abs(mean(p[ix])-mean(t[ix]))/length(t)
  },numeric(1)))}
res<-list();pat<-list();mf<-list()
for(fi in seq_along(sp)){
  split<-sp[[fi]]
  D<-X[split$train_idx,,drop=FALSE];E<-X[split$test_idx,,drop=FALSE]
  yd<-factor(as.character(y[split$train_idx]),levels=levels(y))
  ye<-factor(as.character(y[split$test_idx]),levels=levels(y))
  truth<-as.integer(ye==levels(y)[2L])
  need(min(table(yd))>=4L&&min(table(ye))>=1L,"class support failure")
  pool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
  pairs<-apply_disjointness(pool,K=20L,d=2L)
  original<-rg_parent(D,yd,E,pool,K=20L)
  need(identical(as.character(original$selected$ordered_pair),
                 as.character(pairs$ordered_pair)),"RLR pair mismatch")
  bD<-encode_binary(D,pairs);bE<-encode_binary(E,pairs)
  keep<-which(apply(bD,2L,function(x)length(unique(x))>1L))
  need(length(keep)>=2L,"not enough nonconstant binary features")
  bD<-bD[,keep,drop=FALSE];bE<-bE[,keep,drop=FALSE]
  yD<-as.numeric(yd==levels(yd)[2L])
  nf<-safe_nfolds(yd,5L)
  set.seed(42L)
  en<-glmnet::cv.glmnet(x=bD,y=yD,family="binomial",alpha=.5,
                       nfolds=nf,type.measure="deviance")
  pe<-as.numeric(predict(en,newx=bE,s="lambda.min",type="response"))
  be<-as.numeric(stats::coef(en,s="lambda.min"))[-1L]
  set.seed(42L)
  mc<-ncvreg::cv.ncvreg(X=bD,y=yD,family="binomial",
                       penalty="MCP",alpha=.05,nfolds=nf)
  pm<-as.numeric(predict(mc,X=bE,type="response"))
  bm<-as.numeric(stats::coef(mc))[-1L]
  pr<-as.numeric(original$test)
  br<-as.numeric(stats::coef(original$fit,s="lambda.min"))[-1L]
  preds<-list(RLR20_ORIGINAL=list(p=pr,b=br),
              BINARY_ENET20_PURIST_INSPIRED=list(p=pe,b=be),
              BINARY_MCP20_DECAF_INSPIRED=list(p=pm,b=bm))
  need(identical(names(preds),methods),"method order changed")
  for(mi in seq_along(preds)){
    key<-names(preds)[mi];pp<-preds[[mi]]$p
    need(length(pp)==nrow(E)&&all(is.finite(pp))&&
         all(pp>=0&pp<=1),"bad prediction")
    h<-as.integer(pp>=.5)
    res[[length(res)+1L]]<-data.frame(dataset=ds,fold=fi-1L,
      split_id=split$split_id,method=key,n_train=nrow(D),n_test=nrow(E),
      pair_budget=20L,pairs_selected=nrow(pairs),
      binary_constant_removed=nrow(pairs)-length(keep),
      active_coefficients=sum(abs(preds[[mi]]$b)>1e-9),
      positive_f1=f1(truth,h,1L),
      auc=auc(truth,pp),logloss=ll(truth,pp),
      brier=mean((truth-pp)^2),ece=ece(truth,pp),stringsAsFactors=FALSE)
    pat[[length(pat)+1L]]<-data.frame(dataset=ds,fold=fi-1L,
      split_id=split$split_id,method=key,sample_id=rownames(E),
      truth=truth,prob=pp,pred=h,stringsAsFactors=FALSE)
    cat(sprintf("OUTER5 %s fold=%d method=%s n=%d LL=%.6f F1=%.6f\n",
      ds,fi-1L,key,nrow(E),tail(res,1)[[1]]$logloss,tail(res,1)[[1]]$positive_f1))
  }
  mf[[length(mf)+1L]]<-data.frame(dataset=ds,fold=fi-1L,
       split_id=split$split_id,
       n_train=nrow(D),n_test=nrow(E),
       train_sample_ids=paste(rownames(D),collapse="|"),
       test_sample_ids=paste(rownames(E),collapse="|"),
       same_pair_selection="SWITCHBOX_TOP1000_D2_K20_OUTER_T_ONLY",
       n_pairs=nrow(pairs),n_binary_nonconstant=length(keep),
       cv_folds=nf,source="SCHLIEP_RECONSTRUCTED_ORIGINAL_ARFF_SHA_NOT_VERIFIED",
       stringsAsFactors=FALSE)
}
m<-do.call(rbind,res);p<-do.call(rbind,pat);manifest<-do.call(rbind,mf)
need(nrow(m)==15L&&nrow(manifest)==5L&&nrow(p)==3L*nrow(X) &&
     length(unique(paste(p$method,p$sample_id)))==nrow(p),
     "model-fold-patient expected coverage invalid")
write.csv(m,file.path(out,"fold_metrics.csv"),row.names=FALSE)
write.csv(p,file.path(out,"patients.csv"),row.names=FALSE)
write.csv(manifest,file.path(out,"fold_manifest.csv"),row.names=FALSE)
write.csv(aggregate(cbind(positive_f1,auc,logloss,brier,ece,
    active_coefficients)~dataset+method,m,FUN=mean),
    file.path(out,"dataset_metrics.csv"),row.names=FALSE)
cat("RLR_BINARY_OUTER5_FIVE_FOLDS_THREE_METHODS=PASS\n")
cat("RLR_BINARY_OUTER5_ALL_PATIENTS_EXACTLY_ONCE_PER_METHOD=PASS\n")
