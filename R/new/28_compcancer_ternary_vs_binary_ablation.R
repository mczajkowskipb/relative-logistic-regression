# CompCancer K20 ternary RLR versus two matched binary pair-logistic baselines.
# Requires exact original CompCancer ARFF SHA256 matching paper checksums.
# Run from repository root with RLR_ORIGINAL_DATA_ROOT pointing at data directory.
# Read-only source inputs; writes new results under results/.
source("R/new/07_world_protocol_utils.R")
# Original RLR parent helper copied byte-for-byte from experimental
# rlr_rank_gate_core_v1.R at the frozen original ablation run.
# This avoids any R/experimental dependency in the public release.
rg_need <- function(ok, message) if (!isTRUE(ok)) stop("RANK_GATE_STOP: ", message, call.=FALSE)
rg_parent <- function(D,y,E,pool,K=20L) {
  rg_need(identical(names(D),names(E)),"parent schema mismatch")
  rg_need(nrow(D)==length(y) && length(unique(as.character(y)))==2L,
          "parent training labels invalid")
  pairs <- apply_disjointness(pool,K=K,d=2L)
  rg_need(nrow(pairs)>1L,"original RLR shortlist too short")
  train <- encode_pair_matrix(D,pairs,m=10L)
  test <- encode_pair_matrix(E,pairs,m=10L)
  set.seed(42L)
  cv <- glmnet::cv.glmnet(x=train,y=as.numeric(y==levels(y)[2L]),
         family="binomial",alpha=.5,nfolds=safe_nfolds(y,5L),
         type.measure="deviance")
  list(train=as.numeric(predict(cv,newx=train,s="lambda.min",type="response")),
       test=as.numeric(predict(cv,newx=test,s="lambda.min",type="response")),
       fit=cv, selected=pairs)
}

need<-function(ok,msg)if(!isTRUE(ok))stop("OUTER5_STOP: ",msg,call.=FALSE)
datasets<-c("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
  "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
  "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
ds<-Sys.getenv("RLR_DATASET_ID",unset="")
need(ds%in%datasets,"invalid dataset name")
root<-Sys.getenv("RLR_ORIGINAL_DATA_ROOT",unset="data")
file<-file.path(root,"CompCancer",paste0(ds,".arff"))
out<-file.path("results/rlr_binary_outer5_original_arff_v2",ds)
need(file.exists(file)&&file.exists(file.path(out,"SOURCE_SHA.txt"))&&
       !file.exists(file.path(out,"fold_metrics.csv")),"immutable run inputs")
manifest<-readLines("checksums/DATA_SHA256_EXPECTED.txt",warn=FALSE)
key<-paste0("data/CompCancer/",ds,".arff")
entry<-manifest[endsWith(trimws(manifest),key)]
need(length(entry)==1L,paste("missing exactly one canonical checksum",key))
expected<-strsplit(trimws(entry),"[[:space:]]+")[[1]][1]
actual<-substr(system2("sha256sum",shQuote(file),stdout=TRUE),1L,64L)
need(identical(actual,expected),paste("ORIGINAL_ARFF_SHA256_MISMATCH",ds))
shard<-if(ds%in%c("west-2001","yeoh-2002-v1")) "CompCancer_truth_2x2" else "CompCancer_truth_9x10"
shard_file<-file.path(root,shard,paste0(ds,".arff"))
need(file.exists(shard_file),paste("missing original shard",shard_file))
shard_sha<-substr(system2("sha256sum",shQuote(shard_file),stdout=TRUE),1L,64L)
need(identical(shard_sha,expected),paste("ORIGINAL_SHARD_SHA256_MISMATCH",ds))
xy<-read_xy_file(file)
X<-xy$X;y<-xy$y
need(all(vapply(X,is.numeric,logical(1L)))&&all(is.finite(as.matrix(X))),
     "non-numeric or non-finite original feature matrix")
cat(sprintf("ORIGINAL_ARFF_SHA256_PASS dataset=%s sha=%s n=%d p=%d\n",
            ds,actual,nrow(X),ncol(X)))
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
       cv_folds=nf,source="ORIGINAL_PAPER_ARFF_SHA256_CONFIRMED",
       original_arff_sha256=actual,original_shard_sha256=shard_sha,
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
cat("RLR_BINARY_ORIGINAL_ARFF_SHA256_CONFIRMED=PASS\n")
cat("RLR_BINARY_OUTER5_FIVE_FOLDS_THREE_METHODS=PASS\n")
cat("RLR_BINARY_OUTER5_ALL_PATIENTS_EXACTLY_ONCE_PER_METHOD=PASS\n")
