# RXA-MultiTest margin-deferral cascade: 2 nodes x 3 disjoint gene-pair tests.
# Exploratory v1 locked protocol: docs/RLR_RXA_MULTITEST_DEFER_V1_FROZEN_PROTOCOL.md
# One fixed set of six selected gene pairs shared across all multi-test arms.
# Public Schliep RECONSTRUCTED data; original processed ARFF SHA NOT VERIFIED here.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_rank_gate_core_v1.R")
stop_if<-function(ok,msg)if(!isTRUE(ok))stop("MTDEF_STOP: ",msg,call.=FALSE)
names11<-c("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
ds<-Sys.getenv("RLR_DATASET_ID","")
stop_if(ds%in%names11,"must specify a frozen original-named dataset")
src<-file.path("results/rlr_full11_same_pairs_v1",ds,"source_data",
               paste0(ds,".csv"))
out<-file.path("results/rlr_rxa_multitest_deferral_v1",ds)
stop_if(file.exists(src)&&file.exists(file.path(out,"SOURCE_SHA.txt"))&&
  !file.exists(file.path(out,"fold_metrics.csv")),"missing source/manifest or already completed")
tab<-read.csv(src,check.names=FALSE,stringsAsFactors=FALSE)
stop_if(all(c("sample_id","class")%in%names(tab)),"source schema mismatch")
X<-as.data.frame(tab[,setdiff(names(tab),c("sample_id","class")),drop=FALSE],check.names=FALSE)
rownames(X)<-tab$sample_id
y<-factor(tab$class)
stop_if(nlevels(y)==2L&&!anyDuplicated(rownames(X))&&ncol(X)>10L&&
  all(vapply(X,function(v)is.numeric(v)&&all(is.finite(v)),logical(1))),
  "invalid feature matrix")
splits<-make_repeated_stratified_folds(y,5L,1L,42L)
stop_if(length(splits)==5L && identical(sort(unlist(lapply(splits,function(z)z$test_idx))),
  seq_len(nrow(X))),"invalid disjoint outer split")
methods<-c("ORIGINAL_RLR20","BINARY_MT2_LEAF","TERNARY_MT2_LEAF",
           "BINARY_MT2_RLR_BACKSTOP","TERNARY_MT2_RLR_BACKSTOP",
           "TERNARY_SAME6_ENET")
smooth_weight<-10;min_leaf_n<-5L;vote_threshold<-2L;rank_margin<-10L
stop_if(length(unique(methods))==6L,"arm names duplicated")
clamp<-function(p)pmin(1-1e-15,pmax(1e-15,p))
metric<-function(y,p) {
 q<-clamp(p);pr<-as.integer(p>=.5)
 tp<-sum(y==1L&pr==1L);fp<-sum(y==0L&pr==1L)
 fn<-sum(y==1L&pr==0L)
 n1<-sum(y==1L);n0<-sum(y==0L)
 stop_if(n1>0L&&n0>0L,"missing external test class")
 auc<-(sum(rank(p,ties.method="average")[y==1L])-
       n1*(n1+1)/2)/(n1*n0)
 bins<-pmin(9L,as.integer(floor(p*10)))
 ece<-sum(vapply(0:9,function(b){
  idx<-which(bins==b)
  if(!length(idx))return(0)
  length(idx)/length(y)*abs(mean(y[idx])-mean(p[idx]))
 },numeric(1)))
 c(logloss=mean(-y*log(q)-(1-y)*log1p(-q)),brier=mean((y-p)^2),
   auc=auc,positive_f1=if(2*tp+fp+fn>0)2*tp/(2*tp+fp+fn) else 0,
   ece=ece)
}
# Exactly 3 oriented pair votes per node. Two compatible votes and no
# contradictory vote => decisive. Zero values abstain, conflict defers.
vote3<-function(Z,indices) {
 q<-Z[,indices,drop=FALSE]
 stop_if(ncol(q)==3L&&all(q%in%c(-1,0,1)),"invalid pair states")
 plus<-rowSums(q>0);minus<-rowSums(q<0)
 out<-rep(0L,nrow(Z))
 out[plus>=vote_threshold & minus==0L]<-1L
 out[minus>=vote_threshold & plus==0L]<- -1L
 out
}
route_pair_states<-function(raw_train,raw_test,ytrain) {
 stop_if(ncol(raw_train)==6L&&ncol(raw_test)==6L,
   "expected exactly six independently selected pairs")
 yy<-as.integer(ytrain==levels(ytrain)[2L])
 gap_mean<-colMeans(raw_train[yy==1,,drop=FALSE])-
           colMeans(raw_train[yy==0,,drop=FALSE])
 orientation<-ifelse(gap_mean<0,-1,1)
 stop_if(all(orientation%in%c(-1,1)),"orientation invalid")
 list(T=sweep(raw_train,2L,orientation,"*"),
      E=sweep(raw_test,2L,orientation,"*"),orientation=orientation)
}
fit_terminal_probs<-function(states,ytrain,base_prior){
 stop_if(length(states)==length(ytrain),"train route length")
 values<-c("-1","1")
 fitted<-setNames(rep(base_prior,2L),values)
 counts<-setNames(integer(2L),values)
 for(k in values) {
  ix<-which(as.character(states)==k)
  counts[k]<-length(ix)
  if(length(ix)>=min_leaf_n){
   fitted[k]<-(sum(ytrain[ix])+smooth_weight*base_prior)/
     (length(ix)+smooth_weight)
  }
 }
 list(prob=fitted,count=counts)
}
# D/F routing is identical for both terminal-prediction modes.
train_routing<-function(state_train,y01) {
 a<-vote3(state_train,1:3)
 b<-vote3(state_train,4:6)
 final<-ifelse(a!=0L,paste0("n1_",a),
               ifelse(b!=0L,paste0("n2_",b),"unresolved"))
 prior<-(sum(y01)+1)/(length(y01)+2)
 prob<-setNames(numeric(4L),c("n1_-1","n1_1","n2_-1","n2_1"))
 count<-setNames(integer(4L),names(prob))
 for(k in names(prob)){
  ii<-which(final==k)
  count[k]<-length(ii)
  prob[k]<-if(length(ii)>=min_leaf_n)
    (sum(y01[ii])+smooth_weight*prior)/(length(ii)+smooth_weight) else prior
 }
 list(prob=prob,count=count,prior=prior)
}
predict_cascade<-function(train_state,test_state,y01,global_prob=NULL,n_global=20L){
 fit<-train_routing(train_state,y01)
 a<-vote3(test_state,1:3);b<-vote3(test_state,4:6)
 stage<-ifelse(a!=0L,"NODE1",ifelse(b!=0L,"NODE2","DEFER"))
 terminal<-ifelse(a!=0L,paste0("n1_",a),
                 ifelse(b!=0L,paste0("n2_",b),"unresolved"))
 prob<-rep(fit$prior,length(a));known<-terminal!="unresolved"
 prob[known]<-fit$prob[terminal[known]]
 if(!is.null(global_prob)) {
  stop_if(length(global_prob)==length(a),"backstop length")
  prob[terminal=="unresolved"]<-global_prob[terminal=="unresolved"]
 }
 cost<-ifelse(stage=="NODE1",3L,6L)
 if(!is.null(global_prob))cost[terminal=="unresolved"]<-6L+n_global
 stop_if(all(is.finite(prob))&&all(prob>0&prob<1),"invalid probability")
 list(p=as.numeric(prob),route=stage,checks=cost,
      fallback=if(is.null(global_prob))rep(0L,length(a)) else as.integer(terminal=="unresolved"),
      leaf_count=fit$count,train_unresolved=sum(vote3(train_state,1:3)==0L&
                                            vote3(train_state,4:6)==0L))
}
metric_rows<-list();patient_rows<-list();manifest<-list()
for(fi in seq_along(splits)){
 split<-splits[[fi]]
 D<-X[split$train_idx,,drop=FALSE];E<-X[split$test_idx,,drop=FALSE]
 yd<-factor(as.character(y[split$train_idx]),levels=levels(y))
 ye<-factor(as.character(y[split$test_idx]),levels=levels(y))
 y01<-as.integer(yd==levels(y)[2L])
 truth<-as.integer(ye==levels(y)[2L])
 stop_if(min(table(yd))>=4L&&min(table(ye))>=1L,"class support failure")
 pool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
 parent<-rg_parent(D,yd,E,pool,K=20L)
 pairs<-parent$selected
 stop_if(nrow(pairs)>=6L,"insufficient gene-disjoint relational tests")
 chosen<-pairs[1:6,,drop=FALSE]
 # Orientation determined using binary order on all TRAIN samples only,
 # and is shared across binary/ternary arms. Ternary m10 states are original.
 binary_T<-encode_pair_matrix(D,chosen,m=1L)
 binary_E<-encode_pair_matrix(E,chosen,m=1L)
 ternary_T<-encode_pair_matrix(D,chosen,m=rank_margin)
 ternary_E<-encode_pair_matrix(E,chosen,m=rank_margin)
 O<-route_pair_states(binary_T,binary_E,yd)
 btrain<-O$T;btest<-O$E
 ttrain<-sweep(ternary_T,2L,O$orientation,"*")
 ttest<-sweep(ternary_E,2L,O$orientation,"*")
 nparent<-nrow(pairs)
 pb<-as.numeric(parent$test)
 b_leaf<-predict_cascade(btrain,btest,y01,n_global=nparent)
 t_leaf<-predict_cascade(ttrain,ttest,y01,n_global=nparent)
 b_log<-predict_cascade(btrain,btest,y01,global_prob=pb,n_global=nparent)
 t_log<-predict_cascade(ttrain,ttest,y01,global_prob=pb,n_global=nparent)
 # A pair-budget-matched, no-tree penalized logistic control. No feature
 # search independent of the six fixed original selected pairs.
 cols<-which(apply(ternary_T,2L,function(v)length(unique(v))>1L))
 if(length(cols)>=2L){
  set.seed(42L)
  reg<-glmnet::cv.glmnet(
    x=ternary_T[,cols,drop=FALSE],y=y01,
    family="binomial",alpha=.5,nfolds=safe_nfolds(yd,5L),
    type.measure="deviance")
  p6<-as.numeric(predict(reg,newx=ternary_E[,cols,drop=FALSE],
                        type="response",s="lambda.min"))
 } else p6<-rep((sum(y01)+1)/(length(y01)+2),nrow(E))
 cases<-list(
  ORIGINAL_RLR20=list(p=pb,route=rep("RLR",nrow(E)),
       checks=rep(nparent,nrow(E)),fallback=rep(0L,nrow(E))),
  BINARY_MT2_LEAF=b_leaf,
  TERNARY_MT2_LEAF=t_leaf,
  BINARY_MT2_RLR_BACKSTOP=b_log,
  TERNARY_MT2_RLR_BACKSTOP=t_log,
  TERNARY_SAME6_ENET=list(p=p6,route=rep("ENET",nrow(E)),
       checks=rep(6L,nrow(E)),fallback=rep(0L,nrow(E)))
 )
 stop_if(identical(names(cases),methods),"arm mismatch")
 manifest[[fi]]<-data.frame(dataset=ds,fold=fi-1L,
   split_id=split$split_id,train_ids=paste(rownames(D),collapse="|"),
   test_ids=paste(rownames(E),collapse="|"),
   n_train=nrow(D),n_test=nrow(E),n_parent_pairs=nparent,
   pair_6=paste(chosen$ordered_pair,collapse="|"),
   orientation=paste(O$orientation,collapse="|"),
   source="SCHLIEP_RECONSTRUCTED_NOT_ORIGINAL_ARFF_SHA_VERIFIED",
   strict_gene_disjoint_d2=TRUE,rank_margin=rank_margin,node_pair_count=3L,
   fixed_prior_weight=smooth_weight,min_leaf_support=min_leaf_n,
   stringsAsFactors=FALSE)
 for(key in names(cases)){
  obj<-cases[[key]]
  stop_if(length(obj$p)==nrow(E)&&length(obj$checks)==nrow(E)&&
    length(obj$route)==nrow(E)&&length(obj$fallback)==nrow(E)&&
    all(is.finite(obj$p))&&all(obj$p>=0&obj$p<=1),
    paste("bad model output",key))
  score<-metric(truth,obj$p)
  metric_rows[[length(metric_rows)+1L]]<-data.frame(
   dataset=ds,fold=fi-1L,split_id=split$split_id,method=key,
   n_test=nrow(E),n_train=nrow(D),n_pairs_parent=nparent,
   logloss=score[["logloss"]],brier=score[["brier"]],
   auc=score[["auc"]],positive_f1=score[["positive_f1"]],
   ece=score[["ece"]],mean_pair_checks=mean(obj$checks),
   mean_backstop_rate=mean(obj$fallback),
   node1_exit_rate=mean(obj$route=="NODE1"),
   node2_exit_rate=mean(obj$route=="NODE2"),
   unresolved_rate=mean(obj$route=="DEFER"),
   stringsAsFactors=FALSE)
  patient_rows[[length(patient_rows)+1L]]<-data.frame(
   dataset=ds,fold=fi-1L,split_id=split$split_id,method=key,
   sample_id=rownames(E),truth=truth,prob=obj$p,
   pred=as.integer(obj$p>=.5),exit_route=obj$route,
   pair_checks=obj$checks,backstop=obj$fallback,
   stringsAsFactors=FALSE)
  cat(sprintf("RXA_MTDEF dataset=%s fold=%d method=%s LL=%.6f F1=%.6f pair_checks=%.3f backstop=%.3f\n",
   ds,fi-1L,key,score[["logloss"]],score[["positive_f1"]],
   mean(obj$checks),mean(obj$fallback)))
 }
}
m<-do.call(rbind,metric_rows);p<-do.call(rbind,patient_rows)
f<-do.call(rbind,manifest)
stop_if(nrow(m)==30L&&nrow(f)==5L&&nrow(p)==6L*nrow(X)&&
 length(unique(paste(p$method,p$sample_id)))==nrow(p),
 "incomplete dataset/fold/method coverage")
write.csv(m,file.path(out,"fold_metrics.csv"),row.names=FALSE)
write.csv(p,file.path(out,"patients.csv"),row.names=FALSE)
write.csv(f,file.path(out,"fold_manifest.csv"),row.names=FALSE)
cat("RXA_MTDEF_V1_5FOLDS_6ARMS=PASS\n")
cat("RXA_MTDEF_V1_PATIENTS_ONCE_PER_ARM=PASS\n")
