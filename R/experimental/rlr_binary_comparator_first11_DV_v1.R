# Controlled exploratory pilot: original RLR20 vs PurIST-/DeCAF-inspired binary pair refits.
# Same D-only switchBox top1000, d2 disjoint pairs, K20, same frozen D/V split.
# Paper's original outer E is never read/evaluated. No original ARFF SHA parity.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_rank_gate_core_v1.R")
need<-function(v,s)if(!isTRUE(v))stop(s)
cohorts<-c("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
"chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
"shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
ds<-Sys.getenv("RLR_DATASET_ID",unset="")
need(ds%in%cohorts,"bad dataset")
input<-file.path("results/rlr_full11_same_pairs_v1",ds,"source_data",paste0(ds,".csv"))
out<-file.path("results/rlr_binary_comparator11_DV_v1",ds)
need(file.exists(input)&&file.exists(file.path(out,"SOURCE_SHA.txt"))&&!file.exists(file.path(out,"metrics.csv")),"immutable input/output check")
raw<-read.csv(input,check.names=FALSE,stringsAsFactors=FALSE)
X<-raw[,setdiff(names(raw),c("sample_id","class")),drop=FALSE]
rownames(X)<-raw$sample_id;y<-factor(raw$class)
need(nlevels(y)==2L&&all(vapply(X,function(z)is.numeric(z)&&all(is.finite(z)),logical(1))),"bad schema")
reps<-if(ds%in%tail(cohorts,2))2L else 10L
sp<-make_repeated_stratified_folds(y,5L,reps,42L)[[1L]]
T<-X[sp$train_idx,,drop=FALSE]
yt<-factor(as.character(y[sp$train_idx]),levels=levels(y))
iv<-rg_inner(yt)
D<-T[iv$D,,drop=FALSE];V<-T[iv$V,,drop=FALSE]
yd<-factor(as.character(yt[iv$D]),levels=levels(y))
yv<-factor(as.character(yt[iv$V]),levels=levels(y))
truth<-as.integer(yv==levels(yv)[2L])
need(min(table(yd))>=3L&&min(table(yv))>=1L&&nrow(D)+nrow(V)==nrow(T),"D V split invalid")
pool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
pairs<-apply_disjointness(pool,K=20L,d=2L)
need(nrow(pairs)>=2L,"no shortlist")
original<-rg_parent(D,yd,V,pool,K=20L)
need(identical(as.character(original$selected$ordered_pair),as.character(pairs$ordered_pair)),"pair mismatch")
encode<-function(A){
 z<-vapply(seq_len(nrow(pairs)),function(j)as.numeric(A[[pairs$gene_left[j]]]>A[[pairs$gene_right[j]]]),numeric(nrow(A)))
 if(is.null(dim(z)))z<-matrix(z,ncol=1L)
 z
}
bd<-encode(D);bv<-encode(V)
keep<-which(apply(bd,2L,function(z)length(unique(z))>1L))
need(length(keep)>=2L,"all binary tests constant")
bd<-bd[,keep,drop=FALSE];bv<-bv[,keep,drop=FALSE]
yD<-as.numeric(yd==levels(yd)[2L]);nf<-safe_nfolds(yd,5L)
set.seed(42L)
en<-glmnet::cv.glmnet(x=bd,y=yD,family="binomial",alpha=.5,nfolds=nf,type.measure="deviance")
pe<-as.numeric(predict(en,newx=bv,s="lambda.min",type="response"))
be<-as.numeric(stats::coef(en,s="lambda.min"))[-1L]
set.seed(42L)
mc<-ncvreg::cv.ncvreg(X=bd,y=yD,family="binomial",penalty="MCP",alpha=.05,nfolds=nf)
pm<-as.numeric(predict(mc,X=bv,type="response"))
bm<-as.numeric(stats::coef(mc))[-1L]
pr<-as.numeric(original$test)
br<-as.numeric(stats::coef(original$fit,s="lambda.min"))[-1L]
need(all(vapply(list(pe,pm,pr),function(z)length(z)==nrow(V)&&all(is.finite(z))&&all(z>=0&z<=1),logical(1))),"invalid prediction")
auc<-function(t,p){n1<-sum(t==1);n0<-sum(t==0);(sum(rank(p)[t==1])-n1*(n1+1)/2)/(n1*n0)}
f1<-function(t,h,c){tp<-sum(t==c&h==c);fp<-sum(t!=c&h==c);fn<-sum(t==c&h!=c);if(2*tp+fp+fn==0)0 else 2*tp/(2*tp+fp+fn)}
ece<-function(t,p){bins<-pmin(10,as.integer(floor(p*10))+1L);sum(vapply(1:10,function(b){ix<-which(bins==b);if(!length(ix))0 else length(ix)*abs(mean(p[ix])-mean(t[ix]))/length(t)},numeric(1)))}
arms<-list(RLR20_ORIGINAL=list(p=pr,b=br),BINARY_ENET20_PURIST_INSPIRED=list(p=pe,b=be),BINARY_MCP20_DECAF_INSPIRED=list(p=pm,b=bm))
metrics<-list();patients<-list()
for(i in seq_along(arms)){
 name<-names(arms)[i];z<-arms[[i]]
 p<-pmin(1-1e-15,pmax(1e-15,z$p));h<-as.integer(p>=.5)
 metrics[[i]]<-data.frame(dataset=ds,method=name,n_D=nrow(D),n_V=nrow(V),K_target=20,
 pairs_selected=nrow(pairs),active_coefficients=sum(abs(z$b)>1e-9),
 positive_f1=f1(truth,h,1L),auc=auc(truth,p),
 logloss=mean(-truth*log(p)-(1-truth)*log1p(-p)),
 brier=mean((truth-p)^2),ece=ece(truth,p),stringsAsFactors=FALSE)
 patients[[i]]<-data.frame(dataset=ds,method=name,sample_id=rownames(V),truth=truth,prob=z$p,pred=h)
 cat(sprintf("RLR_BINARY_COMPARE %s %s logloss=%.6f brier=%.6f f1=%.6f auc=%.6f\n",
 ds,name,metrics[[i]]$logloss,metrics[[i]]$brier,metrics[[i]]$positive_f1,metrics[[i]]$auc))
}
write.csv(do.call(rbind,metrics),file.path(out,"metrics.csv"),row.names=FALSE)
write.csv(do.call(rbind,patients),file.path(out,"patients.csv"),row.names=FALSE)
write.csv(data.frame(dataset=ds,n_D=nrow(D),n_V=nrow(V),n_outer_E_UNREAD=length(sp$test_idx),
 D_ids=paste(rownames(D),collapse="|"),V_ids=paste(rownames(V),collapse="|"),
 source="SCHLIEP_RECONSTRUCTED_ARFF_SHA_NOT_VERIFIED",pairs_shared=TRUE,n_pairs=nrow(pairs),
 binary_constant_removed=nrow(pairs)-length(keep),cv_folds=nf),file.path(out,"manifest.csv"),row.names=FALSE)
cat("RLR_BINARY_THREE_METHODS_VALID_DV=PASS\n")
cat("RLR_BINARY_OUTER_E_UNTOUCHED=PASS\n")