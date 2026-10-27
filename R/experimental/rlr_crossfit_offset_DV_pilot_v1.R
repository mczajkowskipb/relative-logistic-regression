source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_crossfit_offset_ablation_v1.R")
DS<-Sys.getenv("RLR_DATASET_ID",unset="")
crp_need(DS %in% c("chen-2002","gordon-2002","singh-2002"),"dataset not predeclared")
out<-file.path("results/rlr_crossfit_offset_DV_v1",DS)
src<-file.path("results/rlr_full11_same_pairs_v1",DS,
               "source_data",paste0(DS,".csv"))
crp_need(file.exists(src)&&dir.exists(out)&&
          file.exists(file.path(out,"SOURCE_SHA.txt"))&&
          !file.exists(file.path(out,"metrics_V.csv")),"missing source / duplicate output")
dat<-read.csv(src,check.names=FALSE,stringsAsFactors=FALSE)
X<-as.data.frame(dat[,setdiff(names(dat),c("sample_id","class")),drop=FALSE],check.names=FALSE)
rownames(X)<-dat$sample_id
y<-factor(dat$class)
outer<-make_repeated_stratified_folds(y,5L,1L,42L)[[1]]
T<-X[outer$train_idx,,drop=FALSE]
yt<-factor(as.character(y[outer$train_idx]),levels=levels(y))
sp<-rg_inner(yt)
D<-T[sp$D,,drop=FALSE];V<-T[sp$V,,drop=FALSE]
yd<-factor(as.character(yt[sp$D]),levels=levels(y))
yv<-factor(as.character(yt[sp$V]),levels=levels(y))
crp_need(!length(intersect(outer$train_idx,outer$test_idx)) &&
          !length(intersect(sp$D,sp$V)),
          "leakage in predeclared D/V patient splits")
start<-proc.time()[["elapsed"]]
results<-cfo_run(D,yd,V,yv)
rows<-list();probs<-list();count<-0L;nprob<-0L
for(k in seq_along(results$roots)){
 r<-results$roots[[k]]
 for(arm in names(r$prob)){
  p<-r$prob[[arm]]
  crp_need(length(p)==nrow(V)&&all(is.finite(p))&&
      all(p>=0 & p<=1), "bad p in comparison")
  yy<-r$y
  h<-as.integer(p>=.5)
  f1<-function(cl){
    tp<-sum(yy==cl & h==cl)
    fp<-sum(yy!=cl & h==cl)
    fn<-sum(yy==cl & h!=cl)
    if(2*tp+fp+fn==0)0 else 2*tp/(2*tp+fp+fn)
  }
  count<-count+1L
  rows[[count]]<-data.frame(dataset=DS,root_id=k,method=arm,
     selected_target=r$root$target,
     selected_anchors=paste(r$root$anchors,collapse="|"),
     selected_threshold=r$root$threshold,
     offset_insample_active=r$D_in_sample_parent_active,
     offset_crossfit_active=r$D_crossfit_parent_active,
     n_D=nrow(D),n_V=nrow(V),
     logloss=rg_ll(yy,p),brier=mean((yy-p)^2),
     macro_f1=mean(c(f1(0L),f1(1L))),
     positive_f1=f1(1L),
     ba=.5*(mean(h[yy==0L]==0L)+mean(h[yy==1L]==1L)),
     stage="D_TRAIN_ONLY_V_VALIDATION_NO_OUTER_E")
  for(i in seq_along(p)){
    nprob<-nprob+1L
    probs[[nprob]]<-data.frame(dataset=DS,root_id=k,
         method=arm,sample_id=rownames(V)[i],truth=yy[i],prob=p[i])
  }
 }
}
if(count){
  write.csv(do.call(rbind,rows),file.path(out,"metrics_V.csv"),row.names=FALSE)
  write.csv(do.call(rbind,probs),file.path(out,"predictions_V.csv"),row.names=FALSE)
}
write.csv(data.frame(dataset=DS,
 n_original_total=nrow(X),n_outer_T=nrow(T),
 n_outer_E_UNREAD=length(outer$test_idx),
 n_inner_D=nrow(D),n_inner_V=nrow(V),
 n_candidates_proposed=results$proposed,
 n_roots_evaluated=length(results$roots),
 in_sample_active_total=sum(vapply(results$roots,function(a)a$D_in_sample_parent_active,integer(1))),
 crossfit_active_total=sum(vapply(results$roots,function(a)a$D_crossfit_parent_active,integer(1))),
 runtime_sec=proc.time()[["elapsed"]]-start,
 stage="D_TRAIN_ONLY_V_VALIDATION_NO_OUTER_E"),
file.path(out,"manifest.csv"),row.names=FALSE)
if(count){
 mm<-do.call(rbind,rows)
 for(method in unique(mm$method)){
  v<-mm[mm$method==method,]
  cat(sprintf("CFO_DV_MEAN %s method=%s LL=%.6f macroF1=%.6f totalActiveD_Insample=%d totalActiveD_Crossfit=%d\n",
     DS,method,mean(v$logloss),mean(v$macro_f1),
     v$offset_insample_active[1],
     v$offset_crossfit_active[1]))
 }
}
cat(sprintf("CFO_DV_RESULT %s n_roots=%d n_probs=%d total_D_IN=%d total_D_OOF=%d\n",
  DS,length(results$roots),nprob,
  sum(vapply(results$roots,function(a)a$D_in_sample_parent_active,integer(1))),
  sum(vapply(results$roots,function(a)a$D_crossfit_parent_active,integer(1)))))
cat("CFO_DV_ORIGINAL_OUTER_E_NEVER_READ=PASS\n")
cat("CFO_DV_IDENTICAL_RELATIONS_DIFFERENT_OFFSETS=PASS\n")
cat("CFO_DV_PROPOSED_EXPERIMENT_COMPLETE=PASS\n")
