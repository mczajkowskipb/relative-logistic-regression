source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_conditional_residual_pairs_v1.R")
DS<-Sys.getenv("RLR_DATASET_ID",unset="")
crp_need(DS %in% c("chen-2002","gordon-2002","singh-2002","yeoh-2002-v1"),
         "unexpected exploratory dataset")
rootpath<-file.path("results/rlr_full11_same_pairs_v1",DS)
file<-file.path(rootpath,"source_data",paste0(DS,".csv"))
out<-file.path("results/rlr_conditional_residual_pairs_DV_v1",DS)
crp_need(dir.exists(out)&&file.exists(file.path(out,"SOURCE_SHA.txt"))&&
         file.exists(file)&&!file.exists(file.path(out,"metrics_V.csv")),
         "source or output manifest invalid")
dat<-read.csv(file,check.names=FALSE)
X<-as.data.frame(dat[,setdiff(names(dat),c("sample_id","class")),
                      drop=FALSE],check.names=FALSE)
rownames(X)<-dat$sample_id
y<-factor(dat$class)
outer<-make_repeated_stratified_folds(y,5L,1L,42L)[[1]]
T<-X[outer$train_idx,,drop=FALSE]
yt<-factor(as.character(y[outer$train_idx]),levels=levels(y))
inner<-rg_inner(yt)
D<-T[inner$D,,drop=FALSE];V<-T[inner$V,,drop=FALSE]
yd<-factor(as.character(yt[inner$D]),levels=levels(y))
yv<-factor(as.character(yt[inner$V]),levels=levels(y))
crp_need(!length(intersect(outer$train_idx,outer$test_idx))&&
         !length(intersect(inner$D,inner$V)),
         "patient split overlap")
start<-proc.time()[["elapsed"]]
result<-crp_eval_V(D,yd,V,yv)
rows<-list();patient<-list();meta<-list()
ri<-0L;pi<-0L
label<-result$yy
f1<-function(y,h,cl=1L){
  tp<-sum(y==cl & h==cl)
  fp<-sum(y!=cl & h==cl)
  fn<-sum(y==cl & h!=cl)
  if(2*tp+fp+fn==0L)0 else 2*tp/(2*tp+fp+fn)
}
mF1<-function(y,h)mean(c(f1(y,h,0L),f1(y,h,1L)))
ba<-function(y,h).5*(mean(h[y==0L]==0L)+mean(h[y==1L]==1L))
if(length(result$roots)){
 for(i in seq_along(result$roots)){
   r<-result$roots[[i]]
   mi<-length(meta)+1L
   meta[[mi]]<-data.frame(dataset=DS,root_id=i,
     root_target=r$root$target,
     root_anchors=paste(r$root$anchors,collapse="|"),
     root_threshold=r$root$threshold,
     D_residual_proxy_score=r$rank_score,
     n_original_local_active=r$old$n_active,
     n_residual_selected_active=r$new$n_active,
     n_original_local_fits=r$old$n_fitted,
     n_residual_local_fits=r$new$n_fitted,
     n_candidate_relations=r$new$n_candidates,
     n_parent_outer_E_UNREAD=length(outer$test_idx),
     n_T=nrow(T),n_D=nrow(D),n_V=nrow(V),
     stage="INNER_DV_DEVELOPMENT_ONLY")
   for(arm in names(r$probs)){
     p<-r$probs[[arm]]
     crp_need(length(p)==nrow(V)&&all(is.finite(p)),
         "bad inner validation probabilities")
     h<-as.integer(p>=.5)
     ri<-ri+1L
     rows[[ri]]<-data.frame(dataset=DS,root_id=i,method=arm,
       n_D=nrow(D),n_V=nrow(V),
       positive_f1=f1(label,h),
       macro_f1=mF1(label,h),ba=ba(label,h),
       logloss=rg_ll(label,p),
       brier=mean((label-p)^2),
       stage="INNER_DV_DEVELOPMENT_ONLY")
     for(k in seq_along(label)){
       pi<-pi+1L
       patient[[pi]]<-data.frame(
        dataset=DS,root_id=i,method=arm,
        sample_id=rownames(V)[k],truth=label[k],prob=p[k])
     }
   }
 }
 utils::write.csv(do.call(rbind,rows),
             file.path(out,"metrics_V.csv"),row.names=FALSE)
 utils::write.csv(do.call(rbind,patient),
             file.path(out,"predictions_V.csv"),row.names=FALSE)
 utils::write.csv(do.call(rbind,meta),
             file.path(out,"root_support_D.csv"),row.names=FALSE)
}
nroot<-length(result$roots)
cat(sprintf("CRP_DV_RESULT dataset=%s D=%d V=%d proposed=%d roots=%d old_active=%d residual_active=%d old_best_deltaLL=%s residual_best_deltaLL=%s elapsed=%.2fsec\n",
  DS,nrow(D),nrow(V),result$n_proposed,nroot,
  sum(vapply(result$roots,function(r)r$old$n_active,integer(1))),
  sum(vapply(result$roots,function(r)r$new$n_active,integer(1))),
  if(nroot)format(max(vapply(result$roots,
      function(r)rg_ll(label,r$probs$ORIGINAL_RLR50)-
         rg_ll(label,r$probs$CLASS_TSP_CHILD_RLR),numeric(1))),
          digits=5) else "NA",
  if(nroot)format(max(vapply(result$roots,
      function(r)rg_ll(label,r$probs$ORIGINAL_RLR50)-
         rg_ll(label,r$probs$RESIDUAL_SCORED_CHILD_RLR),numeric(1))),
          digits=5) else "NA",
  proc.time()[["elapsed"]]-start))
cat("CRP_DV_NO_ORIGINAL_OUTER_E_LABELS_READ=PASS\n")
cat("CRP_DV_FROZEN_COHORT_DIAGNOSTIC_COMPLETE=PASS\n")
