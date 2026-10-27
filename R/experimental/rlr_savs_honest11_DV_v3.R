# Select SAVS3 on D's honest S/H only; report on V (no original E).
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_savs_honest_lookahead_v3.R")
ds<-Sys.getenv("RLR_DATASET_ID",unset="")
cohorts<-c("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
sv3_need<-rg_need
sv3_need(ds %in% cohorts,"bad original manuscript cohort")
src<-file.path("results/rlr_full11_same_pairs_v1",ds,
               "source_data",paste0(ds,".csv"))
out<-file.path("results/rlr_savs_honest11_DV_v3",ds)
sv3_need(file.exists(src)&&file.exists(file.path(out,"SOURCE_SHA.txt"))&&
      !file.exists(file.path(out,"metrics_V.csv")),"data/output manifest")
t<-read.csv(src,check.names=FALSE)
X<-as.data.frame(t[,setdiff(colnames(t),c("sample_id","class")),
                   drop=FALSE],check.names=FALSE)
rownames(X)<-t$sample_id
y<-factor(t$class)
reps<-if(ds %in% c("west-2001","yeoh-2002-v1"))2L else 10L
osplit<-make_repeated_stratified_folds(y,5L,reps,42L)[[1L]]
T<-X[osplit$train_idx,,drop=FALSE]
yt<-factor(as.character(y[osplit$train_idx]),levels=levels(y))
i<-rg_inner(yt)
D<-T[i$D,,drop=FALSE]
V<-T[i$V,,drop=FALSE]
yd<-factor(as.character(yt[i$D]),levels=levels(y))
yv<-factor(as.character(yt[i$V]),levels=levels(y))
sv3_need(!length(intersect(osplit$train_idx,osplit$test_idx))&&
         !length(intersect(i$D,i$V)),"outer or inner overlap")
poolD<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
metrics<-list();predictions<-list();choices<-list();a<-0L;b<-0L;c<-0L
for(mode in c("RLR10","RELATIVE10PCT")){
 sel<-sv3_select(D,yd,mode=mode,max_roots=12L)
 chosen<-sv3_refit(D,yd,poolD,sel)
 root_first<-sv_fit(D,yd,pool=poolD,mode=mode)
 same<-sv_samepair_global(D,yd,root_first)
 global4<-sv_global4(D,yd,poolD,sv_m(ncol(D),mode))
 models<-list(
  HONEST_LOOKAHEAD_V_GATED=chosen,
  V1_ROOT_FIRST_TREE=root_first,
  V1_GLOBAL_SAME_RULES=same,
  GLOBAL_FOUR_RULES=global4)
 c<-c+1L
 choices[[c]]<-data.frame(
  dataset=ds,margin_mode=mode,
  n_outer_E_UNREAD=length(osplit$test_idx),
  n_T=nrow(T),n_D=nrow(D),n_V=nrow(V),
  n_S=sel$n_S,n_H=sel$n_H,
  honest_tree_accepted=as.integer(sel$accepted),
  n_S_proposed_roots=sel$n_proposals,
  H_tree_macroF1=sel$H_best_tree_F1,
  H_global_macroF1=sel$H_global_F1,
  H_tree_BA=sel$H_best_tree_BA,
  H_global_BA=sel$H_global_BA,
  root=if(is.null(sel$root))"" else sel$root$target,
  anchors=if(is.null(sel$root))"" else
                  paste(sel$root$anchors,collapse="|"),
  threshold=if(is.null(sel$root))NA_integer_ else sel$root$threshold,
  reason=sel$reason,
  evidence="S_ONLY_ROOT_POOL_H_SELECTION_D_FINALFIT_V_TEST_NO_ORIGINAL_E")
 for(name in names(models)){
  res<-sv_predict(models[[name]],V)
  met<-sv2_metrics(yv,res$p)
  a<-a+1L
  metrics[[a]]<-data.frame(
    dataset=ds,margin_mode=mode,method=name,
    n_V=nrow(V),positive_f1=met$positiveF1,
    macro_f1=met$macroF1,ba=met$BA,
    logloss=met$LL,
    mean_active=mean(res$n_active),
    median_active=stats::median(res$n_active),
    mean_comparisons=mean(res$n_considered)+res$n_root,
    pct_zero_active=mean(res$n_active==0L),
    pct_le2_active=mean(res$n_active<=2L),
    tree_accepted=as.integer(sel$accepted))
  for(k in seq_len(nrow(V))){
    b<-b+1L
    predictions[[b]]<-data.frame(
      dataset=ds,margin_mode=mode,method=name,
      sample_id=rownames(V)[k],
      truth=sv_truth(yv)[k],
      prob=res$p[k],label=as.integer(res$p[k]>=.5),
      n_active=res$n_active[k],
      n_considered=res$n_considered[k],n_root=res$n_root,
      branch=res$branch[k],vote=res$sum_votes[k])
  }
 }
 cat(sprintf("SAVS3_CASE %s %s S=%d H=%d V=%d gate=%s n_roots=%d F1=%.4f global4F1=%.4f oldTreeF1=%.4f\n",
  ds,mode,sel$n_S,sel$n_H,nrow(V),sel$accepted,
  sel$n_proposals,
  sv2_metrics(yv,sv_predict(chosen,V)$p)$positiveF1,
  sv2_metrics(yv,sv_predict(global4,V)$p)$positiveF1,
  sv2_metrics(yv,sv_predict(root_first,V)$p)$positiveF1))
}
write.csv(do.call(rbind,metrics),file.path(out,"metrics_V.csv"),row.names=FALSE)
write.csv(do.call(rbind,predictions),file.path(out,"patient_V.csv"),row.names=FALSE)
write.csv(do.call(rbind,choices),file.path(out,"honest_selection_D.csv"),row.names=FALSE)
cat("SAVS3_HONEST_ALL_METHODS_V_METRICS=PASS\n")
cat("SAVS3_OUTER_E_NEVER_USED=PASS\n")
