source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_selective_anchor_vote_lookahead_v2.R")
DS<-Sys.getenv("RLR_DATASET_ID",unset="")
sv2_need(DS %in% c("chen-2002","pomeroy-2002-v1",
                   "west-2001","yeoh-2002-v1"),
         "non-predeclared hard cohort")
path<-file.path("results/rlr_full11_same_pairs_v1",DS,
           "source_data",paste0(DS,".csv"))
out<-file.path("results/rlr_savs_lookahead_hard4_DV_v2",DS)
sv2_need(file.exists(path) && file.exists(file.path(out,"SOURCE_SHA.txt")) &&
           !file.exists(file.path(out,"metrics_V.csv")),
           "source or output manifest absent")
dat<-read.csv(path,check.names=FALSE)
X<-as.data.frame(dat[,setdiff(names(dat),c("sample_id","class")),
                    drop=FALSE],check.names=FALSE)
rownames(X)<-dat$sample_id
y<-factor(dat$class)
nrep<-if(DS %in% c("west-2001","yeoh-2002-v1"))2L else 10L
fold<-make_repeated_stratified_folds(y,5L,nrep,42L)[[1L]]
T<-X[fold$train_idx,,drop=FALSE]
yt<-factor(as.character(y[fold$train_idx]),levels=levels(y))
parts<-rg_inner(yt)
D<-T[parts$D,,drop=FALSE]; V<-T[parts$V,,drop=FALSE]
yd<-factor(as.character(yt[parts$D]),levels=levels(y))
yv<-factor(as.character(yt[parts$V]),levels=levels(y))
sv2_need(!length(intersect(parts$D,parts$V)),"train/V subject leakage")
pool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
report<-list()
preds<-list();ri<-0L;pi<-0L
for(mode in c("RLR10","RELATIVE10PCT")){
 select<-sv2_train_only_select(D,yd,pool,mode=mode,max_roots=12L)
 look<-sv2_final_fit(D,yd,pool,select)
 orig<-sv_fit(D,yd,pool,mode=mode)
 same<-sv_samepair_global(D,yd,orig)
 glob<-sv_global4(D,yd,pool,sv_m(ncol(D),mode))
 opts<-list(LOOKAHEAD_TRAIN_ONLY_SELECTED=look,
    V1_ROOT_FIRST_TREE=orig,
    V1_GLOBAL_SAME_RULES=same,
    GLOBAL_FOUR_RULES=glob)
 for(key in names(opts)){
  pred<-sv_predict(opts[[key]],V)
  stat<-sv2_metrics(yv,pred$p)
  ri<-ri+1L
  report[[ri]]<-data.frame(dataset=DS,mode=mode,method=key,
     n_D=nrow(D),n_V=nrow(V),
     accepted_tree_D_CV=as.integer(select$accepted),
     root_candidates_evaluated=select$evaluated,
     root_candidates_proposed=select$proposed,
     D_CV_macroF1=select$cv_macroF1,
     D_CV_global_macroF1=select$cv_global_macroF1,
     positive_f1=stat$positiveF1,
     macro_f1=stat$macroF1,ba=stat$BA,
     logloss=stat$LL,
     mean_active=mean(pred$n_active),
     mean_compared=mean(pred$n_considered)+pred$n_root,
     mean_abstain=mean(pred$n_active==0L),
     source="TRAIN_ONLY_D_NESTED_LOOKAHEAD_V_EXTERNAL_TO_D_NO_E")
  for(i in seq_along(pred$p)){
   pi<-pi+1L
   preds[[pi]]<-data.frame(dataset=DS,mode=mode,method=key,
    sample_id=rownames(V)[i],
    truth=as.integer(yv[i]==levels(yv)[2L]),
    prob=pred$p[i],label=as.integer(pred$p[i]>=.5),
    active=pred$n_active[i],considered=pred$n_considered[i],
    root_tests=pred$n_root)
  }
 }
 cat(sprintf("SAVS_V2_CASE %s mode=%s root_accept=%s roots=%d V_tree_F1=%.5f V_global4_F1=%.5f original_tree_F1=%.5f\n",
   DS,mode,select$accepted,select$evaluated,
   sv2_metrics(yv,sv_predict(look,V)$p)$positiveF1,
   sv2_metrics(yv,sv_predict(glob,V)$p)$positiveF1,
   sv2_metrics(yv,sv_predict(orig,V)$p)$positiveF1))
}
write.csv(do.call(rbind,report),file.path(out,"metrics_V.csv"),row.names=FALSE)
write.csv(do.call(rbind,preds),file.path(out,"patient_V.csv"),row.names=FALSE)
cat("SAVS_V2_FOUR_HARD_COHORT_D_ONLY_LOOKAHEAD_COMPLETE=PASS\n")
cat("SAVS_V2_OUTER_E_UNREAD=PASS\n")
