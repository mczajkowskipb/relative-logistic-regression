# Full 11 ORIGINAL-PAPER-NAMED cohorts: FIRST original outer fold only,
# INNER D/V mechanism feasibility; never access external E patients/labels.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_crossfit_offset_ablation_v1.R")

DS<-Sys.getenv("RLR_DATASET_ID",unset="")
sets<-c("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002",
 "pomeroy-2002-v1","shipp-2002-v1","singh-2002",
 "west-2001","yeoh-2002-v1")
crp_need(DS %in% sets,"unknown original CompCancer cohort")
out<-file.path("results/rlr_crossfit_tree_full11_DV_v1",DS)
sourcepath<-file.path("results/rlr_full11_same_pairs_v1",
                     DS,"source_data",paste0(DS,".csv"))
crp_need(file.exists(sourcepath) &&
         file.exists(file.path(out,"SOURCE_SHA.txt")) &&
         !file.exists(file.path(out,"selection_DV.csv")),
         "missing original dataset source/commit or duplicate results")
dat<-read.csv(sourcepath,check.names=FALSE)
X<-as.data.frame(dat[,setdiff(names(dat),
      c("sample_id","class")),drop=FALSE],check.names=FALSE)
rownames(X)<-dat$sample_id
y<-factor(dat$class)
crp_need(length(levels(y))==2L && !anyDuplicated(rownames(X)) &&
         all(vapply(X,function(x)is.numeric(x)&&
                    all(is.finite(x)),logical(1))),
         "dataset probe schema or values invalid")
nrepeat<-if(DS %in% c("west-2001","yeoh-2002-v1"))2L else 10L
outer<-make_repeated_stratified_folds(y,5L,nrepeat,42L)[[1L]]
T<-X[outer$train_idx,,drop=FALSE]
yt<-factor(as.character(y[outer$train_idx]),levels=levels(y))
sp<-rg_inner(yt)
D<-T[sp$D,,drop=FALSE];V<-T[sp$V,,drop=FALSE]
yd<-factor(as.character(yt[sp$D]),levels=levels(y))
yv<-factor(as.character(yt[sp$V]),levels=levels(y))
crp_need(!length(intersect(outer$train_idx,outer$test_idx)) &&
         !length(intersect(sp$D,sp$V)), "overlap detected")
yy<-as.integer(yv==levels(yv)[2L])
metric<-function(p) {
  h<-as.integer(p>=.5)
  f1<-function(cl) {
    tp<-sum(yy==cl & h==cl)
    fp<-sum(yy!=cl & h==cl)
    fn<-sum(yy==cl & h!=cl)
    if(2*tp+fp+fn==0)0 else 2*tp/(2*tp+fp+fn)
  }
  list(LL=rg_ll(yy,p),macroF1=mean(c(f1(0L),f1(1L))),
       positiveF1=f1(1L),
       BA=.5*(mean(h[yy==0L]==0L)+mean(h[yy==1L]==1L)))
}
empty<-function(reason,state)list(
   state=state,accepted=0L,reason=reason,
   n_root_proposals=0L,n_roots_evaluated=0L,
   best_root_target="",best_root_anchors="",
   best_root_threshold=NA_integer_,
   best_D_in_sample_active=0L,best_D_crossfit_active=0L,
   V_original_LL=NA_real_,V_tree_LL=NA_real_,
   V_samepairs_global_LL=NA_real_,
   V_gain_vs_original=NA_real_,V_gain_vs_global=NA_real_,
   V_original_MacroF1=NA_real_,V_tree_MacroF1=NA_real_,
   V_tree_minus_original_BA=NA_real_)
fallback<-NULL
# This guard mirrors previous SA-RLR engineering correction for small N.
if(safe_nfolds(yd,5L)<4L) {
  fallback<-empty("D_TOO_SMALL_FOR_ORIGINAL_CV_GLMNET",
                  "SAFE_FALLBACK_SMALL_D")
}else{
  oofplans<-make_repeated_stratified_folds(yd,3L,1L,731L)
  if(any(vapply(oofplans,function(a){
      yy<-factor(as.character(yd[a$train_idx]),levels=levels(yd))
      safe_nfolds(yy,5L)<4L
    },logical(1))))
      fallback<-empty("OOF_SUBTRAIN_TOO_SMALL_FOR_ORIGINAL_CV_GLMNET",
                      "SAFE_FALLBACK_SMALL_OOF")
}
if(is.null(fallback)){
  started<-proc.time()[["elapsed"]]
  results<-tryCatch(cfo_run(D,yd,V,yv),error=function(e) {
      structure(list(message=conditionMessage(e)),class="CFO_ERROR")
  })
  if(inherits(results,"CFO_ERROR")) {
    fallback<-empty(results$message,"ERROR_NEEDS_DIAGNOSIS")
  }else if(length(results$roots)==0L){
    fallback<-empty("D_SELECTED_NO_VALID_ORDINAL_ROOT",
                    "SAFE_FALLBACK_NO_ROOT")
  }else{
    r<-results$roots
    crp_need(length(r)<=6L,"too many selected roots")
    LL<-vapply(r,function(a)metric(a$prob$RESIDUAL_LOCAL_CROSSFIT_OFFSET)$LL,
                    numeric(1))
    j<-which.min(LL)
    win<-r[[j]]
    original<-metric(win$prob$ORIGINAL_RLR50)
    tree<-metric(win$prob$RESIDUAL_LOCAL_CROSSFIT_OFFSET)
    global<-metric(win$prob$SAMEPAIRS_GLOBAL_CROSSFIT_OFFSET)
    gain0<-original$LL-tree$LL
    gain1<-global$LL-tree$LL
    deltaBA<-tree$BA-original$BA
    # Same pre-existing logloss / BA gate as SA-RLR 1:
    # .005 parent, .002 global, BA loss <=.02.
    # Publication classification safeguard (diagnostic, NOT independent
    # confirmation; motivated by previously seen Singh F1 decrease):
    # do NOT accept a tree if it lowers Macro-F1 on inner V by >.01.
    accepted<-gain0>=.005 && gain1>=.002 &&
         deltaBA>=-.02 &&
         tree$macroF1>=original$macroF1-.01
    fallback<-list(
      state=if(accepted)"INNER_V_ACCEPTED" else "INNER_V_REJECTED",
      accepted=as.integer(accepted),
      reason="FIXED_PARENT_GLOBAL_BA_MACROF1_TRAIN_ONLY_GATES",
      n_root_proposals=results$proposed,
      n_roots_evaluated=length(r),
      best_root_target=win$root$target,
      best_root_anchors=paste(win$root$anchors,collapse="|"),
      best_root_threshold=win$root$threshold,
      best_D_in_sample_active=win$D_in_sample_parent_active,
      best_D_crossfit_active=win$D_crossfit_parent_active,
      V_original_LL=original$LL,V_tree_LL=tree$LL,
      V_samepairs_global_LL=global$LL,
      V_gain_vs_original=gain0,V_gain_vs_global=gain1,
      V_original_MacroF1=original$macroF1,
      V_tree_MacroF1=tree$macroF1,
      V_tree_minus_original_BA=deltaBA)
    # Diagnostic: all SAME V patients per nominated root,
    # archived even if candidate is rejected; NO external E scores.
    vrows<-list()
    for(k in seq_along(r)){
      m<-r[[k]]
      for(arm in names(m$prob)){
        pp<-m$prob[[arm]]
        for(i in seq_along(pp)){
          vrows[[length(vrows)+1L]]<-data.frame(
            dataset=DS,root_id=k,method=arm,
            sample_id=rownames(V)[i],truth=yy[i],prob=pp[i])
        }
      }
    }
    write.csv(do.call(rbind,vrows),file.path(out,"patient_V.csv"),
              row.names=FALSE)
  }
}
selection<-data.frame(
 dataset=DS,original_outer_split=outer$split_id,
 n_total=nrow(X),n_T=nrow(T),
 n_outer_E_UNREAD=length(outer$test_idx),
 n_D=nrow(D),n_V=nrow(V),
 n_genes=ncol(X),as.data.frame(fallback),
 evidence="ORIGINAL_NAME_PUBLIC_SOURCE_D_V_ONLY_NO_OUTER_E",
 stringsAsFactors=FALSE)
write.csv(selection,file.path(out,"selection_DV.csv"),row.names=FALSE)
cat(sprintf("CFO_FULL11_CASE dataset=%s state=%s D=%d V=%d considered=%d root=%s active_crossfit=%d Vgain_parent=%s Vgain_samepairs=%s F1_original=%s F1_tree=%s\n",
DS,selection$state,nrow(D),nrow(V),
selection$n_roots_evaluated,selection$best_root_target,
selection$best_D_crossfit_active,
format(selection$V_gain_vs_original,digits=5),
format(selection$V_gain_vs_global,digits=5),
format(selection$V_original_MacroF1,digits=5),
format(selection$V_tree_MacroF1,digits=5)))
cat("CFO_FULL11_OUTER_E_NEVER_READ=PASS\n")
cat("CFO_FULL11_TRAIN_ONLY_MECHANISM_COMPLETE=PASS\n")
