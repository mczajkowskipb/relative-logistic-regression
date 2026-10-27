# FIRST original CompCancer fold D/V ONLY, every original manuscript cohort.
# Two predeclared ordinal margins m10 and 10% of assay rank universe.
# Compare true CONDITIONAL small decision stump with SAME active-rule union
# without tree and a standard global four-rule sparse vote.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_selective_anchor_vote_stump_v1.R")

DS<-Sys.getenv("RLR_DATASET_ID",unset="")
datasets<-c("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
sv_need(DS %in% datasets,"unlisted manuscript cohort")
p<-file.path("results/rlr_full11_same_pairs_v1",DS,
             "source_data",paste0(DS,".csv"))
out<-file.path("results/rlr_savs_inner11_v1",DS)
sv_need(file.exists(p)&&
        file.exists(file.path(out,"SOURCE_SHA.txt"))&&
        !file.exists(file.path(out,"metrics_V.csv")),
        "missing public source / source SHA / repeated run")
dat<-read.csv(p,check.names=FALSE)
X<-as.data.frame(dat[,setdiff(names(dat),c("sample_id","class")),
                      drop=FALSE],check.names=FALSE)
rownames(X)<-dat$sample_id
y<-factor(dat$class)
reps<-if(DS %in% c("west-2001","yeoh-2002-v1"))2L else 10L
outer<-make_repeated_stratified_folds(y,5L,reps,42L)[[1L]]
T<-X[outer$train_idx,,drop=FALSE]
yt<-factor(as.character(y[outer$train_idx]),levels=levels(y))
inner<-rg_inner(yt)
D<-T[inner$D,,drop=FALSE]
V<-T[inner$V,,drop=FALSE]
yd<-factor(as.character(yt[inner$D]),levels=levels(y))
yv<-factor(as.character(yt[inner$V]),levels=levels(y))
sv_need(!length(intersect(outer$train_idx,outer$test_idx))&&
         !length(intersect(inner$D,inner$V)),
         "patient overlap in outer or inner splits")
truth<-sv_truth(yv)
metrics<-list();preds<-list();summ<-list()
i<-0L;j<-0L;k<-0L
start<-proc.time()[["elapsed"]]
pool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
for(mode in c("RLR10","RELATIVE10PCT")){
 tree<-sv_fit(D,yd,pool=pool,mode=mode)
 global_same<-sv_samepair_global(D,yd,tree)
 global_four<-sv_global4(D,yd,pool,tree$margin)
 root_only<-sv_root_only(D,yd,tree)
 models<-list(CONDITIONAL_ANCHOR_ABSTAIN=tree,
              GLOBAL_SAME_RULES=global_same,
              GLOBAL_FOUR_RULES=global_four,
              ANCHOR_ROOT_ONLY=root_only)
 for(name in names(models)){
   outp<-sv_predict(models[[name]],V)
   pr<-outp$p;cl<-as.integer(pr>=.5)
   yy<-truth
   f1<-function(z){
     tp<-sum(cl==z & yy==z)
     fp<-sum(cl==z & yy!=z)
     fn<-sum(cl!=z & yy==z)
     if(2*tp+fp+fn==0L)0 else 2*tp/(2*tp+fp+fn)
   }
   auc<-function(pp){
     n1<-sum(yy==1L);n0<-sum(yy==0L)
     rr<-rank(pp,ties.method="average")
     (sum(rr[yy==1L])-n1*(n1+1)/2)/(n0*n1)
   }
   i<-i+1L
   metrics[[i]]<-data.frame(
     dataset=DS,margin_mode=mode,method=name,
     n_D=nrow(D),n_V=nrow(V),
     original_outer_E_unread=length(outer$test_idx),
     median_pairs_active=stats::median(outp$n_active),
     mean_pairs_active=mean(outp$n_active),
     mean_pairs_examined=mean(outp$n_considered),
     root_comparisons=outp$n_root,
     pct_no_vote=mean(outp$abstain==1L),
     pct_le2_active=mean(outp$n_active<=2L),
     positive_f1=f1(1L),
     macro_f1=mean(c(f1(0L),f1(1L))),
     ba=.5*(mean(cl[yy==0L]==0L)+mean(cl[yy==1L]==1L)),
     auc=auc(pr),logloss=rg_ll(yy,pr),
     brier=mean((yy-pr)^2),
     evidence="PUBLIC_SOURCE_ORIGINAL_OUTER_T_ONLY_INNER_V_DEVELOPMENT")
   for(kp in seq_along(yy)){
      j<-j+1L
      preds[[j]]<-data.frame(
        dataset=DS,margin_mode=mode,method=name,
        sample_id=rownames(V)[kp],
        truth=yy[kp],prob=pr[kp],label=cl[kp],
        n_active=outp$n_active[kp],
        n_considered=outp$n_considered[kp],
        n_root=outp$n_root,
        vote_balance=outp$sum_votes[kp],
        branch=outp$branch[kp])
   }
 }
 k<-k+1L
 summ[[k]]<-data.frame(
  dataset=DS,margin_mode=mode,margin=tree$margin,
  root_target=if(is.null(tree$root))"" else tree$root$target,
  root_anchors=if(is.null(tree$root))"" else
            paste(tree$root$anchors,collapse="|"),
  root_threshold=if(is.null(tree$root))NA_integer_ else tree$root$threshold,
  root_n0=tree$training_branches[1L],
  root_n1=tree$training_branches[2L],
  selected_leaf_rules0=tree$n_rules[1L],
  selected_leaf_rules1=tree$n_rules[2L],
  max_per_patient_comparisons=
    max(tree$n_rules)+if(is.null(tree$root))0L else length(tree$root$anchors),
  training_stage="D_ONLY_NOMINATION_V_ONLY_AUDIT")
}
write.csv(do.call(rbind,metrics),file.path(out,"metrics_V.csv"),
           row.names=FALSE)
write.csv(do.call(rbind,preds),file.path(out,"patient_V.csv"),
           row.names=FALSE)
write.csv(do.call(rbind,summ),file.path(out,"structure_D.csv"),
           row.names=FALSE)
for(m in metrics){
 cat(sprintf("SAVS_DV_RESULT %s %s %s F1=%.5f MacroF1=%.5f BA=%.5f AUC=%.5f LL=%.5f mean_active=%.2f pct_abstain=%.3f\n",
 DS,m$margin_mode,m$method,m$positive_f1,m$macro_f1,m$ba,m$auc,
 m$logloss,m$mean_pairs_active,m$pct_no_vote))
}
cat(sprintf("SAVS_DV_DIAGNOSTIC %s D=%d V=%d rows=%d runtime=%.1fsec\n",
 DS,nrow(D),nrow(V),length(preds),
 proc.time()[["elapsed"]]-start))
cat("SAVS_DV_ORIGINAL_E_NEVER_TOUCHED=PASS\n")
cat("SAVS_DV_EIGHT_MATCHED_MODELS_AND_PATIENT_ACTIVITIES=PASS\n")
