# Frozen original 11 CompCancer cohorts: single original protocol E fold per
# dataset, development diagnostic only (not original 470-fold confirmation).
# Actual SWAP.KTSP50 versus original RLR50 and D-OOF decision thresholds.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_conditional_ordinal_threshold_v1.R")

DS<-Sys.getenv("RLR_DATASET_ID",unset="")
datasets<-c("alizadeh-2000-v1","armstrong-2002-v1",
 "bittner-2000","chen-2002","golub-1999-v1","gordon-2002",
 "pomeroy-2002-v1","shipp-2002-v1","singh-2002",
 "west-2001","yeoh-2002-v1")
ct_need(DS %in% datasets,"wrong cohort")
base<-file.path("results/rlr_full11_same_pairs_v1",DS)
path<-file.path(base,"source_data",paste0(DS,".csv"))
out<-file.path("results/rlr_conditional_threshold_11first_v1",DS)
ct_need(file.exists(path)&&
        file.exists(file.path(base,"input_source_manifest.json"))&&
        dir.exists(out)&&
        file.exists(file.path(out,"SOURCE_SHA.txt"))&&
        !file.exists(file.path(out,"method_metrics.csv")),"bad input provenance")
dat<-read.csv(path,check.names=FALSE,stringsAsFactors=FALSE)
X<-as.data.frame(dat[,setdiff(names(dat),
                     c("sample_id","class")),drop=FALSE],check.names=FALSE)
rownames(X)<-dat$sample_id
y<-factor(dat$class)
ct_need(length(levels(y))==2L && !anyDuplicated(rownames(X)) &&
   all(vapply(X,function(z)is.numeric(z)&&all(is.finite(z)),logical(1))),
   "invalid original source")
nreps<-if(DS %in% c("west-2001","yeoh-2002-v1"))2L else 10L
sp<-make_repeated_stratified_folds(y,5L,nreps,42L)[[1L]]
T<-X[sp$train_idx,,drop=FALSE];E<-X[sp$test_idx,,drop=FALSE]
yt<-factor(as.character(y[sp$train_idx]),levels=levels(y))
ye<-factor(as.character(y[sp$test_idx]),levels=levels(y))
truth<-ct_truth(ye)
ct_need(!length(intersect(sp$train_idx,sp$test_idx)),
   "outer T and E patient overlap")
# The continuous parent probabilities are EXACTLY original RLR50 under
# same candidate selection, m10, disjoint d2, alpha=.5, lambda.min.
selected<-ct_select(T,yt)
pool<-compute_top_pairs_switchbox(T,yt,
        top_n=1000L,buffer_mult=5L)
parent<-rg_parent(T,yt,E,pool,K=50L)
p<-parent$test
rE<-aa_ranks(E)
RLR50<-as.integer(p>=.5)
GLOBAL<-as.integer(p>=selected$global_t)
CONDITIONAL<-if(is.null(selected$root))GLOBAL else
  as.integer(p>=ifelse(os_gate(rE,selected$root)==0L,
                       selected$t0,selected$t1))
ACCEPTED<-if(isTRUE(selected$accepted))CONDITIONAL else GLOBAL
# Real switchBox K50, not a proxy, with original native predicted classes.
ktdir<-file.path(out,"kTSP50_work")
fit_ktsp_train_test(T,yt,E,ye,k=50L,out_dir=ktdir,seed=42L)
kt<-read.csv(file.path(ktdir,"predictions_test.csv"),stringsAsFactors=FALSE)
ct_need(identical(as.character(kt$sample_id),rownames(E)),
     "kTSP original predictions E patient misaligned")
kt_label<-as.integer(as.character(kt$y_pred)==levels(ye)[2L])
kt_p<-as.numeric(kt$prob)
ct_need(length(kt_p)==nrow(E)&&all(is.finite(kt_p)),
        "invalid genuine SWAP KTSP50 probabilities")
A<-list(
 KTSP50_NATIVE=list(prob=kt_p,label=kt_label),
 ORIGINAL_RLR50_THRESHOLD05=list(prob=p,label=RLR50),
 RLR50_DOOF_GLOBAL_THRESHOLD=list(prob=p,label=GLOBAL),
 RLR50_ORDINAL_TREE_DIAGNOSTIC=list(prob=p,label=CONDITIONAL),
 RLR50_ORDINAL_TREE_V_GATED=list(prob=p,label=ACCEPTED))
auc<-function(y,p){
  n1<-sum(y==1L);n0<-sum(y==0L)
  a<-rank(p,ties.method="average")
  (sum(a[y==1L])-n1*(n1+1)/2)/(n0*n1)
}
metric<-list();patient<-list();i<-0L;j<-0L
for(method in names(A)){
   prob<-A[[method]]$prob;label<-A[[method]]$label
   ct_need(length(prob)==length(truth)&&
           length(label)==length(truth),"prediction length mismatch")
   i<-i+1L
   metric[[i]]<-data.frame(
      dataset=DS,split_id=sp$split_id,method=method,
      n_train=nrow(T),n_E=nrow(E),
      positive_f1=ct_f1(truth,label),
      macro_f1=ct_macro_f1(truth,label),
      ba=ct_ba(truth,label),auc=auc(truth,prob),
      logloss=rg_ll(truth,prob),
      brier=mean((truth-prob)^2),
      evidence="REUSED_PUBLIC_SOURCE_NOT_ORIGINAL_ARFF_SHA_PARITY")
   for(k in seq_along(truth)){
     j<-j+1L
     patient[[j]]<-data.frame(
         dataset=DS,split_id=sp$split_id,method=method,
         sample_id=rownames(E)[k],truth=truth[k],
         prob=prob[k],label=label[k],
         RLR50_original_prob=p[k])
   }
}
m<-do.call(rbind,metric);pr<-do.call(rbind,patient)
ct_need(nrow(m)==5L && nrow(pr)==5L*nrow(E),
   "missing some metric rows")
write.csv(m,file.path(out,"method_metrics.csv"),row.names=FALSE)
write.csv(pr,file.path(out,"patient_predictions.csv"),row.names=FALSE)
write.csv(data.frame(
  dataset=DS,split_id=sp$split_id,n_train=nrow(T),n_E=nrow(E),
  original_positive_level=levels(y)[2L],
  T_ids=paste(rownames(T),collapse="|"),
  E_ids=paste(rownames(E),collapse="|"),
  original_K=50L,original_m=10L,original_d=2L,
  original_alpha=.5,original_lambda="lambda.min",
  source_stage="EXPOSED_PUBLIC_SOURCE_RECONSTRUCTION_ONLY"),
  file.path(out,"split_manifest.csv"),row.names=FALSE)
write.csv(data.frame(
  dataset=DS,accepted=as.integer(selected$accepted),
  root_found=as.integer(!is.null(selected$root)),
  root_target=if(is.null(selected$root))"" else selected$root$target,
  root_references=if(is.null(selected$root))"" else
                   paste(selected$root$anchors,collapse="|"),
  root_threshold=if(is.null(selected$root))NA_integer_ else
                   selected$root$threshold,
  global_threshold=selected$global_t,
  conditional_t0=selected$t0,conditional_t1=selected$t1,
  D_F1_gain=selected$D_F1_gain,
  V_F1_gain=selected$V_F1_gain,
  V_BA_gain=selected$V_BA_gain,
  roots_proposed=selected$proposed,
  roots_scored=selected$candidates,reason=selected$reason,
  source_stage="EXPOSED_PUBLIC_SOURCE_RECONSTRUCTION_ONLY"),
  file.path(out,"decision_selection.csv"),row.names=FALSE)
# Exact patient-level checks: changing classification cannot change
# original RLR positive probabilities, AUC, logloss or Brier.
r<-m[m$method=="ORIGINAL_RLR50_THRESHOLD05",]
for(arm in c("RLR50_DOOF_GLOBAL_THRESHOLD",
             "RLR50_ORDINAL_TREE_DIAGNOSTIC",
             "RLR50_ORDINAL_TREE_V_GATED")){
  z<-m[m$method==arm,]
  ct_need(identical(z$auc,r$auc)&&
          identical(z$logloss,r$logloss)&&
          identical(z$brier,r$brier),
          "ordinal classification unexpectedly changed original RLR score")
}
message(sprintf(
"COT_CASE %s E=%d root_found=%s accepted=%s global_t=%.2f F1_KTSP=%.4f F1_RLR=%.4f F1_GLOBAL=%.4f F1_TREE_GATED=%.4f BA_RLR=%.4f BA_TREE=%.4f",
DS,nrow(E),!is.null(selected$root),selected$accepted,
selected$global_t,
m$positive_f1[1L],m$positive_f1[2L],
m$positive_f1[3L],m$positive_f1[5L],
m$ba[2L],m$ba[5L]))
cat("COT_5_MATCHED_METHODS_ONE_ORIGINAL_OUTER_FOLD=PASS\n")
cat("COT_SAME_RLR50_PROB_AUC_LOSS_BRIER_ALL_THRESHOLD_ARMS=PASS\n")
cat("COT_E_EXPLICIT_DEVELOPMENT_NOT_FULL470=PASS\n")
