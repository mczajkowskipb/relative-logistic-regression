# ORIGINAL named CompCancer11 same E patients on FIRST original 5-fold split.
# Compare original RLR20/50 to margin-aware selector using exact m10
# ternary rank features; REAL switchBox kTSP nominal K10/15/20/25 and
# train-only matched global-nonzero / train-sample-active RLR budgets.
# No outer E feedback to select any method/hyperparameter.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_margin_aware_budget_core_v1.R")

ds<-Sys.getenv("RLR_DATASET_ID",unset="")
sets<-c("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
mb_need(ds %in% sets,"unknown original dataset")
out<-file.path("results/rlr_margin_budget_first11_v1",ds)
src<-file.path("results/rlr_full11_same_pairs_v1",ds,
               "source_data",paste0(ds,".csv"))
mb_need(file.exists(src)&&file.exists(file.path(out,"SOURCE_SHA.txt"))&&
        !file.exists(file.path(out,"method_metrics.csv")),
        "required source SHA or source data missing")
d<-read.csv(src,stringsAsFactors=FALSE,check.names=FALSE)
X<-as.data.frame(d[,setdiff(names(d),c("sample_id","class")),
                   drop=FALSE],check.names=FALSE)
rownames(X)<-d$sample_id
y<-factor(d$class)
nrep<-if(ds %in% c("west-2001","yeoh-2002-v1"))2L else 10L
fold<-make_repeated_stratified_folds(y,5L,nrep,42L)[[1L]]
T<-X[fold$train_idx,,drop=FALSE]
E<-X[fold$test_idx,,drop=FALSE]
yt<-factor(as.character(y[fold$train_idx]),levels=levels(y))
ye<-factor(as.character(y[fold$test_idx]),levels=levels(y))
truth<-as.integer(ye==levels(ye)[2L])
mb_need(length(intersect(rownames(T),rownames(E)))==0L &&
        min(table(ye))>0L,"outer split overlap or single E class")
start<-proc.time()[["elapsed"]]
pool<-compute_top_pairs_switchbox(T,yt,top_n=1000L,buffer_mult=5L)
newpool<-mb_reorder(T,yt,pool,m=10L)
full<-list()
for(K in c(20L,50L)){
 original<-rg_parent(T,yt,E,pool,K=K)
 modified<-rg_parent(T,yt,E,newpool,K=K)
 nm<-paste0("RLR",K)
 full[[nm]]<-list(fit=original,train=mb_active(original,T),
                  test=mb_active(original,E))
 nm<-paste0("MARGIN_RLR",K)
 full[[nm]]<-list(fit=modified,train=mb_active(modified,T),
                  test=mb_active(modified,E))
}
br20<-full$RLR20$train
br50<-full$RLR50$train
kopts<-c(
 KTSP10=10L,KTSP15=15L,KTSP20=20L,KTSP25=25L,
 KTSP_M20_GLOBAL=mb_budget(br20$globally_nonzero),
 KTSP_M20_LOCAL=mb_budget(mean(br20$locally_active)),
 KTSP_M50_GLOBAL=mb_budget(br50$globally_nonzero),
 KTSP_M50_LOCAL=mb_budget(mean(br50$locally_active)))
# If two matched budgets coincide the same absolute K, still keep
# distinct pre-declared comparisons; no E-controlled re-fitting choices.
for(nm in names(kopts)){
  fit_ktsp_train_test(T,yt,E,ye,k=as.integer(kopts[[nm]]),
    out_dir=file.path(out,"baselines",nm),seed=42L)
}
p_true<-function(z){
   z<-as.numeric(z)
   mb_need(length(z)==length(truth)&&all(is.finite(z)),
           "nonfinite original E probabilities")
   z
}
all<-list()
for(nm in names(full)){
  s<-full[[nm]]
  p<-p_true(s$fit$test)
  h<-as.integer(p>=.5)
  all[[nm]]<-list(p=p,h=h,
    active=s$test$locally_active,
    tested=s$test$tested,
    actual=s$test$globally_nonzero,
    selected=s$test$nominal,
    unique_genes=s$test$n_genes,
    train_active=mean(s$train$locally_active),
    train_global=s$train$globally_nonzero)
}
for(nm in names(kopts)){
  path<-file.path(out,"baselines",nm)
  pred<-read.csv(file.path(path,"predictions_test.csv"),
                 stringsAsFactors=FALSE)
  kk<-nrow(read.csv(file.path(path,"selected_pairs.csv")))
  mb_need(identical(as.character(pred$sample_id),rownames(E)) &&
          kk>=1L && kk<=as.integer(kopts[[nm]]),
          "native switchBox realized K or patient mismatch")
  all[[nm]]<-list(p=p_true(pred$prob),
    h=as.integer(as.character(pred$y_pred)==levels(ye)[2L]),
    active=rep(kk,nrow(E)),tested=rep(kk,nrow(E)),
    actual=kk,selected=kk,
    unique_genes=length(unique(unlist(
      read.csv(file.path(path,"selected_pairs.csv"))[,c("gene_left","gene_right")],
      use.names=FALSE))),
    train_active=kk,train_global=kk)
}
f1<-function(h,cl){
 tp<-sum(truth==cl & h==cl);fp<-sum(truth!=cl & h==cl)
 fn<-sum(truth==cl & h!=cl)
 if(2*tp+fp+fn==0L)0 else 2*tp/(2*tp+fp+fn)
}
auc<-function(p){
 n1<-sum(truth==1L);n0<-sum(truth==0L)
 rr<-rank(p,ties.method="average")
 (sum(rr[truth==1L])-n1*(n1+1L)/2)/(n1*n0)
}
rows<-list();raw<-list();i<-0L;j<-0L
for(nm in names(all)){
  z<-all[[nm]]
  mb_need(length(z$active)==nrow(E)&&length(z$tested)==nrow(E) &&
          all(z$active<=z$tested) && all(z$active>=0L),
          "patient local activity exceeds globally nonzero tests")
  i<-i+1L
  rows[[i]]<-data.frame(dataset=ds,split_id="rep00_fold00",
   method=nm,n_T=nrow(T),n_E=nrow(E),
   nominal_pairs=z$selected,
   global_nonzero_pairs=z$actual,
   test_mean_active_pairs=mean(z$active),
   test_mean_tested_pairs=mean(z$tested),
   train_mean_active_pairs=z$train_active,
   train_global_nonzero_pairs=z$train_global,
   n_genes_used=z$unique_genes,
   positive_f1=f1(z$h,1L),
   macro_f1=mean(c(f1(z$h,0L),f1(z$h,1L))),
   ba=.5*(mean(z$h[truth==0L]==0L)+mean(z$h[truth==1L]==1L)),
   auc=auc(z$p),logloss=rg_ll(truth,z$p),
   brier=mean((truth-z$p)^2),
   evidence="FIRST_OUTER_E_EXPLORATORY_PUBLIC_RECONSTRUCTION")
  for(k in seq_along(truth)){
    j<-j+1L
    raw[[j]]<-data.frame(
      dataset=ds,split_id="rep00_fold00",method=nm,
      sample_id=rownames(E)[k],truth=truth[k],
      prob=z$p[k],label=z$h[k],active=z$active[k],tested=z$tested[k])
  }
}
write.csv(do.call(rbind,rows),file.path(out,"method_metrics.csv"),
          row.names=FALSE)
write.csv(do.call(rbind,raw),file.path(out,"patient_E.csv"),
          row.names=FALSE)
write.csv(data.frame(dataset=ds,split_id="rep00_fold00",
   K20_global=br20$globally_nonzero,
   K20_train_patient=mean(br20$locally_active),
   K50_global=br50$globally_nonzero,
   K50_train_patient=mean(br50$locally_active),
   n_T=nrow(T),n_E=nrow(E),
   T_ids=paste(rownames(T),collapse="|"),
   E_ids=paste(rownames(E),collapse="|"),
   top_switchbox_pool=nrow(pool),
   elapsed_seconds=proc.time()[["elapsed"]]-start,
   provenance="K_MATCH_ALWAYS_CHOSEN_FROM_ORIGINAL_RLR_OUTER_T_ONLY"),
   file.path(out,"train_budget_manifest.csv"),row.names=FALSE)
for(z in rows)cat(sprintf("MB_RLR_CASE %s %s F1=%.5f AUC=%.5f LL=%.5f trained=%d glob=%d train_local=%.2f test_local=%.2f\n",
 ds,z$method,z$positive_f1,z$auc,z$logloss,
 z$nominal_pairs,z$global_nonzero_pairs,
 z$train_mean_active_pairs,z$test_mean_active_pairs))
cat("MB_RLR_ALL_MATCHED_E_METHODS_PREDICTIONS=PASS\n")
cat("MB_RLR_ALL_COMPLEXITY_MATCHING_FROM_T_ONLY=PASS\n")
cat("MB_RLR_ORIGINAL_PROCESSED_ARFF_PARITY_UNVERIFIED\n")
