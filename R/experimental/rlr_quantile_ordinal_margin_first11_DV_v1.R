# QRLR v1 development D/V. Same original 11 paper-named cohorts.
# Strict first outer TRAIN only; split into D (fit all panel thresholds) + V.
# No original E observations, outcomes or repeated outer E tuning.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_quantile_ordinal_margin_core_v1.R")

ds<-Sys.getenv("RLR_DATASET_ID",unset="")
all_ds<-c("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
qr_need(ds %in% all_ds,"outside original manuscript datasets")
source_file<-file.path("results/rlr_full11_same_pairs_v1",
 ds,"source_data",paste0(ds,".csv"))
out<-file.path("results/rlr_quantile_ordinal_margin11_DV_v1",ds)
qr_need(file.exists(source_file) && file.exists(file.path(out,"SOURCE_SHA.txt")) &&
         !file.exists(file.path(out,"metrics_V.csv")),"source/output manifest")
data<-read.csv(source_file,stringsAsFactors=FALSE,check.names=FALSE)
X<-as.data.frame(data[,setdiff(names(data),c("sample_id","class")),
                      drop=FALSE],check.names=FALSE)
rownames(X)<-data$sample_id
y<-factor(data$class)
reps<-if(ds %in% c("west-2001","yeoh-2002-v1"))2L else 10L
outer<-make_repeated_stratified_folds(y,5L,reps,42L)[[1L]]
T<-X[outer$train_idx,,drop=FALSE]
yt<-factor(as.character(y[outer$train_idx]),levels=levels(y))
dvi<-rg_inner(yt)
D<-T[dvi$D,,drop=FALSE];V<-T[dvi$V,,drop=FALSE]
yd<-factor(as.character(yt[dvi$D]),levels=levels(y))
yv<-factor(as.character(yt[dvi$V]),levels=levels(y))
truth<-as.integer(yv==levels(yv)[2L])
qr_need(!length(intersect(outer$train_idx,outer$test_idx))&&
         !length(intersect(dvi$D,dvi$V)) &&
         nrow(D)+nrow(V)==nrow(T),"D/V or outer train overlap")
pool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
f1<-function(h,cl){
 tp<-sum(truth==cl & h==cl)
 fp<-sum(truth!=cl & h==cl)
 fn<-sum(truth==cl & h!=cl)
 if(2*tp+fp+fn==0)0 else 2*tp/(2*tp+fp+fn)
}
auc<-function(p){
 n1<-sum(truth==1L);n0<-sum(truth==0L)
 rr<-rank(p,ties.method="average")
 (sum(rr[truth==1L])-n1*(n1+1L)/2)/(n0*n1)
}
met<-list();pat<-list();schema<-list();i<-0L;j<-0L;k<-0L
for(K in c(10L,20L)){
 original<-rg_parent(D,yd,V,pool,K=K)
 for(q in c(0,.25,.5)){
  fitted<-qr_parent(D,yd,V,pool,K=K,q=q)
  if(q==0){
   qr_need(isTRUE(all.equal(fitted$test,original$test,
                          tolerance=1e-9)),
           paste("new quantile q0 not exact original RLR parity",K))
   qr_need(identical(as.character(fitted$selected$ordered_pair),
                  as.character(original$selected$ordered_pair)),
           "original pair selection mismatch")
  }
  pp<-fitted$test;h<-as.integer(pp>=.5)
  name<-sprintf("QRLR_K%02d_Q%02d",K,as.integer(100*q))
  i<-i+1L
  met[[i]]<-data.frame(dataset=ds,method=name,
    K_requested=K,q=q,n_D=nrow(D),n_V=nrow(V),
    positive_f1=f1(h,1L),macro_f1=mean(c(f1(h,0L),f1(h,1L))),
    ba=.5*(mean(h[truth==0L]==0L)+mean(h[truth==1L]==1L)),
    auc=auc(pp),logloss=rg_ll(truth,pp),
    brier=mean((truth-pp)^2),
    n_global_nonzero=fitted$global_active,
    mean_active_V=mean(fitted$test_active),
    mean_active_D=mean(fitted$train_active),
    n_tested=fitted$n_tested,
    nominal_K=nrow(fitted$selected),
    median_threshold=stats::median(fitted$margin),
    max_threshold=max(fitted$margin),
    evidence="D_FIT_UNSUPERVISED_THRESHOLDS_V_TEST_NO_ORIGINAL_OUTER_E")
  for(z in seq_along(pp)){
   j<-j+1L
   pat[[j]]<-data.frame(dataset=ds,method=name,
      sample_id=rownames(V)[z],truth=truth[z],prob=pp[z],
      label=h[z],n_active=fitted$test_active[z],
      n_tested=fitted$n_tested)
  }
  k<-k+1L
  schema[[k]]<-data.frame(dataset=ds,method=name,
    pair_ids=paste(fitted$selected$ordered_pair,collapse="|"),
    margins=paste(fitted$margin,collapse="|"),
    q=q)
  cat(sprintf("QRLR_DV %s %s F1=%.6f MacroF1=%.6f AUC=%.6f LL=%.6f global=%d local=%.3f medmargin=%.1f\n",
    ds,name,met[[i]]$positive_f1,met[[i]]$macro_f1,
    met[[i]]$auc,met[[i]]$logloss,
    fitted$global_active,mean(fitted$test_active),
    stats::median(fitted$margin)))
 }
}
write.csv(do.call(rbind,met),file.path(out,"metrics_V.csv"),row.names=FALSE)
write.csv(do.call(rbind,pat),file.path(out,"patient_V.csv"),row.names=FALSE)
write.csv(do.call(rbind,schema),file.path(out,"frozen_D_margin_panels.csv"),
          row.names=FALSE)
write.csv(data.frame(dataset=ds,n_D=nrow(D),n_V=nrow(V),n_outer_E_UNREAD=
    length(outer$test_idx),D_ids=paste(rownames(D),collapse="|"),
    V_ids=paste(rownames(V),collapse="|"),
    evidence="ORIGINAL_OUTER_T_ONLY_D_V_SEPARATE"),
    file.path(out,"manifest.csv"),row.names=FALSE)
cat("QRLR_ORIGINAL_Q0_PROBABILITY_PARITY=PASS\n")
cat("QRLR_ALL_SIX_DV_MODELS_AND_PATIENT_ACTIVITY=PASS\n")
cat("QRLR_OUTER_E_NEVER_USED=PASS\n")
