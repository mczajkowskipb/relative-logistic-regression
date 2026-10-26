# Paired real-data test of RLR component changes, NO trees.
# 2 previously exposed CompCancer benchmark cohorts from original RLR paper.
# Training-only candidate pool/shortlist/quantiles/CV. 3-fold E heldout.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_pair_specific_evidence_v1.R")

OUT<-"results/rlr_pair_specific_evidence_v1"
rlr_margin_need(dir.exists(OUT)&&file.exists(file.path(OUT,"SOURCE_SHA.txt"))&&
                 !file.exists(file.path(OUT,"method_metrics.csv")),
                 "missing pinned result directory or previous output")
SETS<-c("chowdary-2006","singh-2002")
MODES<-c("TERNARY_M10","ADAPTIVE_Q25","GRADED_M10_Q75")
KS<-c(10L,20L,50L)
MET<-list()
PRED<-list()
PROTO<-list()
mi<-0L;pi<-0L;si<-0L
for(ds in SETS) {
  path<-file.path(OUT,"source_data",paste0(ds,".csv"))
  rlr_margin_need(file.exists(path),"source dataset unavailable")
  dat<-utils::read.csv(path,check.names=FALSE,stringsAsFactors=FALSE)
  X<-as.data.frame(dat[,setdiff(names(dat),c("sample_id","class")),drop=FALSE],
                   check.names=FALSE)
  rownames(X)<-dat$sample_id
  y<-factor(dat$class)
  rlr_margin_need(length(levels(y))==2L&&
        all(vapply(X,function(z)is.numeric(z)&&all(is.finite(z)),
                   logical(1))) && nrow(X)>=100L,
        "public gene-expression dataset input invalid")
  folds<-make_repeated_stratified_folds(y,n_folds=3L,n_repeats=1L,seed=42L)
  rlr_margin_need(length(folds)==3L,"outer folds unexpected")
  for(i in seq_along(folds)) {
    fold<-folds[[i]]
    tr<-fold$train_idx;te<-fold$test_idx
    rlr_margin_need(!length(intersect(tr,te))&&
      identical(sort(c(tr,te)),seq_len(nrow(X))),
      "outer fold splits overlap")
    T<-X[tr,,drop=FALSE];E<-X[te,,drop=FALSE]
    yt<-y[tr];ye<-as.numeric(y[te]==levels(y)[2L])
    pool<-compute_top_pairs_switchbox(
      X_train=T,y_train=yt,top_n=1000L,buffer_mult=5L,
      out_dir=NULL,save_all=FALSE)
    rlr_margin_need(nrow(pool)>=50L,"publication ranking pool exhausted")
    message(sprintf("RLR_MARGIN_START dataset=%s fold=%d train=%d E=%d",
                    ds,i-1L,length(tr),length(te)))
    for(K in KS) {
      baseline<-NULL
      for(mode in MODES) {
        r<-rlr_margin_one(T,yt,E,pool=pool,K=K,mode=mode,seed=42L)
        if(mode=="TERNARY_M10")baseline<-r$prob
        rlr_margin_need(!is.null(baseline),"missing original RLR")
        p<-pmin(pmax(r$prob,1e-15),1-1e-15)
        po<-pmin(pmax(baseline,1e-15),1-1e-15)
        ll<-function(q)-ye*log(q)-(1-ye)*log1p(-q)
        mi<-mi+1L
        MET[[mi]]<-data.frame(
          dataset=ds,fold=i-1L,K=K,mode=mode,
          n_train=length(tr),n_E=length(te),
          global_active=r$globally_active,
          mean_local_active=r$mean_local_active,
          margin_mean=r$margin_mean,
          margin_min=r$margin_min,margin_max=r$margin_max,
          mean_train_activity=r$train_active_mean,
          mean_test_activity=r$test_active_mean,
          logloss=mean(ll(p)),brier=mean((ye-r$prob)^2),
          delta_logloss=mean(ll(p)-ll(po)),
          scientific_stage="REUSED_COMPCANCER_DEVELOPMENT_ONLY",
          stringsAsFactors=FALSE)
        for(j in seq_along(te)){
          pi<-pi+1L
          PRED[[pi]]<-data.frame(
            dataset=ds,fold=i-1L,K=K,mode=mode,
            sample_id=rownames(E)[[j]],truth=ye[[j]],
            original_prob=baseline[[j]],method_prob=r$prob[[j]],
            stringsAsFactors=FALSE)
        }
        message(sprintf("RLR_MARGIN_CASE %s fold=%d K=%d %s deltaLL=%+.6f active_global=%d margin_mean=%.2f",
                ds,i-1L,K,mode,
                MET[[mi]]$delta_logloss,r$globally_active,r$margin_mean))
      }
    }
    si<-si+1L
    PROTO[[si]]<-data.frame(dataset=ds,fold=i-1L,
       n_train=length(tr),n_E=length(te),
       n_features=ncol(X),root_source="ORIGINAL_RLR_SWITCHBOX_POOL1000",
       test_id_list=paste(rownames(E),collapse="|"),
       stringsAsFactors=FALSE)
  }
}
m<-do.call(rbind,MET)
pr<-do.call(rbind,PRED)
protocol<-do.call(rbind,PROTO)
rlr_margin_need(nrow(m)==2L*3L*3L*3L &&
   nrow(protocol)==6L &&
   all(m$mode %in% MODES) &&
   all(is.finite(m$logloss)),
   "study inventory failed")
write.csv(m,file.path(OUT,"method_metrics.csv"),row.names=FALSE)
write.csv(pr,file.path(OUT,"patient_predictions.csv"),row.names=FALSE)
write.csv(protocol,file.path(OUT,"split_manifest.csv"),row.names=FALSE)
for(K in KS) for(mode in MODES){
  r<-m[m$K==K&m$mode==mode,,drop=FALSE]
  message(sprintf("RLR_MARGIN_SUMMARY K=%d %s mean_deltaLL=%+.7f improved_E=%d/6",
          K,mode,mean(r$delta_logloss),sum(r$delta_logloss<(-1e-9))))
}
cat("RLR_MARGIN_2_COMPCANCER_PUBLICATION_COHORTS=PASS\n")
cat("RLR_MARGIN_6_FROZEN_OUTER_FOLDS=PASS\n")
cat("RLR_MARGIN_54_METHOD_FITS=PASS\n")
cat("RLR_MARGIN_ALL_HELDOUT_PATIENT_PREDICTIONS_ARCHIVED=PASS\n")
cat("RLR_MARGIN_NO_TREE_NO_ORACLE_NO_NEW_COHORT=PASS\n")
