# OAS-DT V1: EXPOSED original-paper development cohorts, no novelty claim.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_ordinal_anchor_score_tree_v1.R")
OUT<-"results/rlr_ordinal_anchor_score_tree_v1"
os_need(dir.exists(OUT)&&file.exists(file.path(OUT,"SOURCE_SHA.txt"))&&
        !file.exists(file.path(OUT,"method_metrics.csv")),
        "isolated output absent or already written")
SETS<-c("chowdary-2006","singh-2002")
ARMS<-c("ORIGINAL_RLR","ANCHOR_SCORE_STUMP","ANCHOR_SCORE_TREE_D2",
        "SAME_ROOT_BEST_SINGLE_PAIR","SAME_COUNT_RIDGE",
        "V_SELECTED_ANCHOR_SCORE_TREE","RLR_V_GATED_ANCHOR_SCORE_TREE")
all_m<-list(); all_p<-list(); all_s<-list(); all_v<-list()
mi<-0L;pi<-0L;si<-0L
metric_f1<-function(y,p){
  pred<-as.integer(p>=.5)
  f1<-function(k) {
    tp<-sum(y==k & pred==k)
    fp<-sum(y!=k & pred==k)
    fn<-sum(y==k & pred!=k)
    if(2*tp+fp+fn==0)return(0)
    2*tp/(2*tp+fp+fn)
  }
  mean(c(f1(0L),f1(1L)))
}
for(ds in SETS) {
  path<-file.path("results/rlr_pair_specific_evidence_v1/source_data",
                   paste0(ds,".csv"))
  os_need(file.exists(path),paste("cohort source not found:",path))
  dat<-read.csv(path,check.names=FALSE,stringsAsFactors=FALSE)
  X<-as.data.frame(dat[,setdiff(names(dat),c("sample_id","class")),
                     drop=FALSE],check.names=FALSE)
  rownames(X)<-dat$sample_id
  y<-factor(dat$class)
  os_need(length(levels(y))==2L&&all(vapply(X,is.numeric,logical(1)))&&
          all(vapply(X,function(z)all(is.finite(z)),logical(1))),
          "input matrix invalid")
  folds<-make_repeated_stratified_folds(y,n_folds=3L,n_repeats=1L,seed=42L)
  for(f in seq_along(folds)) {
    tr<-folds[[f]]$train_idx;te<-folds[[f]]$test_idx
    os_need(!length(intersect(tr,te))&&
       identical(sort(c(tr,te)),seq_len(nrow(X))),"outer E leakage")
    T<-X[tr,,drop=FALSE];E<-X[te,,drop=FALSE]
    yt<-factor(as.character(y[tr]),levels=levels(y))
    ye<-factor(as.character(y[te]),levels=levels(y))
    result<-os_evaluate(T,yt,E,ye)
    choice<-result$selection
    yy<-result$truth
    baseline<-result$probs$ORIGINAL_RLR
    si<-si+1L
    all_s[[si]]<-data.frame(dataset=ds,fold=f-1L,
      chosen=choice$mode,accepted=as.integer(choice$accepted),
      root_gene=choice$root_gene,k=choice$k,
      count_threshold=choice$threshold,
      anchor_genes=choice$anchor_names,
      proposed_children=choice$proposed_children,
      V_gain_vs_RLR=choice$v_gain_rlr,
      V_gain_vs_single=choice$v_gain_single,
      V_gain_vs_ridge=choice$v_gain_global,
      V_BA_diff=choice$v_BA_diff,reason=choice$reason,
      stage="EXPOSED_ORIGINAL_PAPER_DEVELOPMENT_ONLY")
    all_v[[si]]<-data.frame(dataset=ds,fold=f-1L,
      n_train=length(tr),n_E=length(te),n_features=ncol(X),
      training_ids=paste(rownames(T),collapse="|"),
      E_ids=paste(rownames(E),collapse="|"),
      stage="EXPOSED_ORIGINAL_PAPER_DEVELOPMENT_ONLY")
    for(arm in ARMS){
      pp<-result$probs[[arm]]
      os_need(length(pp)==length(te)&&all(is.finite(pp))&&
         all(pp>=0 & pp<=1),"bad predicted probability")
      mi<-mi+1L
      all_m[[mi]]<-data.frame(dataset=ds,fold=f-1L,
        method=arm,n_train=length(tr),n_E=length(te),
        chosen=choice$mode,accepted=as.integer(choice$accepted),
        root_k=choice$k,root_threshold=choice$threshold,
        depth2_children=choice$proposed_children,
        logloss=rg_ll(yy,pp),brier=mean((yy-pp)^2),
        balanced_accuracy=rg_ba(yy,pp),
        macro_f1=metric_f1(yy,pp),
        delta_logloss=rg_ll(yy,pp)-rg_ll(yy,baseline),
        stage="EXPOSED_ORIGINAL_PAPER_DEVELOPMENT_ONLY")
      for(j in seq_along(pp)){
        pi<-pi+1L
        all_p[[pi]]<-data.frame(dataset=ds,fold=f-1L,method=arm,
          sample_id=rownames(E)[j],truth=yy[j],prob=pp[j],
          original_prob=baseline[j])
      }
    }
    message(sprintf("OASDT_CASE %s fold=%d depth=%s accepted=%s children=%d rootK=%d t=%s V_RLR=%s V_SINGLE=%s E_stump=%+.6f E_depth2=%+.6f E_gated=%+.6f",
      ds,f-1L,choice$mode,choice$accepted,
      choice$proposed_children,choice$k,
      format(choice$threshold),format(choice$v_gain_rlr,digits=3),
      format(choice$v_gain_single,digits=3),
      rg_ll(yy,result$probs$ANCHOR_SCORE_STUMP)-rg_ll(yy,baseline),
      rg_ll(yy,result$probs$ANCHOR_SCORE_TREE_D2)-rg_ll(yy,baseline),
      rg_ll(yy,result$probs$RLR_V_GATED_ANCHOR_SCORE_TREE)-rg_ll(yy,baseline)))
  }
}
m<-do.call(rbind,all_m);p<-do.call(rbind,all_p)
s<-do.call(rbind,all_s);v<-do.call(rbind,all_v)
os_need(nrow(m)==42L && nrow(p)==1442L &&
          nrow(s)==6L&&nrow(v)==6L,"expected inventory mismatch")
write.csv(m,file.path(OUT,"method_metrics.csv"),row.names=FALSE)
write.csv(p,file.path(OUT,"patient_predictions.csv"),row.names=FALSE)
write.csv(s,file.path(OUT,"selection_audit.csv"),row.names=FALSE)
write.csv(v,file.path(OUT,"split_manifest.csv"),row.names=FALSE)
for(a in ARMS) {
  q<-m[m$method==a,,drop=FALSE]
  message(sprintf("OASDT_SUMMARY %s mean_deltaLL=%+.7f improved_E=%d/6 mean_delta_BA=%+.7f",
      a,mean(q$delta_logloss),sum(q$delta_logloss<(-1e-9)),
      mean(q$balanced_accuracy-
           m$balanced_accuracy[m$method=="ORIGINAL_RLR"])))
}
cat("OASDT_6_DISJOINT_OUTER_FOLDS=PASS\n")
cat("OASDT_42_FOLD_METHOD_METRICS=PASS\n")
cat("OASDT_1442_E_PATIENT_PROBABILITY_ROWS=PASS\n")
cat("OASDT_NOT_INDEPENDENT_BIOMEDICAL_VALIDATION\n")
