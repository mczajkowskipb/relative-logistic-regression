# Paired, development-only lookahead anchor RLR model tree experiment.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_anchor_lookahead_child_rlr_v1.R")
OUT<-"results/rlr_anchor_lookahead_child_rlr_v1"
la_need(dir.exists(OUT)&&file.exists(file.path(OUT,"SOURCE_SHA.txt"))&&
        !file.exists(file.path(OUT,"method_metrics.csv")),
        "expected isolated output path absent or already used")
SETS<-c("chowdary-2006","singh-2002")
ARMS<-c("ORIGINAL_RLR","GREEDY_ROOT_CHILD_RLR",
        "LOOKAHEAD_CHILD_RLR","LOOKAHEAD_SAME_PAIRS_GLOBAL",
        "LOOKAHEAD_GATE_ONLY","V_GATED_LOOKAHEAD_CHILD_RLR")
all_met<-list();all_predictions<-list();all_selection<-list()
all_manifest<-list()
mi<-0L;pi<-0L;fi<-0L
metric_macro_f1<-function(y,p){
  h<-as.integer(p>=.5)
  f<-function(cl){
    tp<-sum(h==cl & y==cl);fp<-sum(h==cl & y!=cl)
    fn<-sum(h!=cl & y==cl)
    if(2*tp+fp+fn==0)0 else 2*tp/(2*tp+fp+fn)
  }
  mean(c(f(0L),f(1L)))
}
for(ds in SETS){
  path<-file.path("results/rlr_pair_specific_evidence_v1/source_data",
                  paste0(ds,".csv"))
  la_need(file.exists(path),"original cohort source absent")
  tab<-read.csv(path,check.names=FALSE,stringsAsFactors=FALSE)
  X<-as.data.frame(tab[,setdiff(names(tab),c("sample_id","class")),
         drop=FALSE],check.names=FALSE)
  rownames(X)<-tab$sample_id
  y<-factor(tab$class)
  la_need(length(levels(y))==2L && !anyDuplicated(rownames(X)) &&
    all(vapply(X,function(z)is.numeric(z)&&all(is.finite(z)),
             logical(1))),"invalid source cohort")
  folds<-make_repeated_stratified_folds(y,n_folds=3L,
                n_repeats=1L,seed=42L)
  for(f in seq_along(folds)){
    sp<-folds[[f]]
    tr<-sp$train_idx;test<-sp$test_idx
    la_need(!length(intersect(tr,test))&&
        identical(sort(c(tr,test)),seq_len(nrow(X))),
        "outer train/test overlap or lost subjects")
    T<-X[tr,,drop=FALSE];E<-X[test,,drop=FALSE]
    yt<-factor(as.character(y[tr]),levels=levels(y))
    ye<-factor(as.character(y[test]),levels=levels(y))
    res<-la_evaluate(T,yt,E,ye)
    sel<-res$selection
    baseline<-res$probs$ORIGINAL_RLR
    yy<-res$truth
    fi<-fi+1L
    root<-sel$root
    all_selection[[fi]]<-data.frame(
      dataset=ds,fold=f-1L,
      accepted=as.integer(sel$accepted),
      root_target=if(is.null(root))"" else root$target,
      root_anchor_genes=if(is.null(root))"" else
                            paste(root$anchors,collapse="|"),
      root_threshold=if(is.null(root))NA_integer_ else root$threshold,
      n_roots_proposed=sel$proposed,n_roots_screened=sel$considered,
      n_child_RLR_fitted_D=sel$working,
      n_child_active_coefficients_D=sel$active,
      n_child_active_coefficients_T=res$realized_active,
      V_gain_vs_parent=sel$v_gain_parent,
      V_gain_vs_global=sel$v_gain_global,
      V_gain_vs_gate=sel$v_gain_gate,
      V_gain_vs_greedy=sel$v_gain_greedy,
      V_BA_difference=sel$v_ba_delta,
      reason=sel$reason,stage="EXPOSED_RLR_PAPER_DEVELOPMENT_ONLY")
    all_manifest[[fi]]<-data.frame(dataset=ds,fold=f-1L,
      n_T=length(tr),n_E=length(test),n_features=ncol(X),
      T_ids=paste(rownames(T),collapse="|"),
      E_ids=paste(rownames(E),collapse="|"),
      stage="EXPOSED_RLR_PAPER_DEVELOPMENT_ONLY")
    for(arm in ARMS){
      p<-res$probs[[arm]]
      la_need(length(p)==nrow(E)&&all(is.finite(p))&&
         all(p>=0&p<=1),"invalid E predictions")
      mi<-mi+1L
      all_met[[mi]]<-data.frame(dataset=ds,fold=f-1L,
        method=arm,n_T=length(tr),n_E=length(test),
        accepted=as.integer(sel$accepted),
        n_child_active=res$realized_active,
        logloss=rg_ll(yy,p),brier=mean((yy-p)^2),
        ba=rg_ba(yy,p),macro_f1=metric_macro_f1(yy,p),
        delta_logloss=rg_ll(yy,p)-rg_ll(yy,baseline),
        stage="EXPOSED_RLR_PAPER_DEVELOPMENT_ONLY")
      for(i in seq_along(p)){
        pi<-pi+1L
        all_predictions[[pi]]<-data.frame(dataset=ds,fold=f-1L,
          method=arm,sample_id=rownames(E)[i],truth=yy[i],
          prob=p[i],original_prob=baseline[i])
      }
    }
    message(sprintf("LOOKAHEAD_CASE %s fold=%d accepted=%s screened=%d/%d child_models_D=%d active_D=%d active_T=%d Vgain=%s Elook=%+.6f Egreedy=%+.6f Eglobal=%+.6f",
      ds,f-1L,sel$accepted,sel$considered,sel$proposed,
      sel$working,sel$active,res$realized_active,
      format(sel$v_gain_parent,digits=4),
      rg_ll(yy,res$probs$LOOKAHEAD_CHILD_RLR)-rg_ll(yy,baseline),
      rg_ll(yy,res$probs$GREEDY_ROOT_CHILD_RLR)-rg_ll(yy,baseline),
      rg_ll(yy,res$probs$LOOKAHEAD_SAME_PAIRS_GLOBAL)-rg_ll(yy,baseline)))
  }
}
metrics<-do.call(rbind,all_met)
pred<-do.call(rbind,all_predictions)
selected<-do.call(rbind,all_selection)
manifest<-do.call(rbind,all_manifest)
la_need(nrow(metrics)==36L&&nrow(pred)==1236L&&
          nrow(selected)==6L&&nrow(manifest)==6L,
        "incomplete experimental inventory")
write.csv(metrics,file.path(OUT,"method_metrics.csv"),row.names=FALSE)
write.csv(pred,file.path(OUT,"patient_predictions.csv"),row.names=FALSE)
write.csv(selected,file.path(OUT,"selection_audit.csv"),row.names=FALSE)
write.csv(manifest,file.path(OUT,"split_manifest.csv"),row.names=FALSE)
for(arm in ARMS){
  s<-metrics[metrics$method==arm,,drop=FALSE]
  message(sprintf("LOOKAHEAD_SUMMARY %s mean_E_deltaLL=%+.7f E_wins=%d/6 mean_BA=%.5f",
   arm,mean(s$delta_logloss),
   sum(s$delta_logloss<(-1e-9)),mean(s$ba)))
}
cat("LOOKAHEAD_6_E_PARTITIONS=PASS\n")
cat("LOOKAHEAD_36_MODEL_FOLD_METRICS=PASS\n")
cat("LOOKAHEAD_1236_PATIENT_METHOD_PREDICTIONS=PASS\n")
cat("LOOKAHEAD_REUSED_ORIGINAL_PAPER_DEVELOPMENT_ONLY=PASS\n")
