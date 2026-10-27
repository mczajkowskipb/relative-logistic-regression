# Original manuscript-like algorithm comparison on the SAME two reused
# CompCancer source reconstructions and six matched 3-fold outer E partitions.
# Requires R/new/07_world_protocol_utils.R unchanged. Development-only.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_anchor_lookahead_child_rlr_v1.R")
OUT <- "results/rlr_paper_target_tree_v1"
la_need(dir.exists(OUT) && file.exists(file.path(OUT,"SOURCE_SHA.txt")),
        "isolated output not prepared")
ARMS <- c("TSP1","KTSP10","KTSP20","KTSP50",
          "RLR10","RLR20","RLR50","GREEDY_RLR20_TREE",
          "GREEDY_RLR50_TREE","VALUE_LR")
cohorts <- c("chowdary-2006","singh-2002")
met <- list();rows<-list();folds_audit<-list()
mi<-0L;pi<-0L;fi<-0L
ms_f1 <- function(yy,pr) {
  one<-function(c){
    tp<-sum(yy==c & pr==c)
    fp<-sum(yy!=c & pr==c)
    fn<-sum(yy==c & pr!=c)
    if (2*tp+fp+fn==0)0 else 2*tp/(2*tp+fp+fn)
  }
  mean(c(one(0L),one(1L)))
}
ms_f1_pos<-function(yy,pr){
  tp<-sum(yy==1L & pr==1L)
  fp<-sum(yy==0L & pr==1L)
  fn<-sum(yy==1L & pr==0L)
  if(2*tp+fp+fn==0)0 else 2*tp/(2*tp+fp+fn)
}
ms_auc<-function(y,p) {
  one<-which(y==1L);zero<-which(y==0L)
  r<-rank(p,ties.method="average")
  as.numeric((sum(r[one])-length(one)*(length(one)+1)/2)/
                (length(one)*length(zero)))
}
for(ds in cohorts) {
  dat<-read.csv(file.path(
    "results/rlr_pair_specific_evidence_v1/source_data",
    paste0(ds,".csv")),check.names=FALSE)
  X<-as.data.frame(dat[,setdiff(names(dat),
       c("sample_id","class")),drop=FALSE],check.names=FALSE)
  rownames(X)<-dat$sample_id
  y<-factor(dat$class)
  specs<-make_repeated_stratified_folds(y,3L,1L,42L)
  for(f in seq_along(specs)){
    it<-specs[[f]]$train_idx;ie<-specs[[f]]$test_idx
    T<-X[it,,drop=FALSE];E<-X[ie,,drop=FALSE]
    yt<-factor(as.character(y[it]),levels=levels(y))
    ye<-factor(as.character(y[ie]),levels=levels(y))
    yy<-la_truth(ye)
    runfolder<-file.path(OUT,"intermediate",
        ds,paste0("fold",f-1L))
    dir.create(runfolder,recursive=TRUE,showWarnings=FALSE)
    pool<-compute_top_pairs_switchbox(T,yt,top_n=1000L,
                buffer_mult=5L)
    # Genuine historical SWAP.KTSP.* not a hard-vote proxy.
    for(k in c(10L,20L,50L))
      fit_ktsp_train_test(T,yt,E,ye,k=k,
          out_dir=file.path(runfolder,paste0("ktsp",k)),seed=42L)
    fit_tsp1_manual_train_test(T,yt,E,ye,pool,
          out_dir=file.path(runfolder,"tsp1"))
    fit_value_lr_train_test(T,yt,E,ye,
          out_dir=file.path(runfolder,"value_lr"),
          pre_top_features=1000L,alpha=.5,
          lambda_choice="lambda.min",seed=42L)
    parent10<-rg_parent(T,yt,E,pool,K=10L)
    parent20<-rg_parent(T,yt,E,pool,K=20L)
    parent50<-rg_parent(T,yt,E,pool,K=50L)
    # One D-only greedy anchor root/pair set; evaluate conditional RLR
    # under two different RLR offsets, no external E-driven gating.
    chosen<-la_select(T,yt)
    pred20<-parent20$test
    pred50<-parent50$test
    active20<-active50<-0L
    if(!is.null(chosen$greedy_root)){
      fit20<-la_fit_structure(T,yt,parent20$train,
             chosen$greedy_root,pairs_fixed=chosen$greedy_pairs)
      fit50<-la_fit_structure(T,yt,parent50$train,
             chosen$greedy_root,pairs_fixed=chosen$greedy_pairs)
      pred20<-la_predict_tree(E,parent20$test,fit20)
      pred50<-la_predict_tree(E,parent50$test,fit50)
      active20<-fit20$active;active50<-fit50$active
    }
    P<-list(
      TSP1=read.csv(file.path(runfolder,"tsp1","predictions_test.csv")),
      KTSP10=read.csv(file.path(runfolder,"ktsp10","predictions_test.csv")),
      KTSP20=read.csv(file.path(runfolder,"ktsp20","predictions_test.csv")),
      KTSP50=read.csv(file.path(runfolder,"ktsp50","predictions_test.csv")),
      VALUE_LR=read.csv(file.path(runfolder,"value_lr","predictions_test.csv"))
    )
    probs<-list(TSP1=P$TSP1$prob,
      KTSP10=P$KTSP10$prob,
      KTSP20=P$KTSP20$prob,
      KTSP50=P$KTSP50$prob,
      RLR10=parent10$test,
      RLR20=parent20$test,
      RLR50=parent50$test,
      GREEDY_RLR20_TREE=pred20,
      GREEDY_RLR50_TREE=pred50,
      VALUE_LR=P$VALUE_LR$prob)
    fi<-fi+1L
    folds_audit[[fi]]<-data.frame(dataset=ds,fold=f-1L,
      n_train=nrow(T),n_E=nrow(E),
      train_ids=paste(rownames(T),collapse="|"),
      E_ids=paste(rownames(E),collapse="|"),
      RLR10_selected=nrow(parent10$selected),
      RLR20_selected=nrow(parent20$selected),
      RLR50_selected=nrow(parent50$selected),
      greedy_root_found=!is.null(chosen$greedy_root),
      greedy_root=if(is.null(chosen$greedy_root))"" else
                    chosen$greedy_root$target,
      greedy_active20=active20,greedy_active50=active50,
      stage="EXPOSED_COHORT_DEVELOPMENT_ONLY")
    for(a in ARMS) {
      p<-as.numeric(probs[[a]])
      la_need(length(p)==nrow(E)&&
                all(is.finite(p))&&all(p>=0&p<=1),
              paste("invalid method probabilities",a))
      label<-if(a %in% names(P))
        as.integer(as.character(P[[a]]$y_pred)==levels(ye)[2L])
        else as.integer(p>=.5)
      la_need(length(label)==length(p),"bad native predicted class")
      if(a %in% names(P))
        la_need(identical(as.character(P[[a]]$sample_id),rownames(E)),
                "mismatched original baseline predicted patient order")
      mi<-mi+1L
      met[[mi]]<-data.frame(dataset=ds,fold=f-1L,method=a,
        n_train=nrow(T),n_E=nrow(E),
        auc=ms_auc(yy,p),macro_f1=ms_f1(yy,label),
        positive_f1=ms_f1_pos(yy,label),
        ba=.5*(mean(label[yy==0L]==0L)+
                mean(label[yy==1L]==1L)),
        logloss=rg_ll(yy,p),brier=mean((yy-p)^2),
        stage="EXPOSED_COHORT_DEVELOPMENT_ONLY")
      for(i in seq_along(p)){
        pi<-pi+1L
        rows[[pi]]<-data.frame(dataset=ds,fold=f-1L,
           method=a,sample_id=rownames(E)[i],
           truth=yy[i],prob=p[i],pred_label=label[i])
      }
    }
    message(sprintf(
      "PAPER_TARGET_CASE %s fold=%d kTSP50_F1=%.3f RLR50_F1=%.3f Tree50_F1=%.3f RLR50_LL=%.5f Tree50_LL=%.5f",
      ds,f-1L,
      met[[mi-6L]]$macro_f1,met[[mi-3L]]$macro_f1,
      met[[mi-1L]]$macro_f1,
      met[[mi-3L]]$logloss,met[[mi-1L]]$logloss))
  }
}
m<-do.call(rbind,met);p<-do.call(rbind,rows)
audit<-do.call(rbind,folds_audit)
la_need(nrow(m)==60L&&nrow(p)==2060L&&nrow(audit)==6L,
        "benchmark inventory mismatch")
write.csv(m,file.path(OUT,"method_metrics.csv"),row.names=FALSE)
write.csv(p,file.path(OUT,"patient_predictions.csv"),row.names=FALSE)
write.csv(audit,file.path(OUT,"split_manifest.csv"),row.names=FALSE)
for(a in ARMS){
  v<-m[m$method==a,,drop=FALSE]
  message(sprintf("PAPER_TARGET_SUMMARY %-20s AUC=%0.5f MacroF1=%0.5f positiveF1=%0.5f BA=%0.5f LL=%0.5f",
     a,mean(v$auc),mean(v$macro_f1),mean(v$positive_f1),
     mean(v$ba),mean(v$logloss)))
}
cat("PAPER_TARGET_6_REAL_COHORT_E_FOLDS=PASS\n")
cat("PAPER_TARGET_60_ALIGNED_METHOD_METRICS=PASS\n")
cat("PAPER_TARGET_2060_E_PATIENT_PROBABILITY_ROWS=PASS\n")
cat("PAPER_TARGET_NOT_FULL_PAPER_OR_EXTERNAL_VALIDATION\n")
