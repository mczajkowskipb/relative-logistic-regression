# Two original CompCancer paper cohorts, six frozen E outer folds.
# Exploratory, repeatedly exposed cohorts: DO NOT call this fresh confirmation.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_anchor_majority_tree_v1.R")

OUT<-"results/rlr_anchor_majority_tree_v1"
aa_need(dir.exists(OUT) && file.exists(file.path(OUT,"SOURCE_SHA.txt")) &&
        !file.exists(file.path(OUT,"method_metrics.csv")),
        "isolated existing or unprepared results directory")
SETS<-c("chowdary-2006","singh-2002")
ARMS<-c("ORIGINAL_RLR","ANCHOR_STUMP","ANCHOR_TREE_D2",
        "SAME_COMPARISONS_RIDGE","V_SELECTED_ANCHOR_TREE",
        "RLR_GATED_ANCHOR_TREE")
met<-list();prd<-list();sel<-list();man<-list()
ni<-0L;ri<-0L;ki<-0L
for(ds in SETS){
  src<-file.path("results/rlr_pair_specific_evidence_v1/source_data",
                 paste0(ds,".csv"))
  aa_need(file.exists(src),paste("source dataset missing",src))
  dat<-utils::read.csv(src,check.names=FALSE,stringsAsFactors=FALSE)
  X<-as.data.frame(dat[,setdiff(names(dat),c("sample_id","class")),
                         drop=FALSE],check.names=FALSE)
  rownames(X)<-dat$sample_id
  y<-factor(dat$class)
  aa_need(nrow(X)==length(y)&&length(levels(y))==2L&&
          all(vapply(X,is.numeric,logical(1))),"bad original cohort")
  folds<-make_repeated_stratified_folds(y,n_folds=3L,n_repeats=1L,seed=42L)
  aa_need(length(folds)==3L,"outer fold schedule invalid")
  for(f in seq_along(folds)){
    fold<-folds[[f]]
    tr<-fold$train_idx;te<-fold$test_idx
    aa_need(!length(intersect(tr,te)) &&
              identical(sort(c(tr,te)),seq_len(nrow(X))),
            "outer D/E overlap")
    T<-X[tr,,drop=FALSE];E<-X[te,,drop=FALSE]
    yt<-factor(as.character(y[tr]),levels=levels(y))
    ye<-factor(as.character(y[te]),levels=levels(y))
    result<-aa_eval(T,yt,E,ye)
    ps<-result$probs
    choose<-result$selection
    base<-ps$ORIGINAL_RLR
    ki<-ki+1L
    sel[[ki]]<-data.frame(dataset=ds,fold=f-1L,
       inner_mode=choose$mode,
       V_accepted=as.integer(choose$accepted),
       V_logloss_gain=choose$V_delta,
       V_BA_delta=choose$V_ba_delta,
       root_gene=choose$root_gene,reference_genes=choose$anchors,
       realized_anchors=if(nchar(choose$anchors)>0)
                             length(strsplit(choose$anchors,"|",fixed=TRUE)[[1L]])
                        else 0L,
       child_nodes=choose$children,
       reason=choose$reason,
       evidence_stage="EXPOSED_ORIGINAL_RLR_PAPER_DATA_DEVELOPMENT_ONLY")
    man[[ki]]<-data.frame(dataset=ds,fold=f-1L,
      n_train=length(tr),n_E=length(te),n_features=ncol(X),
      training_ids=paste(rownames(T),collapse="|"),
      E_ids=paste(rownames(E),collapse="|"),
      original_parent_K=20L,
      evidence_stage="EXPOSED_ORIGINAL_RLR_PAPER_DATA_DEVELOPMENT_ONLY")
    for(arm in ARMS){
      p<-ps[[arm]]
      aa_need(length(p)==nrow(E)&&all(is.finite(p))&&
                 all(p>=0&p<=1),"nonfinite heldout probability")
      ni<-ni+1L
      met[[ni]]<-data.frame(dataset=ds,fold=f-1L,method=arm,
         n_train=length(tr),n_E=length(te),
         accepted=as.integer(choose$accepted),
         logloss=rg_ll(result$truth,p),
         brier=mean((result$truth-p)^2),
         ba=rg_ba(result$truth,p),
         delta_logloss=rg_ll(result$truth,p)-rg_ll(result$truth,base),
         evidence_stage="EXPOSED_ORIGINAL_RLR_PAPER_DATA_DEVELOPMENT_ONLY")
      for(j in seq_along(p)){
        ri<-ri+1L
        prd[[ri]]<-data.frame(
            dataset=ds,fold=f-1L,method=arm,
            sample_id=rownames(E)[j],
            truth=result$truth[j],prob=p[j],
            original_prob=base[j],stringsAsFactors=FALSE)
      }
    }
    message(sprintf("ANCHOR_TREE_CASE %s fold=%d mode=%s Vaccept=%s children=%d Vgain=%s E_stump=%+.6f E_tree=%+.6f E_selected=%+.6f",
      ds,f-1L,choose$mode,choose$accepted,choose$children,
      format(choose$V_delta,digits=4),
      rg_ll(result$truth,ps$ANCHOR_STUMP)-rg_ll(result$truth,base),
      rg_ll(result$truth,ps$ANCHOR_TREE_D2)-rg_ll(result$truth,base),
      rg_ll(result$truth,ps$V_SELECTED_ANCHOR_TREE)-rg_ll(result$truth,base)))
  }
}
m<-do.call(rbind,met)
p<-do.call(rbind,prd)
s<-do.call(rbind,sel)
v<-do.call(rbind,man)
aa_need(nrow(m)==36L && nrow(p)==1236L && nrow(s)==6L && nrow(v)==6L,
        "study inventory failed")
write.csv(m,file.path(OUT,"method_metrics.csv"),row.names=FALSE)
write.csv(p,file.path(OUT,"patient_predictions.csv"),row.names=FALSE)
write.csv(s,file.path(OUT,"selection_audit.csv"),row.names=FALSE)
write.csv(v,file.path(OUT,"split_manifest.csv"),row.names=FALSE)
for(a in ARMS){
  z<-m[m$method==a,,drop=FALSE]
  message(sprintf("ANCHOR_TREE_SUMMARY %s mean_deltaLL=%+.7f improved=%d/6",
     a,mean(z$delta_logloss),sum(z$delta_logloss<(-1e-9))))
}
cat("ANCHOR_TREE_6_E_FOLDS=PASS\n")
cat("ANCHOR_TREE_36_METHOD_FOLD_METRICS=PASS\n")
cat("ANCHOR_TREE_1236_PATIENT_E_PROBS=PASS\n")
cat("ANCHOR_TREE_PUBLICATION_DATA_DEVELOPMENT_ONLY=PASS\n")
