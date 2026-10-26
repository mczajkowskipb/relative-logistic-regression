# Two existing original-paper CompCancer cohorts; exploratory, NOT fresh validation.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_rank_gate_core_v1.R")

OUT<-"results/rlr_rank_gate_v1"
rg_need(dir.exists(OUT)&&file.exists(file.path(OUT,"SOURCE_SHA.txt")) &&
        !file.exists(file.path(OUT,"method_metrics.csv")),
        "result directory missing or outputs already exist")
SETS<-c("chowdary-2006","singh-2002")
ARMS<-c("ORIGINAL_RLR","GLOBAL_MAIN_CONTROL",
        "RANK_STUMP_UNGATED","RANK_STUMP_V_GATED")
metrics<-list(); predictions<-list(); manifests<-list(); selections<-list()
mi<-0L;pi<-0L;fi<-0L
for(ds in SETS){
  src<-file.path("results/rlr_pair_specific_evidence_v1/source_data",
                  paste0(ds,".csv"))
  rg_need(file.exists(src),paste("dataset missing:",src))
  dat<-utils::read.csv(src,check.names=FALSE,stringsAsFactors=FALSE)
  X<-as.data.frame(dat[,setdiff(names(dat),c("sample_id","class")),
                         drop=FALSE],check.names=FALSE)
  rownames(X)<-dat$sample_id
  y<-factor(dat$class)
  rg_need(all(vapply(X,is.numeric,logical(1)))&&
          length(levels(y))==2L,"invalid source cohort")
  folds<-make_repeated_stratified_folds(y,n_folds=3L,n_repeats=1L,seed=42L)
  rg_need(length(folds)==3L,"unexpected outer protocol")
  for(f in seq_along(folds)){
    spec<-folds[[f]]
    ti<-spec$train_idx;ei<-spec$test_idx
    rg_need(!length(intersect(ti,ei)),"outer leakage")
    T<-X[ti,,drop=FALSE]; E<-X[ei,,drop=FALSE]
    yt<-factor(as.character(y[ti]),levels=levels(y))
    ye<-factor(as.character(y[ei]),levels=levels(y))
    s<-rg_evaluate(T,yt,E,ye)
    orig<-s$probs$ORIGINAL_RLR
    opt<-s$selection
    rootlabel<-if (is.null(opt$pair_table)) "" else
      opt$pair_table$ordered_pair[opt$root]
    partnerlabel<-if (is.null(opt$pair_table)) "" else
      opt$pair_table$ordered_pair[opt$partner]
    fi<-fi+1L
    selections[[fi]]<-data.frame(dataset=ds,fold=f-1L,
      accepted=as.integer(isTRUE(opt$accepted)),
      gate=opt$gate,root_pair=rootlabel,partner_pair=partnerlabel,
      delta_V_vs_RLR=opt$delta_vs_rlr,
      delta_V_vs_global=opt$delta_vs_global,
      candidates=opt$candidate_count,
      realized_original_K=nrow(s$selected_original),
      stage="REUSED_ORIGINAL_PAPER_DATA_DEV_ONLY")
    manifests[[fi]]<-data.frame(dataset=ds,fold=f-1L,
       n_train=length(ti),n_E=length(ei),n_genes=ncol(X),
       train_ids=paste(rownames(T),collapse="|"),
       E_ids=paste(rownames(E),collapse="|"),
       stage="REUSED_ORIGINAL_PAPER_DATA_DEV_ONLY")
    for(a in ARMS){
      prob<-s$probs[[a]]
      rg_need(length(prob)==nrow(E)&&all(is.finite(prob))&&
              all(prob>=0&prob<=1),"invalid probability")
      mi<-mi+1L
      metrics[[mi]]<-data.frame(dataset=ds,fold=f-1L,method=a,
         n_train=length(ti),n_E=length(ei),
         accepted=as.integer(isTRUE(opt$accepted)),
         logloss=rg_ll(s$truth,prob),
         brier=mean((s$truth-prob)^2),
         ba=rg_ba(s$truth,prob),
         delta_logloss=rg_ll(s$truth,prob)-rg_ll(s$truth,orig),
         stage="REUSED_ORIGINAL_PAPER_DATA_DEV_ONLY")
      for(j in seq_len(nrow(E))){
        pi<-pi+1L
        predictions[[pi]]<-data.frame(dataset=ds,fold=f-1L,method=a,
          sample_id=rownames(E)[j],truth=s$truth[j],
          prob=prob[j],baseline_prob=orig[j])
      }
    }
    message(sprintf("RANK_GATE_CASE %s fold=%d accepted=%s K=%d Vgain=%s E_raw_delta=%+.6f",
      ds,f-1L,if(isTRUE(opt$accepted))"YES" else "NO",
      nrow(s$selected_original),
      format(opt$delta_vs_rlr,digits=4),
      rg_ll(s$truth,s$probs$RANK_STUMP_UNGATED)-rg_ll(s$truth,orig)))
  }
}
m<-do.call(rbind,metrics)
p<-do.call(rbind,predictions)
v<-do.call(rbind,manifests)
s<-do.call(rbind,selections)
rg_need(nrow(m)==24L&&nrow(p)==824L&&nrow(v)==6L&&nrow(s)==6L,
        "experiment inventory incomplete")
utils::write.csv(m,file.path(OUT,"method_metrics.csv"),row.names=FALSE)
utils::write.csv(p,file.path(OUT,"patient_predictions.csv"),row.names=FALSE)
utils::write.csv(v,file.path(OUT,"split_manifest.csv"),row.names=FALSE)
utils::write.csv(s,file.path(OUT,"selection_audit.csv"),row.names=FALSE)
cat("RANK_GATE_6_OUTER_PAPER_COHORT_FOLDS=PASS\n")
cat("RANK_GATE_24_METHOD_FOLD_METRICS=PASS\n")
cat("RANK_GATE_824_HELDOUT_PATIENT_PROBABILITY_ROWS=PASS\n")
cat("RANK_GATE_TEST_PREDICTIONS_ARCHIVED=PASS\n")
