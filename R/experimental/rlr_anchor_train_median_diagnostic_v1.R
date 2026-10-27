# Read-only, OUTER-TRAINING-ONLY diagnostic.
# Is the majority of anchors essentially a single median reference gene?
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_anchor_majority_tree_v1.R")
out<-"results/rlr_anchor_median_diagnostic_v1"
aa_need(dir.exists(out)&&file.exists(file.path(out,"SOURCE_SHA.txt")),
        "expected isolated prepared output")
res<-list()
for(ds in c("chowdary-2006","singh-2002")) {
  file<-file.path("results/rlr_pair_specific_evidence_v1/source_data",
                  paste0(ds,".csv"))
  dat<-read.csv(file,check.names=FALSE)
  X<-as.data.frame(dat[,setdiff(names(dat),c("sample_id","class")),
           drop=FALSE],check.names=FALSE)
  rownames(X)<-dat$sample_id
  y<-factor(dat$class)
  folds<-make_repeated_stratified_folds(y,3L,1L,42L)
  for(f in seq_along(folds)){
    ids<-folds[[f]]$train_idx
    T<-X[ids,,drop=FALSE]
    yt<-factor(as.character(y[ids]),levels=levels(y))
    pool<-compute_top_pairs_switchbox(T,yt,top_n=1000L,buffer_mult=5L)
    R<-aa_ranks(T)
    best<-aa_find_node(R,yt,pool,min_leaf=8L,max_hubs=14L)
    if(is.null(best))next
    root<-best$node
    a<-root$anchors
    t<-root$target
    m_idx<-apply(R[,a,drop=FALSE],1L,function(z){
      order(z,names(z))[ (length(z)+1L)%/%2L]
    })
    mid_gene<-a[as.integer(m_idx)]
    med_tab<-sort(table(mid_gene),decreasing=TRUE)
    modal<-names(med_tab)[1L]
    modal_freq<-max(med_tab)/length(mid_gene)
    panel<-aa_call(R,root)
    single<-as.integer(R[,t]>R[,modal])
    agree<-mean(panel==single)
    stats<-aa_rank_stats(R,yt)
    target_drift<-unname(stats$shift[t])
    a_shift<-unname(stats$shift[a])
    stable_half<-sum(a_shift < .5*target_drift)
    order_agreement<-c()
    if(length(a)>1L)for(i in 1:(length(a)-1L))
      for(j in (i+1L):length(a)){
        v<-as.numeric(R[,a[i]]>R[,a[j]])
        order_agreement<-c(order_agreement,max(mean(v),1-mean(v)))
      }
    x<-data.frame(dataset=ds,fold=f-1L,
       n_train=length(ids),target=t,
       anchors=paste(a,collapse="|"),
       anchor_count=length(a),unique_median_anchors=length(med_tab),
       modal_median_gene=modal,
       modal_median_fraction=modal_freq,
       agreement_with_single_modal_reference=agree,
       target_normalized_rank_drift=target_drift,
       median_reference_rank_drift=median(a_shift),
       reference_drift_lt_half_target=stable_half,
       mean_pairwise_anchor_order_stability=mean(order_agreement),
       stage="TRAIN_ONLY_EXPOSED_ORIGINAL_RLR_DATA")
    res[[length(res)+1L]]<-x
    message(sprintf("ANCHOR_MEDIAN_TRAIN %s fold=%d anchors=%d median_unique=%d modal_frac=%.3f stump_agreement=%.3f target_shift=%.4f ref_med_shift=%.4f half_target_refs=%d",
       ds,f-1L,length(a),length(med_tab),modal_freq,agree,
       target_drift,median(a_shift),stable_half))
  }
}
result<-do.call(rbind,res)
aa_need(nrow(result)==6L,"all six training folds required")
write.csv(result,file.path(out,"training_only_anchor_diagnostics.csv"),
          row.names=FALSE)
message(sprintf("ANCHOR_MEDIAN_SUMMARY mean_majority_vs_modal_single_agreement=%.4f",
                mean(result$agreement_with_single_modal_reference)))
message(sprintf("ANCHOR_MEDIAN_SUMMARY mean_distinct_patient_median_genes=%.4f",
                mean(result$unique_median_anchors)))
cat("ANCHOR_MEDIAN_6_ORIGINAL_PAPER_TRAIN_FOLDS_DIAGNOSED=PASS\n")
cat("ANCHOR_MEDIAN_ALL_NEW_E_PATIENT_LABELS_UNREAD=PASS\n")
