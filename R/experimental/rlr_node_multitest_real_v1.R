# Reused public CompCancer biomedical pilot: original RLR vs two leaf sizes.
# Locked split, no outer E labels used for fitting, nomination or V confirmation.
# Development/falsification run, never independent manuscript confirmation.
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v4_rank_repair.R")
source("R/experimental/crit_v7_original_rlr_two_leaf.R")
source("R/experimental/rlr_node_multitest_v1.R")

OUT <- "results/rlr_node_multitest_v1"
node_need(dir.exists(OUT)&&file.exists(file.path(OUT,"SOURCE_SHA.txt"))&&
           !file.exists(file.path(OUT,"method_metrics.csv")),
            "missing destination or previous batch exists")

databases <- c("chowdary-2006","singh-2002")
methods <- c("ORIGINAL_RLR","SINGLE_REL","RLR_NODE_2","RLR_NODE_5")
N_OUTER <- 3L
metric_rows <- list()
pred_rows <- list()
run_rows <- list()
idx <- 0L
pi <- 0L
ri <- 0L
for(nm in databases) {
  path<-file.path(OUT,"source_data",paste0(nm,".csv"))
  sparse_need(file.exists(path),"missing sourced real matrix")
  tab<-utils::read.csv(path,check.names=FALSE,stringsAsFactors=FALSE)
  sparse_need(all(c("sample_id","class") %in% names(tab)),
              "missing phenotype/identifier column")
  X<-as.data.frame(tab[,setdiff(names(tab),c("sample_id","class")),drop=FALSE],
                   check.names=FALSE)
  rownames(X)<-tab$sample_id
  y<-factor(tab[["class"]])
  sparse_need(all(vapply(X,function(z)is.numeric(z)&&all(is.finite(z)),
                          logical(1))) && length(levels(y))==2L &&
              nrow(X)>=100L && ncol(X)>=100L,
              "real-data public matrix schema wrong")
  folds<-make_repeated_stratified_folds(
       y,n_folds=3L,n_repeats=1L,seed=42L)
  sparse_need(length(folds)==3L,
              "outer fold partition count changed")
  for(fi in seq_len(N_OUTER)) {
    fold<-folds[[fi]]
    tr<-fold$train_idx
    te<-fold$test_idx
    sparse_need(!length(intersect(tr,te)) &&
                identical(sort(c(tr,te)),seq_len(nrow(X))),
                "outer E contamination or coverage")
    rid<-paste0(nm,":outer_",fi-1L)
    T<-X[tr,,drop=FALSE]
    E<-X[te,,drop=FALSE]
    message("NODE_REAL_START ",rid,
            " train=",length(tr)," test=",length(te),
            " classes=",paste(table(y[tr]),collapse="/"))
    fit<-tryCatch(node_run(T,y[tr],E),
                  error=function(e) {
       stop("NODE_REAL_FIT_FAILED ",rid,": ",
            conditionMessage(e),call.=FALSE)
    })
    baseline<-fit$preds$ORIGINAL_RLR$prob
    sparse_need(length(baseline)==nrow(E) && 
                all(is.finite(baseline)),
                "original RLR test vector invalid")
    for(method in methods) {
      pred<-fit$preds[[method]]
      sparse_need(!is.null(pred),"missing predeclared arm")
      y01<-as.numeric(y[te]==levels(y)[2L])
      ll<-node_ll(y01,pred$prob)
      oldll<-node_ll(y01,baseline)
      pi<-pi+1L
      metric_rows[[pi]]<-data.frame(
        dataset=nm,outer_fold=fi-1L,method=method,
        n_train=length(tr),n_E=length(te),
        accepted=if(method=="ORIGINAL_RLR") NA_integer_ else
                 as.integer(isTRUE(pred$accepted)),
        E_logloss=mean(ll),E_brier=mean((y01-pred$prob)^2),
        delta_E_logloss=mean(ll-oldll),
        V_mean_gain=pred$v_mean,V_lower_bound=pred$v_lower,
        reason=pred$reason,
        node_relations=pred$node_relations,
        nonzero_node_coeff=pred$nonzero_node_coeff,
        node_gates=pred$gate,
        D_left=pred$D_left,D_right=pred$D_right,
        scientific_stage="EXPOSED_BIOMEDICAL_DEV_ONLY",
        stringsAsFactors=FALSE)
      for(i in seq_along(te)) {
        idx<-idx+1L
        pred_rows[[idx]]<-data.frame(
          dataset=nm,outer_fold=fi-1L,method=method,
          sample_id=rownames(E)[[i]],truth=y01[[i]],
          original_prob=baseline[[i]],method_prob=pred$prob[[i]],
          stringsAsFactors=FALSE)
      }
      if(method!="ORIGINAL_RLR"&&!isTRUE(pred$accepted))
        sparse_need(identical(baseline,pred$prob),
                    "unaccepted model fails exact original fallback")
    }
    ri<-ri+1L
    run_rows[[ri]]<-data.frame(
      dataset=nm,outer_fold=fi-1L,
      n_train=length(tr),n_E=length(te),
      source_gene_count=ncol(X),
      D_n=length(fit$meta$D_ids),
      V_n=length(fit$meta$V_ids),
      original_model_identity="FROZEN_ORIGINAL_RLR_K10_M10_D2",
      n_predeclared_methods=length(methods),
      outer_test_ids_sha256=NA_character_,
      stringsAsFactors=FALSE)
    message("NODE_REAL_COMPLETE ",rid)
  }
}
allm<-do.call(rbind,metric_rows)
allp<-do.call(rbind,pred_rows)
allr<-do.call(rbind,run_rows)
sparse_need(nrow(allm)==length(databases)*N_OUTER*length(methods) &&
            nrow(allr)==length(databases)*N_OUTER &&
            all(allm$method %in% methods),
            "pilot inventory incomplete")
write.csv(allm,file.path(OUT,"method_metrics.csv"),row.names=FALSE)
write.csv(allp,file.path(OUT,"patient_predictions.csv"),row.names=FALSE)
write.csv(allr,file.path(OUT,"run_manifest.csv"),row.names=FALSE)
cat("=== REAL RLR NODE MULTITEST PILOT — REUSED DEVELOPMENT ONLY ===\n")
print(allm[,c("dataset","outer_fold","method","accepted",
              "node_relations","nonzero_node_coeff",
              "V_mean_gain","V_lower_bound",
              "E_logloss","delta_E_logloss","E_brier")],row.names=FALSE)
cat("NODE_REAL_2_DATASETS_X_3_OUTER_FOLDS=PASS\n")
cat("NODE_REAL_4_METHODS_X_6_FOLDS=PASS\n")
cat("NODE_REAL_ORIGINAL_RLR_PREDICTION_VECTORS_ARCHIVED=PASS\n")
cat("NODE_REAL_EXACT_FALLBACK_VERIFIED=PASS\n")
cat("NODE_REAL_SCIENTIFIC_STAGE=EXPOSED_BIOMEDICAL_DEVELOPMENT_ONLY\n")
