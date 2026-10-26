# RLR K50 GATED RESIDUAL TREE V1: frozen original K50 logit + branch corrections.
# Strict original paper cohort names only; no new cohorts, no manuscript edit.
# All candidate ranking D-only, model selection inner V-only, outer E untouched.
# Baseline literal original published RLR K10/K50, TSP1 and kTSP9.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_arlet_evidence_tree_v1.R") # proxy only
source("R/experimental/rlr_k50_gated_residual_tree_v1.R")
DS<-Sys.getenv("RLR_DATASET_ID",unset="")
SETS<-c("chen-2002","gordon-2002","yeoh-2002-v1")
rgt_need(DS %in% SETS,"dataset outside three predeclared published cohorts")
base<-file.path("results/rlr_arlet_v1",DS)
source_file<-file.path(base,"source_data",paste0(DS,".csv"))
rgt_need(file.exists(source_file)&&
   file.exists(file.path(base,"input_source_manifest.json")),
   "no public source provenance manifest")
z<-utils::read.csv(source_file,check.names=FALSE,
                    stringsAsFactors=FALSE)
X<-as.data.frame(z[,setdiff(names(z),c("class","sample_id")),drop=FALSE],
                 check.names=FALSE)
rownames(X)<-z$sample_id
y<-factor(z$class)
rgt_need(length(levels(y))==2L&&nrow(X)>=45L,"wrong cohort class/size")
folds<-make_repeated_stratified_folds(y,n_folds=3L,
                                      n_repeats=1L,seed=42L)
dirr<-file.path("results/rlr_rgt50_v1",DS)
dir.create(dirr,recursive=TRUE,showWarnings=FALSE)
metrics<-list();patient<-list();splits<-list();structure<-list()
nc<-0L;np<-0L
for(f in folds){
 tr<-f$train_idx;te<-f$test_idx
 rgt_need(length(intersect(tr,te))==0L&&
           identical(sort(c(tr,te)),seq_len(nrow(X))),
          "outer training E overlap")
 T<-X[tr,,drop=FALSE];E<-X[te,,drop=FALSE]
 yt<-y[tr];ye<-y[te]
 ey<-as.integer(ye==levels(y)[2L])
 ty<-as.integer(yt==levels(y)[2L])
 splitdir<-file.path(dirr,"original_cache",f$split_id)
 original_pool<-compute_top_pairs_switchbox(T,yt,top_n=1000L,
                                            buffer_mult=5L)
 originals<-list()
 for(K in c(10L,50L)){
   key<-paste0("ORIGINAL_RLR_K",K)
   run<-file.path(splitdir,key)
   fit_rlr_train_test(T,yt,E,ye,pairs_pool=original_pool,
                 K=K,m=10L,d=2L,out_dir=run,
                 alpha=.5,lambda_choice="lambda.min",seed=42L)
   raw<-utils::read.csv(file.path(run,"predictions_test.csv"))
   rgt_need(identical(as.character(raw$sample_id),rownames(E))&&
     identical(as.character(raw$y_true),as.character(ye)),
     paste("original RLR exact E pairing",K))
   originals[[key]]<-as.numeric(raw$prob)
 }
 fit_tsp1_manual_train_test(T,yt,E,ye,original_pool,
           file.path(splitdir,"TSP_1"))
 fit_ktsp_train_test(T,yt,E,ye,k=9L,
           out_dir=file.path(splitdir,"kTSP_9"),seed=42L)
 for(nm in c("TSP_1","kTSP_9")){
   raw<-utils::read.csv(file.path(splitdir,nm,
                                  "predictions_test.csv"))
   rgt_need(identical(as.character(raw$sample_id),rownames(E))&&
      identical(as.character(raw$y_true),as.character(ye)),
      paste("TSP E patient/class mismatch",nm))
   originals[[nm]]<-as.numeric(raw$prob)
 }
 # Tree Gini comparison is a GREEDY PROXY, NOT full CUDA TSPDT.
 orig_tr<-encode_pair_matrix(T,original_pool[seq_len(min(150L,nrow(original_pool))),,drop=FALSE],m=10L)
 orig_te<-encode_pair_matrix(E,original_pool[seq_len(min(150L,nrow(original_pool))),,drop=FALSE],m=10L)
 originals$TSPDT_GREEDY_PROXY<-arlet_tspdt_proxy(
   orig_tr,ty,orig_te,depth=2L,min_leaf=5L)
 # Nested 4fold D/V split strictly WITHIN this outer training split.
 inside<-make_repeated_stratified_folds(yt,n_folds=4L,
                     n_repeats=1L,seed=100L+f$fold)[[1L]]
 D<-T[inside$train_idx,,drop=FALSE]
 V<-T[inside$test_idx,,drop=FALSE]
 d_y<-as.integer(yt[inside$train_idx]==levels(y)[2L])
 v_y<-as.integer(yt[inside$test_idx]==levels(y)[2L])
 rgt_need(!length(intersect(inside$train_idx,inside$test_idx))&&
     all(c(0,1) %in% unique(d_y)),
     "inner D/V leakage or no binary class")
 proposal<-rgt_select(D,yt[inside$train_idx],
   V,yt[inside$test_idx],seed=42L+f$fold,
   path=file.path(splitdir,"inner_D_parent_K50"))
 original_cv<-readRDS(file.path(splitdir,"ORIGINAL_RLR_K50","cvfit.rds"))
 pred_tree<-rgt_refit(T,yt,E,original_pool,proposal,
       originals$ORIGINAL_RLR_K50,original_cv)
 originals$RGT50_GATED_RESIDUAL_TREE<-pred_tree$p
 if(!pred_tree$used_tree)
   rgt_need(max(abs(pred_tree$p-originals$ORIGINAL_RLR_K50))<1e-12,
      "K50 exact original fallback parity failed")
 splits[[length(splits)+1L]]<-data.frame(dataset=DS,
   split_id=f$split_id,n_train=nrow(T),n_E=nrow(E),
   E_ids=paste(rownames(E),collapse="|"))
 structure[[f$split_id]]<-list(
   tree_accepted=pred_tree$used_tree,
   inner_D_n=nrow(D),inner_V_n=nrow(V),
   source_class_counts=as.list(table(yt)),
   V_baseline_logloss=proposal$baseline_ll,
   V_best_logloss=proposal$validation_ll,
   V_improvement=proposal$gain,
   V_gate_accepted=isTRUE(proposal$accept),
   n_root_candidates=proposal$root_nominated,
   main_selected_pairs=apply_disjointness(original_pool,K=50L,d=2L)$ordered_pair,
   root_pair=if(pred_tree$used_tree)
      proposal$candidates$ordered_pair[[proposal$root]] else NULL,
   root_ternary_zero_train=if(pred_tree$used_tree)
      mean(encode_pair_matrix(T,proposal$candidates[proposal$root,,drop=FALSE],m=10L)==0) else NULL,
   conditional_pairs=if(pred_tree$used_tree)
      lapply(proposal$branches,function(js)
         proposal$candidates$ordered_pair[js]) else NULL,
   n_nonzero_correction=pred_tree$n_correction_nonzero,
   path_extra_mean=pred_tree$path_extra_mean,
   max_probability_shift=pred_tree$correction_delta_max,
   V_BA_change=proposal$BAchange)
 for(nm in names(originals)){
   p<-originals[[nm]]
   rgt_need(length(p)==nrow(E)&&all(is.finite(p))&&
            all(p>=0&p<=1),paste("invalid E score",nm))
   nc<-nc+1L
   metrics[[nc]]<-data.frame(
      dataset=DS,split_id=f$split_id,method=nm,
      n_train=nrow(T),n_E=nrow(E),
      logloss=rgt_ll(ey,p),brier=mean((p-ey)^2),
      BA=rgt_ba(ey,p),macro_F1=rgt_f1(ey,p),
      tree_accepted=if(nm=="RGT50_GATED_RESIDUAL_TREE")pred_tree$used_tree else NA,
      global_n_relations=if(nm=="RGT50_GATED_RESIDUAL_TREE")
              pred_tree$global_extra+nrow(apply_disjointness(original_pool,K=50L,d=2L)) else NA_integer_,
      branch_n_extra=if(nm=="RGT50_GATED_RESIDUAL_TREE")
              pred_tree$global_extra else NA_integer_,
      mean_additional_path_evaluated=if(nm=="RGT50_GATED_RESIDUAL_TREE")
              pred_tree$path_extra_mean else NA_real_,
      root_relation_active_rate=if(nm=="RGT50_GATED_RESIDUAL_TREE")
              pred_tree$root_active_rate else NA_real_,
      scientific_stage="RGT50_ORIGINAL_PAPER_RECONSTRUCTED_SOURCE_ARFF_NOT_VERIFIED")
   for(i in seq_len(nrow(E))){
     np<-np+1L
     patient[[np]]<-data.frame(dataset=DS,split_id=f$split_id,
          method=nm,sample_id=rownames(E)[[i]],
          truth=ey[[i]],prob=p[[i]])
   }
   message(sprintf("RGT50_CASE %s %s %-21s BA=%.4f F1=%.4f LL=%.4f",
    DS,f$split_id,nm,rgt_ba(ey,p),rgt_f1(ey,p),rgt_ll(ey,p)))
 }
 message(sprintf("RGT50_OUTER_DONE dataset=%s split=%s tree=%s inner_val_gain=%+.5f proposed_root_n=%d",
  DS,f$split_id,ifelse(pred_tree$used_tree,"YES","FALLBACK"),
  proposal$gain,proposal$root_nominated))
 unlink(splitdir,recursive=TRUE,force=TRUE)
}
df<-do.call(rbind,metrics)
pr<-do.call(rbind,patient)
stopifnot(nrow(df)==18L,nrow(pr)==6L*nrow(X),
  all(df$scientific_stage=="RGT50_ORIGINAL_PAPER_RECONSTRUCTED_SOURCE_ARFF_NOT_VERIFIED"))
write.csv(df,file.path(dirr,"method_metrics.csv"),row.names=FALSE)
write.csv(pr,file.path(dirr,"patient_predictions.csv"),row.names=FALSE)
write.csv(do.call(rbind,splits),file.path(dirr,"split_manifest.csv"),row.names=FALSE)
jsonlite::write_json(structure,file.path(dirr,"tree_structure.json"),
 auto_unbox=TRUE,pretty=TRUE)
for(nm in unique(df$method)){
 w<-df[df$method==nm,,drop=FALSE]
 cat(sprintf("RGT50_SUMMARY %s %s BA=%.5f macroF1=%.5f logloss=%.5f\n",
      DS,nm,mean(w$BA),mean(w$macro_F1),mean(w$logloss)))
}
cat(sprintf("RGT50_3_FOLD_6_ARM_DATASET=%s_PASS\n",DS))
cat("RGT50_FROZEN_ORIGINAL_RLR_TSP_KTSP=PASS\n")
cat("RGT50_OUTER_E_NEVER_USED_FOR_TRAINING_SELECTION=PASS\n")
cat("RGT50_GREEDY_TSPDT_PROXY_NOT_ORIGINAL_CUDA=YES\n")
cat("RGT50_PROCESSED_ARFF_PARITY=NOT_VERIFIED\n")
