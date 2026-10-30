# ORIGINAL frozen DeCAF versus original RLR20 on same CAF subtype labels.
# Authors' 18 classifier genes and independent nontraining PDAC cohorts.
# Frozen before original head-to-head outcome; exploratory; DO NOT confuse with
# original DeCAF generic class/retraining on unrelated CompCancer labels.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_rank_gate_core_v1.R")
require_that<-function(x,msg) if(!isTRUE(x))stop("ORIGINAL_DECAF_HEADTOHEAD_STOP: ",msg,call.=FALSE)
args<-commandArgs(trailingOnly=TRUE)
require_that(length(args)==1L,"usage: Rscript ... <directory containing author inputs>")
root<-normalizePath(args[1],mustWork=TRUE)
data_file<-file.path(root,"combined_classifier_gene_exp_DeCAF.csv")
labels_file<-file.path(root,"combined_subtype_DeCAF.csv")
calls_file<-file.path(root,"combined_classifier_calls_DeCAF.csv")
model_file<-file.path(root,"decaf_classifier.Rdata")
functions_file<-file.path(root,"decaf_functions.R")
require_that(all(file.exists(c(data_file,labels_file,calls_file,model_file,functions_file))),
             "original author inputs missing")
out<-file.path("results","rlr-original-decaf-pdac-headtohead-v1")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
exp<-read.csv(data_file,row.names=1,check.names=FALSE)
meta<-read.csv(labels_file,check.names=FALSE,stringsAsFactors=FALSE)
calls<-read.csv(calls_file,row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
require_that(nrow(exp)>=18&&ncol(exp)>=100&&
             all(vapply(exp,is.numeric,logical(1)))&&
             all(is.finite(as.matrix(exp))),"unexpected classifier gene matrix")
require_that(all(c("ID","dataSet","SCISSORS_CAF","CAF_whitelist") %in% names(meta)),
             "unsupported original subtype metadata schema")
# Author's combined expression and subtype CSV are strictly position-aligned.
# They contain 52 intentionally repeated PACA_AU biological specimen IDs
# (two independent assay platforms), so matching using bare IDs is unsafe.
# Require EXACT column-by-column author IDs rather than deduplicate silently.
norm_id<-function(x)make.names(as.character(x),unique=FALSE)
require_that(ncol(exp)==nrow(meta)&&
             identical(as.character(colnames(exp)),as.character(meta$ID)),
             "source combined expression/subtype CSV row order differs")
require_that(nrow(calls)==nrow(meta),"reference prediction row count mismatch")
aligned_reference_ids<-sum(norm_id(rownames(calls))==norm_id(meta$ID))
require_that(aligned_reference_ids>=0.90*nrow(meta),
             "author reference call and annotated sample order incompatible")
duplicate_original_ids<-sum(duplicated(as.character(meta$ID)))
duplicate_by_assay<-which(meta$dataSet=="PACA_AU_seq"&
            meta$ID %in% meta$ID[meta$dataSet=="PACA_AU_array"])
require_that(duplicate_original_ids==length(duplicate_by_assay),
             "unrecognized repeated biological specimen outside PACA_AU platforms")
cat(sprintf("ORIGINAL_DECAF_REPEATED_BIOLOGICAL_SPECIMENS=%d (PACA_AU_seq/array)\n",
       duplicate_original_ids))
# Official apply_decaf() returns a data.frame indexed by sample names and
# cannot accept duplicate colnames for PACA_AU's same specimen on two platforms.
# Disambiguate names ONLY AFTER verifying exact original author CSV row alignment;
# no sample order, expression value or label is changed.
source_ids<-colnames(exp)
colnames(exp)<-make.unique(source_ids,sep=".platform_repeat_")
require_that(!anyDuplicated(colnames(exp)),"could not disambiguate assay columns")
de<-new.env(parent=globalenv())
load(model_file,envir=de)
sys.source(functions_file,envir=de)
require_that(exists("decaf_classifier",de)&&exists("apply_decaf",de),
             "original official DeCAF object/function missing")
classifier<-de$decaf_classifier$classifier2
require_that(all(unique(c(classifier$TSPs)) %in% rownames(exp)),
             "not all official DeCAF genes available")
pred<-de$apply_decaf(data=exp,classifier=classifier)
require_that(is.data.frame(pred)&&"DeCAF_prob"%in%names(pred)&&
             nrow(pred)==ncol(exp)&&all(is.finite(pred$DeCAF_prob))&&
             all(pred$DeCAF_prob>=0&pred$DeCAF_prob<=1),
             "official classifier produced invalid predictions")
official<-as.numeric(pred$DeCAF_prob)
# Authors' published per-sample predictions provide a model-version integrity check.
# Author reference calls align by original combined table ROW POSITION.
# The duplicated PACA_AU IDs are disambiguated differently in call rownames,
# so never match references with just specimen IDs.
overlap<-seq_len(ncol(exp))
difference<-abs(official-as.numeric(calls$DeCAF_prob))
replay<-data.frame(n_samples=ncol(exp),published_probability_overlap=length(overlap),
 mean_abs_prediction_difference=if(length(overlap))mean(difference[overlap],na.rm=TRUE) else NA_real_,
 max_abs_prediction_difference=if(length(overlap))max(difference[overlap],na.rm=TRUE) else NA_real_,
 exact_probability_match_rate=if(length(overlap))mean(difference[overlap]<1e-7,na.rm=TRUE) else NA_real_)
write.csv(replay,file.path(out,"original_decaf_reference_replay.csv"),row.names=FALSE)
cat(sprintf("ORIGINAL_DECAF_REFERENCE_PROBABILITIES_MATCHED=%d/%d MAE=%.9f max=%.9f\n",
 length(overlap),ncol(exp),replay$mean_abs_prediction_difference,
 replay$max_abs_prediction_difference))
label<-as.character(meta$SCISSORS_CAF)
training_studies<-c("CPTAC","Dijk","Moffitt_GEO_array","TCGA_PAAD")
eligible<-label%in%c("restCAF","proCAF")
train_idx<-which(meta$dataSet%in%training_studies& eligible &
                 !is.na(meta$CAF_whitelist)&meta$CAF_whitelist==TRUE)
require_that(length(train_idx)>=30L&&
             length(unique(label[train_idx]))==2L,
             "inadequate original DeCAF training cohort samples")
train_X<-as.data.frame(t(as.matrix(exp[,train_idx,drop=FALSE])),check.names=FALSE)
rownames(train_X)<-colnames(exp)[train_idx]
train_y<-factor(label[train_idx],levels=c("restCAF","proCAF"))
# Training-fold-only supervised switchBox ranking; same 18 author markers.
pool<-compute_top_pairs_switchbox(train_X,train_y,top_n=1000L,buffer_mult=5L)
test_idx<-setdiff(which(!(meta$dataSet%in%training_studies)&eligible),
                    duplicate_by_assay)
require_that(length(test_idx)>30L&&
             !anyDuplicated(as.character(meta$ID[test_idx])),
             "heldout biological specimen identity overlap after deduplication")
test_X<-as.data.frame(t(as.matrix(exp[,test_idx,drop=FALSE])),check.names=FALSE)
rownames(test_X)<-colnames(exp)[test_idx]
test_y<-factor(label[test_idx],levels=c("restCAF","proCAF"))
rlr<-rg_parent(train_X,train_y,test_X,pool,K=20L)
pR<-as.numeric(rlr$test);pD<-official[test_idx]
require_that(all(is.finite(pR))&&all(pR>=0&pR<=1)&&
             length(pR)==length(test_idx),"RLR test probability invalid")
predictions<-data.frame(sample_id=colnames(exp)[test_idx],
  study=meta$dataSet[test_idx],
  validation_unit=ifelse(meta$dataSet[test_idx]%in%
                       c("PACA_AU_seq","PACA_AU_array"),
                       "PACA_AU",meta$dataSet[test_idx]),
  true_caf_label=as.character(test_y),
  truth=as.integer(test_y=="proCAF"),
  rlr_probability=pR,original_decaf_probability=pD,
  model="ORIGINAL_DECAF_VS_ORIGINAL_RLR20_FIXED_18_AUTHOR_GENES",
  stringsAsFactors=FALSE)
write.csv(predictions,file.path(out,"external_patient_predictions.csv"),row.names=FALSE)
write.csv(data.frame(training_sample_id=colnames(exp)[train_idx],
              training_study=meta$dataSet[train_idx],
              training_truth=label[train_idx]),
     file.path(out,"original_author_training_cohort_manifest.csv"),row.names=FALSE)
write.csv(data.frame(original_combined_rows=nrow(meta),
         repeated_biological_specimen_ids=duplicate_original_ids,
         repeated_PACA_AU_seq_excluded=length(duplicate_by_assay),
         reference_position_aligned_ids=aligned_reference_ids,
         holdout_biological_specimen_unique=!anyDuplicated(meta$ID[test_idx])),
    file.path(out,"source_identity_audit.csv"),row.names=FALSE)
score<-function(y,p){
 q<-pmin(1-1e-15,pmax(1e-15,p))
 h<-as.integer(p>=.5)
 tp<-sum(y==1&h==1);fp<-sum(y==0&h==1);fn<-sum(y==1&h==0)
 n1<-sum(y==1);n0<-sum(y==0)
 auc<-if(n1&&n0)(sum(rank(p,ties.method="average")[y==1])-n1*(n1+1)/2)/(n1*n0) else NA_real_
 c(logloss=mean(-y*log(q)-(1-y)*log1p(-q)),
   brier=mean((y-p)^2),
   auc=auc,
   positive_f1=if(2*tp+fp+fn>0)2*tp/(2*tp+fp+fn) else 0)
}
studies<-sort(unique(predictions$validation_unit))
S<-list();exclusion<-list()
for(ds in studies){
 q<-predictions[predictions$validation_unit==ds,,drop=FALSE]
 yy<-q$truth
 if(length(yy)<6L || length(unique(yy))!=2L){
   exclusion[[length(exclusion)+1L]]<-data.frame(study=ds,n=nrow(q),why="NO_BOTH_CLASSES_OR_TOO_FEW_SAMPLES")
   next
 }
 for(method in c("ORIGINAL_RLR20","ORIGINAL_DeCAF")){
   pp<-if(method=="ORIGINAL_RLR20")q$rlr_probability else q$original_decaf_probability
   sc<-score(yy,pp)
   S[[length(S)+1L]]<-data.frame(study=ds,method=method,
       n=nrow(q),n_proCAF=sum(yy==1),n_restCAF=sum(yy==0),
       logloss=sc["logloss"],brier=sc["brier"],auc=sc["auc"],positive_f1=sc["positive_f1"])
 }
}
require_that(length(S)>0L,"no eligible external evaluation")
metrics<-do.call(rbind,S)
write.csv(metrics,file.path(out,"external_study_metrics.csv"),row.names=FALSE)
if(length(exclusion))write.csv(do.call(rbind,exclusion),file.path(out,"excluded_studies.csv"),row.names=FALSE)
summary<-aggregate(cbind(logloss,brier,auc,positive_f1)~method,metrics,FUN=mean)
write.csv(summary,file.path(out,"dataset_equal_summary.csv"),row.names=FALSE)
cat("ORIGINAL_DECAF_VS_RLR_EXTERNAL_STUDIES_ELIGIBLE=",length(unique(metrics$study)),"\n",sep="")
for(j in seq_len(nrow(summary)))
 cat(sprintf("ORIGINAL_DECAF_HEADTOHEAD method=%s LL=%.6f Brier=%.6f AUC=%.6f F1=%.6f\n",
 summary$method[j],summary$logloss[j],summary$brier[j],summary$auc[j],
 summary$positive_f1[j]))
# One fixed contrast of study means. Never interpret as confirmation.
diff<-vapply(sort(unique(metrics$study)),function(ds){
 p<-metrics[metrics$study==ds,]
 p$logloss[p$method=="ORIGINAL_RLR20"]-p$logloss[p$method=="ORIGINAL_DeCAF"]
},numeric(1))
if(length(diff)<=20){
 grid<-as.matrix(expand.grid(rep(list(c(-1,1)),length(diff))))
 vals<-abs(as.vector(grid%*%diff)/length(diff))
 pval<-mean(vals>=abs(mean(diff))-1e-13)
}else{pval<-NA_real_}
write.csv(data.frame(n_studies=length(diff),mean_RLR_minus_original_DeCAF_logloss=mean(diff),
 exact_two_sided_signflip_p=pval, exploratory=TRUE,
 frozen_original_DeCAF=TRUE,RLR_panel="ORIGINAL_DECAF_18_GENES",
 independent_truth="SCISSORS_CAF_NOT_MODEL_PREDICTION",
 reference_prediction_mae=replay$mean_abs_prediction_difference),
 file.path(out,"paired_exploratory_logloss.csv"),row.names=FALSE)
cat(sprintf("ORIGINAL_DECAF_VS_RLR_EXPLORATORY_SIGNFLIP_P=%.6f\n",pval))
cat("ORIGINAL_DECAF_FROZEN_AUTHORS_PREDICTOR_VS_RLR_MATCHED_TARGET_COMPLETE=PASS\n")
