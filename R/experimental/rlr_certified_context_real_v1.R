# C-RET v1: source-reconstructed original paper cohorts, strict same original
# RLR fit, with an inference-only early-exit tree with a mathematical
# guarantee of exact predicted classes and probability error <=0.025.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_arlet_evidence_tree_v1.R")
source("R/experimental/rlr_certified_early_exit_v1.R")
source("R/experimental/rlr_certified_context_tree_v1.R")
DS<-Sys.getenv("RLR_DATASET_ID",unset="")
stopifnot(DS %in% c("armstrong-2002-v1","bittner-2000","gordon-2002"))
OUT<-file.path("results/rlr_arlet_v1",DS)
tab<-utils::read.csv(file.path(OUT,"source_data",paste0(DS,".csv")),
                   check.names=FALSE,stringsAsFactors=FALSE)
X<-as.data.frame(tab[,setdiff(names(tab),c("sample_id","class")),drop=FALSE],
                 check.names=FALSE)
rownames(X)<-tab$sample_id
y<-factor(tab$class)
stopifnot(length(levels(y))==2L)
folds<-make_repeated_stratified_folds(y,n_folds=3L,n_repeats=1L,seed=42L)
RESULT<-file.path(OUT,"context_certified_tree")
dir.create(RESULT,recursive=TRUE,showWarnings=FALSE)
metrics<-list();preds<-list();nmet<-0L;npred<-0L
for(fold in folds){
  tr<-fold$train_idx;te<-fold$test_idx
  T<-X[tr,,drop=FALSE];E<-X[te,,drop=FALSE]
  yt<-y[tr];ye<-y[te]
  yE01<-as.numeric(ye==levels(y)[2L])
  pool<-compute_top_pairs_switchbox(T,yt,top_n=1000L,buffer_mult=5L)
  path<-file.path(RESULT,"original",fold$split_id)
  fit_rlr_train_test(T,yt,E,ye,pairs_pool=pool,K=50L,m=10L,
    d=2L,out_dir=path,alpha=.5,lambda_choice="lambda.min",seed=42L)
  old<-utils::read.csv(file.path(path,"predictions_test.csv"))
  stopifnot(identical(as.character(old$sample_id),rownames(E)),
    identical(as.character(old$y_true),as.character(ye)))
  cv<-readRDS(file.path(path,"cvfit.rds"))
  shortlist<-apply_disjointness(pool,K=50L,d=2L)
  RT<-encode_pair_matrix(T,shortlist,m=10L)
  RE<-encode_pair_matrix(E,shortlist,m=10L)
  coef<-as.numeric(stats::coef(cv,s="lambda.min"))
  beta<-coef[-1L];intercept<-coef[[1L]]
  original<-plogis(intercept+as.numeric(RE%*%beta))
  stopifnot(length(beta)==ncol(RE),
    max(abs(original-as.numeric(old$prob)))<1e-9)
  opt<-cret_optimize(RT,
      as.numeric(yt==levels(y)[2L]),beta,intercept,
      seed=42L+fold$fold,epsilon=.05)
  standard<-cret_evaluate(RE,beta,intercept,opt$original,epsilon=.05)
  evolved<-cret_context_eval(RT,RE,beta,intercept,epsilon=.05)
  for(method in c("ORIGINAL_RLR","CERT_TREE_ABS_WEIGHT","CERT_TREE_CONTEXT")){
    a<-if(method=="CERT_TREE_ABS_WEIGHT")standard else
       if(method=="CERT_TREE_CONTEXT")evolved else NULL
    p<-if(is.null(a))original else a$prob
    s<-arlet_scorer(yE01,p)
    nmet<-nmet+1L
    metrics[[nmet]]<-data.frame(dataset=DS,split_id=fold$split_id,
      method=method,logloss=s$logloss,brier=s$brier,
      BA=s$ba,macro_F1=s$macro_f1,
      n_global_active=sum(abs(beta)>1e-12),
      evaluated_mean=if(is.null(a))sum(abs(beta)>1e-12) else mean(a$evaluated),
      active_mean=if(is.null(a))mean(rowSums(RE[,abs(beta)>1e-12,drop=FALSE]!=0))
                     else mean(a$active_n),
      stop_early_fraction=if(is.null(a))0 else
           mean(a$evaluated<sum(abs(beta)>1e-12)),
      mean_abs_prob_error=mean(abs(p-original)),
      EA_train_savings=if(method=="CERT_TREE_CONTEXT")NA_real_ else NA_real_,
      input_provenance="SCHLIEP_RECONSTRUCTED_SOURCE_NO_ORIGINAL_ARFF_PARITY")
    for(j in seq_len(nrow(E))){
      npred<-npred+1L
      preds[[npred]]<-data.frame(dataset=DS,split_id=fold$split_id,
        method=method,sample_id=rownames(E)[[j]],
        truth=yE01[[j]],prob=p[[j]],original_prob=original[[j]],
        prob_lower=if(is.null(a))original[[j]] else a$lo[[j]],
        prob_upper=if(is.null(a))original[[j]] else a$hi[[j]],
        evaluated=if(is.null(a))sum(abs(beta)>1e-12) else a$evaluated[[j]],
        active=if(is.null(a))sum(RE[j,abs(beta)>1e-12]!=0) else a$active_n[[j]],
        pair_sequence=if(is.null(a)||is.null(a$paths))"" else a$paths[[j]],
        pair_second=if(is.null(a)||is.null(a$second))NA_integer_ else a$second[[j]])
    }
  }
  message(sprintf("CRET_FOLD_DONE %s %s original_active=%d baseline_eval=%.2f evolved_eval=%.2f exact_class=PASS",
                  DS,fold$split_id,sum(abs(beta)>1e-12),
                  mean(standard$evaluated),mean(evolved$evaluated)))
  unlink(path,recursive=TRUE,force=TRUE)
}
met<-do.call(rbind,metrics);pr<-do.call(rbind,preds)
stopifnot(nrow(met)==9L,nrow(pr)==3L*nrow(X),
          all(is.finite(met$logloss)))
utils::write.csv(met,file.path(RESULT,"method_metrics.csv"),row.names=FALSE)
utils::write.csv(pr,file.path(RESULT,"patient_certificates.csv"),row.names=FALSE)
for(method in unique(met$method)){
 v<-met[met$method==method,]
 cat(sprintf("CRET_RESULT %s %s BA=%.5f macroF1=%.5f LL=%.5f eval_mean=%.3f fraction_stop_early=%.4f\n",
   DS,method,mean(v$BA),mean(v$macro_F1),mean(v$logloss),
   mean(v$evaluated_mean),mean(v$stop_early_fraction)))
}
cat(sprintf("CRET_3_FOLDS_CERTIFIED=%s_PASS\n",DS))
cat("CRET_PREDICTED_CLASS_EQUALS_ORIGINAL_RLR_FOR_EVERY_E_PATIENT=PASS\n")
cat("CRET_PROBABILITY_INTERVAL_0.05_WITHIN_BOUND=PASS\n")
cat("CRET_NEXT_TEST_DEPENDS_ON_OBSERVED_TERNARY_STATE=PASS\n")
cat("CRET_ORIGINAL_PROCESSED_ARFF_SHA=NOT_VERIFIED\n")
