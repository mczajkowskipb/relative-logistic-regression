# C-RET L2 pilot: original RLR vs STATIC, MYOPIC and jointly
# optimized three-way decision tree; binary original paper cohorts ONLY.
# Holdout E never used for policy root/branch selection or model fitting.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_arlet_evidence_tree_v1.R")
source("R/experimental/rlr_certified_early_exit_v1.R")
source("R/experimental/rlr_certified_context_tree_v1.R")
source("R/experimental/rlr_certified_lookahead_tree_v1.R")
DS<-Sys.getenv("RLR_DATASET_ID",unset="")
stopifnot(DS %in% c("alizadeh-2000-v1","pomeroy-2002-v1"))
source_dir<-file.path("results/rlr_arlet_v1",DS)
csvname<-file.path(source_dir,"source_data",paste0(DS,".csv"))
stopifnot(file.exists(csvname),
  file.exists(file.path(source_dir,"input_source_manifest.json")))
z<-utils::read.csv(csvname,check.names=FALSE,
                    stringsAsFactors=FALSE)
X<-as.data.frame(z[,setdiff(names(z),c("sample_id","class")),
  drop=FALSE],check.names=FALSE)
rownames(X)<-z$sample_id
y<-factor(z$class)
stopifnot(nrow(X)>=30L,ncol(X)>=50L,
  length(levels(y))==2L,
  length(unique(z$sample_id))==nrow(X))
folds<-make_repeated_stratified_folds(y,n_folds=3L,n_repeats=1L,seed=42L)
RESULT<-file.path("results/rlr_cl2_v1",DS)
dir.create(RESULT,recursive=TRUE,showWarnings=FALSE)
records<-list();patient<-list();structures<-list()
nr<-0L;np<-0L
all_methods<-c("ORIGINAL_RLR","CERT_STATIC_ABS","CERT_MYOPIC_CONTEXT",
               "CERT_TREE_DEPTH2_LOOKAHEAD")
for(fold in folds){
  tr<-fold$train_idx;te<-fold$test_idx
  stopifnot(!length(intersect(tr,te)),
     identical(sort(c(tr,te)),seq_len(nrow(X))))
  T<-X[tr,,drop=FALSE];E<-X[te,,drop=FALSE]
  yt<-y[tr];ye<-y[te]
  yt01<-as.integer(yt==levels(y)[2L])
  ye01<-as.integer(ye==levels(y)[2L])
  pool<-compute_top_pairs_switchbox(T,yt,top_n=1000L,
                  buffer_mult=5L)
  tmp<-file.path(RESULT,"literal_original_rlr",fold$split_id)
  fit_rlr_train_test(T,yt,E,ye,pairs_pool=pool,
             K=50L,m=10L,d=2L,out_dir=tmp,
             alpha=.5,lambda_choice="lambda.min",seed=42L)
  actual<-utils::read.csv(file.path(tmp,"predictions_test.csv"),
                          stringsAsFactors=FALSE)
  stopifnot(identical(as.character(actual$sample_id),rownames(E)),
    identical(as.character(actual$y_true),as.character(ye)))
  cv<-readRDS(file.path(tmp,"cvfit.rds"))
  pairs<-apply_disjointness(pool,K=50L,d=2L)
  Rtrain<-encode_pair_matrix(T,pairs,m=10L)
  RE<-encode_pair_matrix(E,pairs,m=10L)
  bb<-as.numeric(stats::coef(cv,s="lambda.min"))
  beta<-bb[-1L];b0<-bb[[1L]]
  stopifnot(ncol(Rtrain)==length(beta))
  exact<-plogis(b0+as.numeric(RE%*%beta))
  stopifnot(max(abs(exact-as.numeric(actual$prob)))<1e-9)
  # Original RLR is fitted exclusively on outer TRAIN.
  # The D/V partition is used ONLY for policy optimization/acceptance.
  # For policy selection, CV model is fixed, but V also formed part of
  # its original RLR training; this V gate is not independent model eval.
  inner<-make_repeated_stratified_folds(yt,n_folds=3L,n_repeats=1L,
                           seed=101L+fold$fold)[[1L]]
  RD<-Rtrain[inner$train_idx,,drop=FALSE]
  RV<-Rtrain[inner$test_idx,,drop=FALSE]
  proposal<-cl2_learn(RD,RV,beta,b0,epsilon=.05,
                      max_roots=10L,min_v_gain=.25)
  ord<-cl2_order(beta)
  abs_order<-cret_evaluate(RE,beta,b0,ord,epsilon=.05)
  myopic<-cret_context_eval(RD,RE,beta,b0,epsilon=.05)
  depth2<-cl2_eval(RE,beta,b0,proposal$chosen,epsilon=.05)
  preds<-list(ORIGINAL_RLR=exact,CERT_STATIC_ABS=abs_order$prob,
       CERT_MYOPIC_CONTEXT=myopic$prob,
       CERT_TREE_DEPTH2_LOOKAHEAD=depth2$prob)
  for(method in all_methods){
    a<-if(method=="CERT_STATIC_ABS")abs_order else
       if(method=="CERT_MYOPIC_CONTEXT")myopic else
       if(method=="CERT_TREE_DEPTH2_LOOKAHEAD")depth2 else NULL
    pr<-preds[[method]]
    s<-arlet_scorer(ye01,pr)
    nr<-nr+1L
    records[[nr]]<-data.frame(
       dataset=DS,split_id=fold$split_id,method=method,
       n_train=nrow(T),n_E=nrow(E),
       n_nonzero_RLR_coefficients=sum(abs(beta)>1e-12),
       mean_evaluated=if(is.null(a))sum(abs(beta)>1e-12)
                    else mean(a$evaluated),
       mean_active=if(is.null(a))
            mean(rowSums(RE[,abs(beta)>1e-12,drop=FALSE]!=0))
                    else mean(a$active_n),
       stop_early_fraction=if(is.null(a))0 else
             mean(a$evaluated<sum(abs(beta)>1e-12)),
       BA=s$ba,macro_F1=s$macro_f1,
       logloss=s$logloss,brier=s$brier,
       policy_type=if(method=="CERT_TREE_DEPTH2_LOOKAHEAD")
                         proposal$chosen$type else "",
       inner_V_gain=if(method=="CERT_TREE_DEPTH2_LOOKAHEAD")
                         proposal$V_gain else NA_real_,
       training_provenance="ORIGINAL_MANUSCRIPT_PUBLIC_SOURCE_NOT_ARFF_PARITY")
    for(i in seq_len(nrow(E))){
      np<-np+1L
      patient[[np]]<-data.frame(
        dataset=DS,split_id=fold$split_id,
        method=method,sample_id=rownames(E)[[i]],
        truth=ye01[[i]],prob=pr[[i]],original_prob=exact[[i]],
        lower=if(is.null(a))exact[[i]] else a$lo[[i]],
        upper=if(is.null(a))exact[[i]] else a$hi[[i]],
        evaluated=if(is.null(a))sum(abs(beta)>1e-12)
                  else a$evaluated[[i]],
        active=if(is.null(a))sum(RE[i,abs(beta)>1e-12]!=0)
               else a$active_n[[i]],
        route=if(method=="CERT_TREE_DEPTH2_LOOKAHEAD")
                       depth2$branch[[i]] else "",
        next_test=if(method=="CERT_TREE_DEPTH2_LOOKAHEAD")
                       depth2$second[[i]] else NA_integer_,
        path=if(method=="CERT_TREE_DEPTH2_LOOKAHEAD")
                       depth2$paths[[i]] else "")
    }
    message(sprintf("CL2_CASE %s %s %s BA=%.4f logloss=%.4f eval=%.3f",
       DS,fold$split_id,method,s$ba,s$logloss,
       records[[nr]]$mean_evaluated))
  }
  structures[[fold$split_id]]<-list(
     n_active=sum(abs(beta)>1e-12),
     D_n=nrow(RD),V_n=nrow(RV),
     V_baseline_evaluated=proposal$V_baseline,
     V_candidate_evaluated=proposal$V_best,
     V_gain=proposal$V_gain,
     train_gain=proposal$train_gain,
     candidate_root_count=proposal$roots_considered,
     n_branch_differences=proposal$n_branches_modified,
     chosen_type=proposal$chosen$type,
     candidate_root=proposal$best_train$root,
     candidate_branch_second=proposal$best_train$second,
     chosen_root=proposal$chosen$root,
     chosen_branch_second=proposal$chosen$second,
     n_distinct_E_paths=length(unique(depth2$paths)),
     root_pair_name=if(proposal$chosen$type=="lookahead")
                      pairs$ordered_pair[[proposal$chosen$root]] else NULL,
     alternate_pair_names=if(proposal$chosen$type=="lookahead")
                   pairs$ordered_pair[proposal$chosen$second] else NULL)
  message(sprintf("CL2_FOLD_DONE dataset=%s split=%s policy=%s train_gain=%.3f V_gain=%.3f exact_class=PASS",
        DS,fold$split_id,proposal$chosen$type,proposal$train_gain,
        proposal$V_gain))
  unlink(tmp,recursive=TRUE,force=TRUE)
}
met<-do.call(rbind,records);pa<-do.call(rbind,patient)
stopifnot(nrow(met)==12L,nrow(pa)==4L*nrow(X),
          all(met$training_provenance==
               "ORIGINAL_MANUSCRIPT_PUBLIC_SOURCE_NOT_ARFF_PARITY"))
utils::write.csv(met,file.path(RESULT,"metrics.csv"),row.names=FALSE)
utils::write.csv(pa,file.path(RESULT,"patient_certificates.csv"),row.names=FALSE)
jsonlite::write_json(structures,file.path(RESULT,"tree_structures.json"),
                      auto_unbox=TRUE,pretty=TRUE)
for(nm in all_methods){
  x<-met[met$method==nm,,drop=FALSE]
  cat(sprintf("CL2_SUMMARY %s %s BA=%.5f macroF1=%.5f LL=%.5f EVAL=%.3f\n",
    DS,nm,mean(x$BA),mean(x$macro_F1),
    mean(x$logloss),mean(x$mean_evaluated)))
}
cat(sprintf("CL2_3_FOLDS_4_ARMS_DATASET=%s_PASS\n",DS))
cat("CL2_EXACT_ORIGINAL_RLR_CLASS_AND_POSTERIOR_CERTIFICATE=PASS\n")
cat("CL2_OUTER_E_UNTOUCHED_FOR_TREE_SELECTION=PASS\n")
cat("CL2_ORIGINAL_ARFF_SHA_PARITY=NOT_VERIFIED\n")
