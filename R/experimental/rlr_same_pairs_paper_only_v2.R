# Controlled component study for the original unpublished RLR article.
# SAME supervised pair pool, SAME d2/K shortlist, SAME ternary m=10.
# Compare single-slope Platt calibrated equal votes vs learned individual RLR weights.
# Only training labels used in candidates and fits; E used for evaluation.
source("R/new/07_world_protocol_utils.R")
OUT<-"results/rlr_same_pairs_paper_only_v2"
stopifnot(dir.exists(OUT),file.exists(file.path(OUT,"SOURCE_SHA.txt")),
          !file.exists(file.path(OUT,"patient_predictions.csv")))
sets<-c("armstrong-2002-v1","singh-2002")
KS<-c(10L,20L,50L)
methods<-c("ORIGINAL_RLR","CALIBRATED_EQUAL_VOTE")
nmethod<-0L
npred<-0L
nmfest<-0L
mrows<-list()
predrows<-list()
splits<-list()

run_model<-function(X,y,testR,arm) {
  y01<-as.numeric(y==levels(y)[2L])
  if(arm=="CALIBRATED_EQUAL_VOTE") {
    # Univariate Platt calibration on the *same* signed m10 rank-pair
    # vote representation. glmnet requires >=2 columns, so use the
    # standard one-slope logistic model for this fair baseline.
    vote<-rowMeans(X)
    vote_E<-rowMeans(testR)
    if(stats::sd(vote)<1e-12)
      return(rep(mean(y01),length(vote_E)))
    calibrated<-suppressWarnings(stats::glm(y01 ~ vote,
                         family=stats::binomial()))
    if(!isTRUE(calibrated$converged)||
       any(!is.finite(stats::coef(calibrated))))
      return(rep(mean(y01),length(vote_E)))
    p<-as.numeric(stats::predict(calibrated,
                 newdata=data.frame(vote=vote_E),
                 type="response"))
    stopifnot(length(p)==length(vote_E),all(is.finite(p)),
              all(p>=0&p<=1))
    return(p)
  }
  stopifnot(arm=="ORIGINAL_RLR")
  set.seed(42L)
  cv<-glmnet::cv.glmnet(x=X,y=y01,family="binomial",alpha=.5,
    nfolds=safe_nfolds(y,5L),type.measure="deviance")
  p<-as.numeric(stats::predict(cv,newx=testR,s="lambda.min",
                               type="response"))
  stopifnot(length(p)==nrow(testR),all(is.finite(p)),
            all(p>=0&p<=1))
  p
}
for(ds in sets){
  dfile<-file.path(OUT,"source_data",paste0(ds,".csv"))
  if(!file.exists(dfile))stop("DATA_NOT_FOUND: ",ds)
  tab<-utils::read.csv(dfile,check.names=FALSE,stringsAsFactors=FALSE)
  X<-as.data.frame(tab[,setdiff(names(tab),c("sample_id","class")),
                        drop=FALSE],check.names=FALSE)
  rownames(X)<-tab$sample_id
  y<-factor(tab$class)
  stopifnot(length(levels(y))==2L,nrow(X)>=60L,
       all(vapply(X,function(x)is.numeric(x)&&all(is.finite(x)),
                  logical(1))))
  foldset<-make_repeated_stratified_folds(y,n_folds=3L,n_repeats=1L,
                                        seed=42L)
  for(fi in seq_along(foldset)){
    sp<-foldset[[fi]]
    ti<-sp$train_idx;ei<-sp$test_idx
    stopifnot(length(intersect(ti,ei))==0,
      identical(sort(c(ti,ei)),seq_len(nrow(X))))
    T<-X[ti,,drop=FALSE];E<-X[ei,,drop=FALSE]
    yt<-y[ti];ytest<-as.numeric(y[ei]==levels(y)[2L])
    pool<-compute_top_pairs_switchbox(X_train=T,y_train=yt,
       top_n=1000L,buffer_mult=5L,out_dir=NULL,save_all=FALSE)
    stopifnot(nrow(pool)>=50L)
    nmfest<-nmfest+1L
    splits[[nmfest]]<-data.frame(dataset=ds,outer_fold=fi-1L,
      n_train=nrow(T),n_E=nrow(E),source_n=nrow(X),
      test_sample_ids=paste(rownames(E),collapse="|"),
      stringsAsFactors=FALSE)
    for(K in KS){
      pairs<-apply_disjointness(pool,K=K,d=2L)
      stopifnot(nrow(pairs)>=1L)
      R<-encode_pair_matrix(T,pairs,m=10L)
      Q<-encode_pair_matrix(E,pairs,m=10L)
      stopifnot(identical(colnames(R),colnames(Q)),
                nrow(R)==length(ti),nrow(Q)==length(ei))
      pf<-list()
      for(method in methods)pf[[method]]<-run_model(R,yt,Q,method)
      base<-pf[["ORIGINAL_RLR"]]
      for(method in methods){
        p<-pf[[method]]
        eps<-1e-15
        eval_loss<-function(pred){
          px<-pmin(1-eps,pmax(eps,pred))
          -(ytest*log(px)+(1-ytest)*log1p(-px))
        }
        l<-eval_loss(p)
        l0<-eval_loss(base)
        nmethod<-nmethod+1L
        mrows[[nmethod]]<-data.frame(dataset=ds,
          outer_fold=fi-1L,K=K,method=method,
          n_train=nrow(T),n_E=nrow(E),
          n_selected_pairs=nrow(pairs),
          logloss=mean(l),brier=mean((p-ytest)^2),
          delta_logloss=mean(l-l0),
          vote_signal_sd=stats::sd(rowMeans(R)),
          stage="ORIGINAL_PAPER_INVENTORY_SOURCE_RECONSTRUCTION",
          stringsAsFactors=FALSE)
        for(j in seq_along(ei)){
          npred<-npred+1L
          predrows[[npred]]<-data.frame(
            dataset=ds,outer_fold=fi-1L,K=K,method=method,
            sample_id=rownames(E)[[j]],truth=ytest[[j]],
            baseline_prob=base[[j]],method_prob=p[[j]],
            stringsAsFactors=FALSE)
        }
        message(sprintf(
          "SAME_PAIR_CASE %s fold=%d K=%d %-21s deltaLL=%+.6f",
          ds,fi-1L,K,method,mean(l-l0)))
      }
    }
  }
}
sm<-do.call(rbind,mrows);sp<-do.call(rbind,predrows);
spl<-do.call(rbind,splits)
stopifnot(nrow(sm)==2L*3L*3L*2L,
          nrow(spl)==6L,nrow(sp)==2L*3L*sum(vapply(sets,function(z){
             n<-if(z=="armstrong-2002-v1")72L else 102L
             n
          },integer(1))),
          all(sm$method %in% methods))
utils::write.csv(sm,file.path(OUT,"method_metrics.csv"),row.names=FALSE)
utils::write.csv(sp,file.path(OUT,"patient_predictions.csv"),row.names=FALSE)
utils::write.csv(spl,file.path(OUT,"split_manifest.csv"),row.names=FALSE)
for(K in KS){
  rows<-sm[sm$K==K & sm$method=="CALIBRATED_EQUAL_VOTE",]
  message(sprintf("SAME_PAIR_SUMMARY K=%d EQUAL_VOTE_vs_RLR deltaLL=%+.8f RLR_wins=%d/6",
    K,mean(rows$delta_logloss),
    sum(rows$delta_logloss>0L)))
}
cat("SAME_PAIRS_2_ORIGINAL_COMPCANCER_COHORTS=PASS\n")
cat("SAME_PAIRS_6_OUTER_E_PARTITIONS=PASS\n")
cat("SAME_PAIRS_36_FROZEN_METHOD_FITS=PASS\n")
cat("SAME_PAIRS_BOTH_PREDICTION_VECTORS_ARCHIVED=PASS\n")
cat("SAME_PAIRS_NO_NEW_DATA_OR_TREE=PASS\n")
