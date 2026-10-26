# STRICT ORIGINAL PAPER INVENTORY: same selected ternary m10 relations, separate
# coefficients vs one equally weighted scalar logistic probability. The
# ORIGINAL RLR arm calls fit_rlr_train_test() from the original frozen code.
# On reconstructed Schliep source text; NOT exact ARFF SHA256 reproducibility.
source("R/new/07_world_protocol_utils.R")

dataset<-Sys.getenv("RLR_DATASET_ID",unset="")
max_splits<-as.integer(Sys.getenv("RLR_MAX_SPLITS",unset="1"))
nreps<-if(dataset %in% c("west-2001","yeoh-2002-v1"))2L else 10L
original_11<-c("alizadeh-2000-v1","armstrong-2002-v1",
  "bittner-2000","chen-2002","golub-1999-v1","gordon-2002",
  "pomeroy-2002-v1","shipp-2002-v1","singh-2002",
  "west-2001","yeoh-2002-v1")
stopifnot(dataset %in% original_11, !is.na(max_splits),max_splits>=0L)
dirr<-file.path("results/rlr_full11_same_pairs_v1",dataset)
csvfile<-file.path(dirr,"source_data",paste0(dataset,".csv"))
stopifnot(file.exists(csvfile),file.exists(file.path(dirr,"input_source_manifest.json")))
dat<-utils::read.csv(csvfile,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(all(c("sample_id","class") %in% names(dat)))
X<-as.data.frame(dat[,setdiff(names(dat),c("sample_id","class")),drop=FALSE],
                 check.names=FALSE)
rownames(X)<-dat$sample_id
y<-factor(dat$class)
stopifnot(length(levels(y))==2L,length(unique(dat$sample_id))==nrow(dat),
          all(vapply(X,function(t)is.numeric(t)&&all(is.finite(t)),logical(1))),
          length(unique(names(X)))==ncol(X))
full_splits<-5L*nreps
folds<-make_repeated_stratified_folds(y,n_folds=5L,n_repeats=nreps,seed=42L)
stopifnot(length(folds)==full_splits)
if(max_splits>0L)folds<-head(folds,max_splits)
is_full<-length(folds)==full_splits
outdir<-file.path(dirr,if(is_full)"full" else "smoke")
stopifnot(!dir.exists(outdir))
dir.create(outdir,recursive=TRUE)
metrics<-list(); preds<-list(); evidence<-list();ri<-0L;pi<-0L;ei<-0L
clip<-function(p)pmin(1-1e-15,pmax(1e-15,p))
for(fold in folds){
  tr<-fold$train_idx; te<-fold$test_idx
  stopifnot(!length(intersect(tr,te)),identical(sort(c(tr,te)),
              seq_len(nrow(X))))
  T<-X[tr,,drop=FALSE];E<-X[te,,drop=FALSE]
  yt<-y[tr];ye<-y[te]
  positive<-as.numeric(ye==levels(y)[2L])
  # Original shortlist candidate ranking, once per outer training fold.
  pool<-compute_top_pairs_switchbox(
    X_train=T,y_train=yt,top_n=1000L,buffer_mult=5L,
    out_dir=NULL,save_all=FALSE)
  stopifnot(nrow(pool)>=50L)
  for(K in c(10L,20L,50L)){
    # Call the *literal original publication RLR*, not a copied formula.
    tmp<-file.path(outdir,"original_rlr_work",fold$split_id,paste0("K",K))
    fit_rlr_train_test(X_train=T,y_train=yt,X_test=E,y_test=ye,
      pairs_pool=pool,K=K,m=10L,d=2L,
      out_dir=tmp,alpha=.5,lambda_choice="lambda.min",seed=42L)
    old<-utils::read.csv(file.path(tmp,"predictions_test.csv"),
                         stringsAsFactors=FALSE)
    stopifnot(identical(as.character(old$sample_id),rownames(E)),
              all(as.character(old$y_true)==as.character(ye)))
    p_rlr<-as.numeric(old$prob)
    pairs<-apply_disjointness(pool,K=K,d=2L)
    RTR<-encode_pair_matrix(T,pairs,m=10L)
    ETR<-encode_pair_matrix(E,pairs,m=10L)
    vote<-rowMeans(RTR); vote_E<-rowMeans(ETR)
    y01<-as.numeric(yt==levels(y)[2L])
    if(stats::sd(vote)<1e-12){
      p_equal<-rep(mean(y01),nrow(E))
    }else{
      # Equal relation voting + SAME elastic net alpha and inner CV as RLR.
      # Constant extra column is excluded; only 1 scalar slope may be fitted.
      vtr<-cbind(equal_vote=vote,excluded_constant=rep(0,length(vote)))
      vte<-cbind(equal_vote=vote_E,excluded_constant=rep(0,length(vote_E)))
      set.seed(42L)
      fit<-glmnet::cv.glmnet(x=vtr,y=y01,family="binomial",alpha=.5,
            nfolds=safe_nfolds(yt,5L),type.measure="deviance")
      co<-as.numeric(stats::coef(fit,s="lambda.min"))
      stopifnot(length(co)==3L,abs(co[[3L]])<1e-12)
      p_equal<-as.numeric(stats::predict(fit,newx=vte,s="lambda.min",
                                        type="response"))
    }
    stopifnot(length(p_rlr)==nrow(E),length(p_equal)==nrow(E),
              all(is.finite(p_rlr)),all(is.finite(p_equal)),
              all(p_rlr>=0&p_rlr<=1),all(p_equal>=0&p_equal<=1))
    loss<-function(z){
      p<-clip(z)
      -(positive*log(p)+(1-positive)*log1p(-p))
    }
    la<-loss(p_rlr)
    lb<-loss(p_equal)
    for(method in c("ORIGINAL_RLR","REGULARIZED_EQUAL_VOTE")){
      p<-if(method=="ORIGINAL_RLR")p_rlr else p_equal
      ll<-if(method=="ORIGINAL_RLR")la else lb
      ri<-ri+1L
      metrics[[ri]]<-data.frame(
        dataset=dataset,shard=if(nreps==10L)"9x10" else "2x2",
        split_id=fold$split_id,rep=fold$rep,fold=fold$fold,
        K=K,method=method,n_train=length(tr),n_E=length(te),
        n_selected_pairs=nrow(pairs),
        logloss=mean(ll),brier=mean((p-positive)^2),
        delta_logloss=mean(ll-la),
        scientific_stage="SOURCE_RECONSTRUCTION_NOT_ARFF_PARITY",
        stringsAsFactors=FALSE)
      for(i in seq_along(te)){
        pi<-pi+1L
        preds[[pi]]<-data.frame(
          dataset=dataset,split_id=fold$split_id,K=K,method=method,
          sample_id=rownames(E)[[i]],truth=positive[[i]],
          original_prob=p_rlr[[i]],method_prob=p[[i]])
      }
    }
    # Original per-model output is generated and checked, not used as final
    # evidence package: avoids accidentally shipping unrelated artifacts.
    unlink(tmp,recursive=TRUE,force=TRUE)
  }
  ei<-ei+1L
  evidence[[ei]]<-data.frame(dataset=dataset,split_id=fold$split_id,
     rep=fold$rep,fold=fold$fold,
     n_train=length(tr),n_E=length(te),
     test_ids=paste(rownames(E),collapse="|"))
  cat(sprintf("FULL11_FOLD_DONE dataset=%s split=%s model_fits=6\n",
              dataset,fold$split_id),file=stderr())
}
m<-do.call(rbind,metrics)
pr<-do.call(rbind,preds)
e<-do.call(rbind,evidence)
stopifnot(nrow(m)==6L*length(folds),nrow(e)==length(folds),
          nrow(pr)==6L*sum(e$n_E),
          all(m$method %in% c("ORIGINAL_RLR","REGULARIZED_EQUAL_VOTE")))
utils::write.csv(m,file.path(outdir,"method_metrics.csv"),row.names=FALSE)
utils::write.csv(pr,file.path(outdir,"patient_predictions.csv"),row.names=FALSE)
utils::write.csv(e,file.path(outdir,"split_manifest.csv"),row.names=FALSE)
for(K in c(10L,20L,50L)){
 mm<-m[m$K==K & m$method=="REGULARIZED_EQUAL_VOTE",]
 cat(sprintf("FULL11_SAME_PAIRS dataset=%s K=%d EQUAL_minus_RLR_logloss=%+.7f RLR_better=%d/%d\n",
       dataset,K,mean(mm$delta_logloss),
       sum(mm$delta_logloss>0),nrow(mm)))
}
cat(sprintf("FULL11_DATASET_FINISHED=%s\n",dataset))
cat(sprintf("FULL11_SPLITS_ACTUAL=%d EXPECTED=%d\n",length(folds),full_splits))
cat(sprintf("FULL11_MODE=%s\n",if(is_full)"FULL" else "SMOKE"))
cat("FULL11_RLR_BASELINE=ORIGINAL_FIT_RLR_TRAIN_TEST\n")
cat("FULL11_SOURCE_RECONSTRUCTION_ARFF_SHA256=NOT_VERIFIED\n")
