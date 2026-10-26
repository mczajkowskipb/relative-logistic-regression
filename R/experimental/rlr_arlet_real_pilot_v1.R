# RLR ARLET v1: finite real publication-dataset pilot.
# ONLY chen-2002, singh-2002, yeoh-2002-v1 from locked original paper.
# Original RLR TSP kTSP baseline functions used as-is from main.
# Greedy TSPDT comparison is transparently a proxy, not CUDA implementation.
# Outer E untouched by pair ranking, EA evaluation and hyperparameter choice.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_arlet_evidence_tree_v1.R")
ID<-Sys.getenv("RLR_DATASET_ID",unset="")
SETS<-c("chen-2002","singh-2002","yeoh-2002-v1")
arlet_need(ID %in% SETS,"pilot dataset outside paper inventory")
OUT<-file.path("results/rlr_arlet_v1",ID)
input<-file.path(OUT,"source_data",paste0(ID,".csv"))
arlet_need(file.exists(input)&&file.exists(file.path(OUT,"input_source_manifest.json")),
           "source input provenance unavailable")
dat<-utils::read.csv(input,check.names=FALSE,stringsAsFactors=FALSE)
X<-as.data.frame(dat[,setdiff(names(dat),c("class","sample_id")),
  drop=FALSE],check.names=FALSE)
rownames(X)<-dat$sample_id
y<-factor(dat$class)
arlet_need(length(levels(y))==2L,"dataset not binary")
arlet_need(nrow(X)>=90&&ncol(X)>=50,"source dimensions too narrow")
folds<-make_repeated_stratified_folds(y,n_folds=3L,n_repeats=1L,seed=42L)
arlet_need(length(folds)==3L,"pilot exactly three folds")
dir.create(file.path(OUT,"pilot"),recursive=TRUE,showWarnings=FALSE)
rows<-list();predicts<-list();splits<-list();models<-list()
mi<-0L;pi<-0L
for(fold in folds){
  tr<-fold$train_idx;te<-fold$test_idx
  arlet_need(!length(intersect(tr,te))&&
     identical(sort(c(tr,te)),seq_len(nrow(X))),"outer E overlap")
  T<-X[tr,,drop=FALSE];E<-X[te,,drop=FALSE]
  yt<-y[tr];ye<-y[te]
  ytrain<-as.numeric(yt==levels(y)[2L])
  ye01<-as.numeric(ye==levels(y)[2L])
  splits[[length(splits)+1L]]<-data.frame(
    dataset=ID,split_id=fold$split_id,
    n_train=nrow(T),n_E=nrow(E),
    E_patient_ids=paste(rownames(E),collapse="|"))
  # All original controls use exactly original supervised pair nomination
  # only from OUTER TRAIN data, and original manuscript code for learners.
  original_pool<-compute_top_pairs_switchbox(T,yt,top_n=1000L,
                             buffer_mult=5L)
  wdir<-file.path(OUT,"pilot","original_wrapper",fold$split_id)
  fit_rlr_train_test(T,yt,E,ye,pairs_pool=original_pool,K=50L,
      m=10L,d=2L,out_dir=file.path(wdir,"rlr"),
      alpha=.5,lambda_choice="lambda.min",seed=42L)
  fit_tsp1_manual_train_test(T,yt,E,ye,pairs_pool=original_pool,
      out_dir=file.path(wdir,"tsp"))
  fit_ktsp_train_test(T,yt,E,ye,k=9L,
      out_dir=file.path(wdir,"ktsp"),seed=42L)
  original_pred<-list()
  for(nm in c("rlr","tsp","ktsp")){
     tab<-utils::read.csv(file.path(wdir,nm,"predictions_test.csv"))
     arlet_need(identical(as.character(tab$sample_id),rownames(E)),
                 paste("original benchmark E order",nm))
     arlet_need(identical(as.character(tab$y_true),as.character(ye)),
                 paste("original benchmark class labels",nm))
     original_pred[[nm]]<-as.numeric(tab$prob)
  }
  # Inner D/V split lives fully inside OUTER TRAIN and is independent of E.
  inn<-make_repeated_stratified_folds(yt,n_folds=4L,n_repeats=1L,
                                     seed=42L+fold$fold+100L)[[1L]]
  D<-T[inn$train_idx,,drop=FALSE]
  V<-T[inn$test_idx,,drop=FALSE]
  dy<-as.numeric(yt[inn$train_idx]==levels(y)[2L])
  vy<-as.numeric(yt[inn$test_idx]==levels(y)[2L])
  arlet_need(length(intersect(inn$train_idx,inn$test_idx))==0L,
           "D/V overlap")
  # Critically: proposal pair pool is nominated ONLY from D, not V.
  candidate<-arlet_pair_pool(D,yt[inn$train_idx],
          seed=42L+fold$fold)
  zv<-arlet_encode_both(D,V,candidate)
  search<-arlet_search(zv$train,dy,zv$test,vy,candidate,
                      seed=42L+fold$fold)
  # Refit *selected structures* once on ALL outer training data.
  # This is standard nested selection/refitting; E remains unseen.
  z<-arlet_encode_both(T,E,candidate)
  detfit<-arlet_fit(z$train,ytrain,candidate,
    search$deterministic$proposal,search$deterministic$s,
    search$deterministic$conf)
  eafit<-arlet_fit(z$train,ytrain,candidate,
    search$evolutionary$proposal,search$evolutionary$s,
    search$evolutionary$conf)
  dp<-arlet_predict(detfit,z$test)
  ep<-arlet_predict(eafit,z$test)
  tspp<-arlet_tspdt_proxy(z$train,ytrain,z$test,
                        depth=2L,min_leaf=7L)
  predictions<-list(ORIGINAL_RLR=original_pred$rlr,
    TSP_1=original_pred$tsp,kTSP_9=original_pred$ktsp,
    TSPDT_GREEDY_PROXY=tspp,ARLET_DETERMINISTIC=dp$prob,
    ARLET_EA=ep$prob)
  for(method in names(predictions)){
     p<-predictions[[method]]
     arlet_need(length(p)==nrow(E)&&all(is.finite(p))&&
                   all(p>=0&p<=1),paste("invalid E score",method))
     score<-arlet_scorer(ye01,p)
     mi<-mi+1L
     thisfit<-if(method=="ARLET_DETERMINISTIC")detfit else
              if(method=="ARLET_EA")eafit else NULL
     rels<-if(!is.null(thisfit))length(unique(c(
            thisfit$gate$top,thisfit$gate$second))) else NA_integer_
     used<-if(method=="ARLET_DETERMINISTIC")dp else
           if(method=="ARLET_EA")ep else NULL
     rows[[mi]]<-data.frame(
       dataset=ID,split_id=fold$split_id,method=method,
       n_train=nrow(T),n_E=nrow(E),
       logloss=score$logloss,brier=score$brier,BA=score$ba,
       macro_F1=score$macro_f1,
       total_tree_relations=rels,
       root_relations=if(!is.null(thisfit))thisfit$gate$s else NA_integer_,
       n_deferred=if(!is.null(used))sum(used$deferred) else NA_integer_,
       mean_active_relations=if(!is.null(used))
                            mean(used$active_evidence) else NA_real_,
       mean_evaluated_relations=if(!is.null(used))
                            mean(used$evaluated) else NA_real_,
       EA_val_gain=if(method=="ARLET_EA")search$ea_val_gain else NA_real_,
       evidence_stage="EXPOSED_PUBLICATION_SOURCE_RECONSTRUCTION_NOT_ARFF_SHA",
       stringsAsFactors=FALSE)
     for(j in seq_len(nrow(E))){
       pi<-pi+1L
       predicts[[pi]]<-data.frame(
         dataset=ID,split_id=fold$split_id,
         sample_id=rownames(E)[[j]],
         method=method,truth=ye01[[j]],prob=p[[j]],
         route=if(!is.null(used))used$route[[j]] else "",
         active_relations=if(!is.null(used))
                       used$active_evidence[[j]] else NA_integer_)
     }
     message(sprintf("ARLET_CASE %s %s %-21s BA=%.4f F1=%.4f LL=%.4f",
               ID,fold$split_id,method,
               score$ba,score$macro_f1,score$logloss))
  }
  models[[fold$split_id]]<-list(
     inner_D_ids=rownames(D),inner_V_ids=rownames(V),
     outer_E_ids=rownames(E),pool_n=nrow(candidate),
     pool_first_n=1000L,
     unique_genes=length(unique(c(candidate$gene_left,
                                   candidate$gene_right))),
     deterministic=list(s=detfit$gate$s,
               confidence=detfit$gate$confidence,
               top_pairs=candidate$ordered_pair[detfit$gate$top],
               deferred_pairs=candidate$ordered_pair[detfit$gate$second]),
     EA=list(s=eafit$gate$s,confidence=eafit$gate$confidence,
          val_gain=search$ea_val_gain,
          top_pairs=candidate$ordered_pair[eafit$gate$top],
          deferred_pairs=candidate$ordered_pair[eafit$gate$second]))
  unlink(wdir,recursive=TRUE,force=TRUE)
  message(sprintf("ARLET_OUTER_FINISHED dataset=%s fold=%s pool=%d",
           ID,fold$split_id,nrow(candidate)))
}
met<-do.call(rbind,rows)
pred<-do.call(rbind,predicts)
stopifnot(nrow(met)==18L, nrow(pred)==6L*nrow(X),
          all(met$evidence_stage=="EXPOSED_PUBLICATION_SOURCE_RECONSTRUCTION_NOT_ARFF_SHA"))
utils::write.csv(met,file.path(OUT,"pilot","method_metrics.csv"),
                 row.names=FALSE)
utils::write.csv(pred,file.path(OUT,"pilot","patient_predictions.csv"),
                 row.names=FALSE)
utils::write.csv(do.call(rbind,splits),
   file.path(OUT,"pilot","split_manifest.csv"),row.names=FALSE)
jsonlite::write_json(models,file.path(OUT,"pilot","tree_models.json"),
             pretty=TRUE,auto_unbox=TRUE)
for(method in unique(met$method)){
 x<-met[met$method==method,]
 cat(sprintf("ARLET_SUMMARY %s %s mean_BA=%.5f mean_macroF1=%.5f mean_LL=%.5f\n",
      ID,method,mean(x$BA),mean(x$macro_F1),mean(x$logloss)))
}
cat(sprintf("ARLET_3_FOLD_6_METHODS_DATASET=%s_PASS\n",ID))
cat("ARLET_NO_E_LABELS_USED_FOR_MODEL_SELECTION=PASS\n")
cat("ARLET_ORIGINAL_RLR_TSP_KTSP_USED_UNMODIFIED=PASS\n")
cat("ARLET_TSPDT_CONTROL=GREEDY_PROXY_NOT_CUDA_REPRODUCTION\n")
cat("ARLET_PROCESSED_ORIGINAL_ARFF_SHA=NOT_VERIFIED\n")
