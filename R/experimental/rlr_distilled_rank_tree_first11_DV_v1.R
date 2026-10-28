# RLR-TREE v1: exploratory FIRST original outer T split, inner D/V only.
# 11 original paper-named CompCancer cohorts. No original E ever read.
# Compare RLR OOF-teacher distilled autonomous tiny rank tree to:
# a same-candidate label-only tree, a TOP20 TSP-like label-only tree,
# original RLR20 teacher and true native switchBox kTSP5.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_distilled_rank_tree_core_v1.R")

ds<-Sys.getenv("RLR_DATASET_ID",unset="")
cohorts<-c("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
dt_need(ds %in% cohorts,"unknown original manuscript dataset")
f<-file.path("results/rlr_full11_same_pairs_v1",ds,
             "source_data",paste0(ds,".csv"))
out<-file.path("results/rlr_distilled_compact_tree11_DV_v1",ds)
dt_need(file.exists(f)&&file.exists(file.path(out,"SOURCE_SHA.txt"))&&
        !file.exists(file.path(out,"metrics_V.csv")),
        "missing public source/SHA or duplicated run")
tab<-read.csv(f,stringsAsFactors=FALSE,check.names=FALSE)
X<-as.data.frame(tab[,setdiff(names(tab),c("sample_id","class")),
                     drop=FALSE],check.names=FALSE)
rownames(X)<-tab$sample_id
y<-factor(tab$class)
reps<-if(ds %in% c("west-2001","yeoh-2002-v1"))2L else 10L
outer<-make_repeated_stratified_folds(y,5L,reps,42L)[[1L]]
T<-X[outer$train_idx,,drop=FALSE]
yt<-factor(as.character(y[outer$train_idx]),levels=levels(y))
dv<-rg_inner(yt)
D<-T[dv$D,,drop=FALSE]
V<-T[dv$V,,drop=FALSE]
yd<-factor(as.character(yt[dv$D]),levels=levels(y))
yv<-factor(as.character(yt[dv$V]),levels=levels(y))
truth<-as.integer(yv==levels(yv)[2L])
dt_need(!length(intersect(outer$train_idx,outer$test_idx)) &&
        !length(intersect(dv$D,dv$V)) &&
        nrow(D)+nrow(V)==nrow(T),"non-disjoint outer or D/V split")
started<-proc.time()[["elapsed"]]
# Original unpublished RLR exactly preserved: fit on D, independently
# calculate TRUE 3fold teacher predictions to avoid train-self imitation.
poolD<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
parent<-rg_parent(D,yd,V,poolD,K=20L)
oof<-dt_teacher_oof(D,yd,K=20L,folds=3L,seed=119L)
dt_need(length(oof)==nrow(D)&&all(oof>0&oof<1),
        "cross-fitted teacher not patient-aligned")
rlrpairs<-dt_parent_candidates(parent,max_candidates=20L)
tsp_pairs<-poolD[seq_len(min(20L,nrow(poolD))),
                 c("gene_left","gene_right","ordered_pair"),drop=FALSE]
minleaf<-max(4L,as.integer(floor(.07*nrow(D))))
alg<-list(
  LABEL_BINARY_SAME_RLR20=dt_train(D,yd,rlrpairs,oof,alpha=0,
           mode="BINARY",max_depth=2L,max_internal=3L,
           min_leaf=minleaf,penalty=.003),
  LABEL_BINARY_TSP_TOP20=dt_train(D,yd,tsp_pairs,oof,alpha=0,
           mode="BINARY",max_depth=2L,max_internal=3L,
           min_leaf=minleaf,penalty=.003),
  DISTILL_BINARY_RLR20=dt_train(D,yd,rlrpairs,oof,alpha=.55,
           mode="BINARY",max_depth=2L,max_internal=3L,
           min_leaf=minleaf,penalty=.003),
  DISTILL_TERNARY_RLR20=dt_train(D,yd,rlrpairs,oof,alpha=.55,
           mode="TERNARY",max_depth=2L,max_internal=3L,
           min_leaf=minleaf,penalty=.003))
# EXACT SAME V population, native kTSP labels and original RLR20 probs.
fit_ktsp_train_test(D,yd,V,yv,k=5L,
    out_dir=file.path(out,"baselines","ktsp5"),seed=42L)
kp<-read.csv(file.path(out,"baselines","ktsp5","predictions_test.csv"),
             stringsAsFactors=FALSE,check.names=FALSE)
k_pairs<-read.csv(file.path(out,"baselines","ktsp5","selected_pairs.csv"),
             check.names=FALSE)
dt_need(identical(as.character(kp$sample_id),rownames(V)) &&
        nrow(k_pairs)>=1L&&nrow(k_pairs)<=5L,
        "real kTSP5 V predictions/realized K invalid")
results<-lapply(alg,function(model)dt_predict(model,V))
results$RLR20_TEACHER<-list(prob=parent$test,
  label=as.integer(parent$test>=.5),
  depth=rep(sum(abs(as.numeric(stats::coef(parent$fit,
      s="lambda.min"))[-1L])>1e-9),nrow(V)),
  leaf=rep(NA_integer_,nrow(V)),
  path=rep("GLOBAL_RLR20",nrow(V)))
results$KTSP5_NATIVE<-list(prob=as.numeric(kp$prob),
  label=as.integer(as.character(kp$y_pred)==levels(yv)[2L]),
  depth=rep(nrow(k_pairs),nrow(V)),
  leaf=rep(NA_integer_,nrow(V)),
  path=rep("NATIVE_SWITCHBOX_KTSP",nrow(V)))
f1<-function(h,cl){
 tp<-sum(truth==cl & h==cl)
 fp<-sum(truth!=cl & h==cl)
 fn<-sum(truth==cl & h!=cl)
 if(2*tp+fp+fn==0)0 else 2*tp/(2*tp+fp+fn)
}
auc<-function(p){
 n1<-sum(truth==1L);n0<-sum(truth==0L)
 rr<-rank(p,ties.method="average")
 (sum(rr[truth==1L])-n1*(n1+1)/2)/(n0*n1)
}
metrics<-list();preds<-list();struc<-list();i<-0L;j<-0L;k<-0L
for(method in names(results)){
 pr<-results[[method]]
 student<-alg[[method]]
 q<-as.numeric(pr$prob);h<-as.integer(pr$label)
 dt_need(length(q)==nrow(V) && all(is.finite(q)) &&
        all(q>=0 & q<=1),"invalid E-independent V model outputs")
 st<-dt_fidelity(parent$test,q)
 i<-i+1L
 metrics[[i]]<-data.frame(
  dataset=ds,method=method,n_D=nrow(D),n_V=nrow(V),
  positive_f1=f1(h,1L),
  macro_f1=mean(c(f1(h,0L),f1(h,1L))),
  ba=.5*(mean(h[truth==0L]==0L)+mean(h[truth==1L]==1L)),
  auc=auc(q),logloss=rg_ll(truth,q),
  brier=mean((truth-q)^2),
  mean_evaluated_pairs=mean(pr$depth),
  max_evaluated_pairs=max(pr$depth),
  n_internal_nodes=if(is.null(student))NA_integer_ else student$internal,
  n_leaves=if(is.null(student))NA_integer_ else student$n_leaves,
  n_unique_model_genes=if(is.null(student))NA_integer_ else student$n_unique_genes,
  n_tree_pairs=if(is.null(student))NA_integer_ else length(student$used_pairs),
  student_teacher_agreement=st$agreement,
  teacher_soft_cross_entropy=st$soft_logloss,
  teacher_fidelity_brier=st$teacher_brier,
  evidence="ORIGINAL_OUTER_T_D_ONLY_TRAIN_V_ONLY_COMPARISON_NO_E")
 for(r in seq_len(nrow(V))){
  j<-j+1L
  preds[[j]]<-data.frame(dataset=ds,method=method,
   sample_id=rownames(V)[r],truth=truth[r],prob=q[r],
   label=h[r],teacher_prob=parent$test[r],
   evaluated_pairs=pr$depth[r],
   terminal_leaf=pr$leaf[r],
   actual_path=pr$path[r],stringsAsFactors=FALSE)
 }
 if(!is.null(student)){
  for(nd in student$nodes){
   k<-k+1L
   pair<-if(is.na(nd$split))"" else
     as.character(student$pairs$ordered_pair[nd$split])
   struc[[k]]<-data.frame(dataset=ds,method=method,
    node_id=nd$id,depth=nd$depth,training_n=length(nd$idx),
    split_pair=pair,children=paste(nd$children,collapse="|"),
    leaf_prob=nd$p,
    mode=student$mode,
    source="D_ONLY_NO_OUTER_E",stringsAsFactors=FALSE)
  }
 }
 cat(sprintf("RLRTREE_DV_RESULT %s %s F1=%.5f MacroF1=%.5f AUC=%.5f LL=%.5f path=%.2f nodes=%s leaves=%s fidelity=%.4f\n",
   ds,method,metrics[[i]]$positive_f1,metrics[[i]]$macro_f1,
   metrics[[i]]$auc,metrics[[i]]$logloss,
   metrics[[i]]$mean_evaluated_pairs,
   if(is.null(student))"NA" else as.character(student$internal),
   if(is.null(student))"NA" else as.character(student$n_leaves),
   st$agreement))
}
write.csv(do.call(rbind,metrics),file.path(out,"metrics_V.csv"),row.names=FALSE)
write.csv(do.call(rbind,preds),file.path(out,"patient_V.csv"),row.names=FALSE)
write.csv(do.call(rbind,struc),file.path(out,"tree_nodes_D.csv"),row.names=FALSE)
write.csv(data.frame(dataset=ds,n_total=nrow(X),n_outer_T=nrow(T),
 n_D=nrow(D),n_V=nrow(V),n_outer_E_unread=length(outer$test_idx),
 D_ids=paste(rownames(D),collapse="|"),
 V_ids=paste(rownames(V),collapse="|"),
 teacher_K=20L,teacher_OOF_folds=3L,teacher_alpha=.55,
 max_depth=2L,max_internal=3L,min_leaf=minleaf,
 soft_ce_node_penalty=.003,
 source="PUBLIC_SCHLIEP_RECONSTRUCTION_ARFF_PARITY_UNVERIFIED",
 elapsed_seconds=proc.time()[["elapsed"]]-started),
 file.path(out,"manifest.csv"),row.names=FALSE)
cat("RLRTREE_TEACHER_ACTUALLY_OUT_OF_FOLD_ON_D=PASS\n")
cat("RLRTREE_INFERENCE_AUTONOMOUS_MAX_DEPTH_TWO=PASS\n")
cat("RLRTREE_COMPLETE_SIX_V_METHODS_ORIGINAL_E_UNREAD=PASS\n")
