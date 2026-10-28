# RLR-WT v4: 11 original named cohorts, FIRST outer T -> strict D/V.
# ZERO outer E use. RLR-beta-weighted two-pair node tree vs sign ablation,
# label-only vs 25% independent OOF teacher, real kTSP5 and frozen RLR20.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_weighted_pairpair_tree_core_v4.R")

ds<-Sys.getenv("RLR_DATASET_ID",unset="")
cohorts<-c("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
wt_need(ds %in% cohorts,"unknown original publication named cohort")
src_path<-file.path("results/rlr_full11_same_pairs_v1",
                    ds,"source_data",paste0(ds,".csv"))
out<-file.path("results/rlr_weighted_pairpair11_DV_v4",ds)
wt_need(file.exists(src_path) && file.exists(file.path(out,"SOURCE_SHA.txt")) &&
        !file.exists(file.path(out,"metrics_V.csv")),
        "missing public source or existing immutable experiment output")
raw<-read.csv(src_path,stringsAsFactors=FALSE,check.names=FALSE)
X<-as.data.frame(raw[,setdiff(names(raw),c("sample_id","class")),
                       drop=FALSE],check.names=FALSE)
rownames(X)<-raw$sample_id
y<-factor(raw$class)
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
wt_need(nrow(D)>0&&nrow(V)>0&&
        !length(intersect(outer$train_idx,outer$test_idx))&&
        !length(intersect(dv$D,dv$V))&&
        nrow(D)+nrow(V)==nrow(T),"outer or inner leakage")
timer<-proc.time()[["elapsed"]]
pool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
parent<-rg_parent(D,yd,V,pool,K=20L)
oof<-dt_teacher_oof(D,yd,K=20L,folds=3L,seed=119L)
wt_need(length(oof)==nrow(D) && all(oof>0&oof<1),"OOF RLR teacher invalid")
n_fallback<-as.integer(attr(oof,"ridge_fallback_folds"))
wt_need(length(n_fallback)==1L&&is.finite(n_fallback)&&n_fallback>=0L,
        "OOF original/rare-class surrogate folds not disclosed")
minleaf<-max(4L,as.integer(floor(.07*nrow(D))))
# ONE model architecture and ONE dataset. Only two experimental factors
# vary, fixed before observing V: coef beta vs sign and OOF vs labels.
algs<-list(
 WT_BETA_LABEL=wt_fit(D,yd,parent,oof,coeff_mode="BETA",alpha=0,
   leaf_alpha=0,max_depth=2L,max_internal=3L,min_leaf=minleaf,
   penalty=.003),
 WT_BETA_OOF=wt_fit(D,yd,parent,oof,coeff_mode="BETA",alpha=.25,
   leaf_alpha=0,max_depth=2L,max_internal=3L,min_leaf=minleaf,
   penalty=.003),
 WT_SIGN_LABEL=wt_fit(D,yd,parent,oof,coeff_mode="SIGN",alpha=0,
   leaf_alpha=0,max_depth=2L,max_internal=3L,min_leaf=minleaf,
   penalty=.003),
 WT_SIGN_OOF=wt_fit(D,yd,parent,oof,coeff_mode="SIGN",alpha=.25,
   leaf_alpha=0,max_depth=2L,max_internal=3L,min_leaf=minleaf,
   penalty=.003))
tsp<-pool[seq_len(min(20L,nrow(pool))),
          c("gene_left","gene_right","ordered_pair"),drop=FALSE]
algs$LABEL_BINARY_TSP_TOP20<-dt_train(D,yd,tsp,oof,
   alpha=0,leaf_alpha=0,mode="BINARY",
   max_depth=2L,max_internal=3L,min_leaf=minleaf,penalty=.003)
fit_ktsp_train_test(D,yd,V,yv,k=5L,
  out_dir=file.path(out,"baseline_ktsp5"),seed=42L)
ktpred<-read.csv(file.path(out,"baseline_ktsp5","predictions_test.csv"),
                 stringsAsFactors=FALSE,check.names=FALSE)
ktpairs<-read.csv(file.path(out,"baseline_ktsp5","selected_pairs.csv"),
                  stringsAsFactors=FALSE,check.names=FALSE)
wt_need(identical(as.character(ktpred$sample_id),rownames(V))&&
  nrow(ktpairs)>=1L&&nrow(ktpairs)<=5L,
  "native switchBox K or V patient ID mismatch")
results<-lapply(names(algs),function(name){
  model<-algs[[name]]
  if(startsWith(name,"WT_")){
    p<-wt_predict(model,V)
    list(prob=p$prob,label=p$label,checks=p$n_evaluated,
         leaf=p$leaf,path=p$path)
  } else {
    p<-dt_predict(model,V)
    list(prob=p$prob,label=p$label,checks=p$depth,
         leaf=p$leaf,path=p$path)
  }
})
names(results)<-names(algs)
coef<-as.numeric(stats::coef(parent$fit,s="lambda.min"))[-1L]
results$RLR20_TEACHER<-list(
  prob=parent$test,label=as.integer(parent$test>=.5),
  checks=rep(sum(abs(coef)>1e-9),nrow(V)),
  leaf=rep(NA_integer_,nrow(V)),
  path=rep("ORIGINAL_GLOBAL_RLR20",nrow(V)))
results$KTSP5_NATIVE<-list(
  prob=as.numeric(ktpred$prob),
  label=as.integer(as.character(ktpred$y_pred)==levels(yv)[2L]),
  checks=rep(nrow(ktpairs),nrow(V)),
  leaf=rep(NA_integer_,nrow(V)),
  path=rep("REAL_NATIVE_SWITCHBOX_KTSP5",nrow(V)))
wt_need(length(results)==7L,"wrong fixed controlled method count")
f1<-function(h,cls) {
 tp<-sum(truth==cls & h==cls)
 fp<-sum(truth!=cls & h==cls)
 fn<-sum(truth==cls & h!=cls)
 if(2*tp+fp+fn==0)0 else 2*tp/(2*tp+fp+fn)
}
auc<-function(p) {
 n1<-sum(truth==1L);n0<-sum(truth==0L)
 wt_need(min(n1,n0)>0L,"V missing class for AUC")
 (sum(rank(p,ties.method="average")[truth==1L])-n1*(n1+1)/2)/(n1*n0)
}
tab_metrics<-list();tab_patients<-list();tab_nodes<-list()
i<-0L;j<-0L;k<-0L
for(method in names(results)){
 z<-results[[method]]
 fit<-algs[[method]]
 is_tree<-!is.null(fit)
 weighted<-is_tree && startsWith(method,"WT_")
 p<-as.numeric(z$prob);h<-as.integer(z$label)
 checks<-as.integer(z$checks)
 wt_need(length(p)==nrow(V)&&all(is.finite(p))&&
         all(p>=0&p<=1)&&length(checks)==nrow(V),
         "invalid heldout V predictions")
 if(weighted)wt_need(all(checks %in% c(0L,2L,4L)),
                     "actual two-pair nodes not audited")
 if(is_tree && !weighted)wt_need(all(checks>=0L&checks<=2L),
                               "TSP-like tree >2 path comparisons")
 nf<-if(is_tree)if(weighted)fit$n_internal else fit$internal else NA_integer_
 nl<-if(is_tree)fit$n_leaves else NA_integer_
 np<-if(is_tree)if(weighted)fit$n_unique_pairs else
       length(fit$used_pairs) else NA_integer_
 ng<-if(is_tree)fit$n_unique_genes else NA_integer_
 soft<-dt_fidelity(parent$test,p)
 i<-i+1L
 tab_metrics[[i]]<-data.frame(
  dataset=ds,method=method,n_D=nrow(D),n_V=nrow(V),
  positive_f1=f1(h,1L),
  macro_f1=mean(c(f1(h,0L),f1(h,1L))),
  ba=.5*(mean(h[truth==0L]==0L)+mean(h[truth==1L]==1L)),
  auc=auc(p),logloss=rg_ll(truth,p),
  brier=mean((truth-p)^2),
  mean_evaluated_pairs=mean(checks),
  max_evaluated_pairs=max(checks),
  n_internal_nodes=nf,n_leaves=nl,
  n_unique_relations=np,n_unique_model_genes=ng,
  teacher_decision_agreement=soft$agreement,
  teacher_soft_cross_entropy=soft$soft_logloss,
  teacher_fidelity_brier=soft$teacher_brier,
  training_scope="ORIGINAL_OUTER_T_TO_D_V_NO_E")
 for(idx in seq_len(nrow(V))){
   j<-j+1L
   tab_patients[[j]]<-data.frame(dataset=ds,method=method,
    sample_id=rownames(V)[idx],truth=truth[idx],prob=p[idx],
    label=h[idx],teacher_prob=parent$test[idx],
    n_evaluated=checks[idx],leaf=z$leaf[idx],
    actual_path=z$path[idx],stringsAsFactors=FALSE)
 }
 if(is_tree){
   for(nd in fit$nodes){
     if(weighted && !is.na(nd$split)){
       c<-fit$selector$candidates[nd$split,,drop=FALSE]
       desc<-paste0(fit$selector$pairs$ordered_pair[c$gene1[1]],
         "*",sprintf("%.8f",c$weight1[1])," + ",
         fit$selector$pairs$ordered_pair[c$gene2[1]],
         "*",sprintf("%.8f",c$weight2[1])," > 0")
     } else if(!weighted && !is.na(nd$split)){
       desc<-as.character(fit$pairs$ordered_pair[nd$split])
     } else desc<-""
     k<-k+1L
     tab_nodes[[k]]<-data.frame(dataset=ds,method=method,
       node_id=nd$id,depth=nd$depth,training_n=length(nd$idx),
       node_predicate=desc,children=paste(nd$children,collapse="|"),
       leaf_probability=nd$p,source="D_ONLY",stringsAsFactors=FALSE)
   }
 }
 cat(sprintf("RLRWT_DV_RESULT %s %s F1=%.6f MacroF1=%.6f AUC=%.6f LL=%.6f pairs=%.3f internal=%s\n",
   ds,method,tab_metrics[[i]]$positive_f1,tab_metrics[[i]]$macro_f1,
   tab_metrics[[i]]$auc,tab_metrics[[i]]$logloss,
   tab_metrics[[i]]$mean_evaluated_pairs,
   if(is_tree)as.character(nf)else"NA"))
}
write.csv(do.call(rbind,tab_metrics),file.path(out,"metrics_V.csv"),row.names=FALSE)
write.csv(do.call(rbind,tab_patients),file.path(out,"patient_V.csv"),row.names=FALSE)
write.csv(do.call(rbind,tab_nodes),file.path(out,"nodes_D.csv"),row.names=FALSE)
write.csv(data.frame(dataset=ds,n_total=nrow(X),
 n_outer_T=nrow(T),n_outer_E_UNREAD=length(outer$test_idx),
 n_D=nrow(D),n_V=nrow(V),
 D_ids=paste(rownames(D),collapse="|"),
 V_ids=paste(rownames(V),collapse="|"),
 RLR20_K=20L,teacher_oof_folds=3L,
 teacher_ridge_fallback_folds=n_fallback,
 max_depth=2L,max_internal=3L,
 min_leaf=minleaf,penalty=.003,
 coefficient_candidates_top8=8L,
 teacher_alpha_fixed=.25,
 source="SCHLIEP_RECONSTRUCTED_ORIGINAL_ARFF_SHA_UNVERIFIED",
 elapsed_seconds=proc.time()[["elapsed"]]-timer),
 file.path(out,"manifest.csv"),row.names=FALSE)
cat("RLRWT_ALL_SEVEN_METHODS_REAL_DV=PASS\n")
cat("RLRWT_RLR_BETA_VS_SIGN_AND_OOF_ABLATION=PASS\n")
cat("RLRWT_ORIGINAL_OUTER_E_UNREAD=PASS\n")
