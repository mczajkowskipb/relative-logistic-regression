# Original CompCancer 11-cohort, original 5fold x 10/2-repetition E protocol.
# One immutable (dataset, rep, fold) per independent CI artifact/job.
# Full T-only internal D/V structural choice, refit on ALL T with T-only
# crossfitted original RLR50 offsets; never use E outcomes for selection.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_crossfit_offset_ablation_v1.R")

ds<-Sys.getenv("RLR_DATASET_ID",unset="")
rep_no<-as.integer(Sys.getenv("RLR_REP",unset="0"))
fold_no<-as.integer(Sys.getenv("RLR_FOLD",unset="0"))
sets<-c("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
        "chen-2002","golub-1999-v1","gordon-2002",
        "pomeroy-2002-v1","shipp-2002-v1","singh-2002",
        "west-2001","yeoh-2002-v1")
crp_need(ds %in% sets && is.finite(rep_no) &&
         is.finite(fold_no) && fold_no %in% 0:4,
         "unknown cohort/repetition/fold")
nreps<-if(ds %in% c("west-2001","yeoh-2002-v1"))2L else 10L
crp_need(rep_no %in% 0:(nreps-1L),"repetition outside original 470 protocol")
split_id<-sprintf("rep%02d_fold%02d",rep_no,fold_no)
base<-file.path("results/rlr_full11_same_pairs_v1",ds)
out<-file.path("results/rlr_original11_outer_v1",ds,split_id)
src<-file.path(base,"source_data",paste0(ds,".csv"))
crp_need(file.exists(src)&&file.exists(file.path(base,"input_source_manifest.json")) &&
         file.exists(file.path(out,"SOURCE_SHA.txt")) &&
         !file.exists(file.path(out,"method_metrics.csv")),
         "missing public-source manifest/output SHA or duplicate execution")
df<-read.csv(src,check.names=FALSE,stringsAsFactors=FALSE)
X<-as.data.frame(df[,setdiff(names(df),c("sample_id","class")),
                    drop=FALSE],check.names=FALSE)
rownames(X)<-df$sample_id
y<-factor(df$class)
crp_need(length(levels(y))==2L &&
         !anyDuplicated(rownames(X)) &&
         !anyDuplicated(names(X)) &&
         all(vapply(X,function(z)is.numeric(z) && all(is.finite(z)),
                    logical(1))),"binary cohort/schema invalid")
specs<-make_repeated_stratified_folds(y,5L,nreps,42L)
slot<-which(vapply(specs,function(q)identical(q$split_id,split_id),logical(1)))
crp_need(length(slot)==1L,"original split missing")
sp<-specs[[slot]]
T<-X[sp$train_idx,,drop=FALSE]
E<-X[sp$test_idx,,drop=FALSE]
yt<-factor(as.character(y[sp$train_idx]),levels=levels(y))
ye<-factor(as.character(y[sp$test_idx]),levels=levels(y))
truth<-as.integer(ye==levels(ye)[2L])
crp_need(length(intersect(rownames(T),rownames(E)))==0 &&
         nrow(T)+nrow(E)==nrow(X) &&
         all(table(ye)>0L),"outer patient partition invalid")

start<-proc.time()[["elapsed"]]
# FIT all original comparators on exactly T; never reuse another method's
# learned normalizer, feature shortlist, kTSP pairs or class predictions.
poolT<-compute_top_pairs_switchbox(T,yt,top_n=1000L,buffer_mult=5L)
orig<-list(
  RLR10=rg_parent(T,yt,E,poolT,K=10L),
  RLR20=rg_parent(T,yt,E,poolT,K=20L),
  RLR50=rg_parent(T,yt,E,poolT,K=50L))
outbase<-file.path(out,"baselines")
dir.create(outbase,showWarnings=FALSE,recursive=TRUE)
for(k in c(10L,20L,50L)) {
   fit_ktsp_train_test(T,yt,E,ye,k=k,
      out_dir=file.path(outbase,paste0("ktsp",k)),seed=42L)
}
fit_tsp1_manual_train_test(T,yt,E,ye,poolT,
      out_dir=file.path(outbase,"tsp1"))
fit_value_lr_train_test(T,yt,E,ye,
      out_dir=file.path(outbase,"value_lr"),
      pre_top_features=1000L,alpha=.5,
      lambda_choice="lambda.min",seed=42L)

# T ONLY: propose and select model structure on its internal D/V.
iv<-rg_inner(yt)
D<-T[iv$D,,drop=FALSE]
V<-T[iv$V,,drop=FALSE]
yd<-factor(as.character(yt[iv$D]),levels=levels(yt))
yv<-factor(as.character(yt[iv$V]),levels=levels(yt))
state<-"EXACT_ORIGINAL_RLR_FALLBACK"
reason<-"NESTED_TRAIN_SUPPORT"
accepted<-FALSE
root<-NULL
nomination<-NULL
gain_parent<-gain_global<-deltaF1<-deltaBA<-NA_real_
roots_proposed<-roots_evaluated<-0L

if(safe_nfolds(yd,5L)>=4L) {
  oofplans<-make_repeated_stratified_folds(yd,3L,1L,731L)
  enough<-all(vapply(oofplans,function(p) {
     z<-factor(as.character(yd[p$train_idx]),levels=levels(yd))
     safe_nfolds(z,5L)>=4L
  },logical(1)))
  if(enough) {
    result<-cfo_run(D,yd,V,yv)
    roots_proposed<-result$proposed
    roots_evaluated<-length(result$roots)
    if(roots_evaluated>0L){
      yyv<-as.integer(yv==levels(yv)[2L])
      LL<-vapply(result$roots,function(a)
           rg_ll(yyv,a$prob$RESIDUAL_LOCAL_CROSSFIT_OFFSET),
           numeric(1))
      k<-which.min(LL)
      best<-result$roots[[k]]
      p0<-best$prob$ORIGINAL_RLR50
      pt<-best$prob$RESIDUAL_LOCAL_CROSSFIT_OFFSET
      pg<-best$prob$SAMEPAIRS_GLOBAL_CROSSFIT_OFFSET
      h0<-as.integer(p0>=.5)
      ht<-as.integer(pt>=.5)
      macro<-function(h){
        f1<-function(cl){
          tp<-sum(yyv==cl & h==cl)
          fp<-sum(yyv!=cl & h==cl)
          fn<-sum(yyv==cl & h!=cl)
          if(2*tp+fp+fn==0L)0 else 2*tp/(2*tp+fp+fn)
        }
        mean(c(f1(0L),f1(1L)))
      }
      ba<-function(h){
        .5*(mean(h[yyv==0L]==0L)+mean(h[yyv==1L]==1L))
      }
      gain_parent<-rg_ll(yyv,p0)-rg_ll(yyv,pt)
      gain_global<-rg_ll(yyv,pg)-rg_ll(yyv,pt)
      deltaF1<-macro(ht)-macro(h0)
      deltaBA<-ba(ht)-ba(h0)
      accepted<-isTRUE(gain_parent>=.005 &&
                   gain_global>=.002 &&
                   deltaBA>=-.02 &&
                   deltaF1>=-.01)
      # Root, threshold, and local gene PAIR IDs remain D-only forever.
      root<-best$root
      ppD<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
      pD<-sa_parent_oof(D,yd,K=50L,seed=731L)
      nomination<-cfo_nominate(D,yd,pD,root,ppD)
      state<-if(accepted)"TREE_ACCEPTED_ON_INTERNAL_V" else
             "TREE_REJECTED_ON_INTERNAL_V"
      reason<-"FIXED_D_OOF_V_LOGLOSS_SAMEPAIRS_BA_MACROF1_GATE"
    }else{
      reason<-"NO_ELIGIBLE_D_ONLY_ROOT"
    }
  }else{reason<-"CROSSFIT_D_SUBTRAIN_TOO_SMALL"}
}
# Refitting the accepted tree uses full-T OOF parent logits for its
# conditional coefficients. NO branch pair nomination on T is allowed.
tree_raw<-global_same<-orig$RLR50$test
n_active<-n_fitted<-0L
if(!is.null(root) && accepted) {
   if(safe_nfolds(yt,5L)>=4L) {
      pT_oof<-sa_parent_oof(T,yt,K=50L,seed=731L)
      final<-cfo_fit_frozen(T,yt,pT_oof,root,nomination)
      n_active<-final$n_active
      n_fitted<-final$n_fitted
      tree_raw<-sa_predict(E,orig$RLR50$test,final)
      model_all<-final
      model_all$pairs<-nomination
      global_same<-la_global_control(T,yt,pT_oof,
                   E,orig$RLR50$test,model_all)
   }else{
      accepted<-FALSE
      state<-"EXACT_ORIGINAL_RLR_FALLBACK"
      reason<-"OUTER_T_OOF_TOO_SMALL"
   }
}
tree_gated<-if(accepted)tree_raw else orig$RLR50$test
global_gated<-if(accepted)global_same else orig$RLR50$test

# Original switchBox kTSP native class decisions, not an unvalidated
# 0.5 reinterpretation of its vote score.
read_baseline<-function(sub) {
  a<-read.csv(file.path(outbase,sub,"predictions_test.csv"),
             stringsAsFactors=FALSE,check.names=FALSE)
  crp_need(identical(as.character(a$sample_id),rownames(E)),
           paste("real comparator patient order mismatch",sub))
  crp_need(all(as.character(a$y_pred) %in% levels(ye)),
           paste("invalid real comparator labels",sub))
  list(prob=as.numeric(a$prob),
       label=as.integer(as.character(a$y_pred)==levels(ye)[2L]))
}
P<-list(TSP1=read_baseline("tsp1"),
      KTSP10=read_baseline("ktsp10"),
      KTSP20=read_baseline("ktsp20"),
      KTSP50=read_baseline("ktsp50"),
      VALUE_LR=read_baseline("value_lr"))
for(nm in c("RLR10","RLR20","RLR50")){
  p<-orig[[nm]]$test
  P[[nm]]<-list(prob=p,label=as.integer(p>=.5))
}
P$CFO_TREE_V_GATED<-list(prob=tree_gated,
                        label=as.integer(tree_gated>=.5))
P$CFO_SAMEPAIRS_GLOBAL_V_GATED<-list(prob=global_gated,
                        label=as.integer(global_gated>=.5))
P$CFO_TREE_UNGATED_DIAGNOSTIC<-list(prob=tree_raw,
                        label=as.integer(tree_raw>=.5))
arms<-c("TSP1","KTSP10","KTSP20","KTSP50",
        "RLR10","RLR20","RLR50","VALUE_LR",
        "CFO_TREE_V_GATED","CFO_SAMEPAIRS_GLOBAL_V_GATED",
        "CFO_TREE_UNGATED_DIAGNOSTIC")
calc_f1<-function(h,cl){
  tp<-sum(truth==cl & h==cl)
  fp<-sum(truth!=cl & h==cl)
  fn<-sum(truth==cl & h!=cl)
  if(2*tp+fp+fn==0L)0 else 2*tp/(2*tp+fp+fn)
}
auc<-function(p){
  n1<-sum(truth==1L);n0<-sum(truth==0L)
  ranks<-rank(p,ties.method="average")
  (sum(ranks[truth==1L])-n1*(n1+1L)/2)/(n0*n1)
}
M<-list();R<-list();i<-0L;j<-0L
for(a in arms){
  p<-P[[a]]$prob;h<-P[[a]]$label
  crp_need(length(p)==nrow(E) && length(h)==nrow(E) &&
           all(is.finite(p)) && all(p>=0&p<=1) &&
           all(h %in% 0:1), paste("invalid E prediction",a))
  i<-i+1L
  M[[i]]<-data.frame(dataset=ds,split_id=split_id,
    method=a,n_train=nrow(T),n_E=nrow(E),
    positive_f1=calc_f1(h,1L),
    macro_f1=mean(c(calc_f1(h,0L),calc_f1(h,1L))),
    ba=.5*(mean(h[truth==0L]==0L)+mean(h[truth==1L]==1L)),
    auc=auc(p),logloss=rg_ll(truth,p),brier=mean((truth-p)^2),
    stage="PUBLIC_RECONSTRUCTION_ORIGINAL_ARFF_PARITY_UNVERIFIED")
  for(k in seq_along(truth)){
    j<-j+1L
    R[[j]]<-data.frame(dataset=ds,split_id=split_id,
        method=a,sample_id=rownames(E)[k],
        truth=truth[k],prob=p[k],label=h[k])
  }
}
write.csv(do.call(rbind,M),file.path(out,"method_metrics.csv"),
          row.names=FALSE)
write.csv(do.call(rbind,R),file.path(out,"patient_predictions.csv"),
          row.names=FALSE)
write.csv(data.frame(dataset=ds,split_id=split_id,
  n_total=nrow(X),n_outer_T=nrow(T),n_outer_E=nrow(E),
  n_inner_D=nrow(D),n_inner_V=nrow(V),
  T_ids=paste(rownames(T),collapse="|"),
  E_ids=paste(rownames(E),collapse="|"),
  original_positive_class=levels(ye)[2L],
  RLR50_selected=nrow(orig$RLR50$selected),
  state=state,reason=reason,accepted=as.integer(accepted),
  roots_proposed=roots_proposed,roots_evaluated=roots_evaluated,
  root_gene=if(is.null(root))"" else root$target,
  root_anchors=if(is.null(root))"" else paste(root$anchors,collapse="|"),
  root_threshold=if(is.null(root))NA_integer_ else root$threshold,
  D_pair_ids=if(is.null(nomination))"" else paste(unique(unlist(
    lapply(nomination,function(z)if(is.null(z))character()
      else z$ordered_pair),use.names=FALSE)),collapse="|"),
  refit_active_coefficients=n_active,refit_active_branches=n_fitted,
  V_parent_logloss_gain=gain_parent,V_global_logloss_gain=gain_global,
  V_macro_f1_change=deltaF1,V_ba_change=deltaBA,
  elapsed_seconds=proc.time()[["elapsed"]]-start,
  evidence="ORIGINAL_PAPER_NAMED_PUBLIC_SOURCE_RECONSTRUCTION_DEVELOPMENT"),
  file.path(out,"split_manifest.csv"),row.names=FALSE)
crp_need(i==length(arms)&&j==length(arms)*nrow(E),
         "incomplete original paper comparison")
cat(sprintf("CFO_OUTER_CASE %s %s accepted=%s state=%s real_kTSP50_F1=%.4f RLR50_F1=%.4f CFO_F1=%.4f RLR_LL=%.5f CFO_LL=%.5f elapsed=%.1fsec\n",
 ds,split_id,accepted,state,
 M[[4L]]$positive_f1,M[[7L]]$positive_f1,M[[9L]]$positive_f1,
 M[[7L]]$logloss,M[[9L]]$logloss,
 proc.time()[["elapsed"]]-start))
cat("CFO_OUTER_11_METHODS_TRUE_MATCHED_E_EVALUATION=PASS\n")
cat("CFO_OUTER_ALL_CLASSIFIER_HYPERPARAMETERS_TRAIN_ONLY=PASS\n")
cat("CFO_OUTER_ORIGINAL_PAPER_ARFF_PARITY_NOT_VERIFIED\n")
