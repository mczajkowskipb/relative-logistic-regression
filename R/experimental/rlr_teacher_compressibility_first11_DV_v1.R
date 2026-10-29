# Exact, D-only RLR20 compressibility frontier on original paper named
# Schliep reconstructed cohort. NO original outer E accessed at evaluation.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_distilled_rank_tree_core_v1.R")
source("R/experimental/rlr_teacher_compressibility_core_v1.R")
ds<-Sys.getenv("RLR_DATASET_ID",unset="")
cohorts<-c("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
rc_need(ds %in% cohorts,"invalid paper cohort")
input<-file.path("results/rlr_full11_same_pairs_v1",ds,
                 "source_data",paste0(ds,".csv"))
out<-file.path("results/rlr_teacher_capacity11_DV_v1",ds)
rc_need(file.exists(input)&&file.exists(file.path(out,"SOURCE_SHA.txt")) &&
        !file.exists(file.path(out,"metrics_V.csv")),
        "missing verified public reconstruction or immutable output")
df<-read.csv(input,stringsAsFactors=FALSE,check.names=FALSE)
X<-as.data.frame(df[,setdiff(colnames(df),c("sample_id","class")),
                    drop=FALSE],check.names=FALSE)
rownames(X)<-df$sample_id
y<-factor(df$class)
rc_need(nlevels(y)==2L&&!anyDuplicated(rownames(X))&&
        all(vapply(X,function(z)is.numeric(z)&&all(is.finite(z)),
                   logical(1))),"original binary cohort schema")
reps<-if(ds %in% c("west-2001","yeoh-2002-v1"))2L else 10L
outer<-make_repeated_stratified_folds(y,5L,reps,42L)[[1L]]
T<-X[outer$train_idx,,drop=FALSE]
yt<-factor(as.character(y[outer$train_idx]),levels=levels(y))
iv<-rg_inner(yt)
D<-T[iv$D,,drop=FALSE]
V<-T[iv$V,,drop=FALSE]
yd<-factor(as.character(yt[iv$D]),levels=levels(y))
yv<-factor(as.character(yt[iv$V]),levels=levels(y))
truth<-as.integer(yv==levels(yv)[2L])
rc_need(nrow(D)+nrow(V)==nrow(T) &&
  !length(intersect(rownames(D),rownames(V))) &&
  !length(intersect(outer$train_idx,outer$test_idx)),"D V or E overlap")
pool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
tsp<-pool[seq_len(min(20L,nrow(pool))),
          c("gene_left","gene_right","ordered_pair"),drop=FALSE]
original<-rg_parent(D,yd,V,pool,K=20L)
oof<-dt_teacher_oof(D,yd,K=20L,folds=3L,seed=119L)
fallback<-as.integer(attr(oof,"ridge_fallback_folds"))
rc_need(length(oof)==nrow(D)&&all(is.finite(oof))&&
        all(oof>0&oof<1)&&length(fallback)==1L&&
        fallback>=0L && fallback<=3L,"OOF error")
ll<-function(t,p){
  pp<-pmax(1e-15,pmin(1-1e-15,p))
  mean(-t*log(pp)-(1-t)*log1p(-pp))
}
f1<-function(truth,h,cls){
 tp<-sum(truth==cls & h==cls);fp<-sum(truth!=cls & h==cls)
 fn<-sum(truth==cls & h!=cls)
 if(2*tp+fp+fn==0L)0 else 2*tp/(2*tp+fp+fn)
}
auc<-function(truth,p){
 n1<-sum(truth==1L);n0<-sum(truth==0L)
 rc_need(min(n1,n0)>0L,"invalid validation class support")
 (sum(rank(p,ties.method="average")[truth==1L])-n1*(n1+1)/2)/(n1*n0)
}
teacherV<-as.numeric(original$test)
rc_need(length(teacherV)==nrow(V)&&all(teacherV>0&teacherV<1),
        "original RLR20 V predictions invalid")
ell<-qlogis(pmin(.99999,pmax(.00001,oof)))
base_sse<-sum((ell-mean(ell))^2)
xD<-seq_len(nrow(D))
prior<-(sum(as.integer(yd==levels(yd)[2L]))+1)/(nrow(D)+2)
root_prob<-(sum(as.integer(yd==levels(yd)[2L]))+4*prior)/(nrow(D)+4)
base_label_ce<-sum(-as.numeric(yd==levels(yd)[2L])*log(root_prob)-
 (1-as.numeric(yd==levels(yd)[2L]))*log1p(-root_prob))
models<-list()
for(goal in c("TEACHER","LABEL"))for(cap in 0:3){
  name<-paste0(goal,"_CAP",cap)
  models[[name]]<-rc_fit(D,yd,tsp,oof,cap=cap,
    min_leaf=max(4L,as.integer(floor(.07*nrow(D)))),goal=goal)
}
rc_need(length(models)==8L,"fixed 8 arm frontier invalid")
M<-list();P<-list();i<-0L;j<-0L
for(name in names(models)){
 fit<-models[[name]]
 predV<-rc_predict(fit,V)
 predD<-rc_predict(fit,D)
 pl<-predV$prob_label;ps<-predV$prob_teacher
 h<-as.integer(pl>=.5)
 # Teacher projection loss on D, identically applied to BOTH topology
 # families, even for LABEL-only induced tree.
 oofSSE<-sum((ell-qlogis(predD$prob_teacher))^2)
 yd01<-as.integer(yd==levels(yd)[2L])
 labelCE<-sum(-yd01*log(predD$prob_label)-
            (1-yd01)*log1p(-predD$prob_label))
 i<-i+1L
 M[[i]]<-data.frame(dataset=ds,method=name,goal=fit$goal,
  cap=fit$cap,n_D=nrow(D),n_V=nrow(V),
  positive_f1=f1(truth,h,1L),
  macro_f1=mean(c(f1(truth,h,0L),f1(truth,h,1L))),
  ba=.5*(mean(h[truth==0L]==0L)+mean(h[truth==1L]==1L)),
  auc=auc(truth,pl),logloss=ll(truth,pl),
  brier=mean((truth-pl)^2),
  teacher_fidelity_brier=mean((teacherV-ps)^2),
  teacher_soft_cross_entropy=ll(teacherV,ps),
  teacher_agreement=mean(as.integer(teacherV>=.5)==as.integer(ps>=.5)),
  mean_checks=mean(predV$checks),max_checks=max(predV$checks),
  actual_nodes=fit$internal,n_leaves=fit$leaves,
  n_unique_pairs=fit$n_unique_pairs,n_unique_genes=fit$n_unique_genes,
  train_optimized_loss=fit$train_loss,
  train_optimized_fraction=fit$train_fraction,
  train_teacher_sse_fraction=if(base_sse>1e-12)oofSSE/base_sse else 0,
  train_true_label_ce_fraction=if(base_label_ce>1e-12)
    labelCE/base_label_ce else 0,
  scope="ORIGINAL_OUTER_FIRST_T_ONLY_D_V_NO_E",
  stringsAsFactors=FALSE)
 for(k in seq_len(nrow(V))){
  j<-j+1L
  P[[j]]<-data.frame(dataset=ds,method=name,
    sample_id=rownames(V)[k],truth=truth[k],
    prob_label=pl[k],prob_teacher=ps[k],
    teacherV=teacherV[k],label=h[k],
    checks=predV$checks[k],leaf=predV$leaf[k],
    path=predV$path[k],stringsAsFactors=FALSE)
 }
 cat(sprintf("RLR_CAPACITY_COHORT %s %s cap=%d used=%d train_teacher_fraction=%.6f heldout_teacher_Brier=%.6f F1=%.6f checks=%.4f\n",
 ds,name,fit$cap,fit$internal,M[[i]]$train_teacher_sse_fraction,
 M[[i]]$teacher_fidelity_brier,M[[i]]$positive_f1,
 M[[i]]$mean_checks))
}
# Parent comparator: teacher fidelity reference is exactly zero; not a tree
coef<-as.numeric(stats::coef(original$fit,s="lambda.min"))[-1L]
parent_checks<-sum(abs(coef)>1e-9)
ph<-as.integer(teacherV>=.5)
M[[length(M)+1L]]<-data.frame(dataset=ds,method="RLR20_REFERENCE",
 goal="REFERENCE",cap=-1L,n_D=nrow(D),n_V=nrow(V),
 positive_f1=f1(truth,ph,1L),
 macro_f1=mean(c(f1(truth,ph,0L),f1(truth,ph,1L))),
 ba=.5*(mean(ph[truth==0L]==0L)+mean(ph[truth==1L]==1L)),
 auc=auc(truth,teacherV),logloss=ll(truth,teacherV),
 brier=mean((truth-teacherV)^2),teacher_fidelity_brier=0,
 teacher_soft_cross_entropy=ll(teacherV,teacherV),
 teacher_agreement=1,mean_checks=parent_checks,max_checks=parent_checks,
 actual_nodes=NA_integer_,n_leaves=NA_integer_,
 n_unique_pairs=NA_integer_,n_unique_genes=NA_integer_,
 train_optimized_loss=NA_real_,train_optimized_fraction=NA_real_,
 train_teacher_sse_fraction=NA_real_,
 train_true_label_ce_fraction=NA_real_,
 scope="ORIGINAL_OUTER_FIRST_T_ONLY_D_V_NO_E")
for(k in seq_len(nrow(V)))P[[length(P)+1L]]<-data.frame(
 dataset=ds,method="RLR20_REFERENCE",sample_id=rownames(V)[k],
 truth=truth[k],prob_label=teacherV[k],prob_teacher=teacherV[k],
 teacherV=teacherV[k],label=ph[k],checks=parent_checks,
 leaf=NA_integer_,path="ORIGINAL_RLR20_REFERENCE",
 stringsAsFactors=FALSE)
metrics<-do.call(rbind,M);patient<-do.call(rbind,P)
rc_need(nrow(metrics)==9L&&nrow(patient)==9L*nrow(V),
 "fixed nine arm patient coverage failure")
write.csv(metrics,file.path(out,"metrics_V.csv"),row.names=FALSE)
write.csv(patient,file.path(out,"patient_V.csv"),row.names=FALSE)
write.csv(data.frame(dataset=ds,n_D=nrow(D),n_V=nrow(V),
 n_outer_T=nrow(T),n_outer_E_UNREAD=length(outer$test_idx),
 D_ids=paste(rownames(D),collapse="|"),
 V_ids=paste(rownames(V),collapse="|"),
 source="SCHLIEP_RECONSTRUCTED_ARFF_SHA_UNVERIFIED",
 positive_class=levels(yd)[2L],
 candidate_pool_D_only="SWITCHBOX_TOP20",
 teacher="ORIGINAL_RLR20_D_OOF_3FOLD",
 surrogate_OOF_folds=fallback,
 parent_nonzero_beta=parent_checks,
 max_depth=2L,max_internal=3L,
 min_leaf=max(4L,as.integer(floor(.07*nrow(D))))),
 file.path(out,"manifest.csv"),row.names=FALSE)
cat("RLR_CAPACITY_9_FIXED_ARMS_DV=PASS\n")
cat("RLR_CAPACITY_0_1_2_3_TRAIN_ONLY=PASS\n")
cat("RLR_CAPACITY_ORIGINAL_OUTER_E_UNREAD=PASS\n")
