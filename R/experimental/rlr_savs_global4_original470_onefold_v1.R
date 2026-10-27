# 11 paper-named CompCancer original 470 outer folds (9*10*5+2*2*5).
# Parsimonious NON-TREE selective signed abstaining rank votes (4 rules).
# Two PRE-DECLARED margin variants vs exactly matched actual SWAP.KTSP
# K4/10/20/50, RLR10/20/50, original TSP1 and value LR.
# All training, panel selection, fitted priors inside outer T only.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_selective_anchor_vote_stump_v1.R")

ds<-Sys.getenv("RLR_DATASET_ID",unset="")
repeat_id<-as.integer(Sys.getenv("RLR_REP",unset="0"))
fold_id<-as.integer(Sys.getenv("RLR_FOLD",unset="0"))
cohorts<-c("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
sv_need(ds %in% cohorts && is.finite(repeat_id) &&
         is.finite(fold_id) && fold_id %in% 0:4,"bad original cohort/split")
nrep<-if(ds %in% c("west-2001","yeoh-2002-v1"))2L else 10L
sv_need(repeat_id %in% 0:(nrep-1L),"outside original repeated CV")
sid<-sprintf("rep%02d_fold%02d",repeat_id,fold_id)
base<-file.path("results/rlr_full11_same_pairs_v1",ds)
src<-file.path(base,"source_data",paste0(ds,".csv"))
out<-file.path("results/savs_global4_full470_v1",ds,sid)
sv_need(file.exists(src)&&file.exists(file.path(base,"input_source_manifest.json")) &&
         file.exists(file.path(out,"SOURCE_SHA.txt"))&&
         !file.exists(file.path(out,"method_metrics.csv")),
         "source manifest or output mismatch")
dat<-read.csv(src,stringsAsFactors=FALSE,check.names=FALSE)
X<-as.data.frame(dat[,setdiff(names(dat),c("sample_id","class")),
                      drop=FALSE],check.names=FALSE)
rownames(X)<-dat$sample_id
y<-factor(dat$class)
sv_need(length(levels(y))==2L && !anyDuplicated(rownames(X)) &&
       all(vapply(X,function(z)is.numeric(z)&&all(is.finite(z)),logical(1))),
       "invalid original cohort input")
spec<-make_repeated_stratified_folds(y,5L,nrep,42L)
sel<-which(vapply(spec,function(z)identical(z$split_id,sid),logical(1)))
sv_need(length(sel)==1L,"bad split ID")
sp<-spec[[sel]]
T<-X[sp$train_idx,,drop=FALSE];E<-X[sp$test_idx,,drop=FALSE]
yt<-factor(as.character(y[sp$train_idx]),levels=levels(y))
ye<-factor(as.character(y[sp$test_idx]),levels=levels(y))
yy<-as.integer(ye==levels(ye)[2L])
sv_need(!length(intersect(rownames(T),rownames(E)))&&
         nrow(T)+nrow(E)==nrow(X)&&min(table(ye))>0L,
         "patient leakage or binary cohort no support")
started<-proc.time()[["elapsed"]]
pool<-compute_top_pairs_switchbox(T,yt,top_n=1000L,buffer_mult=5L)
rmodels<-lapply(c(10L,20L,50L),function(k)rg_parent(T,yt,E,pool,K=k))
names(rmodels)<-paste0("RLR",c(10L,20L,50L))
models<-list()
for(mode in c("RLR10","RELATIVE10PCT")){
 fit<-sv_global4(T,yt,pool,sv_m(ncol(T),mode))
 q<-sv_predict(fit,E)
 nm<-if(mode=="RLR10")"SELECTIVE4_M10" else "SELECTIVE4_M10PCT"
 models[[nm]]<-list(prob=q$p,label=as.integer(q$p>=.5),
                    n_active=q$n_active,n_considered=q$n_considered,
                    panel_size=nrow(fit$rule_panels[[1L]]),
                    gene_count=length(unique(c(fit$rule_panels[[1L]]$gene_left,
                                               fit$rule_panels[[1L]]$gene_right))),
                    rules=fit$rule_panels[[1L]],
                    margin=fit$margin)
}
for(k in c(4L,10L,20L,50L))
 fit_ktsp_train_test(T,yt,E,ye,k=k,
     out_dir=file.path(out,"baselines",paste0("ktsp",k)),seed=42L)
fit_tsp1_manual_train_test(T,yt,E,ye,pool,
     out_dir=file.path(out,"baselines","tsp1"))
fit_value_lr_train_test(T,yt,E,ye,
     out_dir=file.path(out,"baselines","value_lr"),
     pre_top_features=1000L,alpha=.5,
     lambda_choice="lambda.min",seed=42L)
read_comparator<-function(name) {
 path<-file.path(out,"baselines",name,"predictions_test.csv")
 z<-read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
 sv_need(identical(as.character(z$sample_id),rownames(E)),
         paste("baseline E patient mismatch",name))
 sv_need(all(z$y_pred %in% levels(ye)),paste("invalid label",name))
 list(prob=as.numeric(z$prob),
      label=as.integer(as.character(z$y_pred)==levels(ye)[2L]))
}
for(k in c(4L,10L,20L,50L)){
 nm<-paste0("KTSP",k)
 x<-read_comparator(paste0("ktsp",k))
 realized<-nrow(read.csv(file.path(out,"baselines",
                paste0("ktsp",k),"selected_pairs.csv")))
 sv_need(realized>=1L && realized<=k,
        "real switchBox kTSP panel outside nominal K")
 x$n_active<-rep(realized,nrow(E))
 x$n_considered<-rep(realized,nrow(E))
 x$panel_size<-realized
 models[[nm]]<-x
}
tsp<-read_comparator("tsp1")
tsp$n_active<-rep(1L,nrow(E));tsp$n_considered<-rep(1L,nrow(E))
tsp$panel_size<-1L;models$TSP1<-tsp
x<-read_comparator("value_lr")
x$n_active<-rep(NA_integer_,nrow(E))
x$n_considered<-rep(NA_integer_,nrow(E))
x$panel_size<-NA_integer_;models$VALUE_LR<-x
for(nm in names(rmodels)){
  z<-rmodels[[nm]]
  pp<-z$test
  coeff<-as.numeric(stats::coef(z$fit,s="lambda.min"))[-1L]
  R<-encode_pair_matrix(E,z$selected,m=10L)
  sv_need(length(coeff)==ncol(R),"RLR coefficient/pair mismatch")
  live<-abs(coeff)>1e-9
  n_active<-if(any(live))
     rowSums(R[,live,drop=FALSE]!=0) else rep(0L,nrow(E))
  models[[nm]]<-list(prob=pp,label=as.integer(pp>=.5),
     n_active=as.integer(n_active),
     n_considered=rep(sum(live),nrow(E)),
     panel_size=nrow(z$selected),
     n_global_active=sum(live))
}
arms<-c("TSP1","KTSP4","KTSP10","KTSP20","KTSP50",
        "RLR10","RLR20","RLR50","VALUE_LR",
        "SELECTIVE4_M10","SELECTIVE4_M10PCT")
positive_f1<-function(h,cl){
 tp<-sum(yy==cl&h==cl);fp<-sum(yy!=cl&h==cl)
 fn<-sum(yy==cl&h!=cl)
 if(2*tp+fp+fn==0L)0 else 2*tp/(2*tp+fp+fn)
}
auc<-function(prob){
 n1<-sum(yy==1L);n0<-sum(yy==0L)
 rr<-rank(prob,ties.method="average")
 (sum(rr[yy==1L])-n1*(n1+1L)/2)/(n0*n1)
}
rows<-list();predictions<-list();idx<-0L;idxP<-0L
for(name in arms){
  x<-models[[name]]
  pr<-as.numeric(x$prob);h<-as.integer(x$label)
  sv_need(length(pr)==nrow(E)&&length(h)==nrow(E) &&
       all(is.finite(pr)) && all(pr>=0&pr<=1) &&
       all(h %in% 0:1),paste("prediction invalid",name))
  idx<-idx+1L
  rows[[idx]]<-data.frame(
    dataset=ds,split_id=sid,method=name,
    n_train=nrow(T),n_test=nrow(E),
    positive_f1=positive_f1(h,1L),
    macro_f1=mean(c(positive_f1(h,0L),positive_f1(h,1L))),
    ba=.5*(mean(h[yy==0L]==0L)+mean(h[yy==1L]==1L)),
    auc=auc(pr),logloss=rg_ll(yy,pr),
    brier=mean((yy-pr)^2),
    nominal_panel=x$panel_size,
    mean_active=if(all(is.na(x$n_active)))NA_real_ else mean(x$n_active),
    mean_tested=if(all(is.na(x$n_considered)))NA_real_ else mean(x$n_considered),
    frac_active_le2=if(all(is.na(x$n_active)))NA_real_ else mean(x$n_active<=2L),
    class_source="REAL_NATIVE_KTSP_WHEN_APPLICABLE",
    stage="EXPLORATORY_PUBLIC_RECONSTRUCTION_ORIGINAL_ARFF_UNVERIFIED")
  for(i in seq_along(pr)){
    idxP<-idxP+1L
    predictions[[idxP]]<-data.frame(
      dataset=ds,split_id=sid,method=name,
      sample_id=rownames(E)[i],truth=yy[i],
      prob=pr[i],label=h[i],
      n_active=x$n_active[i],n_tested=x$n_considered[i])
  }
}
write.csv(do.call(rbind,rows),file.path(out,"method_metrics.csv"),
          row.names=FALSE)
write.csv(do.call(rbind,predictions),
          file.path(out,"patient_predictions.csv"),row.names=FALSE)
met<-data.frame(dataset=ds,split_id=sid,
 n_total=nrow(X),n_T=nrow(T),n_E=nrow(E),
 T_ids=paste(rownames(T),collapse="|"),
 E_ids=paste(rownames(E),collapse="|"),
 margin_m10=sv_m(ncol(T),"RLR10"),
 margin_m10pct=sv_m(ncol(T),"RELATIVE10PCT"),
 gene_features=ncol(T),
 trained_rules_m10=nrow(models$SELECTIVE4_M10$rules),
 trained_rules_m10pct=nrow(models$SELECTIVE4_M10PCT$rules),
 unique_genes_m10=models$SELECTIVE4_M10$gene_count,
 unique_genes_m10pct=models$SELECTIVE4_M10PCT$gene_count,
 elapsed_seconds=proc.time()[["elapsed"]]-started,
 evidence="ALL_FEATURE_SELECTION_OUTER_T_ONLY_ORIGINAL_E_TEST_ONLY")
write.csv(met,file.path(out,"manifest.csv"),row.names=FALSE)
sv_need(idx==11L&&idxP==11L*nrow(E),"not exactly 11 original benchmark arms")
cat(sprintf("SAVS470_CASE %s %s kTSP4_F1=%.4f kTSP50_F1=%.4f RLR50_F1=%.4f four_M10_F1=%.4f four_M10pct_F1=%.4f mean_active_M10pct=%.3f\n",
 ds,sid,rows[[2L]]$positive_f1,rows[[5L]]$positive_f1,
 rows[[8L]]$positive_f1,rows[[10L]]$positive_f1,
 rows[[11L]]$positive_f1,rows[[11L]]$mean_active))
cat("SAVS470_MATCHED_E_TRUE_KTSP4_KTSP50_RLR50_GLOBAL4=PASS\n")
cat("SAVS470_RANK_ONLY_TRAIN_SELECTED_CLASSIFIER_NO_OUTER_LEAK=PASS\n")
