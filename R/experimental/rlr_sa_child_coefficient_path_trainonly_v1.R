# Diagnostics ONLY. No outer E and no inner V outcomes ever used.
# Does SA-RLR fail because child models cannot be fitted or because their
# penalized coefficients are zero? Test D-only coefficient paths, no tuning.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_selective_anchor_gated_v1.R")
DS<-Sys.getenv("RLR_DATASET_ID",unset="")
sa_need(DS %in% c("chen-2002","gordon-2002","singh-2002"),
 "diagnostic cohort outside predeclared list")
csv<-file.path("results/rlr_full11_same_pairs_v1",DS,
        "source_data",paste0(DS,".csv"))
sa_need(file.exists(csv),"missing reconstructed source")
dat<-read.csv(csv,check.names=FALSE)
X<-as.data.frame(dat[,setdiff(names(dat),c("sample_id","class")),
                      drop=FALSE],check.names=FALSE)
rownames(X)<-dat$sample_id
y<-factor(dat$class)
sp<-make_repeated_stratified_folds(y,5L,1L,42L)[[1L]]
T<-X[sp$train_idx,,drop=FALSE]
yt<-factor(as.character(y[sp$train_idx]),levels=levels(y))
inner<-rg_inner(yt)
D<-T[inner$D,,drop=FALSE]
yd<-factor(as.character(yt[inner$D]),levels=levels(y))
pool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
orig<-rg_parent(D,yd,D,pool,K=50L)
oof<-sa_parent_oof(D,yd,K=50L)
candidates<-sa_screen_roots(D,yd,oof,pool,max_top=6L)
rows<-list()
n<-0L
for(j in seq_along(candidates$roots)){
  entry<-candidates$roots[[j]]
  g<-os_gate(aa_ranks(D),entry$root)
  for(branch in 0:1){
    ix<-which(g==branch)
    xx<-D[ix,,drop=FALSE]
    yy<-factor(as.character(yd[ix]),levels=levels(y))
    if(length(ix)<12L||min(table(yy))<3L){
      n<-n+1L
      rows[[n]]<-data.frame(dataset=DS,proposal=j,branch=branch,
           n=length(ix),class0=sum(yy==levels(yy)[1L]),
           class1=sum(yy==levels(yy)[2L]),
           D_only_root_score=candidates$scores[j],
           pair_count=0L,n_ternary_nonconstant=0L,
           n_active_lambda_010=0L,n_active_lambda_003=0L,
           n_active_lambda_001=0L,
           state="BRANCH_TOO_SMALL_FOR_RLR")
      next
    }
    base<-la_fit_child(xx,yy,orig$train[ix],
                select_pairs=TRUE,lambda=.10)
    P<-base$pairs
    valid<-!is.null(P)&&nrow(P)>=2L
    if(valid){
      z<-encode_pair_matrix(xx,P,m=10L)
      nonconstant<-sum(apply(z,2L,stats::var)>0)
      alt03<-la_fit_child(xx,yy,orig$train[ix],pairs=P,
                         select_pairs=FALSE,lambda=.03)
      alt01<-la_fit_child(xx,yy,orig$train[ix],pairs=P,
                         select_pairs=FALSE,lambda=.01)
    }else {
      nonconstant<-0L
      alt03<-list(active=0L);alt01<-list(active=0L)
    }
    n<-n+1L
    rows[[n]]<-data.frame(dataset=DS,proposal=j,branch=branch,
      n=length(ix),class0=sum(yy==levels(yy)[1L]),
      class1=sum(yy==levels(yy)[2L]),
      D_only_root_score=candidates$scores[j],
      pair_count=if(valid)nrow(P) else 0L,
      n_ternary_nonconstant=nonconstant,
      n_active_lambda_010=base$active,
      n_active_lambda_003=alt03$active,
      n_active_lambda_001=alt01$active,
      state=base$state)
  }
}
out<-file.path("results/sa_rlr_local_diagnostic_v1",DS)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(length(rows)) {
  d<-do.call(rbind,rows)
  utils::write.csv(d,file.path(out,"D_only_child_path.csv"),row.names=FALSE)
  cat(sprintf("SA_LOCAL_DIAG %s roots=%d branches=%d viable=%d coefficients_lambda010=%d coefficients_lambda003=%d coefficients_lambda001=%d\n",
     DS,length(candidates$roots),nrow(d),
     sum(d$state!="BRANCH_TOO_SMALL_FOR_RLR"),
     sum(d$n_active_lambda_010),
     sum(d$n_active_lambda_003),
     sum(d$n_active_lambda_001)))
}else {
  write.csv(data.frame(dataset=DS,reason="NO_SCREENED_ROOTS"),
            file.path(out,"D_only_child_path.csv"),row.names=FALSE)
  cat(sprintf("SA_LOCAL_DIAG %s roots=0\n",DS))
}
cat("SA_LOCAL_DIAG_TRAINING_ONLY_NO_V_OR_E_USED=PASS\n")
