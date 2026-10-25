# CRIT Stage G: 32-gene exhaustive-pool intervention on 48 independent train-test studies.
# 24 exact reused v9 studies + 24 freshly seeded DEVELOPMENT simulations.
# Method is deliberately NOT scalable to multi-omics.
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v4_rank_repair.R")
source("R/experimental/crit_v7_original_rlr_two_leaf.R")
source("R/experimental/crit_v9_d_crossfit_one_shot.R")
source("R/experimental/critG_exhaustive32_one_shot.R")

args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("STAGEG_STOP: provide CRIT9_WORKTREE",call.=FALSE)
V9<-normalizePath(args[[1L]],mustWork=TRUE)
OLD<-file.path(V9,"results/crit_v9_d_crossfit_2x3_v1")
OUT<-"results/crit_v9_stageG_exhaustive32_v1"
need<-function(ok,msg)if(!isTRUE(ok))stop("STAGEG_STOP: ",msg,call.=FALSE)
need(dir.exists(OUT)&&file.exists(file.path(OUT,"CODE_SHA.txt")) &&
     !file.exists(file.path(OUT,"study_summary.csv")),
     "result root absent or would overwrite")
need(identical(trimws(readLines(file.path(OLD,"status.txt"))[[1L]]),
               "CRIT9_BATCH_EXIT=0"),"v9 archive incomplete")
arch<-read.csv(file.path(OLD,"synthetic_panel_patients.csv"),
               stringsAsFactors=FALSE)
arch_summary<-read.csv(file.path(OLD,"synthetic_panel.csv"),
                       stringsAsFactors=FALSE)
need(nrow(arch)==17280L && nrow(arch_summary)==144L,
     "v9 synthetic archive incomplete")

generator<-function(n,seed,prefix,scenario) {
  set.seed(as.integer(seed))
  A<-matrix(stats::rnorm(n*32L),nrow=n,ncol=32L)
  colnames(A)<-sprintf("gene%02d",1:32)
  gate<-ifelse(A[,1L]>A[,2L],1,-1)
  rel<-ifelse(A[,3L]>A[,4L],1,-1)
  eta<-switch(scenario,null=rep(0,n),global_pair=1.6*rel,
         gate_only=1.6*gate,xor_interaction=2.4*gate*rel,
         stop("STAGEG_STOP: invalid scenario"))
  yy<-stats::runif(n)<stats::plogis(eta)
  y<-factor(ifelse(yy,"P","N"),levels=c("N","P"))
  X<-as.data.frame(A,check.names=FALSE)
  rownames(X)<-paste0(prefix,seq_len(n))
  list(X=X,y=y)
}
has_pair<-function(pool,a,b)
  any(pool$gene_a==a & pool$gene_b==b)
pair_is<-function(pair,a,b)
  !is.null(pair) && nrow(pair)==1L &&
    identical(as.character(pair$gene_a[[1L]]),a) &&
    identical(as.character(pair$gene_b[[1L]]),b)
scenarios<-c("null","global_pair","gate_only","xor_interaction")
rows<-list();patient_rows<-list();candidate_rows<-list()
si<-0L;pi<-0L;ci<-0L

for(s in seq_along(scenarios))for(rep in 1:12) {
  scenario<-scenarios[[s]]
  reuse<-rep<=6L
  seed<-if(reuse) as.integer(20261100L+s*1000L+rep*23L) else
    as.integer(20271100L+s*1000L+(rep-6L)*103L)
  group<-if(reuse) "REUSED_DEV" else "FRESH_DEV"
  train<-generator(180L,seed,
            sprintf("%s_r%02d_t_",scenario,rep),scenario)
  test<-generator(120L,seed+9L,
            sprintf("%s_r%02d_e_",scenario,rep),scenario)
  need(min(table(train$y))>=35L && min(table(test$y))>=12L,
       "inadequate independent study class support")
  fit<-critG_fit(train$X,train$y)
  pred<-critG_predict(fit,test$X)
  need(identical(names(pred$probs),critG_modes) &&
       identical(fit$model_identity,
         "CRITG_EXHAUSTIVE32_D_RANK_SINGLE_V_FIVE_ARM_ABLATION"),
       "ablation identity/model modes wrong")
  need(identical(fit$pool_size,496L) && fit$nomination_n<=6L,
       "full 496-pair pool or candidate shortlist violated")
  root<-as.numeric(pred$probs$ORIGINAL_RLR)
  if(reuse) {
    src<-arch[arch$scenario==scenario & arch$replicate==rep &
              arch$method=="ORIGINAL_RLR",,drop=FALSE]
    need(nrow(src)==120L &&
         identical(as.character(src$sample_id),rownames(test$X)) &&
         identical(as.character(src$y_true),as.character(test$y)) &&
         max(abs(root-as.numeric(src$prob)))<1e-10,
         "exact archived original RLR patient parity failed")
  }

  idsD<-fit$partition$discovery
  idsV<-fit$partition$confirmation
  d<-match(idsD,rownames(train$X))
  v<-match(idsV,rownames(train$X))
  need(all(!is.na(c(d,v))) &&
       !length(intersect(d,v)) &&
       identical(sort(c(d,v)),seq_len(nrow(train$X))),
       "D/V sample leakage")
  D<-train$X[d,,drop=FALSE]
  yd<-train$y[d]
  pd<-rlmt4_parent(D,yd)
  oldpool<-crit7_pool(D,pd$pairs,crit9_cfg())
  complete<-critG_complete32_pool(D,pd$pairs,critG_cfg())
  need(nrow(complete)==496L &&
       has_pair(complete,"gene01","gene02") &&
       has_pair(complete,"gene03","gene04") &&
       nrow(oldpool)<496L,
       "coverage intervention not isolated")

  # Frozen v9 shortlist proposals are recovered from the selected v9 source
  # candidate provenance only for the ABLATED model; no test-E information enters.
  cands<-fit$candidate_provenance
  need(length(cands)==fit$nomination_n &&
       nrow(fit$ranking)==fit$nomination_n,
       "candidate rank provenance missing")
  truth_gate<-vapply(cands,function(c)
           pair_is(c$gate,"gene01","gene02"),logical(1))
  truth_joint<-vapply(cands,function(c)
           pair_is(c$gate,"gene01","gene02") &&
           pair_is(c$interaction_pair,"gene03","gene04"),logical(1))
  truth_leaf<-vapply(cands,function(c)
           pair_is(c$gate,"gene01","gene02") &&
           pair_is(c$pair_false,"gene03","gene04") &&
           pair_is(c$pair_true,"gene03","gene04"),logical(1))
  selected<-fit$selected_D_rank
  need(length(selected)==1L && selected>=1L &&
       selected<=fit$nomination_n &&
       identical(as.integer(selected),
                 as.integer(which.max(fit$ranking$D_cv_lower))),
       "selected D candidate differs from D-only CV ranking")
  for(j in seq_len(fit$nomination_n)) {
     c<-cands[[j]]
     ci<-ci+1L
     candidate_rows[[ci]]<-data.frame(
       scenario=scenario,replicate=rep,group=group,D_rank=j,
       gate=paste(c$gate$gene_a,c$gate$gene_b,sep=">"),
       false_pair=paste(c$pair_false$gene_a,c$pair_false$gene_b,sep=">"),
       true_pair=paste(c$pair_true$gene_a,c$pair_true$gene_b,sep=">"),
       interaction_pair=paste(c$interaction_pair$gene_a,
                              c$interaction_pair$gene_b,sep=">"),
       planted_gate=as.integer(truth_gate[[j]]),
       planted_joint=as.integer(truth_joint[[j]]),
       planted_both_leaves=as.integer(truth_leaf[[j]]),
       chosen_by_D=as.integer(j==selected),
       D_cv_lower=fit$ranking$D_cv_lower[[j]],
       stringsAsFactors=FALSE)
  }

  truth_p<-as.integer(test$y=="P")
  logloss<-function(p) {
    p<-pmax(1e-15,pmin(1-1e-15,p))
    mean(-(truth_p*log(p)+(1-truth_p)*log1p(-p)))
  }
  base<-logloss(root)
  for(method in critG_modes) {
    p<-as.numeric(pred$probs[[method]])
    accepted<-if(method=="ORIGINAL_RLR")FALSE else
              isTRUE(fit$accepted[[method]])
    need(length(p)==length(root) && all(is.finite(p)) &&
         all(p>=0&p<=1),"invalid test probabilities")
    if(!accepted)
      need(identical(p,root),
           "exact baseline fallback failed")
    si<-si+1L
    rows[[si]]<-data.frame(
       scenario=scenario,replicate=rep,group=group,seed=seed,
       method=method,n_train=nrow(train$X),n_D=length(d),
       n_V=length(v),n_test=nrow(test$X),
       exhaustive_pair_count=nrow(complete),
       original_pool_size=nrow(oldpool),
       original_gate_present=as.integer(
          has_pair(oldpool,"gene01","gene02")),
       original_relation_present=as.integer(
          has_pair(oldpool,"gene03","gene04")),
       exhaustive_true_gate_present=1L,
       exhaustive_true_relation_present=1L,
       shortlisted_n=fit$nomination_n,
       shortlisted_planted_gate=as.integer(any(truth_gate)),
       shortlisted_planted_joint=as.integer(any(truth_joint)),
       shortlisted_planted_both_leaves=as.integer(any(truth_leaf)),
       D_selected_rank=selected,
       D_selected_planted_gate=as.integer(truth_gate[[selected]]),
       D_selected_planted_joint=as.integer(truth_joint[[selected]]),
       accepted=as.integer(accepted),
       original_logloss=base,
       logloss=logloss(p),
       delta_vs_original=logloss(p)-base,
       brier=mean((p-truth_p)^2),
       scientific_stage="SYNTHETIC_DEVELOPMENT_ONLY",
       stringsAsFactors=FALSE)
    pi<-pi+1L
    patient_rows[[pi]]<-data.frame(
       scenario=scenario,replicate=rep,group=group,
       method=method,sample_id=rownames(test$X),
       y_true=as.character(test$y),positive=truth_p,
       original_prob=root,prob=p,
       accepted=as.integer(accepted),stringsAsFactors=FALSE)
  }
  message("STAGEG_STUDY ",scenario," rep=",rep," ",group,
     " original_gate=",has_pair(oldpool,"gene01","gene02"),
     " full_pool=496 shortlist=",fit$nomination_n,
     " shortlist_truegate=",as.integer(any(truth_gate)),
     " shortlist_truejoint=",as.integer(any(truth_joint)),
     " D_choice=",selected,
     " selected_truejoint=",as.integer(truth_joint[[selected]]),
     " accepted_tree=",as.integer(
         isTRUE(fit$accepted[["CONDITIONAL_TREE"]])),
     " tree_deltaLL=",
     sprintf("%+.7f",logloss(pred$probs$CONDITIONAL_TREE)-base))
}
summary<-do.call(rbind,rows)
predictions<-do.call(rbind,patient_rows)
candidates<-do.call(rbind,candidate_rows)
need(nrow(summary)==48L*length(critG_modes) &&
     nrow(predictions)==48L*length(critG_modes)*120L &&
     nrow(candidates)==sum(summary$shortlisted_n[
                                 summary$method=="ORIGINAL_RLR"]) &&
     all(table(paste(summary$scenario,summary$group,
                     summary$method))==6L),
     "48-study model and patient inventory incorrect")
write.csv(summary,file.path(OUT,"study_summary.csv"),row.names=FALSE)
write.csv(predictions,file.path(OUT,"test_patients.csv"),row.names=FALSE)
write.csv(candidates,file.path(OUT,"candidate_provenance.csv"),row.names=FALSE)
message("STAGEG_48_STUDIES_288_METHODS=PASS")
message("STAGEG_34560_TEST_PATIENT_VECTORS=PASS")
message("STAGEG_EXHAUSTIVE_496_RELATIONS=PASS")
message("STAGEG_ORIGINAL_RLR_REUSED24_PATIENT_PARITY=PASS")
message("STAGEG_SCIENTIFIC_STAGE=DEVELOPMENT_ONLY")
