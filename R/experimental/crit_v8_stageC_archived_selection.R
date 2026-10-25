# CRIT-v8 Stage C: purely retrospective **saved-RDS** validation-selection diagnostic.
# NO RLR training. NO use of outer-test labels for model/threshold selection.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop("STAGEC_STOP: ORIGINAL_REF CRIT7_REF CRIT8_REF required")
REF<-normalizePath(args[[1L]],mustWork=TRUE)
V7<-normalizePath(args[[2L]],mustWork=TRUE)
V8<-normalizePath(args[[3L]],mustWork=TRUE)
OUT<-"results/crit_v8_stageC_selection_v1"
d7<-file.path(V7,"results/crit_v7_original_rlr_2x3_v1")
d8<-file.path(V8,"results/crit_v8_one_shot_2x3_v1")
old<-file.path(REF,"results/rlmt_dc_published_data_5fold_v1")
need<-function(ok,msg)if(!isTRUE(ok))stop("STAGEC_STOP: ",msg,call.=FALSE)
need(dir.exists(OUT) && file.exists(file.path(OUT,"CODE_SHA.txt")) &&
     !file.exists(file.path(OUT,"selection_folds.csv")),
     "new output root missing or already used")
need(identical(trimws(readLines(file.path(old,"status.txt"))[[1L]]),
               "RUNNER_EXIT=0") &&
     identical(trimws(readLines(file.path(d7,"status.txt"))[[1L]]),
               "CRIT7_PILOT_BATCH_EXIT=0") &&
     identical(trimws(readLines(file.path(d8,"status.txt"))[[1L]]),
               "CRIT8_PILOT_BATCH_EXIT=0"),"source experiment statuses not PASS")
m7<-read.csv(file.path(d7,"run_manifest.csv"),stringsAsFactors=FALSE)
m8<-read.csv(file.path(d8,"run_manifest.csv"),stringsAsFactors=FALSE)
need(nrow(m7)==24L && nrow(m8)==36L,"unexpected previous manifest lengths")
getrow<-function(m,ds,fold,method) {
 r<-m[m$dataset==ds & m$split_id==fold & m$method==method,,drop=FALSE]
 need(nrow(r)==1L,paste("missing manifest",ds,fold,method))
 r
}
samepair<-function(a,b) {
 !is.null(a)&&!is.null(b)&&nrow(a)==1L&&nrow(b)==1L&&
 identical(as.character(a$gene_a[[1L]]),as.character(b$gene_a[[1L]]))&&
 identical(as.character(a$gene_b[[1L]]),as.character(b$gene_b[[1L]]))
}
candidate_rows<-list()
fold_rows<-list()
validation_rows<-list()
ci<-0L
fi<-0L
vi<-0L
for(ds in c("chen-2002","gordon-2002"))for(fold in sprintf("rep00_fold%02d",0:2)) {
 base7<-file.path(d7,"models",ds,fold)
 base8<-file.path(d8,"models",ds,fold)
 old7<-readRDS(file.path(base7,"complete_crit7.rds"))
 new8<-readRDS(file.path(base8,"complete_crit8.rds"))
 need(identical(old7$original$model_identity,"FROZEN_ORIGINAL_RLR_K10_M10_D2") &&
      identical(new8$original$model_identity,"FROZEN_ORIGINAL_RLR_K10_M10_D2"),
      "not genuine original RLR parent models")
 need(!is.null(old7$partition)&&!is.null(new8$partition)&&
      identical(as.character(old7$partition$train),
                as.character(new8$partition$discovery))&&
      identical(as.character(old7$partition$confirmation),
                as.character(new8$partition$confirmation)),
      "v7-v8 D/V inner patient partition mismatch")
 need(!length(intersect(new8$partition$discovery,
                        new8$partition$confirmation)),
      "D/V patients overlap")
 c7<-old7$confirm
 need(!is.null(c7)&&is.data.frame(c7$confirmation)&&
      nrow(c7$confirmation)>=1L && nrow(c7$confirmation)<=6L,
      "v7 confirmation has no candidates or exceeds expected six")
 need(identical(new8$nomination_n,1L) &&
      is.data.frame(new8$validation) &&
      nrow(new8$validation)==5L &&
      identical(as.character(new8$validation$method),
                c("GLOBAL_CAL","GATE_ONLY","ADDITIVE_PAIRS",
                  "GLOBAL_INTERACTION","CONDITIONAL_TREE")),
      "v8 one nominee/five controls violated")
 pool<-c7$candidate_pool
 need(is.data.frame(pool) && nrow(pool)>0L,"v7 pool not archived")
 cs<-c7$confirmation
 cols<-c("gate_index","rel_false_index","rel_true_index",
         "non_tree_interaction_index","confirmation_mean_gain",
         "confirmation_se","confirmation_lower")
 need(all(cols %in% names(cs)),"v7 confirmation candidates missing fields")
 candidate_pair<-function(i)pool[as.integer(i),c("gene_a","gene_b"),drop=FALSE]
 first<-cs[1L,,drop=FALSE]
 need(samepair(candidate_pair(first$gate_index),new8$gate)&&
      samepair(candidate_pair(first$rel_false_index),new8$pair_false)&&
      samepair(candidate_pair(first$rel_true_index),new8$pair_true)&&
      samepair(candidate_pair(first$non_tree_interaction_index),
               new8$interaction_pair),
      "v8 D-only nominee is not v7 D-rank first candidate")
 sel8<-new8$validation[new8$validation$method=="CONDITIONAL_TREE",,drop=FALSE]
 need(nrow(sel8)==1L && isTRUE(sel8$available[[1L]]),
      "v8 tree validation not available")
 eq<-abs(as.numeric(first$confirmation_lower)-
         as.numeric(sel8$validation_lower))
 need(is.finite(eq) && eq<1e-8,
      paste("v7 candidate rank-one and v8 lower score mismatch",ds,fold,eq))
 eqmean<-abs(as.numeric(first$confirmation_mean_gain)-
             as.numeric(sel8$validation_mean))
 eqse<-abs(as.numeric(first$confirmation_se)-
           as.numeric(sel8$validation_se))
 need(eqmean<1e-8 && eqse<1e-8,
      "candidate-1 v7/v8 validation mean/SE mismatch")
 need(all(is.finite(cs$confirmation_lower)),
      "nonfinite archived v7 candidate results")
 best<-which.max(cs$confirmation_lower)
 need(length(best)==1L && abs(max(cs$confirmation_lower)-
               c7$best_confirmed_lower)<1e-8,
      "v7 selected-max confirmation lower does not match archive")
 expect7<-max(cs$confirmation_lower)>old7$cfg$confirmation_min_gain
 accept7<-!is.null(old7$tree_coef)
 need(identical(accept7,expect7),
      "v7 tree acceptance differs from best-V statistic")
 expect8<-sel8$validation_lower[[1L]]>
                       new8$cfg$confirmation_min_gain
 accept8<-isTRUE(new8$accepted[["CONDITIONAL_TREE"]])
 need(identical(accept8,expect8),
      "v8 tree acceptance differs from first-D V-statistic")
 for(j in seq_len(nrow(cs))) {
   x<-cs[j,,drop=FALSE]
   g<-candidate_pair(x$gate_index)
   lf<-candidate_pair(x$rel_false_index)
   lt<-candidate_pair(x$rel_true_index)
   ci<-ci+1L
   candidate_rows[[ci]]<-data.frame(
      dataset=ds,split_id=fold,D_rank=j,
      candidate_gate=paste(g$gene_a,g$gene_b,sep=">"),
      candidate_false=paste(lf$gene_a,lf$gene_b,sep=">"),
      candidate_true=paste(lt$gene_a,lt$gene_b,sep=">"),
      discovery_score=as.numeric(x$discovery_score),
      confirmation_mean=as.numeric(x$confirmation_mean_gain),
      confirmation_se=as.numeric(x$confirmation_se),
      confirmation_lower=as.numeric(x$confirmation_lower),
      selected_by_v7=as.integer(j==best),
      first_D_nominee=as.integer(j==1L),
      stringsAsFactors=FALSE)
 }
 for(j in seq_len(nrow(new8$validation))) {
   x<-new8$validation[j,,drop=FALSE]
   method<-as.character(x$method[[1L]])
   v8row<-getrow(m8,ds,fold,method)
   need(as.integer(x$accepted[[1L]])==as.integer(v8row$accepted),
       "v8 RDS vs manifest acceptance mismatch")
   vi<-vi+1L
   validation_rows[[vi]]<-data.frame(
     dataset=ds,split_id=fold,method=method,
     available=as.integer(x$available),
     accepted=as.integer(x$accepted),
     validation_mean=as.numeric(x$validation_mean),
     validation_se=as.numeric(x$validation_se),
     confirmation_lower=as.numeric(x$validation_lower),
     outer_deltaLL=as.numeric(v8row$delta_vs_original),
     stringsAsFactors=FALSE)
 }
 v7row<-getrow(m7,ds,fold,"CRIT_TREE")
 v8row<-getrow(m8,ds,fold,"CONDITIONAL_TREE")
 need(as.integer(v7row$tree_accepted)==as.integer(accept7)&&
      as.integer(v8row$accepted)==as.integer(accept8),
      "archive manifests disagree with saved RDS")
 fi<-fi+1L
 fold_rows[[fi]]<-data.frame(
     dataset=ds,split_id=fold,n_D=length(new8$partition$discovery),
     n_V=length(new8$partition$confirmation),
     candidate_count=nrow(cs),v7_best_D_rank=best,
     first_D_equals_v7_best=as.integer(best==1L),
     v7_best_confirm_lower=as.numeric(c7$best_confirmed_lower),
     first_D_confirm_lower=as.numeric(first$confirmation_lower),
     v8_confirm_lower=as.numeric(sel8$validation_lower),
     v7_accepted=as.integer(accept7),v8_accepted=as.integer(accept8),
     v7_outer_deltaLL=as.numeric(v7row$delta_vs_original),
     v8_outer_deltaLL=as.numeric(v8row$delta_vs_original),
     v7_v8_first_candidate_parity_err=eq,
     stringsAsFactors=FALSE)
 cat("STAGEC_VERIFIED",ds,fold,
     "v7_best_D_rank=",best,"num_V_candidates=",nrow(cs),
     "v7_accept=",as.integer(accept7),
     "v8_accept=",as.integer(accept8),
     "v7_v8_rank1_err=",format(eq,digits=3),"\n")
}
a<-do.call(rbind,candidate_rows)
b<-do.call(rbind,fold_rows)
v<-do.call(rbind,validation_rows)
need(nrow(b)==6L && nrow(v)==30L &&
     nrow(a)>=6L && nrow(a)<=36L,
     "not complete 6fold x5-arm validation audit")
write.csv(a,file.path(OUT,"v7_candidate_rankings.csv"),row.names=FALSE)
write.csv(b,file.path(OUT,"selection_folds.csv"),row.names=FALSE)
write.csv(v,file.path(OUT,"v8_confirmation_arms.csv"),row.names=FALSE)
cat("STAGEC_V7_V8_D_FIRST_CANDIDATE_EXACT_PARITY=PASS\n")
cat("STAGEC_RDS_SELECTION_AUDIT=PASS\n")
cat("STAGEC_FOLDS=6 V8_VALIDATION_DECISIONS=30\n")
cat("STAGEC_EXPOSED_RETROSPECTIVE_DEV_ONLY\n")
