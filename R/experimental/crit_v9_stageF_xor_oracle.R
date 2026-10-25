# CRIT-v9 Stage F: mechanism diagnostics, D discovery vs known XOR oracle.
# 6 XOR + 6 null fresh synthetic datasets reused from frozen v9 panel.
# This is a retrospective DEVELOPMENT mechanism audit, not external validation.
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v4_rank_repair.R")
source("R/experimental/crit_v7_original_rlr_two_leaf.R")
source("R/experimental/crit_v9_d_crossfit_one_shot.R")
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("STAGEF_STOP: provide pinned CRIT9_CHECKOUT",call.=FALSE)
V9<-normalizePath(args[[1L]],mustWork=TRUE)
ARCH<-file.path(V9,"results/crit_v9_d_crossfit_2x3_v1")
OUT<-"results/crit_v9_stageF_xor_oracle_v1"
need<-function(x,msg)if(!isTRUE(x))stop("STAGEF_STOP: ",msg,call.=FALSE)
need(dir.exists(OUT)&&file.exists(file.path(OUT,"CODE_SHA.txt")) &&
     !file.exists(file.path(OUT,"mechanism_summary.csv")),
     "output absent or previously used")
need(identical(trimws(readLines(file.path(ARCH,"status.txt"))[[1L]]),
               "CRIT9_BATCH_EXIT=0"),"v9 archive was not completed")
archsum<-read.csv(file.path(ARCH,"synthetic_panel.csv"),stringsAsFactors=FALSE)
archpatients<-read.csv(file.path(ARCH,"synthetic_panel_patients.csv"),
                       stringsAsFactors=FALSE,
                       colClasses=c(sample_id="character",y_true="character"))
need(nrow(archsum)==144L && nrow(archpatients)==17280L,
     "frozen v9 synthetic archive incomplete")

# EXACT copied original v9 synthetic design: 32 genes, 180+120 patients.
make_xy<-function(n,seed,prefix,scenario) {
 set.seed(as.integer(seed))
 p<-32L
 A<-matrix(stats::rnorm(n*p),nrow=n,ncol=p)
 colnames(A)<-sprintf("gene%02d",seq_len(p))
 gate<-ifelse(A[,1L]>A[,2L],1,-1)
 pair<-ifelse(A[,3L]>A[,4L],1,-1)
 eta<-switch(scenario,
     null=rep(0,n),global_pair=1.6*pair,
     gate_only=1.6*gate,xor_interaction=2.4*gate*pair,
     stop("STAGEF_STOP: unknown scenario"))
 yy<-stats::runif(n)<stats::plogis(eta)
 y<-factor(ifelse(yy,"P","N"),levels=c("N","P"))
 X<-as.data.frame(A,check.names=FALSE)
 rownames(X)<-paste0(prefix,seq_len(n))
 list(X=X,y=y)
}
oracle_design<-function(X,method) {
 g<-ifelse(X[["gene01"]]>X[["gene02"]],1,-1)
 r<-ifelse(X[["gene03"]]>X[["gene04"]],1,-1)
 if(method=="oracle_tree")return(crit7_design(g,r,r))
 if(method=="oracle_global")return(crit7_interaction_design(g,r))
 stop("STAGEF_STOP: unknown oracle design")
}
ll<-function(y,p)crit7_logloss(y,p)
same_pair<-function(tbl,one,two) {
 if(is.null(tbl)||!nrow(tbl))return(FALSE)
 all(c(one,two) %in% c(as.character(tbl$gene_a[[1L]]),
                         as.character(tbl$gene_b[[1L]])))
}
p_rows<-list()
s_rows<-list()
cand_rows<-list()
n_record<-0L
n_summary<-0L
n_candidate<-0L
for(s in c(1L,4L)) {
 scenario<-c("null","global_pair","gate_only","xor_interaction")[[s]]
 for(rep in seq_len(6L)) {
   seed<-as.integer(20261100L+s*1000L+rep*23L)
   train<-make_xy(180L,seed,sprintf("%s_r%02d_t_",scenario,rep),
                  scenario)
   test<-make_xy(120L,seed+9L,sprintf("%s_r%02d_e_",scenario,rep),
                 scenario)
   X<-train$X
   y<-train$y
   parentcfg<-rlmt4_config()
   cfg<-crit9_cfg()
   model<-crit9_fit(X,y,cfg,parentcfg)
   learned<-crit9_predict(model,test$X)
   old<-archpatients[archpatients$scenario==scenario &
           archpatients$replicate==rep &
           archpatients$method=="ORIGINAL_RLR",,drop=FALSE]
   need(nrow(old)==120L &&
        identical(as.character(old$sample_id),rownames(test$X)) &&
        identical(as.character(old$y_true),as.character(test$y)) &&
        max(abs(as.numeric(old$prob)-
                    learned$probs$ORIGINAL_RLR))<1e-10,
        "synthetic replay vs pinned v9 RLR patient probabilities differs")
   archivedtree<-archsum[archsum$scenario==scenario &
                    archsum$replicate==rep &
                    archsum$method=="CONDITIONAL_TREE",,drop=FALSE]
   need(nrow(archivedtree)==1L &&
        max(abs(ll(test$y,learned$probs$CONDITIONAL_TREE)-
                archivedtree$logloss))<1e-10,
        "synthetic replay v9 conditional test metric differs")
   need(!is.null(model$partition),"no v9 D/V split available")
   d<-match(model$partition$discovery,rownames(X))
   v<-match(model$partition$confirmation,rownames(X))
   need(all(!is.na(c(d,v))) &&
        !length(intersect(d,v)) &&
        identical(sort(c(d,v)),seq_len(nrow(X))),
        "D/V split does not cover T exactly")
   D<-X[d,,drop=FALSE]; V<-X[v,,drop=FALSE]
   yd<-y[d]; yv<-y[v]
   pd<-rlmt4_parent(D,yd,parentcfg)
   dood<-rlmt4_strict_oof(D,yd,parentcfg)
   need(!is.null(dood),"D-only RLR OOF unavailable")
   pV<-rlmt4_parent_predict(pd,V)$prob
   selection<-crit9_rank_D_candidates(D,yd,dood$prob,pd,cfg)
   need(length(selection$candidates)==model$nomination_n &&
        identical(selection$ranking,model$ranking) &&
        identical(selection$selected_rank,model$selected_D_rank),
        "frozen candidate ranking did not replay deterministically")
   pool<-selection$pool
   truegate<-which(pool$gene_a=="gene01" & pool$gene_b=="gene02")
   truerel<-which(pool$gene_a=="gene03" & pool$gene_b=="gene04")
   n_pool<-nrow(pool)
   k<-length(selection$candidates)
   need(k<=6L,"D candidate budget exceeded 6")
   need(length(truegate)<=1L&&length(truerel)<=1L,
        "duplicate true feature in candidate pool")
   gate_hits<-integer()
   interaction_hits<-integer()
   bothleaf_hits<-integer()
   selected_gate<-FALSE
   selected_interaction<-FALSE
   for(j in seq_len(k)) {
      c<-selection$candidates[[j]]
      gate_match<-length(truegate)==1L &&
                    c$gate==truegate[[1L]]
      rel_match<-length(truerel)==1L &&
                   c$interaction_rel==truerel[[1L]]
      leaf_match<-length(truerel)==1L &&
                  c$rel_false==truerel[[1L]] &&
                  c$rel_true==truerel[[1L]]
      gate_hits<-c(gate_hits,as.integer(gate_match))
      interaction_hits<-c(interaction_hits,
                          as.integer(gate_match&&rel_match))
      bothleaf_hits<-c(bothleaf_hits,
                       as.integer(gate_match&&leaf_match))
      n_candidate<-n_candidate+1L
      cand_rows[[n_candidate]]<-data.frame(
        scenario=scenario,replicate=rep,D_rank=j,
        gate_gene_a=pool$gene_a[[c$gate]],
        gate_gene_b=pool$gene_b[[c$gate]],
        false_gene_a=pool$gene_a[[c$rel_false]],
        false_gene_b=pool$gene_b[[c$rel_false]],
        true_gene_a=pool$gene_a[[c$rel_true]],
        true_gene_b=pool$gene_b[[c$rel_true]],
        interaction_gene_a=pool$gene_a[[c$interaction_rel]],
        interaction_gene_b=pool$gene_b[[c$interaction_rel]],
        gate_is_planted=as.integer(gate_match),
        gate_interaction_is_planted=as.integer(gate_match&&rel_match),
        both_leaf_relations_planted=as.integer(gate_match&&leaf_match),
        selected_by_D=as.integer(j==selection$selected_rank),
        D_cv_lower=model$ranking$D_cv_lower[[j]],
        stringsAsFactors=FALSE)
   }
   if(!is.null(selection$selected_rank)) {
      selected_gate<-gate_hits[[selection$selected_rank]]==1L
      selected_interaction<-
         interaction_hits[[selection$selected_rank]]==1L
   }
   variants<-list(
      ORIGINAL_RLR=list(p=learned$probs$ORIGINAL_RLR,accepted=FALSE),
      DISCOVERED_CRIT9=list(p=learned$probs$CONDITIONAL_TREE,
         accepted=isTRUE(model$accepted[["CONDITIONAL_TREE"]])))
   parentT_oof<-rlmt4_strict_oof(X,y,parentcfg)
   need(!is.null(parentT_oof),"full T RLR strict OOF unavailable")
   oracle_stats<-list()
   for(method in c("oracle_tree","oracle_global")){
      designD<-oracle_design(D,method)
      designV<-oracle_design(V,method)
      designT<-oracle_design(X,method)
      designE<-oracle_design(test$X,method)
      betaD<-crit7_fit_offset(yd,dood$prob,designD,cfg)
      valp<-crit7_correction_prob(pV,designV,betaD,cfg)
      gains<-crit7_fold_gain(yv,pV,valp)
      lower<-mean(gains)-cfg$confirmation_SE_mult*
                    stats::sd(gains)/sqrt(length(gains))-
                       cfg$confirmation_cost
      accept<-is.finite(lower)&&lower>cfg$confirmation_min_gain
      betaT<-crit7_fit_offset(y,parentT_oof$prob,designT,cfg)
      allp<-crit7_correction_prob(learned$probs$ORIGINAL_RLR,
                                  designE,betaT,cfg)
      variants[[paste0(method,"_unconditional")]]<-
           list(p=allp,accepted=TRUE)
      variants[[paste0(method,"_V_accepted")]]<-
           list(p=if(accept)allp else learned$probs$ORIGINAL_RLR,
                accepted=accept)
      oracle_stats[[method]]<-list(lower=lower,accept=accept)
   }
   truth_cross_percentile<-NA_real_
   truth_cross_rank<-NA_integer_
   if(length(truegate)==1L && length(truerel)==1L) {
      relD<-crit7_relations(D,pool)
      residual<-rlmt_y01(yd)-dood$prob
      cross<-abs(crossprod(relD,relD*as.numeric(residual)))
      observed<-cross[truegate[[1L]],truerel[[1L]]]
      truth_cross_percentile<-mean(cross<=observed)
      truth_cross_rank<-as.integer(1L+sum(cross>observed))
   }
   for(method in names(variants)) {
      x<-variants[[method]]
      p<-as.numeric(x$p)
      need(length(p)==nrow(test$X)&&
           all(is.finite(p)) && all(p>=0&p<=1),
           "invalid heldout predictions")
      if(method %in% c("oracle_tree_V_accepted",
                       "oracle_global_V_accepted") && !x$accepted)
        need(identical(p,as.numeric(learned$probs$ORIGINAL_RLR)),
             "rejected oracle did not fall back bit-exact to original RLR")
      n_record<-n_record+1L
      p_rows[[n_record]]<-data.frame(
        scenario=scenario,replicate=rep,method=method,
        sample_id=rownames(test$X),
        y_true=as.character(test$y),
        positive=as.integer(test$y=="P"),
        prob=p,original_prob=learned$probs$ORIGINAL_RLR,
        accepted=as.integer(x$accepted),
        stringsAsFactors=FALSE)
   }
   base<-ll(test$y,learned$probs$ORIGINAL_RLR)
   for(method in names(variants)) {
      x<-variants[[method]]
      n_summary<-n_summary+1L
      s_rows[[n_summary]]<-data.frame(
        scenario=scenario,replicate=rep,seed=seed,
        method=method,n_train=nrow(X),n_D=nrow(D),n_V=nrow(V),
        n_test=nrow(test$X),pool_size=n_pool,
        true_gate_in_pool=as.integer(length(truegate)==1L),
        true_relation_in_pool=as.integer(length(truerel)==1L),
        true_cross_percentile=truth_cross_percentile,
        true_cross_rank=truth_cross_rank,
        shortlisted_candidates=k,
        shortlisted_true_gate=as.integer(any(gate_hits==1L)),
        shortlisted_true_gate_relation=
            as.integer(any(interaction_hits==1L)),
        shortlisted_true_both_leaves=
            as.integer(any(bothleaf_hits==1L)),
        selected_D_rank=if(is.null(model$selected_D_rank))NA_integer_ else
                        model$selected_D_rank,
        selected_true_gate=as.integer(selected_gate),
        selected_true_gate_relation=as.integer(selected_interaction),
        oracle_tree_V_lower=oracle_stats$oracle_tree$lower,
        oracle_tree_V_accept=as.integer(oracle_stats$oracle_tree$accept),
        oracle_global_V_lower=oracle_stats$oracle_global$lower,
        oracle_global_V_accept=as.integer(oracle_stats$oracle_global$accept),
        accepted=as.integer(x$accepted),
        original_logloss=base,
        logloss=ll(test$y,x$p),
        delta_vs_original=ll(test$y,x$p)-base,
        brier=mean((x$p-as.integer(test$y=="P"))^2),
        scientific_stage="REUSED_SYNTHETIC_DEV_MECHANISM_ONLY",
        stringsAsFactors=FALSE)
   }
   message("STAGEF_REPLAY ",scenario," rep=",rep,
      " pool=",n_pool," truepool=",length(truegate),",",length(truerel),
      " shortlist_gate=",as.integer(any(gate_hits==1L)),
      " D_rank=",model$selected_D_rank,
      " D_true=",as.integer(selected_interaction),
      " oracle_tree_accept=",as.integer(oracle_stats$oracle_tree$accept),
      " oracle_XOR_deltaLL=",
      sprintf("%+.6f",ll(test$y,
       variants$oracle_tree_unconditional$p)-base))
 }
}
p<-do.call(rbind,p_rows)
summary<-do.call(rbind,s_rows)
candidate<-do.call(rbind,cand_rows)
need(nrow(summary)==12L*6L &&
     nrow(p)==12L*6L*120L &&
     all(table(paste(summary$scenario,summary$replicate))==6L),
     "Stage F 12 by 6 inventory incomplete")
write.csv(candidate,file.path(OUT,"D_candidate_provenance.csv"),row.names=FALSE)
write.csv(summary,file.path(OUT,"mechanism_summary.csv"),row.names=FALSE)
write.csv(p,file.path(OUT,"mechanism_test_patients.csv"),row.names=FALSE)
message("STAGEF_12_PAIRED_SYNTHETIC_STUDIES=PASS")
message("STAGEF_8640_PATIENT_PREDICTIONS=PASS")
message("STAGEF_GENUINE_ORIGINAL_RLR_REPLAY=PASS")
message("STAGEF_SCIENTIFIC_STAGE=REUSED_SYNTHETIC_DEV_ONLY")
