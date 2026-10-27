# SA-RLR original-paper 11 cohort TRAINING-ONLY feasibility gate, one
# frozen original-paper 5-fold split per cohort; no outer E fit or assessment.
source("R/new/07_world_protocol_utils.R")
source("R/experimental/rlr_selective_anchor_gated_v1.R")

DS<-Sys.getenv("RLR_DATASET_ID",unset="")
EXPECTED<-c("alizadeh-2000-v1","armstrong-2002-v1",
 "bittner-2000","chen-2002","golub-1999-v1","gordon-2002",
 "pomeroy-2002-v1","shipp-2002-v1","singh-2002",
 "west-2001","yeoh-2002-v1")
sa_need(DS %in% EXPECTED,"dataset outside frozen 11")
base<-file.path("results/rlr_full11_same_pairs_v1",DS)
src<-file.path(base,"source_data",paste0(DS,".csv"))
sa_need(file.exists(src)&&file.exists(file.path(base,"input_source_manifest.json")),
        "expected source matrix/manifest absent")
out<-file.path("results/sa_rlr_trainonly_v1",DS)
sa_need(dir.exists(out) &&
    file.exists(file.path(out,"SOURCE_SHA.txt")) &&
    !file.exists(file.path(out,"train_only_selection.csv")),
    "isolated output directory not ready")
df<-utils::read.csv(src,check.names=FALSE,stringsAsFactors=FALSE)
X<-as.data.frame(df[,setdiff(names(df),c("sample_id","class")),
                      drop=FALSE],check.names=FALSE)
rownames(X)<-df$sample_id
y<-factor(df$class)
sa_need(length(levels(y))==2L && length(unique(rownames(X)))==nrow(X) &&
   length(unique(colnames(X)))==ncol(X) &&
   all(vapply(X,function(q)is.numeric(q)&&all(is.finite(q)),
                   logical(1))),"source not fit for rank model")
nrep<-if(DS %in% c("west-2001","yeoh-2002-v1"))2L else 10L
spec<-make_repeated_stratified_folds(y,5L,nrep,42L)[[1L]]
T<-X[spec$train_idx,,drop=FALSE]
yt<-factor(as.character(y[spec$train_idx]),levels=levels(y))
sa_need(nrow(T)+length(spec$test_idx)==nrow(X) &&
        !length(intersect(spec$train_idx,spec$test_idx)),
        "train/test patient overlap")
# The E rows are not passed to sa_select and no E endpoint is calculated.
inner<-rg_inner(yt)
D<-T[inner$D,,drop=FALSE]
yd<-factor(as.character(yt[inner$D]),levels=levels(yt))
pool<-compute_top_pairs_switchbox(D,yd,top_n=1000L,buffer_mult=5L)
proposed<-sa_candidate_roots(D,yd,pool,max_proposals=240L)
keys<-vapply(proposed,function(z)paste(z$g,collapse=""),character(1))
sa_need(!anyDuplicated(keys),"duplicate ordinal partition in root proposal")
swapped<-vapply(proposed,function(z) {
  conventional<-tryCatch(aa_anchors(aa_ranks(D),yd,pool,
                   z$root$target,length(z$root$anchors)),
                   error=function(e)NULL)
  if(is.null(conventional))FALSE else
        !setequal(conventional,z$root$anchors)
},logical(1))
supports<-vapply(proposed,function(z)sum(z$valid_branches),integer(1))
cat(sprintf("SA_TRAIN_PROGRESS dataset=%s T=%d D=%d V=%d p=%d root_partitions=%d root_swapped=%d one_branch=%d two_branch=%d\n",
  DS,nrow(T),length(inner$D),length(inner$V),ncol(T),
  length(proposed),sum(swapped),sum(supports==1L),sum(supports==2L)))
start<-proc.time()[["elapsed"]]
result<-tryCatch(sa_select(T,yt),
                 error=function(e) structure(list(message=conditionMessage(e)),
                                           class="sa_training_error"))
elapsed<-proc.time()[["elapsed"]]-start
if(inherits(result,"sa_training_error")){
   state<-"ERROR"
   accepted<-FALSE
   root<-NULL
   considered<-NA_integer_
   nfit<-NA_integer_
   active<-NA_integer_
   proposed_final<-NA_integer_
   vparent<-NA_real_;vglobal<-NA_real_;vba<-NA_real_
   reason<-result$message
}else{
   state<-if(is.null(result$root))"FALLBACK_NO_VALID_TREE"
          else if(isTRUE(result$accepted))"ACCEPTED_INNER_V"
          else "REJECTED_INNER_V"
   accepted<-isTRUE(result$accepted)
   root<-result$root
   considered<-result$considered
   nfit<-result$working
   active<-result$active
   proposed_final<-result$proposed
   vparent<-result$v_gain_parent
   vglobal<-result$v_gain_samepairs
   vba<-result$v_BA
   reason<-result$reason
}
record<-data.frame(dataset=DS,
  original_protocol_split=spec$split_id,
  n_original_total=nrow(X),n_outer_T=nrow(T),
  n_outer_E_UNREAD=length(spec$test_idx),
  n_inner_D=length(inner$D),n_inner_V=length(inner$V),
  n_genes=ncol(T),n_unique_proposed_ordinal_roots=length(proposed),
  n_one_swapped_panels=sum(swapped),
  n_roots_with_one_viable_branch=sum(supports==1L),
  n_roots_with_two_viable_branches=sum(supports==2L),
  selection_state=state,accepted_V=as.integer(accepted),
  n_roots_final_screen=considered,
  n_proposals_screen_recorded=proposed_final,
  selected_target=if(is.null(root))"" else root$target,
  selected_anchors=if(is.null(root))"" else
                          paste(root$anchors,collapse="|"),
  selected_t=if(is.null(root))NA_integer_ else root$threshold,
  n_working_experts_D=nfit,n_active_coefficients_D=active,
  V_logloss_gain_vs_RLR50=vparent,
  V_logloss_gain_vs_samepair_nontree=vglobal,
  V_BA_change=vba,
  elapsed_selection_seconds=elapsed,
  reason=reason,
  evidence_level="INNER_TRAINING_ONLY_NOT_EXTERNAL_VALIDATION")
utils::write.csv(record,file.path(out,"train_only_selection.csv"),
                 row.names=FALSE)
cat(sprintf("SA_TRAIN_RESULT %s state=%s V_accept=%s fits=%s active=%s V_parent=%s V_vs_global=%s elapsed=%.2fsec\n",
 DS,state,as.character(accepted),format(nfit),
 format(active),format(vparent,digits=5),
 format(vglobal,digits=5),elapsed))
cat("SA_TRAIN_ORIGINAL_5FOLD_SPLIT_ZERO_T_ONLY=PASS\n")
cat("SA_TRAIN_EXTERNAL_E_LABELS_NOT_USED=PASS\n")
# A scientific negative/error is still an auditable measured outcome. The
# aggregation process will flag any ERROR rather than label it a model win.
cat("SA_TRAIN_COHORT_DIAGNOSTIC_COMPLETE=PASS\n")
