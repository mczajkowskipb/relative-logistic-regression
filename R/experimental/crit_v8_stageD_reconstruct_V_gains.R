# Stage D: reconstruct actual V patient-specific candidate gains from frozen v7 RDS.
# Six candidate corrections per fold. No outer test labels or predictions used.
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v4_rank_repair.R")
source("R/experimental/crit_v7_original_rlr_two_leaf.R")
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)
  stop("STAGED_STOP: provide ORIGINAL_REF and CRIT7_REF",call.=FALSE)
REF<-normalizePath(args[[1L]],mustWork=TRUE)
V7<-normalizePath(args[[2L]],mustWork=TRUE)
OLD<-file.path(REF,"results/rlmt_dc_published_data_5fold_v1")
PRIOR<-file.path(V7,"results/crit_v7_original_rlr_2x3_v1")
OUT<-"results/crit_v8_stageD_null_simulation_v1"
need<-function(ok,msg)if(!isTRUE(ok))stop("STAGED_STOP: ",msg,call.=FALSE)
need(dir.exists(OUT)&&file.exists(file.path(OUT,"CODE_SHA.txt")) &&
     !file.exists(file.path(OUT,"candidate_V_gains.csv")),
     "missing result root or would overwrite previous reconstruction")
need(identical(trimws(readLines(file.path(OLD,"status.txt"))[[1L]]),
               "RUNNER_EXIT=0") &&
     identical(trimws(readLines(file.path(PRIOR,"status.txt"))[[1L]]),
               "CRIT7_PILOT_BATCH_EXIT=0"),
     "archived status not completed")
orig<-read.csv(file.path(OLD,"run_manifest.csv"),stringsAsFactors=FALSE)
inp<-read.csv(file.path(OLD,"input_manifest.csv"),stringsAsFactors=FALSE)
need(nrow(orig)==285L && nrow(inp)==19L,"original archive manifests incomplete")
verify_sha<-function(p,hash){
 line<-suppressWarnings(system2("sha256sum",
          shQuote(normalizePath(p,mustWork=TRUE)),stdout=TRUE))
 need(length(line)==1L &&
      grepl("^[a-f0-9]{64}[[:space:]]",line[[1L]]) &&
      identical(tolower(substr(line[[1L]],1L,64L)),tolower(hash)),
      "original dataset SHA mismatch")
}
results<-list()
meta<-list()
rr<-0L
mm<-0L
for(dataset in c("chen-2002","gordon-2002")){
  src<-inp[inp$dataset==dataset,,drop=FALSE]
  need(nrow(src)==1L,"missing original source dataset")
  verify_sha(src$source_path[[1L]],src$sha256[[1L]])
  xy<-read_xy_file(src$source_path[[1L]])
  folds<-make_repeated_stratified_folds(
               xy$y,n_folds=5L,n_repeats=1L,seed=42L)
  need(length(folds)==5L,"original split inventory unexpected")
  for(f in folds[1:3]){
    fit<-readRDS(file.path(PRIOR,"models",dataset,
                          f$split_id,"complete_crit7.rds"))
    need(!is.null(fit$confirm) &&
         is.data.frame(fit$confirm$confirmation) &&
         nrow(fit$confirm$confirmation)==6L,
         "expected six archived candidates per fold")
    T<-xy$X[f$train_idx,,drop=FALSE]
    Y<-xy$y[f$train_idx]
    idsD<-fit$partition$train
    idsV<-fit$partition$confirmation
    d<-match(idsD,rownames(T))
    v<-match(idsV,rownames(T))
    need(all(!is.na(d))&&all(!is.na(v))&&
         length(intersect(d,v))==0L &&
         identical(sort(c(d,v)),seq_len(nrow(T))) &&
         identical(rownames(T)[d],idsD) &&
         identical(rownames(T)[v],idsV),
         "archived D/V patient identifiers do not match original X")
    D<-T[d,,drop=FALSE]
    V<-T[v,,drop=FALSE]
    yd<-Y[d]
    yv<-Y[v]
    # Deterministic D-only original RLR fits exactly as archived v7.
    parentD<-rlmt4_parent(D,yd)
    dood<-rlmt4_strict_oof(D,yd)
    need(!is.null(dood) && length(dood$prob)==length(d),
         "genuine original D-only RLR OOF unavailable")
    pV<-rlmt4_parent_predict(parentD,V)$prob
    pool<-fit$confirm$candidate_pool
    cands<-fit$confirm$confirmation
    gD<-crit7_relations(D,pool)
    gV<-crit7_relations(V,pool)
    for(k in seq_len(nrow(cands))){
      row<-cands[k,,drop=FALSE]
      gate<-as.integer(row$gate_index)
      left<-as.integer(row$rel_false_index)
      right<-as.integer(row$rel_true_index)
      need(all(c(gate,left,right)>=1L) &&
           all(c(gate,left,right)<=nrow(pool)),
           "archived candidate relation index out of range")
      xd<-crit7_design(gD[,gate],gD[,left],gD[,right])
      xv<-crit7_design(gV[,gate],gV[,left],gV[,right])
      theta<-crit7_fit_offset(yd,dood$prob,xd,fit$cfg)
      p<-crit7_correction_prob(pV,xv,theta,fit$cfg)
      gains<-crit7_fold_gain(yv,pV,p)
      meanv<-mean(gains)
      sev<-stats::sd(gains)/sqrt(length(gains))
      lower<-meanv-fit$cfg$confirmation_SE_mult*sev-
                      fit$cfg$confirmation_cost
      err<-max(abs(c(
        meanv-as.numeric(row$confirmation_mean_gain),
        sev-as.numeric(row$confirmation_se),
        lower-as.numeric(row$confirmation_lower))))
      need(is.finite(err)&&err<1e-7,
           sprintf("saved v7 V score cannot be exactly reconstructed %s %s rank %d err=%g",
                   dataset,f$split_id,k,err))
      rr<-rr+1L
      results[[rr]]<-data.frame(
        dataset=dataset,split_id=f$split_id,
        D_rank=k,sample_id=rownames(V),
        y_true=as.character(yv),
        baseline_prob=pV,candidate_prob=p,
        patient_gain=gains,stringsAsFactors=FALSE)
      mm<-mm+1L
      meta[[mm]]<-data.frame(
        dataset=dataset,split_id=f$split_id,D_rank=k,
        n_D=length(d),n_V=length(v),
        confirmation_mean=meanv,confirmation_se=sev,
        confirmation_lower=lower,
        saved_score_parity_max_abs=err,stringsAsFactors=FALSE)
    }
    message("STAGED_RECONSTRUCTED ",dataset," ",f$split_id,
            " candidates=",nrow(cands)," n_V=",length(v))
  }
}
p<-do.call(rbind,results)
m<-do.call(rbind,meta)
need(nrow(m)==36L && nrow(p)==sum(m$n_V) &&
     max(m$saved_score_parity_max_abs)<1e-7,
     "V gain reconstruction inventory wrong")
write.csv(p,file.path(OUT,"candidate_V_gains.csv"),row.names=FALSE)
write.csv(m,file.path(OUT,"candidate_score_manifest.csv"),row.names=FALSE)
message("STAGED_ALL_36_CANDIDATE_PATIENT_GAIN_VECTORS=PASS")
message("STAGED_ARCHIVED_V7_SCORE_RECONSTRUCTION=PASS")
message("STAGED_BIOMEDICAL_EXPOSED_DEV_ONLY")
