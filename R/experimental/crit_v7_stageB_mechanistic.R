# CRIT-v7 Stage B: one-batch MECHANISM controls using archived 2x3 original folds.
# DEVELOPMENT only. New result root; no changes to CRIT-v7 source/results.
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v4_rank_repair.R")
source("R/experimental/crit_v7_original_rlr_two_leaf.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop("STAGEB_STOP: usage Rscript stageB.R ORIGINAL_RLR FIX3")
REF <- normalizePath(args[[1L]],mustWork=TRUE)
FIX3 <- normalizePath(args[[2L]],mustWork=TRUE)
OUT <- "results/crit_v7_stageB_mechanistic_v1"
OLD <- file.path(REF,"results/rlmt_dc_published_data_5fold_v1")
PRIOR <- file.path(FIX3,"results/crit_v7_original_rlr_2x3_v1")

need <- function(ok,msg) {
 if(!isTRUE(ok))stop("STAGEB_STOP: ",msg,call.=FALSE)
}
need(dir.exists(OUT) && file.exists(file.path(OUT,"CODE_SHA.txt")) &&
     !file.exists(file.path(OUT,"comparison_manifest.csv")),
     "output is missing or used")
need(identical(trimws(readLines(file.path(OLD,"status.txt"))[[1L]]),"RUNNER_EXIT=0"),
     "original reference status")
need(identical(trimws(readLines(file.path(PRIOR,"status.txt"))[[1L]]),
               "CRIT7_PILOT_BATCH_EXIT=0"),"fix3 pilot status")
orig <- read.csv(file.path(OLD,"run_manifest.csv"),stringsAsFactors=FALSE)
old_inputs <- read.csv(file.path(OLD,"input_manifest.csv"),stringsAsFactors=FALSE)
prior <- read.csv(file.path(PRIOR,"run_manifest.csv"),stringsAsFactors=FALSE)
need(nrow(orig)==285L && nrow(old_inputs)==19L && nrow(prior)==24L,
     "incorrect original or pilot inventories")

methods <- c("ORIGINAL_RLR","GLOBAL_INTERACTION","CRIT_TREE",
             "GATE_ONLY","ADDITIVE_PAIRS","UNCONDITIONAL_TREE")
all_rows <- list()
all_patient <- list()
ii <- 0L
jj <- 0L
for(dataset in c("chen-2002","gordon-2002")) {
 inp <- old_inputs[old_inputs$dataset==dataset,,drop=FALSE]
 need(nrow(inp)==1L,"missing dataset input")
 # Recheck the source-file byte identity, not just its cohort label.
 sha_line <- suppressWarnings(system2("sha256sum",
     shQuote(normalizePath(inp$source_path[[1L]],mustWork=TRUE)),
     stdout=TRUE))
 need(length(sha_line)==1L && grepl("^[a-f0-9]{64}[[:space:]]",sha_line[[1L]]) &&
      identical(tolower(substr(sha_line[[1L]],1L,64L)),
                tolower(inp$sha256[[1L]])),"original dataset SHA256 mismatch")
 # Original exact input data and original patient IDs are checked on every fold.
 xy <- read_xy_file(inp$source_path[[1L]])
 folds <- make_repeated_stratified_folds(
   xy$y,n_folds=5L,n_repeats=1L,seed=42L)
 need(length(folds)==5L,"original folds missing")
 for(f in folds[1:3]) {
   sid <- f$split_id
   tr <- f$train_idx
   te <- f$test_idx
   Xt <- xy$X[tr,,drop=FALSE]
   Yt <- xy$y[tr]
   Xv <- xy$X[te,,drop=FALSE]
   Yv <- xy$y[te]
   match <- prior[prior$dataset==dataset &
                  prior$split_id==sid &
                  prior$method=="CRIT_TREE",,drop=FALSE]
   om <- orig[orig$dataset==dataset & orig$split_id==sid &
              orig$method=="RLR",,drop=FALSE]
   need(nrow(match)==1L && nrow(om)==1L,"missing exact frozen manifest key")
   odir <- normalizePath(file.path(REF,om$model_dir[[1L]]),mustWork=TRUE)
   prdir <- file.path(PRIOR,"models",dataset,sid)
   need(file.exists(file.path(prdir,"complete_crit7.rds")),
        "missing exact stageB source model")
   original <- read.csv(file.path(odir,"predictions_test.csv"),
       stringsAsFactors=FALSE,colClasses=c(sample_id="character",
             y_true="character",y_pred="character"))
   model <- readRDS(file.path(prdir,"complete_crit7.rds"))
   pred <- crit7_predict(model,Xv)
   root <- as.numeric(pred$probs$ORIGINAL_RLR)
   need(identical(as.character(original$sample_id),rownames(Xv)) &&
        identical(as.character(original$y_true),as.character(Yv)),
        "original patient IDs/labels do not match")
   need(max(abs(root-original$prob))<1e-9,
        "archived original RLR is not reproduced")
   hasgate <- !is.null(model$gate) && !is.null(model$pair_false) &&
              !is.null(model$pair_true)
   need(hasgate,"StageB expects one nominated gate and two leaf relations in each fold")
   gate <- crit7_relations(Xt,model$gate)[,1L]
   gv <- crit7_relations(Xv,model$gate)[,1L]
   rf <- crit7_relations(Xt,model$pair_false)[,1L]
   rt <- crit7_relations(Xt,model$pair_true)[,1L]
   rfv <- crit7_relations(Xv,model$pair_false)[,1L]
   rtv <- crit7_relations(Xv,model$pair_true)[,1L]

   # Strict OUTER-TRAIN OOF reference; no outer-test outcomes enter estimation.
   oo <- rlmt4_strict_oof(Xt,Yt)
   need(!is.null(oo),"original parent strict train OOF unavailable")
   need(identical(names(oo$prob),NULL) || length(oo$prob)==length(Yt),
        "OOF shape invalid")

   leaf_train <- crit7_design(gate,rf,rt)
   leaf_test <- crit7_design(gv,rfv,rtv)
   gate_train <- leaf_train[,c("leaf_false","leaf_true"),drop=FALSE]
   gate_test <- leaf_test[,c("leaf_false","leaf_true"),drop=FALSE]
   additive_train <- cbind(gate_train,rel_false=rf,rel_true=rt)
   additive_test <- cbind(gate_test,rel_false=rfv,rel_true=rtv)

   theta_gate <- crit7_fit_offset(Yt,oo$prob,gate_train)
   theta_add <- crit7_fit_offset(Yt,oo$prob,additive_train)
   theta_tree <- crit7_fit_offset(Yt,oo$prob,leaf_train)
   source_tree <- pred$probs$CRIT_TREE
   probs <- list(
     ORIGINAL_RLR=root,
     GLOBAL_INTERACTION=as.numeric(pred$probs$GLOBAL_INTERACTION),
     CRIT_TREE=as.numeric(source_tree),
     GATE_ONLY=crit7_correction_prob(root,gate_test,theta_gate),
     ADDITIVE_PAIRS=crit7_correction_prob(root,additive_test,theta_add),
     UNCONDITIONAL_TREE=crit7_correction_prob(root,leaf_test,theta_tree)
   )
   need(identical(names(probs),methods),"incorrect method order")
   need(all(vapply(probs,function(x)
     length(x)==length(Yv)&&all(is.finite(x))&&
       all(x>=0&x<=1),logical(1))),"invalid probabilities")

   # A leaf-conditioned stump is exactly a single global feature
   # expansion with interactions; this is a numerical identity,
   # not evidence for a better model class.
   z <- as.numeric(gv>0)
   expanded <- cbind(
     intercept=rep(1,length(z)),
     gate=z,
     false_pair=rfv,
     gate_false_pair=z*rfv,
     gate_true_pair=z*rtv)
   t <- as.numeric(theta_tree)
   expanded_coef <- c(t[[1L]], t[[2L]]-t[[1L]],
                      t[[3L]],-t[[3L]],t[[4L]])
   equivalent <- crit7_correction_prob(root,expanded,expanded_coef)
   parity <- max(abs(equivalent-probs$UNCONDITIONAL_TREE))
   need(parity<1e-12,"tree vs explicit global GLM representation not equal")
   if(!is.null(model$tree_coef)) {
     need(max(abs(theta_tree-model$tree_coef))<1e-7 &&
          max(abs(probs$UNCONDITIONAL_TREE-source_tree))<1e-8,
          "accepted CRIT source and recomputed OOF coefficients disagree")
   } else {
     need(max(abs(source_tree-root))<1e-12,
          "rejected source tree should be original RLR")
   }

   positive <- levels(Yv)[[2L]]
   yy <- as.integer(Yv==positive)
   logloss <- function(p) {
     p <- pmax(1e-15,pmin(1-1e-15,p))
     mean(-(yy*log(p)+(1-yy)*log1p(-p)))
   }
   ref_ll <- logloss(root)
   for(method in methods) {
     p <- probs[[method]]
     ii <- ii+1L
     all_rows[[ii]] <- data.frame(
       dataset=dataset,split_id=sid,method=method,
       status="OK",n_train=length(tr),n_test=length(te),
       original_logloss=ref_ll,logloss=logloss(p),
       delta_vs_original=logloss(p)-ref_ll,
       brier=mean((p-yy)^2),
       tree_accepted=as.integer(!is.null(model$tree_coef)),
       gate_gene_a=as.character(model$gate$gene_a[[1L]]),
       gate_gene_b=as.character(model$gate$gene_b[[1L]]),
       leaf_false_gene_a=as.character(model$pair_false$gene_a[[1L]]),
       leaf_false_gene_b=as.character(model$pair_false$gene_b[[1L]]),
       leaf_true_gene_a=as.character(model$pair_true$gene_a[[1L]]),
       leaf_true_gene_b=as.character(model$pair_true$gene_b[[1L]]),
       validation_lower=model$confirm$best_confirmed_lower,
       matched_glm_parity=parity,stage="EXPOSED_RETROSPECTIVE_DEV_ONLY",
       stringsAsFactors=FALSE)
     for(k in seq_along(p)) {
       jj <- jj+1L
       all_patient[[jj]] <- data.frame(
         dataset=dataset,split_id=sid,method=method,
         sample_id=rownames(Xv)[[k]],
         y_true=as.character(Yv[[k]]),
         y_positive=yy[[k]],
         prob=p[[k]],original_prob=root[[k]],
         leaf_id=ifelse(gv[[k]]>0,2L,1L),
         tree_accepted=as.integer(!is.null(model$tree_coef)),
         stringsAsFactors=FALSE)
     }
   }
   cat("STAGEB_COMPLETED",dataset,sid,
       "tree_accepted=",as.integer(!is.null(model$tree_coef)),
       "equivalent_GLM_max_abs=",format(parity,digits=4),"\n")
 }
}
m <- do.call(rbind,all_rows)
p <- do.call(rbind,all_patient)
need(nrow(m)==36L && nrow(p)==sum(m$n_test) &&
     length(unique(paste(m$dataset,m$split_id,m$method)))==36L,
     "incomplete 6-fold by 6-method diagnostic")
write.csv(m,file.path(OUT,"comparison_manifest.csv"),row.names=FALSE)
write.csv(p,file.path(OUT,"patient_predictions.csv"),row.names=FALSE)
cat("STAGEB_36_METHOD_FOLDS=PASS\n")
cat("STAGEB_GLOBAL_GLM_EXACT_REPRESENTATION=PASS\n")
cat("STAGEB_SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY\n")
