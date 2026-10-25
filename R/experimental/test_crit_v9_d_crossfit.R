# Strict synthetic structural tests for CRIT-v8. No biomedical samples.
for(pkg in c("switchBox","glmnet","pROC","jsonlite"))
 if(!requireNamespace(pkg,quietly=TRUE))
   stop("CRIT9_SYNTH_STOP: missing package ",pkg)
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v4_rank_repair.R")
source("R/experimental/crit_v7_original_rlr_two_leaf.R")
source("R/experimental/crit_v9_d_crossfit_one_shot.R")
need <- function(ok,msg) {
 if(!isTRUE(ok))stop("CRIT9_SYNTH_STOP: ",msg,call.=FALSE)
 message("PASS ",msg)
}
make_xy <- function(n,seed,prefix,mode="interaction") {
 set.seed(seed)
 p <- 32L
 m <- matrix(stats::rnorm(n*p),nrow=n,ncol=p)
 colnames(m) <- sprintf("gene%02d",seq_len(p))
 gate <- ifelse(m[,1]>m[,2],1,-1)
 relation <- ifelse(m[,3]>m[,4],1,-1)
 eta <- if(mode=="interaction")2.7*gate*relation else
         if(mode=="global")1.4*relation else
            rep(0,n)
 y <- factor(ifelse(runif(n)<plogis(eta),"P","N"),
             levels=c("N","P"))
 X <- as.data.frame(m,check.names=FALSE)
 rownames(X) <- paste0(prefix,seq_len(n))
 list(X=X,y=y)
}
train<-make_xy(240L,4181L,"train_")
test<-make_xy(90L,4182L,"test_")
need(min(table(train$y))>=40L && min(table(test$y))>=15L,
     "synthetic class support")
fit<-crit9_fit(train$X,train$y)
pred<-crit9_predict(fit,test$X)
need(identical(names(pred$probs),crit9_modes) &&
     identical(fit$original$model_identity,
       "FROZEN_ORIGINAL_RLR_K10_M10_D2"),
     "six model variants with exact original RLR parent")
need(!is.null(fit$partition) &&
     !length(intersect(fit$partition$discovery,
                       fit$partition$confirmation)) &&
     length(fit$partition$discovery)+
       length(fit$partition$confirmation)==nrow(train$X),
     "D/V partition disjoint inside outer training")
need(fit$nomination_n<=6L &&
     identical(fit$cfg$shortlist,6L),
     "up to six D-only candidates, exactly one routed to V")
need((fit$nomination_n==0L && is.null(fit$selected_D_rank)) ||
     (fit$nomination_n>0L &&
      nrow(fit$ranking)==fit$nomination_n &&
      fit$selected_D_rank>=1L &&
      fit$selected_D_rank<=fit$nomination_n &&
      identical(as.integer(fit$selected_D_rank),
                as.integer(which.max(fit$ranking$D_cv_lower)))),
     "D-only candidate ranking selects maximal crossfit gain")
need(nrow(fit$validation)==5L &&
     identical(as.character(fit$validation$method),
               crit9_modes[-1L]),
     "five predeclared, separately-confirmed controls")
for(i in seq_len(nrow(fit$validation))) {
 mode<-fit$validation$method[[i]]
 ok<-isTRUE(fit$validation$accepted[[i]])
 need(identical(ok,isTRUE(fit$accepted[[mode]])) &&
      (ok==!is.null(fit$coefficients[[mode]])),
      paste("acceptance/fit consistency",mode))
 if(!ok)
   need(max(abs(pred$probs[[mode]]-
        pred$probs$ORIGINAL_RLR))<1e-12,
        paste("exact RLR fallback",mode))
}
need(all(vapply(pred$probs,function(v)length(v)==nrow(test$X) &&
                  all(is.finite(v))&&all(v>=0&v<=1),logical(1))),
     "complete valid heldout probabilities")
X2 <- as.data.frame(lapply(test$X,function(v)v*2+13),
                    check.names=FALSE)
names(X2)<-names(test$X)
rownames(X2)<-rownames(test$X)
other<-crit9_predict(fit,X2)
need(all(vapply(crit9_modes,function(mode)
       isTRUE(all.equal(pred$probs[[mode]],other$probs[[mode]],
                        tolerance=1e-12)),logical(1))),
     "affine within-patient rank invariance")
f<-tempfile(fileext=".rds")
saveRDS(fit,f)
replay<-crit9_predict(readRDS(f),test$X)
unlink(f)
need(identical(pred,replay),"bit-exact RDS replay")
if(isTRUE(fit$accepted[["CONDITIONAL_TREE"]])) {
 ds<-crit9_designs(test$X,fit)$CONDITIONAL_TREE
 g<-crit7_relations(test$X,fit$gate)[,1L]
 rf<-crit7_relations(test$X,fit$pair_false)[,1L]
 rt<-crit7_relations(test$X,fit$pair_true)[,1L]
 z<-as.numeric(g>0)
 x<-cbind(1,z,rf,z*rf,z*rt)
 theta<-fit$coefficients$CONDITIONAL_TREE
 global_coef<-c(theta[[1L]],theta[[2L]]-theta[[1L]],
                theta[[3L]],-theta[[3L]],theta[[4L]])
 exact<-crit7_correction_prob(
   pred$probs$ORIGINAL_RLR,x,global_coef,fit$cfg)
 need(max(abs(exact-pred$probs$CONDITIONAL_TREE))<1e-12,
      "tree vs explicit global-GLM algebraic parity")
}
# Mandatory candidate-absence fallback test: cannot nominate
# a split when each requested leaf size exceeds the entire D cohort.
cfg0<-crit9_cfg()
cfg0$min_leaf_n<-100000L
zero<-crit9_fit(train$X,train$y,cfg0)
pzero<-crit9_predict(zero,test$X)
need(zero$nomination_n==0L &&
     !any(zero$accepted[c("GATE_ONLY","ADDITIVE_PAIRS",
       "GLOBAL_INTERACTION","CONDITIONAL_TREE")]) &&
     all(vapply(c("GATE_ONLY","ADDITIVE_PAIRS",
       "GLOBAL_INTERACTION","CONDITIONAL_TREE"),function(mode)
       identical(pzero$probs[[mode]],pzero$probs$ORIGINAL_RLR),
       logical(1))),
     "zero candidates produce exact RLR fallback in all four arms")
# Idempotent D-only selection from the same seed and patient rows.
fit2<-crit9_fit(train$X,train$y)
need(identical(fit$selected_D_rank,fit2$selected_D_rank) &&
     identical(fit$ranking,fit2$ranking) &&
     identical(fit$accepted,fit2$accepted),
     "deterministic D-only crossfit nominee and V decisions")
message("CRIT9_SYNTHETIC_ALL_GATES=PASS")
message("BIOMEDICAL_DATA_USED=FALSE")
