# CRIT-v9 predeclared fresh synthetic 4 mechanism x 6 seeds train/test.
# This is a DEVELOPMENT power/null feasibility panel, NOT calibrated FWER.
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")
source("R/new/40_rlmt_dc_utils.R")
source("R/new/42_rlmt_dc_run.R")
source("R/experimental/rlmt_v4_rank_repair.R")
source("R/experimental/crit_v7_original_rlr_two_leaf.R")
source("R/experimental/crit_v9_d_crossfit_one_shot.R")
OUT<-"results/crit_v9_d_crossfit_2x3_v1"
need<-function(ok,msg)if(!isTRUE(ok))stop("CRIT9_PANEL_STOP: ",msg,call.=FALSE)
need(dir.exists(OUT) && file.exists(file.path(OUT,"CODE_SHA.txt")) &&
     !file.exists(file.path(OUT,"synthetic_panel.csv")),
     "panel output missing or previously saved")
scenarios<-c("null","global_pair","gate_only","xor_interaction")
seed_data<-function(n,seed,prefix,scenario) {
 set.seed(as.integer(seed))
 p<-32L
 A<-matrix(stats::rnorm(n*p),nrow=n,ncol=p)
 colnames(A)<-sprintf("gene%02d",seq_len(p))
 gate<-ifelse(A[,1L]>A[,2L],1,-1)
 pair<-ifelse(A[,3L]>A[,4L],1,-1)
 eta<-switch(scenario,
     null=rep(0,n),
     global_pair=1.6*pair,
     gate_only=1.6*gate,
     xor_interaction=2.4*gate*pair,
     stop("CRIT9_PANEL_STOP: unknown scenario"))
 yy<-stats::runif(n)<stats::plogis(eta)
 y<-factor(ifelse(yy,"P","N"),levels=c("N","P"))
 X<-as.data.frame(A,check.names=FALSE)
 rownames(X)<-paste0(prefix,seq_len(n))
 list(X=X,y=y)
}
rows<-list()
j<-0L
for(s in seq_along(scenarios)){
 scenario<-scenarios[[s]]
 for(rep in seq_len(6L)){
   seed<-as.integer(20261100L+s*1000L+rep*23L)
   train<-seed_data(180L,seed,
                  sprintf("%s_r%02d_t_",scenario,rep),scenario)
   test<-seed_data(120L,seed+9L,
                  sprintf("%s_r%02d_e_",scenario,rep),scenario)
   need(min(table(train$y))>=35L &&
        min(table(test$y))>=12L,
        "synthetic train/test class support insufficient")
   model<-crit9_fit(train$X,train$y)
   pred<-crit9_predict(model,test$X)
   need(model$nomination_n<=6L &&
        identical(names(pred$probs),crit9_modes),
        "v9 nomination or method inventory invalid")
   y<-as.numeric(test$y==levels(test$y)[[2L]])
   logloss<-function(p){
      p<-pmax(1e-15,pmin(1-1e-15,p))
      mean(-(y*log(p)+(1-y)*log1p(-p)))
   }
   base<-logloss(pred$probs$ORIGINAL_RLR)
   for(method in crit9_modes) {
     p<-pred$probs[[method]]
     accepted<-if(method=="ORIGINAL_RLR")FALSE else
               isTRUE(model$accepted[[method]])
     if(!accepted)
       need(max(abs(p-pred$probs$ORIGINAL_RLR))<1e-12,
            "synthetic rejected correction must match original")
     j<-j+1L
     rows[[j]]<-data.frame(
       scenario=scenario,replicate=rep,
       seed=seed,n_train=nrow(train$X),n_test=nrow(test$X),
       method=method,
       accepted=as.integer(accepted),
       n_D_candidates=model$nomination_n,
       selected_D_rank=if(is.null(model$selected_D_rank))
             NA_integer_ else model$selected_D_rank,
       original_logloss=base,logloss=logloss(p),
       delta_vs_original=logloss(p)-base,
       brier=mean((p-y)^2),
       scientific_stage="NEW_INDEPENDENT_SYNTHETIC_DEV_ONLY",
       stringsAsFactors=FALSE)
   }
   message("CRIT9_PANEL_COMPLETED ",scenario,":",rep,
      " nominated=",model$nomination_n,
      " D_choice=",if(is.null(model$selected_D_rank))
                       "none" else model$selected_D_rank,
      " tree_accept=",
       as.integer(isTRUE(model$accepted[["CONDITIONAL_TREE"]])),
      " tree_deltaLL=",
       sprintf("%+.7f",logloss(pred$probs$CONDITIONAL_TREE)-base))
 }
}
all<-do.call(rbind,rows)
need(nrow(all)==4L*6L*length(crit9_modes),"panel 24x6 incomplete")
write.csv(all,file.path(OUT,"synthetic_panel.csv"),row.names=FALSE)
message("CRIT9_24_FRESH_SYNTHETIC_DATASETS_SIX_ARMS=PASS")
message("CRIT9_SYNTHETIC_PANEL_DEVELOPMENT_ONLY")
