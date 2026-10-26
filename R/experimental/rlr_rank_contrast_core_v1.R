# Rank-contrast stump: compare two WITHIN-PATIENT pair rank separations.
# Example test: [rank(B)-rank(A)] > [rank(D)-rank(C)].
# No raw expression, gene ratios, cohort-global cutoffs, or absolute levels.
source("R/experimental/rlr_rank_gate_core_v1.R")

rc_candidates <- function(X,y,prob,pool,max_roots=40L,
                          max_contrasts=240L,max_return=8L){
  pairs<-rg_pairs(X,pool,max_pool=650L,max_extra=700L)
  D<-rg_gaps(rg_ranks(X),pairs)
  Z<-rg_ternary(D)
  yy<-as.integer(y==levels(y)[2L])
  e<-yy-prob
  w<-pmax(1e-5,prob*(1-prob))
  roots<-unique(pmax(1L,as.integer(round(seq(1,nrow(pairs),length.out=max_roots)))))
  combos<-t(utils::combn(roots,2L))
  if(nrow(combos)>max_contrasts){
    set.seed(414L)
    combos<-combos[sort(sample.int(nrow(combos),max_contrasts)),,drop=FALSE]
  }
  out<-list(); k<-0L
  for(ii in seq_len(nrow(combos))) {
    r1<-combos[ii,1L];r2<-combos[ii,2L]
    g<-as.numeric(D[,r1]>D[,r2])
    a<-which(g==0);b<-which(g==1)
    if(min(length(a),length(b))<12L)next
    if(min(table(factor(yy[a],levels=0:1)),
           table(factor(yy[b],levels=0:1)))<3L)next
    # A scored contextual interaction; later compared with SAME pairs
    # without interaction on independent inner V.
    B<-Z*(g-mean(g))
    score<-as.numeric(crossprod(e,B))^2/
           (as.numeric(crossprod(w,B^2))+1)
    score[c(r1,r2)]<--Inf
    if(!any(is.finite(score)))next
    q<-which.max(score)
    k<-k+1L
    out[[k]]<-data.frame(root_left=r1,root_right=r2,
                          partner=q,score=score[q],
                          root_n_left=length(a),root_n_right=length(b))
  }
  rg_need(k>0L,"no balanced rank contrast; fallback")
  choices<-do.call(rbind,out)
  choices<-head(choices[order(-choices$score,
                      choices$root_left,choices$root_right),,drop=FALSE],
                max_return)
  rownames(choices)<-NULL
  list(candidates=choices,pairs=pairs)
}

rc_design<-function(X,pairs,r1,r2,partner){
  rk<-rg_ranks(X)
  d<-rg_gaps(rk,pairs[c(r1,r2,partner),,drop=FALSE])
  g<-as.numeric(d[,1L]>d[,2L])
  z<-rg_ternary(d[,3L])
  cbind(root=g,partner=z,interaction=g*z)
}

rc_select<-function(T,yT){
  split<-rg_inner(yT)
  D<-T[split$D,,drop=FALSE];V<-T[split$V,,drop=FALSE]
  yD<-factor(as.character(yT[split$D]),levels=levels(yT))
  yV<-factor(as.character(yT[split$V]),levels=levels(yT))
  pool<-compute_top_pairs_switchbox(D,yD,top_n=1000L,buffer_mult=5L)
  base<-rg_parent(D,yD,V,pool,K=20L)
  yy<-as.integer(yV==levels(yV)[2L])
  fallback<-list(accepted=FALSE,root_left=NA_integer_,
      root_right=NA_integer_,partner=NA_integer_,pair_table=NULL,
      delta_vs_rlr=NA_real_,delta_vs_global=NA_real_,
      candidate_count=0L,reason="no_balanced_pair_contrast")
  candidates<-tryCatch(rc_candidates(D,yD,base$train,pool),
       error=function(e){message(conditionMessage(e));NULL})
  if(is.null(candidates))return(fallback)
  best<-NULL
  for(i in seq_len(nrow(candidates$candidates))){
    a<-candidates$candidates[i,,drop=FALSE]
    xd<-rc_design(D,candidates$pairs,
                   a$root_left,a$root_right,a$partner)
    xv<-rc_design(V,candidates$pairs,
                   a$root_left,a$root_right,a$partner)
    pt<-tryCatch(rg_correct(xd,yD,base$train,xv,base$test,TRUE),
                 error=function(e)NULL)
    pg<-tryCatch(rg_correct(xd,yD,base$train,xv,base$test,FALSE),
                 error=function(e)NULL)
    if(is.null(pt)||is.null(pg))next
    gain<-rg_ll(yy,base$test)-rg_ll(yy,pt)
    global_gain<-rg_ll(yy,pg)-rg_ll(yy,pt)
    ba_change<-rg_ba(yy,pt)-rg_ba(yy,base$test)
    can<-list(accepted=gain>=.005 && global_gain>=.002 && ba_change>=-.02,
      root_left=a$root_left,root_right=a$root_right,
      partner=a$partner,pair_table=candidates$pairs,
      delta_vs_rlr=gain,delta_vs_global=global_gain,
      candidate_count=nrow(candidates$candidates),
      reason="V_GATE",ba_delta=ba_change)
    if(is.null(best)||gain>best$delta_vs_rlr)best<-can
  }
  if(is.null(best))return(fallback)
  best
}

rc_evaluate<-function(T,yT,E,yE){
  chosen<-rc_select(T,yT)
  poolT<-compute_top_pairs_switchbox(T,yT,top_n=1000L,buffer_mult=5L)
  base<-rg_parent(T,yT,E,poolT,K=20L)
  global<-base$test;raw<-base$test
  if(!is.null(chosen$pair_table)){
    xd<-rc_design(T,chosen$pair_table,chosen$root_left,
                  chosen$root_right,chosen$partner)
    xe<-rc_design(E,chosen$pair_table,chosen$root_left,
                  chosen$root_right,chosen$partner)
    global<-rg_correct(xd,yT,base$train,xe,base$test,FALSE)
    raw<-rg_correct(xd,yT,base$train,xe,base$test,TRUE)
  }
  gated<-if(isTRUE(chosen$accepted))raw else base$test
  list(probs=list(ORIGINAL_RLR=base$test,
                  GLOBAL_MAIN_CONTROL=global,
                  RANK_STUMP_UNGATED=raw,
                  RANK_STUMP_V_GATED=gated),
       selection=chosen,
       selected_original=base$selected,
       truth=as.integer(yE==levels(yE)[2L]))
}
