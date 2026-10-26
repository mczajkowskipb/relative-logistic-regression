# CRIT Stage H1: isolated D-only pair-pool retrieval and interaction ranking.
# This pilot DOES NOT fit original RLR, trees or V/E confirmation.
# It uses the known Bernoulli-null 0.5 residual, only for the synthetic D screen.
# The goal is retrieval-versus-budget scaling, not biomarker/model performance.
source("R/experimental/crit_stageH1_pool_guard_fix1.R")
OUT <- "results/crit_stageH1_bounded_screening_fix1_v1"
need <- function(ok,msg)if(!isTRUE(ok))
    stop("STAGEH1_STOP: ",msg,call.=FALSE)
need(dir.exists(OUT)&&file.exists(file.path(OUT,"CODE_SHA.txt")) &&
     !file.exists(file.path(OUT,"screening_summary.csv")),
     "missing result preflight or would overwrite")
p_values <- c(32L,64L,128L,256L,512L)
scenarios <- c("null","xor")
repetitions <- 8L
nD <- 120L
budget <- 320L
top_edges <- 128L
top_genes <- 24L
top_joint <- 6L

makeD <- function(p,seed,scenario) {
  set.seed(as.integer(seed))
  x<-matrix(rnorm(nD*p),nD,p)
  colnames(x)<-sprintf("gene%03d",seq_len(p))
  z<-ifelse(x[,1L]>x[,2L],1,-1)
  r<-ifelse(x[,3L]>x[,4L],1,-1)
  eta<-if(scenario=="null")rep(0,nD) else 2.4*z*r
  y<-as.integer(runif(nD)<plogis(eta))
  list(X=x,y=y,residual=y-.5,seed=seed)
}
allpairs <- function(p) {
  mat<-utils::combn(seq_len(p),2L)
  matrix(as.integer(mat),nrow=2L)
}
pair_key <- function(p)paste(p[1L,],p[2L,],sep=":")
contains <- function(pool,a,b)
  any(pool[1L,]==a & pool[2L,]==b)
random_pool <- function(pairs,target,seed,exclude=character()) {
  if(length(exclude)) {
    keep<-which(!(pair_key(pairs) %in% exclude))
  } else keep<-seq_len(ncol(pairs))
  set.seed(as.integer(seed))
  chosen<-sample(keep,min(target,length(keep)),replace=FALSE)
  pairs[,chosen,drop=FALSE]
}
aware_pool <- function(X,resid,universe,seed) {
  p<-ncol(X)
  # Feature ranks are computed within patient; no V or test E.
  q<-t(apply(X,1L,rank,ties.method="average"))
  centered<-sweep(q,2L,colMeans(q),"-")
  scale<-sqrt(colSums(centered^2)+1e-8)
  q<-sweep(centered,2L,scale,"/")
  weighted<-abs(crossprod(q*as.numeric(resid),q))
  diag(weighted)<-0
  ij<-which(upper.tri(weighted),arr.ind=TRUE)
  weights<-weighted[ij]
  ord<-order(-weights,ij[,1L],ij[,2L])
  take<-ord[seq_len(min(top_edges,length(ord)))]
  degree<-numeric(p)
  for(j in take) {
    a<-ij[j,1L]; b<-ij[j,2L]; w<-weights[[j]]
    degree[[a]]<-degree[[a]]+w
    degree[[b]]<-degree[[b]]+w
  }
  nodes<-sort(order(-degree,seq_len(p))[seq_len(min(top_genes,p))])
  induced<-utils::combn(nodes,2L)
  need(ncol(induced)<=budget,
       "induced full node pair set exceeds candidate budget")
  extra<-random_pool(universe,budget-ncol(induced),seed,
                     exclude=pair_key(induced))
  pool<-cbind(induced,extra)
  need(ncol(pool)==budget && !anyDuplicated(pair_key(pool)),
       "aware pool duplicate or budget mismatch")
  pool
}
interaction_rank <- function(X,residual,pool) {
  n<-nrow(X)
  k<-ncol(pool)
  # Pair signs R_i in {-1,+1}; no marginal association requirement.
  R<-vapply(seq_len(k),function(j)
       ifelse(X[,pool[1L,j]]>X[,pool[2L,j]],1,-1),
       numeric(n))
  if(is.null(dim(R))) R<-matrix(R,nrow=n)
  C<-abs(crossprod(R,R*as.numeric(residual)))
  # Undirected pair-by-pair screen. Only gene-disjoint pairs permitted.
  idx<-which(upper.tri(C),arr.ind=TRUE)
  # Vectorized exact equivalent of gene-disjoint pair test.
  a<-idx[,1L]; b<-idx[,2L]
  disjoint<-pool[1L,a]!=pool[1L,b] &
            pool[1L,a]!=pool[2L,b] &
            pool[2L,a]!=pool[1L,b] &
            pool[2L,a]!=pool[2L,b]
  idx<-idx[disjoint,,drop=FALSE]
  scores<-C[idx]
  ord<-order(-scores,idx[,1L],idx[,2L])
  idx<-idx[ord[seq_len(min(top_joint,length(ord)))],,drop=FALSE]
  list(idx=idx,scores=scores[ord[seq_len(nrow(idx))]])
}
summary_rows<-list()
pool_rows<-list()
screen_rows<-list()
si<-0L;pi<-0L;ji<-0L
for(sc in seq_along(scenarios))for(p in p_values)for(rep in seq_len(repetitions)) {
  scenario<-scenarios[[sc]]
  seed<-as.integer(20261026L+sc*100000L+p*300L+rep*19L)
  dat<-makeD(p,seed,scenario)
  universe<-allpairs(p)
  strategies<-list(
    FIXED_RANDOM=random_pool(universe,min(budget,ncol(universe)),659L),
    SEED_RANDOM=random_pool(universe,min(budget,ncol(universe)),seed+617L)
  )
  # Time the new D-only interaction-aware second-order graph.
  t0<-proc.time()[["elapsed"]]
  strategies$AWARE_GRAPH<-aware_pool(dat$X,dat$residual,universe,seed+7319L)
  graph_seconds<-proc.time()[["elapsed"]]-t0
  if(p==32L)strategies$EXHAUSTIVE32<-universe
  for(method in names(strategies)) {
    pool<-strategies[[method]]
    K<-ncol(pool)
    critH1_need_pool(pool,method,budget=budget)
    start<-proc.time()[["elapsed"]]
    ranked<-interaction_rank(dat$X,dat$residual,pool)
    screening_seconds<-proc.time()[["elapsed"]]-start
    z<-which(pool[1L,]==1L & pool[2L,]==2L)
    r<-which(pool[1L,]==3L & pool[2L,]==4L)
    containsZ<-length(z)==1L
    containsR<-length(r)==1L
    match_joint<-containsZ && containsR &&
         any(vapply(seq_len(nrow(ranked$idx)),function(j)
           setequal(ranked$idx[j,],c(z,r)),logical(1)))
    si<-si+1L
    summary_rows[[si]]<-data.frame(
       scenario=scenario,p=p,replicate=rep,seed=seed,
       method=method,n_D=nD,pool_budget=K,
       reference_gate_in_pool=as.integer(containsZ),
       reference_pair_in_pool=as.integer(containsR),
       both_reference_pairs_in_pool=as.integer(containsZ&&containsR),
       reference_joint_in_top6=as.integer(match_joint),
       aware_graph_seconds=if(method=="AWARE_GRAPH")graph_seconds else NA_real_,
       screening_seconds=screening_seconds,
       scientific_stage="D_ONLY_SYNTHETIC_SCREENING_DEV",
       stringsAsFactors=FALSE)
    for(j in seq_len(K)) {
      pi<-pi+1L
      pool_rows[[pi]]<-data.frame(
        scenario=scenario,p=p,replicate=rep,method=method,
        candidate_index=j,gene_a=pool[1L,j],gene_b=pool[2L,j],
        stringsAsFactors=FALSE)
    }
    for(j in seq_len(nrow(ranked$idx))) {
      a<-ranked$idx[j,1L];b<-ranked$idx[j,2L]
      ji<-ji+1L
      screen_rows[[ji]]<-data.frame(
        scenario=scenario,p=p,replicate=rep,method=method,
        joint_rank=j,
        pairA_index=a,pairB_index=b,
        score=ranked$scores[[j]],
        stringsAsFactors=FALSE)
    }
    message(sprintf(
       "STAGEH1_CASE %s p=%d r=%d %s K=%d reference_pairs=%d/%d joint_top6=%d graph_s=%.3f screen_s=%.3f",
       scenario,p,rep,method,K,
       as.integer(containsZ),as.integer(containsR),as.integer(match_joint),
       if(method=="AWARE_GRAPH")graph_seconds else 0,screening_seconds))
  }
}
s<-do.call(rbind,summary_rows)
pools<-do.call(rbind,pool_rows)
screen<-do.call(rbind,screen_rows)
need(nrow(s)==2L*5L*8L*3L+2L*8L &&
     nrow(screen)==6L*nrow(s) &&
     !anyNA(s$reference_joint_in_top6),
     "method count or top-six result inventory incorrect")
write.csv(s,file.path(OUT,"screening_summary.csv"),row.names=FALSE)
write.csv(pools,file.path(OUT,"candidate_pools.csv"),row.names=FALSE)
write.csv(screen,file.path(OUT,"top6_pair_interactions.csv"),row.names=FALSE)
message("STAGEH1_80_INDEPENDENT_D_ONLY_SYNTHETIC_STUDIES=PASS")
message("STAGEH1_256_SCREENING_RUNS_COMPLETE=PASS")
message("STAGEH1_SOURCE_POOL_AND_TOP6_PROVENANCE=PASS")
message("STAGEH1_SCIENTIFIC_STAGE=D_ONLY_SYNTHETIC_SCREENING_DEVELOPMENT_ONLY")
