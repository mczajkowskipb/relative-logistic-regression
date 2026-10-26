# CRIT Stage H2: oracle-anchor conditional-rank-pair search at growing p.
# A discovery-mechanism diagnostic, NOT a feasible blind/scalable CRIT model.
# All data generated independently; only D labels used for second-pair scores.
# Known anchor is privileged PLANTED TRUTH and is not learned from outcomes.
OUT <- "results/crit_stageH2_anchored_completion_v1"
stopifnot(dir.exists(OUT),file.exists(file.path(OUT,"SOURCE_SHA.txt")),
          !file.exists(file.path(OUT,"study_summary.csv")))
need <- function(ok,msg) if(!isTRUE(ok))
  stop("STAGEH2_STOP: ",msg,call.=FALSE)

ps <- c(32L,64L,128L,256L,512L,1024L)
scenarios <- c("null","xor")
replications <- 8L
nD <- 120L
seed_base <- 20270117L
block_n <- 2048L
top_k <- 6L
reference_K <- 320L

makeD <- function(p,seed,scenario) {
  set.seed(as.integer(seed))
  x <- matrix(stats::rnorm(nD*p),nrow=nD,ncol=p)
  colnames(x) <- sprintf("gene%04d",seq_len(p))
  z <- ifelse(x[,1L]>x[,2L],1,-1)
  r <- ifelse(x[,3L]>x[,4L],1,-1)
  eta <- if(identical(scenario,"xor"))2.4*z*r else rep(0,nD)
  yy <- as.integer(stats::runif(nD)<stats::plogis(eta))
  list(X=x,y=yy)
}
rank_one_anchor <- function(X,y,anchor,target,block_n,top_k) {
  p <- ncol(X)
  pairs <- utils::combn(seq_len(p),2L)
  good <- pairs[1L,]!=anchor[[1L]] &
          pairs[2L,]!=anchor[[1L]] &
          pairs[1L,]!=anchor[[2L]] &
          pairs[2L,]!=anchor[[2L]]
  pairs <- pairs[,good,drop=FALSE]
  K <- ncol(pairs)
  need(K==(p-2L)*(p-3L)/2L,"invalid disjoint candidate count")
  anchor_sign <- ifelse(X[,anchor[[1L]]]>X[,anchor[[2L]]],1,-1)
  weights <- (as.numeric(y)-.5)*anchor_sign
  values <- numeric(K)
  then <- proc.time()[["elapsed"]]
  for(start in seq.int(1L,K,by=block_n)) {
    end <- min(start+block_n-1L,K)
    a <- pairs[1L,start:end]
    b <- pairs[2L,start:end]
    # Every operation is a D-only patient-wise rank comparison.
    signatures <- ifelse(X[,a,drop=FALSE]>X[,b,drop=FALSE],1,-1)
    values[start:end] <- abs(drop(crossprod(weights,signatures)))
  }
  elapsed <- proc.time()[["elapsed"]]-then
  order_ix <- order(-values,seq_along(values))
  chosen <- order_ix[seq_len(top_k)]
  trueix <- which(pairs[1L,]==target[[1L]] &
                  pairs[2L,]==target[[2L]])
  need(length(trueix)==1L,"planted target excluded from oracle scan")
  trank <- which(order_ix==trueix)
  need(length(trank)==1L,"missing target candidate rank")
  list(K=K,pairs=pairs[,chosen,drop=FALSE],
       scores=values[chosen],true_rank=as.integer(trank),
       true_score=as.numeric(values[trueix]),scan_seconds=elapsed)
}

rows <- list()
tops <- list()
si <- 0L
ti <- 0L
for(sidx in seq_along(scenarios))for(p in ps)for(rep in seq_len(replications)) {
  scenario <- scenarios[[sidx]]
  seed <- as.integer(seed_base + sidx*100000L + p*311L + rep*37L)
  D <- makeD(p,seed,scenario)
  arms <- list(
    ORACLE_Z_TO_R=list(anchor=c(1L,2L),target=c(3L,4L)),
    ORACLE_R_TO_Z=list(anchor=c(3L,4L),target=c(1L,2L)),
    SHAM_56_TO_R=list(anchor=c(5L,6L),target=c(3L,4L))
  )
  for(method in names(arms)) {
    arm<-arms[[method]]
    got<-rank_one_anchor(D$X,D$y,arm$anchor,arm$target,
                         block_n,top_k)
    si<-si+1L
    rows[[si]]<-data.frame(
      scenario=scenario,p=p,replicate=rep,seed=seed,n_D=nD,
      method=method,anchor_gene_a=arm$anchor[[1L]],
      anchor_gene_b=arm$anchor[[2L]],target_gene_a=arm$target[[1L]],
      target_gene_b=arm$target[[2L]],candidate_pair_count=got$K,
      true_target_rank=got$true_rank,
      target_in_top1=as.integer(got$true_rank==1L),
      target_in_top6=as.integer(got$true_rank<=top_k),
      true_target_score=got$true_score,
      max_score=got$scores[[1L]],
      screening_seconds=got$scan_seconds,
      random_two_pair_coverage_K320=if(p>=32L){
        M<-p*(p-1)/2
        (reference_K/M)*((reference_K-1)/(M-1))
      } else NA_real_,
      scientific_stage="ORACLE_ANCHOR_D_ONLY_SYNTHETIC_DEVELOPMENT",
      stringsAsFactors=FALSE)
    for(j in seq_len(top_k)) {
      ti<-ti+1L
      tops[[ti]]<-data.frame(
        scenario=scenario,p=p,replicate=rep,method=method,
        joint_rank=j,found_pair_a=got$pairs[1L,j],
        found_pair_b=got$pairs[2L,j],
        score=got$scores[[j]],
        true_target=as.integer(got$pairs[1L,j]==arm$target[[1L]] &&
                               got$pairs[2L,j]==arm$target[[2L]]),
        stringsAsFactors=FALSE)
    }
    message(sprintf(
      "H2_CASE %s p=%d rep=%d %s candidates=%d true_rank=%d hit1=%d hit6=%d scan_s=%.3f",
      scenario,p,rep,method,got$K,got$true_rank,
      as.integer(got$true_rank==1L),
      as.integer(got$true_rank<=top_k),got$scan_seconds))
  }
}
s<-do.call(rbind,rows)
t<-do.call(rbind,tops)
need(nrow(s)==2L*length(ps)*replications*3L &&
     nrow(t)==6L*nrow(s) && !anyNA(s$true_target_rank),
     "incomplete H2 case inventory")
need(all(s$candidate_pair_count==(s$p-2L)*(s$p-3L)/2L),
     "candidate cardinality violation")
write.csv(s,file.path(OUT,"study_summary.csv"),row.names=FALSE)
write.csv(t,file.path(OUT,"top6_candidates.csv"),row.names=FALSE)
message("STAGEH2_96_D_ONLY_SYNTHETIC_DATASETS=PASS")
message("STAGEH2_288_ANCHORED_SCANS=PASS")
message("STAGEH2_TRUE_TARGET_RANK_AND_TOP6_PROVENANCE=PASS")
message("STAGEH2_SCIENTIFIC_STAGE=ORACLE_ANCHOR_DEVELOPMENT_ONLY")
