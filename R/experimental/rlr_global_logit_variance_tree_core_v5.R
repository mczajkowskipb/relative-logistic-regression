# RLR-GLVT v5: globally optimal bounded rank tree for a FIXED candidate
# pool, honest training-only RLR OOF logit conditional-variance objective.
# NOT a publication or novelty claim. No original RLR source is modified.
# All model leaf predictions use TRUE training labels, never the teacher.
rrv5_need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call.=FALSE)

rrv5_ranks <- function(X) {
  m <- as.matrix(X)
  rrv5_need(is.numeric(m) && nrow(m)>0L && ncol(m)>1L &&
             all(is.finite(m)), "numeric finite expression matrix required")
  z <- t(apply(m, 1L, rank, ties.method="average"))
  dimnames(z) <- dimnames(m)
  z
}

rrv5_matrix <- function(X, pairs) {
  rrv5_need(all(c("gene_left","gene_right") %in% names(pairs)) &&
             nrow(pairs)>0L, "missing ranked-pair candidate panel")
  r <- rrv5_ranks(X)
  rrv5_need(all(c(pairs$gene_left,pairs$gene_right) %in% colnames(r)),
             "unavailable gene/probe (fixed feature universe required)")
  z <- vapply(seq_len(nrow(pairs)), function(j) {
    as.integer(r[,pairs$gene_right[j]]-r[,pairs$gene_left[j]] > 0)
  }, integer(nrow(r)))
  if (is.null(dim(z))) z <- matrix(z,ncol=1L)
  storage.mode(z) <- "integer"
  z
}

rrv5_fit <- function(D, y, pairs, p_oof=NULL, weight=0.25,
                     max_depth=2L, max_internal=3L,
                     min_leaf=4L, node_penalty=0.005,
                     path_penalty=0.005) {
  rrv5_need(is.factor(y) && nlevels(y)==2L && length(y)==nrow(D),
             "binary training labels required")
  rrv5_need(weight>=0 && weight<=1 && is.finite(weight) &&
             max_depth==2L && max_internal==3L &&
             min_leaf>=2L && node_penalty>=0 && path_penalty>=0,
             "v5 finite depth/budget/weights contract")
  if (weight>0) rrv5_need(length(p_oof)==nrow(D) &&
      all(is.finite(p_oof)) && all(p_oof>0 & p_oof<1),
      "genuinely heldout training-only OOF RLR probabilities required")
  yy <- as.integer(y==levels(y)[2L])
  a <- rrv5_matrix(D,pairs)
  n <- length(yy)
  parent_p <- (sum(yy)+1)/(n+2)
  leaf_prob <- function(idx) (sum(yy[idx])+4*parent_p)/(length(idx)+4)
  ce <- function(idx) {
    p <- leaf_prob(idx)
    -sum(yy[idx]*log(p)+(1-yy[idx])*log1p(-p))
  }
  all_idx <- seq_len(n)
  ref_ce <- max(ce(all_idx), 1e-10)
  # OOF logit only enters the *structure* objective; neither leaf
  # probability nor autonomous inference depends on RLR.
  if (weight>0) {
    ell <- qlogis(pmin(1-1e-4,pmax(1e-4,as.numeric(p_oof))))
    ref_sse <- sum((ell-mean(ell))^2)
    if (ref_sse < 1e-10) {
      weight <- 0
      ref_sse <- 1
    }
  } else {
    ell <- rep(0,n)
    ref_sse <- 1
  }
  loss <- function(idx) {
    sse <- if (weight==0) 0 else sum((ell[idx]-mean(ell[idx]))^2)
    (1-weight)*ce(idx)/ref_ce + weight*sse/ref_sse
  }
  # At depth 2 every child of the root can have at most one more split.
  # The best child split is independent of its sibling. Consequently,
  # enumerating all possible roots and both best children IS exhaustive
  # global optimization over 0..3 internal nodes in this restricted class.
  best_child <- function(idx, exclude) {
    ans <- list(feature=NA_integer_, gain=0)
    if (length(idx)<2L*min_leaf) return(ans)
    old <- loss(idx)
    for (j in seq_len(ncol(a))) {
      if (j==exclude) next
      left <- idx[a[idx,j]==0L]
      right <- idx[a[idx,j]==1L]
      if (min(length(left),length(right))<min_leaf) next
      gain <- old-loss(left)-loss(right)-node_penalty-
              path_penalty*length(idx)/n
      if (gain > ans$gain + 1e-12)
        ans <- list(feature=j,gain=gain,left=left,right=right)
    }
    ans
  }
  winner <- list(feature=NA_integer_, score=loss(all_idx))
  for (j in seq_len(ncol(a))) {
    left <- all_idx[a[,j]==0L]
    right <- all_idx[a[,j]==1L]
    if (min(length(left),length(right))<min_leaf) next
    cl <- best_child(left,j)
    cr <- best_child(right,j)
    score <- loss(left)+loss(right)+node_penalty+path_penalty-
             cl$gain-cr$gain
    if (score < winner$score - 1e-12) {
      winner <- list(feature=j,score=score,left=left,right=right,
                     left_split=cl,right_split=cr)
    }
  }
  nodes <- list()
  insert_node <- function(idx,depth,plan=NULL) {
    id <- length(nodes)+1L
    feature <- if (is.null(plan)) NA_integer_ else plan$feature
    nodes[[id]] <<- list(id=id,n=length(idx),depth=depth,
                         prob=leaf_prob(idx),feature=feature,
                         children=integer())
    if (!is.na(feature)) {
      subleft <- if(depth==0L) plan$left_split else NULL
      subright <- if(depth==0L) plan$right_split else NULL
      id0 <- insert_node(plan$left,depth+1L,subleft)
      id1 <- insert_node(plan$right,depth+1L,subright)
      nodes[[id]]$children <<- c(id0,id1)
    }
    id
  }
  insert_node(all_idx,0L,winner)
  used <- sort(unique(vapply(nodes,function(nd)nd$feature,integer(1))),
               na.last=NA)
  compact_pairs <- pairs[used,c("gene_left","gene_right"),drop=FALSE]
  rownames(compact_pairs) <- NULL
  for (i in seq_along(nodes))
    if (!is.na(nodes[[i]]$feature))
      nodes[[i]]$feature <- match(nodes[[i]]$feature,used)
  nn <- sum(vapply(nodes,function(nd)!is.na(nd$feature),logical(1)))
  rrv5_need(nn<=3L && length(nodes)==2L*nn+1L,
             "invalid exact-search tree construction")
  structure(list(nodes=nodes,pairs=compact_pairs,
    n_internal=nn,n_leaves=nn+1L,n_unique_pairs=length(used),
    n_unique_genes=length(unique(c(compact_pairs$gene_left,
                                   compact_pairs$gene_right))),
    objective=winner$score,teacher_weight_effective=weight,
    max_depth=2L,max_evaluated_pairs=2L,
    n_training=n,training_mode="D_ONLY_OOF_RLR_LOGIT_VARIANCE"),
    class="rlr_global_logit_variance_tree")
}

rrv5_predict <- function(model,E) {
  ranks <- rrv5_ranks(E)
  rrv5_need(all(c(model$pairs$gene_left,model$pairs$gene_right) %in%
                 colnames(ranks)), "missing trained probes at inference")
  n <- nrow(ranks)
  p <- numeric(n);checks <- integer(n)
  ends <- integer(n);paths <- character(n)
  for (i in seq_len(n)) {
    node_id <- 1L
    route <- character()
    repeat {
      nd <- model$nodes[[node_id]]
      if (is.na(nd$feature)) {
        p[i] <- nd$prob
        ends[i] <- node_id
        break
      }
      f <- nd$feature
      gap <- ranks[i,model$pairs$gene_right[f]]-
             ranks[i,model$pairs$gene_left[f]]
      b <- as.integer(gap>0)+1L
      checks[i] <- checks[i]+1L
      route <- c(route,paste0(model$pairs$gene_left[f],">",
                   model$pairs$gene_right[f],":",b))
      node_id <- nd$children[b]
    }
    paths[i] <- paste(route,collapse=" / ")
  }
  rrv5_need(all(is.finite(p)) && all(p>0 & p<1) &&
              all(checks<=2L), "invalid autonomous tree inference")
  list(prob=p,label=as.integer(p>=.5),n_evaluated=checks,
       leaf=ends,path=paths)
}
