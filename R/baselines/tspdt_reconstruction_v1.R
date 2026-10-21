# TSPDT 2011-inspired, independently reconstructed greedy node-local TSP tree.
# NOT the original authors' code and NOT BigTSP::tsp.tree.
# This self-contained module does not consume or modify RLMT results.
# Deliberate departures: training-only MAD f100 (paper: ReliefF f1000),
# min_leaf=5, max_depth=3, deterministic gene-order ties,
# optional approximate Wilson pessimistic pruning, Laplace leaf score.

tspdt_recon_check_X <- function(X,training_names=NULL) {
  if (!is.matrix(X) && !is.data.frame(X))
    stop("TSPDT_RECON: X must be a matrix/data.frame")
  X <- as.matrix(X)
  if (!is.numeric(X) || nrow(X)<1L || ncol(X)<2L ||
      is.null(colnames(X)) || anyNA(colnames(X)) ||
      anyDuplicated(colnames(X)) || any(!nzchar(colnames(X))) ||
      anyNA(X) || any(!is.finite(X)))
    stop("TSPDT_RECON: numeric finite matrix with unique gene names required")
  if (!is.null(training_names) && !identical(colnames(X),training_names))
    stop("TSPDT_RECON: predictor schema differs from training")
  storage.mode(X) <- "double"
  X
}

tspdt_recon_select <- function(X,top_n=100L) {
  if (length(top_n)!=1L || is.na(top_n) ||
      !is.finite(top_n) || top_n!=as.integer(top_n) ||
      top_n<2L || top_n>1000L)
    stop("TSPDT_RECON: invalid feature cap")
  mad <- apply(X,2L,stats::mad,constant=1)
  indices <- order(-mad,colnames(X),method="radix")
  indices <- indices[seq_len(min(as.integer(top_n),ncol(X)))]
  list(names=colnames(X)[indices],
       mad=unname(mad[indices]),indices=indices)
}

tspdt_recon_best_pair <- function(X,y,min_leaf=5L) {
  if (!is.matrix(X) || !is.factor(y) ||
      length(levels(y))!=2L || nrow(X)!=length(y) ||
      anyNA(y))
    stop("TSPDT_RECON: invalid node search input")
  n <- nrow(X)
  p <- ncol(X)
  class0 <- which(y==levels(y)[1L])
  class1 <- which(y==levels(y)[2L])
  if (!length(class0) || !length(class1) ||
      p<2L || n<2L*min_leaf) return(NULL)
  best <- 0
  winner <- NULL
  # O(n*p^2) work with O(n*p) temporary storage, NOT a dense all-pairs X.
  for (a in seq_len(p-1L)) {
    remaining <- (a+1L):p
    comparisons <- sweep(X[,remaining,drop=FALSE],1L,
                         X[,a],FUN=">")
    n_left <- colSums(comparisons)
    feasible <- which(n_left>=min_leaf & n-n_left>=min_leaf)
    if (!length(feasible)) next
    q0 <- colMeans(comparisons[class0,feasible,drop=FALSE])
    q1 <- colMeans(comparisons[class1,feasible,drop=FALSE])
    scores <- abs(q0-q1)
    top <- which.max(scores)
    score <- as.numeric(scores[[top]])
    # First encountered feature index wins scores tied within 1e-12.
    if (is.finite(score) && score>best+1e-12) {
      best <- score
      b <- remaining[feasible[top]]
      winner <- list(gene_a=colnames(X)[a],gene_b=colnames(X)[b],
                     a=as.integer(a),b=as.integer(b),score=score,
                     p_class0=as.numeric(q0[top]),
                     p_class1=as.numeric(q1[top]),
                     n_left=as.integer(n_left[feasible[top]]),
                     n_right=as.integer(n-n_left[feasible[top]]))
    }
  }
  winner
}

tspdt_recon_wilson_error <- function(errors,n,cf=0.25) {
  if (n<1L || errors<0 || errors>n ||
      cf<=0 || cf>=0.5) stop("TSPDT_RECON: Wilson input invalid")
  z <- stats::qnorm(1-cf)
  r <- errors/n
  upper <- (r+z^2/(2*n)+z*sqrt(r*(1-r)/n+z^2/(4*n^2)))/
           (1+z^2/n)
  as.numeric(n*upper)
}

tspdt_recon_size <- function(node) {
  if (is.null(node$split)) return(c(internal=0L,leaves=1L))
  c(internal=1L,leaves=0L) +
    tspdt_recon_size(node$left)+tspdt_recon_size(node$right)
}

tspdt_recon_fit <- function(X_train,y_train,top_n=100L,
                            min_leaf=5L,max_depth=3L,
                            pruning=c("wilson","none"),
                            confidence=0.25,laplace_alpha=1) {
  X <- tspdt_recon_check_X(X_train)
  if (!is.factor(y_train) || length(levels(y_train))!=2L ||
      length(y_train)!=nrow(X) || anyNA(y_train) ||
      any(table(y_train)==0L))
    stop("TSPDT_RECON: both classes required")
  pruning <- match.arg(pruning)
  if (length(min_leaf)!=1L || is.na(min_leaf) ||
      min_leaf!=as.integer(min_leaf) || min_leaf<1L ||
      length(max_depth)!=1L || is.na(max_depth) ||
      max_depth!=as.integer(max_depth) || max_depth<0L ||
      max_depth>8L || length(laplace_alpha)!=1L ||
      !is.finite(laplace_alpha) || laplace_alpha<0 ||
      length(confidence)!=1L || !is.finite(confidence) ||
      confidence<=0 || confidence>=0.5)
    stop("TSPDT_RECON: invalid tree options")
  selected <- tspdt_recon_select(X,top_n)
  Xt <- X[,selected$names,drop=FALSE]
  y <- factor(y_train,levels=levels(y_train))
  node_id <- 0L
  grow <- function(indices,depth) {
    node_id <<- node_id+1L
    local_y <- y[indices]
    counts <- as.integer(table(local_y))
    n <- length(indices)
    errors <- n-max(counts)
    leaf_err <- tspdt_recon_wilson_error(errors,n,confidence)
    node <- list(
      id=as.integer(node_id),depth=as.integer(depth),
      n=as.integer(n),class0_n=counts[1L],class1_n=counts[2L],
      majority=if(counts[2L]>counts[1L]) levels(y)[2L] else levels(y)[1L],
      prob_class1=(counts[2L]+laplace_alpha)/(n+2*laplace_alpha),
      estimated_leaf_error=leaf_err,
      estimated_subtree_error=leaf_err,
      split=NULL,left=NULL,right=NULL)
    if (depth>=max_depth || any(counts==0L) ||
        n<2L*min_leaf) return(node)
    best <- tspdt_recon_best_pair(Xt[indices,,drop=FALSE],
                                   local_y,min_leaf=min_leaf)
    if (is.null(best)) return(node)
    route <- Xt[indices,best$a] < Xt[indices,best$b]
    i_left <- indices[route]
    i_right <- indices[!route]
    if (length(i_left)<min_leaf || length(i_right)<min_leaf)
      stop("TSPDT_RECON: child support invariant failed")
    left <- grow(i_left,depth+1L)
    right <- grow(i_right,depth+1L)
    child_err <- left$estimated_subtree_error+
                 right$estimated_subtree_error
    if (pruning=="wilson" && leaf_err<=child_err+1e-12)
      return(node)
    node$split <- best
    node$left <- left
    node$right <- right
    node$estimated_subtree_error <- child_err
    node
  }
  tree <- grow(seq_len(nrow(Xt)),0L)
  sz <- tspdt_recon_size(tree)
  structure(list(
    method="TSPDT2011_INSPIRED_RECONSTRUCTION_V1",
    source_claim="not_original_TSPDT_2011",
    classes=levels(y),
    training_genes=colnames(X),selected_genes=selected$names,
    train_feature_mad=selected$mad,
    hyperparameters=list(top_n=top_n,min_leaf=min_leaf,
                         max_depth=max_depth,pruning=pruning,
                         confidence=confidence,
                         laplace_alpha=laplace_alpha),
    tree=tree,n_train=nrow(X),
    internal_nodes=as.integer(sz[["internal"]]),
    leaves=as.integer(sz[["leaves"]])),
    class="tspdt_recon_v1")
}

tspdt_recon_predict <- function(model,X_test) {
  if (!inherits(model,"tspdt_recon_v1"))
    stop("TSPDT_RECON: wrong model class")
  X <- tspdt_recon_check_X(X_test,model$training_genes)
  X <- X[,model$selected_genes,drop=FALSE]
  prob <- rep(NA_real_,nrow(X))
  assigned <- rep(NA_integer_,nrow(X))
  route <- function(node,ix) {
    if (!length(ix)) return(invisible(NULL))
    if (is.null(node$split)) {
      prob[ix] <<- as.numeric(node$prob_class1)
      assigned[ix] <<- as.integer(node$id)
    } else {
      left <- X[ix,node$split$gene_a] < X[ix,node$split$gene_b]
      route(node$left,ix[left])
      route(node$right,ix[!left])
    }
    invisible(NULL)
  }
  route(model$tree,seq_len(nrow(X)))
  if (anyNA(prob) || anyNA(assigned) ||
      any(!is.finite(prob)) || any(prob<0|prob>1))
    stop("TSPDT_RECON: invalid prediction route")
  klass <- factor(ifelse(prob>0.5,model$classes[2L],model$classes[1L]),
                  levels=model$classes)
  list(class=klass,prob_class1=prob,leaf_id=assigned,
       probability_note="Laplace_training_leaf_NOT_original_2011")
}

tspdt_recon_nodes <- function(model) {
  if (!inherits(model,"tspdt_recon_v1"))
    stop("TSPDT_RECON: invalid model")
  rows <- list()
  walk <- function(node) {
    rows[[length(rows)+1L]] <<- data.frame(
      node_id=node$id,depth=node$depth,n=node$n,
      class0_n=node$class0_n,class1_n=node$class1_n,
      split_a=if (is.null(node$split)) "" else node$split$gene_a,
      split_b=if (is.null(node$split)) "" else node$split$gene_b,
      score=if (is.null(node$split)) NA_real_ else node$split$score,
      is_leaf=is.null(node$split),stringsAsFactors=FALSE)
    if (!is.null(node$split)) {
      walk(node$left)
      walk(node$right)
    }
  }
  walk(model$tree)
  do.call(rbind,rows)
}
