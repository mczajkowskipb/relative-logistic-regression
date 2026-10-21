# Non-equivalent comparator: BigTSP 1.0 global pair-feature CART.
# Based on reviewed CRAN/BigTSP R/LDCA.R (2012); NOT TSPDT 2011.
# Reuses the independently selected train-only MAD feature subset.
# Requires archived BigTSP 1.0 and its tree/glmnet/randomForest/gbm deps.
# This adapter does not attempt installation or modify the user R library.

bigtsp_pair_cart_fit <- function(X_train,y_train,top_n=100L,
                                 mincut=2L,minsize=5L,mindev=0.001) {
  if (!exists("tspdt_recon_check_X",mode="function") ||
      !exists("tspdt_recon_select",mode="function"))
    stop("BIGTSP_PAIR_CART: source tspdt_reconstruction_v1.R first")
  X <- tspdt_recon_check_X(X_train)
  if (!is.factor(y_train) || length(levels(y_train))!=2L ||
      any(table(y_train)==0L) || anyNA(y_train) ||
      length(y_train)!=nrow(X))
    stop("BIGTSP_PAIR_CART: invalid binary response")
  if (!requireNamespace("BigTSP",quietly=TRUE) ||
      !requireNamespace("tree",quietly=TRUE))
    stop("BIGTSP_PAIR_CART: archived BigTSP 1.0/tree is not installed")
  if (as.character(utils::packageVersion("BigTSP"))!="1.0")
    stop("BIGTSP_PAIR_CART: requires independently reviewed BigTSP 1.0")
  cap <- tspdt_recon_select(X,top_n)
  # 100 feature cap prevents an uncontrolled p(p-1)/2 pair expansion.
  if (length(cap$names)>100L)
    stop("BIGTSP_PAIR_CART: safety cap exceeded (100 genes)")
  if (length(mincut)!=1L || mincut!=as.integer(mincut) ||
      mincut<1L || length(minsize)!=1L ||
      minsize!=as.integer(minsize) || minsize<2L ||
      length(mindev)!=1L || !is.finite(mindev) ||
      mindev<0 || mindev>1)
    stop("BIGTSP_PAIR_CART: invalid tree.control parameters")
  ctrl <- tree::tree.control(
    nobs=nrow(X),mincut=mincut,minsize=minsize,mindev=mindev)
  fit <- BigTSP::tsp.tree(
    X=X[,cap$names,drop=FALSE],
    response=y_train,control=ctrl,split="deviance")
  n_nodes <- nrow(fit$frame)
  n_leaves <- sum(fit$frame$var=="<leaf>")
  structure(list(
    method="BIGTSP_1_0_PAIR_CART_F100",
    source_claim="not_original_TSPDT_2011",
    tree=fit,training_genes=colnames(X),
    selected_genes=cap$names,train_feature_mad=cap$mad,
    classes=levels(y_train),
    tree_nodes=as.integer(n_nodes),
    internal_nodes=as.integer(n_nodes-n_leaves),
    leaves=as.integer(n_leaves),
    hyperparameters=list(top_n=top_n,mincut=mincut,
                         minsize=minsize,mindev=mindev)),
    class="bigtsp_pair_cart_v1")
}

bigtsp_pair_cart_predict <- function(model,X_test) {
  if (!inherits(model,"bigtsp_pair_cart_v1"))
    stop("BIGTSP_PAIR_CART: wrong model class")
  X <- tspdt_recon_check_X(X_test,model$training_genes)
  X <- X[,model$selected_genes,drop=FALSE]
  # Legacy BigTSP::predict.tsp.tree strips the FIRST newdata column
  # unconditionally. Supply a dummy label, NEVER y_test.
  dummy <- factor(rep(model$classes[1L],nrow(X)),
                  levels=model$classes)
  newdata <- data.frame(response=dummy,X,check.names=FALSE)
  pred <- stats::predict(model$tree,newdata=newdata,type="class")
  pred <- factor(as.character(pred),levels=model$classes)
  if (anyNA(pred) || length(pred)!=nrow(X))
    stop("BIGTSP_PAIR_CART: missing or incorrectly ordered test predictions")
  list(class=pred,
       prob_class1=rep(NA_real_,nrow(X)),
       probability_note="No_calibrated_probabilities_original_BigTSP_hard_class")
}
