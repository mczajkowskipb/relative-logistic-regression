# Usage: Rscript R/new/29_rlmt_gate_check.R configs/rlmt/compcancer_gate.json
# Fails closed: PASS is printed only after checking the written outer-test artifacts.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 29_rlmt_gate_check.R <config.json>")
suppressPackageStartupMessages(library(jsonlite))
source("R/new/07_world_protocol_utils.R")

cfg <- jsonlite::read_json(args[[1]], simplifyVector = TRUE)
root <- cfg$out_root
world <- if (is.data.frame(cfg$worlds)) as.list(cfg$worlds[1, ]) else cfg$worlds[[1]]
dataset <- as.character(world$datasets_include[[1]])
split <- "rep00_fold00"
base <- file.path(root, "worlds", world$name, "datasets", dataset, "splits", split)
r_name <- "RLR_K10_m10_d2_pool1000"
l1_name <- "RLMT_K10_m10_d2_pool1000_L1"
l2_name <- "RLMT_K10_m10_d2_pool1000_L2"
dirs <- list(rlr = file.path(base, "rlr", r_name),
             l1 = file.path(base, "rlmt", l1_name),
             l2 = file.path(base, "rlmt", l2_name))
needed <- c("predictions_test.csv", "metrics_test.json")
for (k in names(dirs)) {
  for (f in needed) if (!file.exists(file.path(dirs[[k]], f)))
    stop("Missing ", k, "/", f)
  if (file.exists(file.path(dirs[[k]], "error.json")))
    stop("Error artifact found in ", k)
}
pred <- lapply(dirs, function(d) read.csv(file.path(d, "predictions_test.csv"),
                                          stringsAsFactors = FALSE))
met <- lapply(dirs, function(d) jsonlite::read_json(
  file.path(d, "metrics_test.json"), simplifyVector = TRUE))
close_num <- function(a, b, tag, tol = 1e-12) {
  if (length(a) != length(b) || anyNA(a) || anyNA(b) ||
      any(!is.finite(a)) || any(!is.finite(b)) ||
      !isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))) {
    stop("PARITY MISMATCH: ", tag)
  }
}
for (field in c("sample_id", "y_true", "y_pred")) {
  if (!identical(pred$rlr[[field]], pred$l1[[field]]))
    stop("GATE 1 parity mismatch in ", field)
}
close_num(pred$rlr$prob, pred$l1$prob, "prob")
close_num(pred$rlr$score, pred$l1$score, "score")
for (field in c("auc", "f1", "logloss", "brier", "ece",
                "global_active_pairs", "unique_genes")) {
  close_num(met$rlr[[field]], met$l1[[field]], field)
}
p1 <- read.csv(file.path(dirs$rlr, "active_pairs.csv"))
p2 <- read.csv(file.path(dirs$l1, "active_pairs.csv"))
if (!identical(as.character(p1$ordered_pair), as.character(p2$ordered_pair)))
  stop("GATE 1 selected active pairs mismatch")
co1 <- read.csv(file.path(dirs$rlr, "coefficients.csv"))
co2 <- read.csv(file.path(dirs$l1, "coefficients.csv"))
if (!identical(co1$ordered_pair, co2$ordered_pair))
  stop("GATE 1 coefficient names mismatch")
close_num(co1$beta, co2$beta, "coefficients")
if (met$l1$method_family != "RLMT" || met$l1$realized_leaves != 1L)
  stop("Incorrect L1 metadata")
cat("GATE_1=PASS (RLR parity: probabilities, classes, metrics, active pairs, coefficients)\n")

for (k in c("l1", "l2")) {
  for (f in c("tree.json", "nodes.csv", "test_routes.csv"))
    if (!file.exists(file.path(dirs[[k]], f))) stop("Missing ", k, "/", f)
}
tree <- jsonlite::read_json(file.path(dirs$l2, "tree.json"),
                            simplifyVector = TRUE)
nodes <- read.csv(file.path(dirs$l2, "nodes.csv"))
routes <- read.csv(file.path(dirs$l2, "test_routes.csv"),
                   stringsAsFactors = FALSE)
p <- pred$l2
if (anyDuplicated(p$sample_id) || anyDuplicated(routes$sample_id) ||
    nrow(routes) != nrow(p) || !identical(routes$sample_id, p$sample_id))
  stop("GATE 2 test samples must have one route, ordered as predictions")
if (!identical(p$y_pred, ifelse(p$prob >= 0.5,
                                levels(factor(p$y_true))[2],
                                levels(factor(p$y_true))[1])))
  stop("GATE 2 class threshold disagreement")
if (anyNA(p$prob) || any(!is.finite(p$prob)) ||
    any(p$prob < 0 | p$prob > 1))
  stop("GATE 2 invalid probabilities")
if (met$l2$method_family != "RLMT" ||
    met$l2$realized_leaves != tree$size ||
    tree$size != nrow(nodes[which(nodes$is_leaf), , drop = FALSE]))
  stop("GATE 2 tree/metric nodes mismatch")

# Reconstruct the ONE outer split from the configured raw CompCancer dataset.
arff <- file.path(world$comp_dir, paste0(dataset, ".arff"))
if (!file.exists(arff)) stop("Missing CompCancer ARFF: ", arff)
xy <- read_xy_file(arff)
spec <- make_repeated_stratified_folds(
  xy$y, n_folds = as.integer(world$n_folds), n_repeats = 1L,
  seed = as.integer(cfg$models$seed))[[1]]
X_test <- xy$X[spec$test_idx, , drop = FALSE]
if (!identical(rownames(X_test), p$sample_id))
  stop("GATE 2 raw outer-test sample ids do not match")

if (isTRUE(tree$accepted_split)) {
  if (tree$size != 2L || nrow(nodes) != 3L ||
      !setequal(routes$leaf_id, c(2L, 3L)) &&
        any(!routes$leaf_id %in% c(2L, 3L)))
    stop("GATE 2 accepted split structure mismatch")
  a <- tree$split$gene_a; b <- tree$split$gene_b
  actual <- ifelse(X_test[[a]] > X_test[[b]], 2L, 3L)
  if (!identical(as.integer(routes$leaf_id), as.integer(actual)))
    stop("GATE 2 relational routing mismatch")
  if (!is.finite(tree$split$mean_gain - tree$split$se_gain) ||
      tree$split$mean_gain - tree$split$se_gain <= 0)
    stop("GATE 2 nonconservative split accepted")
  for (id in c(2L, 3L)) {
    leaf <- file.path(dirs$l2, "leaf_models", paste0("node_", id))
    for (f in c("cvfit.rds", "coefficients.csv", "active_pairs.csv"))
      if (!file.exists(file.path(leaf, f))) stop("Missing leaf artifact: ", leaf, "/", f)
    co <- read.csv(file.path(leaf, "coefficients.csv"))
    if (!nrow(co)) stop("Empty child coefficients in leaf ", id)
  }
  cat("GATE_2_SPLIT=ACCEPTED\n")
} else {
  if (tree$size != 1L || nrow(nodes) != 1L ||
      !all(routes$leaf_id == 1L))
    stop("GATE 2 rejected split must be single root leaf")
  close_num(p$prob, pred$rlr$prob, "L2 stopped at root")
  cat("GATE_2_SPLIT=REJECTED (valid pre-pruning)\n")
}
# Recheck the exported evaluation against untouched outer-test labels.
y01 <- as.numeric(p$y_true == levels(xy$y)[2])
eps <- 1e-15
pc <- pmin(pmax(p$prob, eps), 1 - eps)
ll <- -mean(y01 * log(pc) + (1 - y01) * log(1 - pc))
close_num(ll, met$l2$logloss, "L2 logloss")
close_num(mean((p$prob - y01)^2), met$l2$brier, "L2 brier")
cat("GATE_2=PASS (routes, tree, probabilities, untouched outer-test metrics)\n")
