# Test-only forced L2 acceptance. NEVER use generated metrics as research results.
# Rscript --vanilla R/new/30_rlmt_forced_split_smoke.R results/rlmt_gate2_chen/config.json
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Specify existing Chen config.json")
suppressPackageStartupMessages(library(jsonlite))
source("R/new/07_world_protocol_utils.R")
source("R/new/28_rlmt_utils.R")

cfg <- jsonlite::read_json(args[1], simplifyVector = TRUE)
n_worlds <- if (is.data.frame(cfg$worlds)) nrow(cfg$worlds) else length(cfg$worlds)
if (n_worlds != 1L) stop("Expected exactly one CompCancer world")
world <- if (is.data.frame(cfg$worlds)) as.list(cfg$worlds[1, , drop = FALSE]) else cfg$worlds[[1]]
if (
    !identical(as.character(world$datasets_include[[1]]), "chen-2002") ||
    as.integer(world$n_folds) != 5L ||
    as.integer(world$max_splits_per_dataset) != 1L ||
    as.integer(cfg$models$seed) != 42L) stop("Expected single-split Chen seed-42 config")
dataset <- "chen-2002"
base <- file.path(cfg$out_root, "worlds", world$name, "datasets",
                  dataset, "splits", "rep00_fold00")
pool_file <- file.path(base, "pairs_train", "top1000", "top_1000_pairs.csv")
l2 <- file.path(base, "rlmt", "RLMT_K10_m10_d2_pool1000_L2")
audit_file <- file.path(l2, "lookahead.csv")
screen_file <- file.path(l2, "screen_top.csv")
for (p in c(pool_file, audit_file, screen_file))
  if (!file.exists(p)) stop("Missing prior artifact: ", p)
pool <- read.csv(pool_file, stringsAsFactors = FALSE)
audit <- read.csv(audit_file, stringsAsFactors = FALSE)
screen_old <- read.csv(screen_file, stringsAsFactors = FALSE)
if (nrow(pool) < 10L || nrow(audit) != 20L || nrow(screen_old) != 20L)
  stop("Unexpected saved Chen screening artifact sizes")
valid <- which(!is.na(audit$valid) & audit$valid &
               is.finite(audit$mean_gain) & audit$mean_gain > 0)
if (!length(valid)) stop("No previously valid positive-gain candidate")
chosen <- audit[valid[which.max(audit$mean_gain[valid])], , drop = FALSE]
if (chosen$mean_gain - chosen$se_gain > 0)
  stop("Candidate passes real 1-SE; forced smoke would be misleading")
if (!any(screen_old$gene_a == chosen$gene_a &
         screen_old$gene_b == chosen$gene_b))
  stop("Chosen candidate was not originally screened")

xy <- read_xy_file(file.path(world$comp_dir, "chen-2002.arff"))
sp <- make_repeated_stratified_folds(xy$y, n_folds = 5L,
                                     n_repeats = 1L, seed = 42L)[[1]]
Xtr <- xy$X[sp$train_idx, , drop = FALSE]
ytr <- xy$y[sp$train_idx]
Xte <- xy$X[sp$test_idx, , drop = FALSE]
yte <- xy$y[sp$test_idx]
conf <- rlmt_defaults(cfg$models$rlmt)
zt <- Xtr[[chosen$gene_a]] > Xtr[[chosen$gene_b]]
if (!rlmt_children_ok(zt, ytr, conf))
  stop("Selected split violates outer-training child support")

# Separate output, intentionally outside worlds/ used by the result collector.
out <- file.path(cfg$out_root, "functional_smoke_FORCED_NOT_FOR_PAPER")
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("Output already exists; no overwrite: ", out)
ensure_dir(out)
writeLines(c("TEST ONLY: SPLIT FORCED IN MEMORY; NO SCIENTIFIC VALIDITY.",
             "Do not collect, compare, publish or interpret these metrics.",
             paste("Forced relation:", chosen$gene_a, ">", chosen$gene_b),
             paste("Original mean_gain - SE:", chosen$mean_gain - chosen$se_gain)),
           file.path(out, "NOTICE_TEST_ONLY.txt"))

# Override ONLY this function in this test's R process. Core Git files unchanged.
rlmt_lookahead <- function(X, y, pool, parent_oof, folds, screen, conf) {
  if (!nrow(screen) ||
      !any(screen$gene_a == chosen$gene_a & screen$gene_b == chosen$gene_b))
    stop("Chosen comparison did not survive the reproduced OOF screen")
  if (!rlmt_children_ok(X[[chosen$gene_a]] > X[[chosen$gene_b]], y, conf))
    stop("Candidate no longer admissible")
  list(winner = chosen, audit = audit)
}
cat("TEST_ONLY forced split:", chosen$gene_a, ">", chosen$gene_b, "\n")
rlmt_run(Xtr, ytr, Xte, yte, pool, max_leaves = 2L,
         config = cfg$models$rlmt, out_dir = out,
         extra_metrics = list(world_name = "FUNCTIONAL_SMOKE_ONLY",
                              dataset_id = dataset, split_id = "rep00_fold00"))

pred <- read.csv(file.path(out, "predictions_test.csv"), stringsAsFactors = FALSE)
routes <- read.csv(file.path(out, "test_routes.csv"), stringsAsFactors = FALSE)
nodes <- read.csv(file.path(out, "nodes.csv"))
tree_path <- file.path(out, "tree.json")
met_path <- file.path(out, "metrics_test.json")
tree <- jsonlite::read_json(tree_path, simplifyVector = TRUE)
met <- jsonlite::read_json(met_path, simplifyVector = TRUE)
check <- function(cond, msg) if (!isTRUE(cond)) stop("SMOKE FAIL: ", msg)
check(isTRUE(tree$accepted_split) && tree$size == 2L &&
      met$realized_leaves == 2L && nrow(nodes) == 3L &&
      identical(as.integer(nodes$node_id), c(1L,2L,3L)) &&
      identical(as.logical(nodes$is_leaf), c(FALSE,TRUE,TRUE)),
      "two-leaf tree structure")
check(nrow(pred) == nrow(Xte) && nrow(routes) == nrow(Xte) &&
      !anyDuplicated(pred$sample_id) && !anyDuplicated(routes$sample_id) &&
      identical(pred$sample_id, rownames(Xte)) &&
      identical(routes$sample_id, pred$sample_id) &&
      identical(pred$y_true, as.character(yte)),
      "outer-test IDs, labels, exactly one route per sample")
expected <- ifelse(Xte[[chosen$gene_a]] > Xte[[chosen$gene_b]], 2L, 3L)
check(identical(as.integer(routes$leaf_id), as.integer(expected)),
      "binary routing")
check(!anyNA(pred$prob) && all(is.finite(pred$prob)) &&
      all(pred$prob >= 0 & pred$prob <= 1) &&
      identical(pred$y_pred, ifelse(pred$prob >= 0.5,
                levels(yte)[2], levels(yte)[1])), "probabilities and class labels")
check(identical(as.integer(nodes$n[2:3]), as.integer(c(sum(zt),sum(!zt)))),
      "training node membership and sizes")
check(identical(as.character(tree$split$gene_a), as.character(chosen$gene_a)) &&
      identical(as.character(tree$split$gene_b), as.character(chosen$gene_b)),
      "serialized comparison")

for (id in c(2L, 3L)) {
  fd <- file.path(out, "leaf_models", paste0("node_", id))
  for (fn in c("cvfit.rds", "coefficients.csv", "active_pairs.csv"))
    check(file.exists(file.path(fd, fn)), paste("persisted leaf", id, fn))
  co <- read.csv(file.path(fd, "coefficients.csv"), stringsAsFactors = FALSE)
  active <- read.csv(file.path(fd, "active_pairs.csv"), stringsAsFactors = FALSE)
  fit <- readRDS(file.path(fd, "cvfit.rds"))
  check(nrow(co) > 0L && nrow(active) == sum(co$beta != 0) &&
        nodes$leaf_globally_active_relation_count[id] == nrow(active),
        paste("coefficient / active pair counts in leaf", id))
  pieces <- strsplit(co$ordered_pair, " > ", fixed = TRUE)
  check(all(lengths(pieces) == 2L), "stored ordered pair names")
  pairs <- data.frame(
    gene_left = vapply(pieces, function(x) x[[1L]], character(1)),
    gene_right = vapply(pieces, function(x) x[[2L]], character(1)),
    ordered_pair = co$ordered_pair, stringsAsFactors = FALSE)
  test_rows <- which(routes$leaf_id == id)
  if (length(test_rows)) {
    design <- encode_pair_matrix(Xte[test_rows, , drop = FALSE],
                                 pairs, m = conf$m)
    reconstructed <- as.numeric(predict(
      fit, newx = design, s = conf$lambda_choice, type = "response"))
    check(isTRUE(all.equal(reconstructed, pred$prob[test_rows],
                           tolerance = 1e-12)),
          paste("restored fitted leaf predictions", id))
  }
}
y01 <- rlmt_y01(yte)
check(abs(rlmt_logloss(y01, pred$prob) - met$logloss) <= 5.1e-5,
      "recalculated outer-test logloss")
check(abs(mean((pred$prob - y01)^2) - met$brier) <= 5.1e-5,
      "recalculated outer-test Brier")

# Mark both machine-readable files as scientifically invalid test fixtures.
tree$forced_split_test_only <- TRUE
tree$scientific_valid <- FALSE
tree$original_lower_gain <- chosen$mean_gain - chosen$se_gain
met$forced_split_test_only <- TRUE
met$scientific_valid <- FALSE
met$method <- "RLMT_FORCED_SMOKE_NOT_SCIENTIFIC"
jsonlite::write_json(tree, tree_path, auto_unbox = TRUE, pretty = TRUE)
jsonlite::write_json(met, met_path, auto_unbox = TRUE, pretty = TRUE)
cat("FORCED_SPLIT_FUNCTIONAL_SMOKE=PASS\n")
cat("TRAIN_CHILD_SIZES=", sum(zt), ",", sum(!zt), "\n", sep = "")
cat("OUTER_TEST_LEAF_COUNTS=", sum(routes$leaf_id == 2L), ",",
    sum(routes$leaf_id == 3L), "\n", sep = "")
cat("SCIENTIFIC_RESULT=NO\n")
