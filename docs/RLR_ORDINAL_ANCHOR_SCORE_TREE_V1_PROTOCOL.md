# Ordinal Anchor-Score Decision Tree V1 — prespecified protocol

**Science objective:** Extend the already finished, unpublished original RLR with a genuinely induced top-down relational tree while preserving *strict within-patient invariance*. Original RLR core and manuscript lock must remain unchanged.

## New mechanism

For one nominated target gene g and 3 or 5 reference genes A from the same sample, compute

S(g,A,x) = sum[a in A] 1(x_g > x_a), giving an integer score 0..|A|.

Split an internal node via 1[S(g,A,x) >= t], with t chosen from 1..|A| on TRAINING data. This includes majority as a special case but does NOT restrict the node to majority voting.

The test uses only ordinal comparisons within the individual patient. No original expression scale, sample-normalization ranges, expression ratios, or between-patient abundance thresholds enter a test. Freeze exact gene/probe identity and universe; independently increasing patient-wise transforms preserve every node score. Individual gene-specific transforms, missing genes or ties introduced by discretization can change the representation.

Mathematical caveat: fixed t corresponds to comparing target g with a particular anchor-set order statistic, whose *identity* can change between patients. Even when it does, the entire panel might agree with one well-chosen fixed reference gene in most patients; hence **SAME_ROOT_BEST_SINGLE_PAIR** is a mandatory control. This is not assumed to be novel versus TSP/Top-Scoring-Group variants.

## Explicit model construction (no E labels during design)

1. On an outer training partition T, make an inner D/V stratified split using seed 420. The outer E partition is never accessed until after model selection.
2. On D, construct original supervised switchBox top1000, then nominate at most 14 hub genes based on occurrences among top candidate pairs BEFORE gene-disjointness.
3. For each candidate target, construct ONE reference panel of size 3 and ONE of size 5 via the previously frozen anchor-selection heuristic: neighbor TSP pairs plus high-stability normalized within-sample rank candidates, ranked by directional separation and penalties for across-class rank drift, spread, proximity and near-duplicate comparisons. The panel does not use raw-abundance stability.
4. Evaluate each integer t=1..k by smoothed D leaf logloss, minimum 10 training samples and minimum 2 patients from each class in both root branches. Add the fixed 0.0025*(k-1) penalty to prefer short signatures.
5. Induce a genuine second-level candidate separately in each eligible D branch with its own local supervised pair pool and newly selected target, panel and t. Eligibility: at least 22 D patients, 5 of each class, 6 in each new leaf and split gain >=0.004 D logloss. These are proposals only, never forced into the final tree.
6. On inner V, select depth2 only if it exceeds the D-trained stump by >=0.005 logloss. If valid tree model also beats original RLR K20 (m10,d2, alpha .5, lambda.min) by >=.005 V logloss and the best single reference from its own root panel by >=.002, without balanced accuracy decreasing more than .02, accept; otherwise return exact full-T original RLR on E.
7. For heldout E scores, keep the chosen D structure fixed; estimate leaf class probabilities with fixed empirical-Bayes smoothing strength 8 from full outer T. The original RLR parent is freshly fitted by its original code on outer T. No E-based model choice.

## Seven predeclared, mutually patient-aligned arms

- ORIGINAL_RLR: unchanged K20 original code.
- ANCHOR_SCORE_STUMP: original root ordinal score threshold, estimated leaf probabilities.
- ANCHOR_SCORE_TREE_D2: original root plus any D-supported independently induced second-level nodes.
- SAME_ROOT_BEST_SINGLE_PAIR: best one target-versus-reference pair among exactly the selected root panel, fitted on T with the same leaf smoothing.
- SAME_COUNT_RIDGE: globally fitted ridge logistic regression on the chosen node count scores, without tree routing, fixed lambda 0.10; a permanently inert column is present only when glmnet requires >=2 predictor columns.
- V_SELECTED_ANCHOR_SCORE_TREE: stump or depth2 as picked on V regardless of RLR fallback, diagnostic.
- RLR_V_GATED_ANCHOR_SCORE_TREE: the accepted tree if V supports it, otherwise bit-exact original RLR.

All seven arms evaluate precisely matched E patients. The matched no-tree ridge control uses the same count scores as the selected D-proposed structural nodes, but its model degrees of freedom are not exactly the same as tree-leaf smoothing; do not claim a perfectly matched-capacity superiority test.

## Frozen exploratory data and audit plan

- Reused original-paper CompCancer cohort names: Chowdary-2006 (104 patients) and Singh-2002 (102 patients). Public Schliep-source *reconstructions* not yet SHA equivalent to original processed ARFF.
- Fixed 3-fold stratified outer CV × 1 repetition, seed42, 6 E folds.
- Seven methods * 6 folds = 42 metrics and 206 * 7 = 1442 patient-method E probabilities.
- Primary diagnostic: heldout E logloss; also raw Brier, BA, Macro-F1; per-fold accepted trees, proposed child nodes, selected k,t, V delta versus original RLR and single-pair control. Independent Python auditor recomputes every metric, E coverage, class identity and exact original-RLR fallback.
- Mechanistic tests: strictly increasing per-patient transforms of 130 synthetic patient vectors; reference-mutation count bound; score threshold t=3/5 exactly reduces to majority/median; explicit depth2 routing with four synthetic leaves; original data schema enforced. Technical PASS does not establish scientific advantage.

## Stop and interpretation rules

This is the **single** frozen ordinal score mechanism V1, not license to retrospectively retune lambda, candidate counts, margins, thresholds, depth or acceptance using exposed E folds. If no reproducible E benefit, record negative even if isolated fold wins. If branches cannot be learned, record infeasibility instead of claiming a successful depth2 tree. Future confirmation must use genuinely untouched datasets and the actual author's TSPDT CUDA implementation rather than proxy. Original completed RLR manuscript is the fallback.
