# Original frozen CompCancer ARFF: RLR20 versus controlled binary pair-model comparators

Status: **local Ubuntu run successfully completed** (owner-provided terminal output dated 2026-10-30); 11/11 source and shard ARFF SHA256 matched original publication manifest, 11×5 outer folds×3 methods, independent all-patient metric audit PASS, exit code 0. **This report is grounded in the owner's supplied local stdout/CSV, not in a GitHub-hosted original-ARFF artifact.** Local result path: \`/opt/2026/RLR/rlr-binary-original-arff-v2/results/rlr_binary_outer5_original_arff_audit_v2/\`.

Original data filenames \`/opt/2026/RLR/data/CompCancer/*.arff\`, sharded original files at \`/opt/2026/RLR/data/CompCancer_truth_9x10\` and \`/opt/2026/RLR/data/CompCancer_truth_2x2\`. Frozen checksum manifest: \`checksums/DATA_SHA256_EXPECTED.txt\`. Owner also independently ran \`sha256sum -c\`: **45/45 OK**. Original R runner reported \`ORIGINAL_ARFF_SHA256_PASS\`, source and shard equality on each cohort, \`RLR_ORIGINAL_ARFF_OUTER5_ALL11_RUN_AND_INDEPENDENT_AUDIT=PASS\`, \`BENCHMARK_EXIT_CODE=0\`.

## Fixed method variants

All models are trained on the same original five outer folds, with the **same training-fold-only** switchBox top1000 pairs restricted to K20 gene-disjoint pairs (d=2). Original RLR20 uses signed ternary within-sample rank margins (m=10) and glmnet elastic net alpha=.5; binary elastic net uses simple within-sample pair ordering (PurIST-inspired); binary MCP uses ordering and ncvreg MCP alpha=.05 (DeCAF-inspired). **Only controlled methodological adaptations, not full original clinical PurIST/DeCAF algorithms.**

## 11-dataset equally weighted means; full 55 heldout folds

| Method | F1 ↑ | AUC ↑ | Logloss ↓ | Brier ↓ | ECE ↓ | Active pairs |
|---|---:|---:|---:|---:|---:|---:|
| **RLR20 original** | **0.876179704348** | **0.961673478529** | **0.237943894013** | **0.062637675982** | **0.092092417592** | 16.818182 |
| Binary ENET20 (PurIST-inspired) | 0.823985216632 | 0.948486849809 | 0.290442282242 | 0.077383726834 | 0.110473770917 | 16.745455 |
| Binary MCP20 (DeCAF-inspired) | 0.813640709539 | 0.947842990405 | 0.291233444479 | 0.085454757775 | 0.134419941655 | 11.527273 |

The original ARFF summary values are **numerically identical** (up to printed floating precision) to public Schliep-reconstruction benchmark [#41064335985](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/41064335985), a useful cross-execution confirmation of the preprocessing equivalence for this pipeline. The source SHA256 gate on original ARFF is independently stronger provenance than the reconstructed pilot.

## Paired exact 2^11 sign-flips (11 dataset means; eight exploratory comparisons)

| Comparator | Metric | Mean RLR minus comparator | Raw p | Exploratory BH q |
|---|---|---:|---:|---:|
| Binary ENET | logloss | -0.052498388229 | 0.0341796875 | **0.0390625** |
| Binary ENET | brier | -0.014746050852 | 0.0087890625 | **0.0234375** |
| Binary ENET | auc | +0.013186628720 | 0.0781250000 | 0.0781250 |
| Binary ENET | positive_f1 | +0.052194487716 | 0.0312500000 | **0.0390625** |
| Binary MCP | logloss | -0.053289550466 | 0.0263671875 | **0.0390625** |
| Binary MCP | brier | -0.022817081793 | 0.0019531250 | **0.0156250** |
| Binary MCP | auc | +0.013830488125 | 0.0312500000 | **0.0390625** |
| Binary MCP | positive_f1 | +0.062538994809 | 0.0078125000 | **0.0234375** |

Seven of eight exploratory BH-adjusted differences were below 0.05, except ENET AUC. The statistical unit was 11 datasets, **not 55 independent folds**. For ENET, relative RLR logloss decrease was ~18.1%, Brier decrease ~19.1%.

## Scientific caveats and method integrity

1. Full original publication used an additional repeated outer-CV schedule (470 outer folds/model), **not reproduced here**; this comparison uses the first five original outer splits per dataset.
2. These models do **NOT** reproduce the pretrained/published PurIST and DeCAF clinical classifiers. Claims of 'RLR significantly better than PurIST or DeCAF' are unsupported; use 'RLR versus matched binary-pair model ablations inspired by prior methods'.
3. The experiments and BH family were planned after earlier analyses involving the same CompCancer datasets and an initial D/V pilot. Statistical findings are **exploratory** and should not be reported as a priori confirmatory evidence.
4. The MCP model emitted \`ncvreg: Model saturated; exiting...\` warnings on at least the Shipp cohort (visible in owner's logs). All three method predictions passed replay, but a saturation sensitivity check remains advisable before publication.
5. The owner-provided CSV \`original_paper_ARFF_SHA256_parity\` currently incorrectly says \`NOT_VERIFIED\` because the original pilot's standalone auditor hardcoded that field. This is **an audit metadata bug**, not a data-integrity failure. The local R runner independently verified SHA256 and all 11 sources. A revised audit script supports \`--original-data-root\` to rehash original and sharded ARFF and then write \`PASS\` into the inferential CSV **without re-fitting any models**.
6. No main manuscript data/results or original historical Git commits were changed. Maintain complete positive and negative comparison results.

## Suggested accurate manuscript wording (Supplement / optional ablation)

> In a controlled exploratory ablation on the 11 original CompCancer data matrices (verified against the frozen SHA256 provenance manifest), we evaluated the original RLR-20 model against two binary gene-pair logistic-regression variants: elastic net and MCP, inspired by previous pair-based probabilistic classifiers. Identical training-only pair selection (K=20, gene-disjoint) and five matched outer folds per dataset were used. RLR showed lower mean logloss (0.2379) than binary elastic net (0.2904) and binary MCP (0.2912). The corresponding Brier scores were 0.0626, 0.0774, and 0.0855. At the level of 11 dataset-wise fold means, exploratory exact paired sign-flip tests yielded BH-adjusted q=0.0391 for both logloss comparisons; the matched RLR-vs-binary-elastic-net Brier comparison yielded q=0.0234. These findings indicate benefits from the margin-aware ternary relation representation over binary pair encoding in this controlled setting, but do not establish superiority to complete published clinical PurIST or DeCAF implementations. 

## Rehash + rerun **audit only**, not models

From the separate Ubuntu checkout:

\`\`\`bash
cd /opt/2026/RLR/rlr-binary-original-arff-v2
git pull --ff-only
python3 scripts/audit_rlr_binary_outer5_v2.py \
  --input results/rlr_binary_outer5_original_arff_v2 \
  --out results/rlr_binary_outer5_original_arff_audit_v2 \
  --original-data-root /opt/2026/RLR/data
cat results/rlr_binary_outer5_original_arff_audit_v2/original_arff_provenance_attestation.txt
\`\`\`

No new training, data modifications, or extra dependencies. If the audit hash gate fails, preserve the original results and report the failure; do not force PASS.
