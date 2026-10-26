# Full original-paper 11-cohort matched-pairs audit — complete, provenance-limited

**Execution date:** 2026-10-26. **Status: ALL 11/11 complete; 470/470 outer folds; 2820 original/control fits + patient-level audits PASS.**

- Full GitHub Actions: https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40338028247
- Overall evidence ZIP artifact: `rlr-full11-overall-470-evidence-40338028247` (artifact ID `12661375940`, SHA256 of archive reported by GitHub `530491fd951c23ec3a3c09bfd1942a0fbcfb1b92b45f89c1ecaf669a10ea6073`).
- Original RLR code on `main` and manuscript lock untouched.
- These data are **public reconstructed Schliep filtered source matrices** from the **11 exact names in the original paper**, NOT independently SHA256-verified originals of the processed ARFF files. The Golub source has 1868 probes vs 1877 catalogue; source header annotations in Golub/West required positional sample IDs. Therefore **DO NOT call these results an exact reproduction of the original manuscript benchmark**.
- Original `fit_rlr_train_test()` used directly, original supervised switchBox pool, signed ternary rank-gap m10, strict gene-disjoint d2, alpha .5, 5-fold inner elastic-net CV/lambda.min. Same selected feature vector is averaged for REGULARIZED_EQUAL_VOTE, whose single scalar coefficient is fitted with identical glmnet alpha .5, five-fold inner CV and lambda.min; dummy all-zero column is checked to be inactive.
- Nine datasets: 10 repetitions × 5 CV folds, West/Yeoh: 2 × 5 CV folds; seed 42. **470 outer E folds for each K and both methods**.
- The 11 study-level differences, not individual repeated folds, are used for signed-rank analysis. Exact two-sided 2^11 sign enumeration used; descriptive Holm adjustment across K=10/20/50.

## Primary comparison: equal-vote logloss minus original RLR

**Positive** means RLR has lower (better) logloss. K is a *requested maximum*, not guaranteed number of actually selected disjoint gene pairs.

| Dataset | Folds | K10 ΔLL | K20 ΔLL | K50 ΔLL |
|---|---:|---:|---:|---:|
| alizadeh-2000-v1 | 50 | -0.04686 | -0.06316 | -0.06046 |
| armstrong-2002-v1 | 50 | +0.00340 | +0.00308 | +0.00278 |
| bittner-2000 | 50 | +0.10406 | +0.03697 | +0.07350 |
| chen-2002 | 50 | +0.03274 | +0.03572 | +0.03034 |
| golub-1999-v1 | 50 | +0.05436 | +0.04568 | +0.04012 |
| gordon-2002 | 50 | +0.00043 | -0.00341 | -0.00584 |
| pomeroy-2002-v1 | 50 | +0.22256 | +0.17330 | +0.18135 |
| shipp-2002-v1 | 50 | +0.04472 | +0.09318 | +0.17907 |
| singh-2002 | 50 | +0.02461 | +0.07163 | +0.07943 |
| west-2001 | 10 | -0.11776 | +0.06144 | +0.09369 |
| yeoh-2002-v1 | 10 | +0.05638 | +0.05638 | +0.05638 |

| K | Dataset-mean ΔLL | Dataset median ΔLL | RLR better datasets | Exact Wilcoxon 2-sided p | Holm across 3 K |
|---:|---:|---:|---:|---:|---:|
| 10 | +0.034422 | +0.032737 | 9/11 | 0.147461 | 0.147461 |
| 20 | +0.046438 | +0.045677 | 9/11 | 0.041992 | 0.083984 |
| 50 | +0.060943 | +0.056383 | 9/11 | 0.024414 | 0.073242 |

**Caution:** none of the exploratory family of 3 comparisons passes Holm-adjusted 0.05. A simpler two-sided exact sign test for 9/11 directions also gives p=0.06543. Avoid significance claims and post-hoc cherry-picking of K50. Original paper's existing evidence and its primary endpoint are distinct from this exploratory ablation.

## Secondary metric: Brier

The dataset-mean Brier difference (equal vote minus RLR, positive favors RLR) is **+0.000944 (6/11 RLR wins) at K10, +0.000833 (6/11) at K20, +0.005960 (7/11) at K50**. Unlike logloss, Brier direction is not consistently in favor of RLR across a strong majority of datasets. Do not claim universally improved calibration from this component comparison; assessing calibration properly also requires ECE/reliability analysis.

## Effective disjointness budget audit

All 470 folds are complete, but `apply_disjointness` can return fewer than K relations from its supervised pool. Examples (actual selected pairs across relevant folds):

- Yeoh-2002-v1: **7–10** selected for K10, K20, *and* K50; all three K give exactly identical predictions and summary statistics. This is a budget truncation, **not three independent supporting configurations**.
- Armstrong: K50 gives **25–35** pairs.
- Chen: K50 gives **24–29** pairs.
- Golub: K50 gives **30–45** pairs.
- Gordon: K50 gives **31–44** pairs.
- West: K50 gives **30–44** pairs.
- Remaining cohorts may also fall short of 50; full selected-pair counts are archived per fold in `method_metrics.csv`.

This truncation is shared exactly by original RLR and the equal-vote control within a fold, so it **does not violate the paired-comparison fairness**, but limits claims that the K50 result is based on 50 distinct gene-pair features.

## Final interpretation and finish line

**Established by the audited exploratory reconstruction:** original RLR with independently penalized relation coefficients has better *average logloss* than a one-coefficient, equally penalized vote score on 9 of 11 named paper cohorts for all K. This is useful descriptive support for the original manuscript's probabilistic relational model. However the per-dataset effect size varies, and Alizadeh is consistently negative.

**Not established:** superior predictive performance under Holm-adjusted significance, better Brier across most cohorts, proof that individual weights alone constitute a novel algorithm, exact parity with processed original ARFF, independent prospective biomolecular validation or a positive decision-tree component.

**Next required publication audit:** obtain original 11 processed ARFF files on the local Ubuntu research environment and independently compare the SHA256 values in `checksums/DATA_SHA256_EXPECTED.txt` with `scripts/rlr_original_arff_sha_gate_v1.py`. If all match, rerun/compare in an isolated environment using the genuine locked ARFF and manuscript's preprocessing pipeline before editing any results. If they do not match, label this analysis reconstructed-source exploratory and keep the locked original paper results separate. Either way, a journal submission can be evaluated on the merits of the *original integrated RLR method* without fabricating further extensions.

**Research decision:** no new trees, local logistic leaf models, ad-hoc per-pair thresholds, XOR-LSH branches or post-hoc parameter sweeping on exposed data. Focus on publication-ready scope and originality/novelty review for *Computer Methods and Programs in Biomedicine*.

**Evidence generated:** GitHub-hosted per-cohort archives for 11 datasets; independently recomputed E patient-level logloss/Brier; full `dataset_level.csv`, `summary_exploratory.csv`, `manifest.json`, and `full11_stats.log`. All run steps SUCCESS.
